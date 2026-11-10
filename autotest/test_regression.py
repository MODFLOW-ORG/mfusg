"""
Regression tests of the models in examples/, test-reg/, and test-cmp/.

Each model is run with each available executable in bin/ and with the same
precision of the reference release, which is downloaded from the mfusg GitHub
releases once per session.

Cases:
  - test_regression : both runs terminate normally, the list file budgets agree
                      within 0.1 percent of the total flow, and the binary head
                      and drawdown arrays agree within 0.01.
"""

import shutil
import stat
import subprocess
import sys
import urllib.request
import warnings
import zipfile
from pathlib import Path

import flopy
import numpy as np
import pytest
from filelock import FileLock
from modflow_devtools.ostags import get_ostag

ROOT = Path(__file__).resolve().parents[1]
EXT = ".exe" if sys.platform == "win32" else ""
TARGETS = ["mfusg", "mfusgdbl"]
GROUPS = ["examples", "test-reg", "test-cmp"]
ARRAY_SUFFIXES = {".hds", ".hed", ".ddn", ".clnhds"}
HEAD_ATOL = 1e-2
REF_VERSION = "1.5.2"
REF_URL = "https://github.com/MODFLOW-ORG/mfusg/releases/download/v{}/{}.zip"

# budget terms may differ by this fraction of the total flow
BUDGET_RTOL = 1e-3

# initial heads are read from a single precision binary file
SKIP = {("mfusgdbl", "test-reg/test019_complex_problem")}


def get_namefiles():
    namefiles = []
    for group in GROUPS:
        for fpth in sorted((ROOT / group).rglob("*")):
            if fpth.suffix.lower() == ".nam" and ".cmp" not in fpth.parent.name:
                namefiles.append(fpth)
    return namefiles


NAMEFILES = get_namefiles()


@pytest.fixture(scope="session")
def ref_bin(tmp_path_factory, worker_id):
    root = tmp_path_factory.getbasetemp()
    if worker_id != "master":
        root = root.parent
    bindir = root / "ref_bin"
    with FileLock(str(root / "ref_bin.lock")):
        if not (bindir / f"mfusg{EXT}").is_file():
            bindir.mkdir(exist_ok=True)
            url = REF_URL.format(REF_VERSION, get_ostag())
            fpth = root / "ref_bin.zip"
            try:
                urllib.request.urlretrieve(url, fpth)
            except OSError as err:
                pytest.skip(f"cannot download {url}: {err}")
            with zipfile.ZipFile(fpth) as zf:
                zf.extractall(bindir)
            for exe in bindir.iterdir():
                exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return bindir


def read_namefile(fpth):
    files = {}
    for line in fpth.read_text().splitlines():
        tokens = line.split()
        if len(tokens) >= 3 and not tokens[0].startswith("#"):
            files.setdefault(tokens[0].upper(), []).append(tokens[2])
    return files


def run_model(exe, namefile, ws):
    shutil.copytree(
        namefile.parent, ws, ignore=shutil.ignore_patterns("*.cmp", "output")
    )
    proc = subprocess.run(
        [str(exe), namefile.name], cwd=ws, capture_output=True, text=True
    )
    return "Normal termination" in proc.stdout, proc.stdout + proc.stderr


def read_header(buf, pos, nbyte):
    ints = np.frombuffer(buf, "<i4", 2, pos)
    pos += 8
    times = np.frombuffer(buf, f"<f{nbyte}", 2, pos)
    pos += 2 * nbyte
    text = buf[pos : pos + 16]
    dims = np.frombuffer(buf, "<i4", 3, pos + 16)
    return ints, times, text, dims, pos + 28


def is_header(buf, pos, nbyte):
    if pos + 36 + 2 * nbyte > len(buf):
        return False
    ints, _, text, _, _ = read_header(buf, pos, nbyte)
    return (
        ints[0] > 0
        and ints[1] > 0
        and text.isascii()
        and text.strip() != b""
        and text.decode().isprintable()
    )


def read_arrays(fpth):
    """Read the records of a structured or unstructured binary array file."""
    buf = fpth.read_bytes()
    for nbyte in (4, 8):
        if is_header(buf, 0, nbyte):
            break
    else:
        raise ValueError(f"{fpth.name} is not a binary array file")
    records = []
    pos = 0
    while pos < len(buf):
        ints, times, text, dims, pos = read_header(buf, pos, nbyte)
        # structured records hold ncol*nrow values, unstructured nodes n1 to n2
        for size in (int(dims[0]) * int(dims[1]), int(dims[1]) - int(dims[0]) + 1):
            end = pos + size * nbyte
            if size > 0 and (end == len(buf) or is_header(buf, end, nbyte)):
                break
        else:
            raise ValueError(f"cannot size record {len(records) + 1} in {fpth.name}")
        key = (text.decode().strip(), int(ints[0]), int(ints[1]), int(dims[2]))
        records.append((key, np.frombuffer(buf, f"<f{nbyte}", size, pos)))
        pos = end
    return records


def compare_budgets(fpth, fpth_ref):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lst = flopy.utils.MfusgListBudget(fpth)
        lst_ref = flopy.utils.MfusgListBudget(fpth_ref)
        budgets = [
            (lst.get_cumulative(), lst_ref.get_cumulative()),
            (lst.get_incremental(), lst_ref.get_incremental()),
        ]
    for bud, bud_ref in budgets:
        assert bud.shape == bud_ref.shape, f"{fpth.name} budget lengths differ"
        # percent discrepancy is derived from IN-OUT, which is compared
        names = [
            n
            for n in bud_ref.dtype.names
            if n.isupper() and n != "PERCENT_DISCREPANCY"
        ]
        for i in range(bud_ref.shape[0]):
            total = max(abs(bud_ref["TOTAL_IN"][i]), abs(bud_ref["TOTAL_OUT"][i]))
            for name in names:
                diff = abs(bud[name][i] - bud_ref[name][i])
                assert diff <= BUDGET_RTOL * max(total, 1.0), (
                    f"{fpth.name} {name} differs by {diff} in entry {i + 1}"
                )


def compare_arrays(fpth, fpth_ref):
    if fpth_ref.stat().st_size == 0:
        assert fpth.stat().st_size == 0, f"{fpth.name} is not empty"
        return
    records = read_arrays(fpth)
    records_ref = read_arrays(fpth_ref)
    assert len(records) == len(records_ref), f"{fpth.name} record counts differ"
    for (key, arr), (key_ref, arr_ref) in zip(records, records_ref):
        assert key == key_ref, f"{fpth.name} record {key} is not {key_ref}"
        diff = np.max(np.abs(arr - arr_ref))
        assert diff <= HEAD_ATOL, f"{fpth.name} {key} differs by {diff}"


@pytest.mark.parametrize(
    "namefile",
    NAMEFILES,
    ids=[p.parent.relative_to(ROOT).as_posix() for p in NAMEFILES],
)
@pytest.mark.parametrize("target", TARGETS)
def test_regression(tmp_path, ref_bin, target, namefile):
    exe = ROOT / "bin" / f"{target}{EXT}"
    if not exe.is_file():
        pytest.skip(f"{exe.name} is not built")
    if (target, namefile.parent.relative_to(ROOT).as_posix()) in SKIP:
        pytest.skip(f"{namefile.name} cannot be run with {target}")

    ws, ws_ref = tmp_path / "dev", tmp_path / "ref"
    success, output = run_model(exe, namefile, ws)
    assert success, output

    success, output = run_model(ref_bin / f"{target}{EXT}", namefile, ws_ref)
    assert success, output

    files = read_namefile(namefile)
    lst = files["LIST"][0]
    compare_budgets(ws / lst, ws_ref / lst)
    for fname in files.get("DATA(BINARY)", []):
        if Path(fname).suffix.lower() in ARRAY_SUFFIXES:
            compare_arrays(ws / fname, ws_ref / fname)
