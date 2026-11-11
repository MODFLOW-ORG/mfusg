"""
Performance tests of the PCGU solver against the reference release.

Each case builds a 3-layer, 200 by 200 model with FloPy, runs it with
bin/mfusg and the reference release three times each, and compares the
fastest runs. Run with `pixi run perf`; the file is not collected by
`pixi run test` because timings need a quiet machine.

Cases:
  - bcgs_milu     : BiCGSTAB with MILU0 and a convertible top layer that dries,
                    no more than 10 percent slower than the reference.
  - cg_ilu_rcm    : CG with ILU0, symmetric scaling, and RCM ordering, no more
                    than 10 percent slower than the reference.
  - bcgs_inactive : half of each layer inactive, solved with and without the
                    active node option, which must not be slower than the
                    reference.
"""

import subprocess
import sys
import time
from pathlib import Path

import flopy
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from test_regression import EXT, ROOT, ref_bin  # noqa: E402, F401

NLAY, NROW, NCOL = 3, 200, 200
NRUN = 3
SLOWDOWN = 1.10
HEAD_ATOL = 1e-2

CASES = {
    "bcgs_milu": dict(
        laytyp=[1, 0, 0], inactive=False, clin="bcgs", ipc=3, iscl=0, iord=0
    ),
    "cg_ilu_rcm": dict(
        laytyp=[0, 0, 0], inactive=False, clin="cg", ipc=2, iscl=1, iord=1
    ),
    "bcgs_inactive": dict(
        laytyp=[1, 0, 0], inactive=True, clin="bcgs", ipc=3, iscl=0, iord=0
    ),
}

RESULTS = []


def build(ws, case, iord):
    """Write a case to ws, with iord overriding the case ordering."""
    rng = np.random.default_rng(20261111)
    name = "perf"
    m = flopy.mfusg.MfUsg(modelname=name, model_ws=str(ws), structured=True)
    flopy.modflow.ModflowDis(
        m,
        nlay=NLAY,
        nrow=NROW,
        ncol=NCOL,
        delr=50.0,
        delc=50.0,
        top=100.0,
        botm=[60.0, 20.0, -20.0],
        nper=1,
        steady=True,
    )
    ibound = np.ones((NLAY, NROW, NCOL), dtype=int)
    if case["inactive"]:
        i, j = np.meshgrid(np.arange(NROW), np.arange(NCOL), indexing="ij")
        outside = (i - NROW / 2) ** 2 + (j - NCOL / 2) ** 2 > (0.4 * NROW) ** 2
        ibound[:, outside] = 0
    ibound[:, :, 0] = np.where(ibound[:, :, 0] != 0, -1, 0)
    if case["inactive"]:
        ibound[:, NROW // 2, int(0.1 * NCOL) + 1] = -1
    flopy.mfusg.MfUsgBas(m, ibound=ibound, strt=90.0)
    hk = 10.0 ** rng.normal(1.0, 0.5, (NLAY, NROW, NCOL))
    flopy.mfusg.MfUsgLpf(m, laytyp=case["laytyp"], hk=hk, vka=0.1 * hk)
    flopy.modflow.ModflowRch(m, rech=1.0e-3)
    wells = []
    for k in range(10):
        r, c = rng.integers(20, NROW - 20), rng.integers(20, NCOL - 20)
        if ibound[2, r, c] > 0:
            wells.append([2, r, c, -2000.0])
    flopy.mfusg.MfUsgWel(m, stress_period_data={0: wells})
    flopy.mfusg.MfUsgSms(
        m,
        hclose=1e-4,
        hiclose=1e-5,
        mxiter=200,
        iter1=500,
        iprsms=1,
        linmeth=2,
        clin=case["clin"],
        ipc=case["ipc"],
        iscl=case["iscl"],
        iord=iord,
        rclosepcgu=1e-3,
        relaxpcgu=0.97,
    )
    flopy.mfusg.MfUsgOc(m, stress_period_data={(0, 0): ["save head"]})
    m.write_input()
    return ws / f"{name}.nam"


def run(exe, nam):
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(exe), nam.name], cwd=nam.parent, capture_output=True, text=True
    )
    elapsed = time.perf_counter() - t0
    assert "Normal termination" in proc.stdout, proc.stdout[-2000:]
    return elapsed


def fastest(exes, nams):
    """Alternate the runs so machine load affects each executable alike."""
    times = {k: [] for k in exes}
    for _ in range(NRUN):
        for k, exe in exes.items():
            times[k].append(run(exe, nams[k]))
    return {k: min(v) for k, v in times.items()}


def heads(nam):
    return flopy.utils.HeadFile(nam.with_suffix(".hds")).get_data()


@pytest.mark.parametrize("name", list(CASES))
def test_performance(name, ref_bin, tmp_path):
    exe = ROOT / "bin" / f"mfusg{EXT}"
    if not exe.is_file():
        pytest.skip(f"{exe} not built")
    case = CASES[name]
    nams = {
        "ref": build(tmp_path / "ref", case, case["iord"]),
        "new": build(tmp_path / "new", case, case["iord"]),
    }
    exes = {"ref": ref_bin / f"mfusg{EXT}", "new": exe}
    if case["inactive"]:
        nams["active"] = build(tmp_path / "active", case, case["iord"] + 10)
        exes["active"] = exe
    t = fastest(exes, nams)

    href = heads(nams["ref"])
    for k in nams:
        if k == "ref":
            continue
        h = heads(nams[k])
        wet = (np.abs(href) < 1e20) & (np.abs(h) < 1e20)
        assert np.allclose(h[wet], href[wet], atol=HEAD_ATOL), f"{k} heads"

    RESULTS.append((name, t))
    for k in t:
        if k != "ref":
            print(f"\n{name} {k}: {t[k]:.2f} s, reference {t['ref']:.2f} s")
    assert t["new"] <= SLOWDOWN * t["ref"], f"{name} slower than reference"
    if "active" in t:
        assert t["active"] <= SLOWDOWN * t["ref"], f"{name} active slower"


def teardown_module():
    if not RESULTS:
        return
    print("\n\n| Case | Reference (s) | This build (s) | Ratio | Active nodes (s) |")
    print("|---|---|---|---|---|")
    for name, t in RESULTS:
        active = f"{t['active']:.2f}" if "active" in t else ""
        print(
            f"| {name} | {t['ref']:.2f} | {t['new']:.2f} "
            f"| {t['new'] / t['ref']:.2f} | {active} |"
        )
