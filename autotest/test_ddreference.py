"""
Tests CLN drawdown after the drawdown reference is reset with DDREFERENCE.

test-reg/test053_clnmodel is run with CLN drawdown saved and the drawdown
reference reset at the end of stress period 1 (one time step).

Cases:
  - test_cln_ddreference : CLN drawdown in later stress periods is the CLN head
                           at the end of stress period 1 minus the current CLN
                           head.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from test_regression import read_arrays

ROOT = Path(__file__).resolve().parents[1]
EXT = ".exe" if sys.platform == "win32" else ""
TARGETS = ["mfusg", "mfusgdbl"]
MODEL = ROOT / "test-reg" / "test053_clnmodel"
NAM = "DEN4Tc.nam"
CLN_HEADS = "DEN4Tc_cln.hds"
CLN_DDN = "DEN4Tc_cln.ddn"
ATOL = 1e-3

# values at or above this are inactive or dry
HNOFLO = 1e20


def setup_model(ws):
    shutil.copytree(MODEL, ws)

    # save CLN drawdown (ICLNDD, the fifth value) on unit 35
    fpth = ws / "DEN4Tc.cln"
    lines = fpth.read_text().splitlines()
    assert lines[0].split() == ["66", "-378", "24", "34", "0", "0", "378", "66"]
    lines[0] = " 66 -378 24 34 35 0 378 66"
    fpth.write_text("\n".join(lines) + "\n")

    fpth = ws / NAM
    fpth.write_text(fpth.read_text() + f"DATA(BINARY) 35 {CLN_DDN}\n")

    # reset the drawdown reference at the end of stress period 1
    fpth = ws / "DEN4Tc.oc"
    text = fpth.read_text()
    old = "PERIOD 1  STEP 1\n"
    assert text.count(old) == 1
    fpth.write_text(text.replace(old, "PERIOD 1  STEP 1 DDREFERENCE\n"))


@pytest.mark.parametrize("target", TARGETS)
def test_cln_ddreference(tmp_path, target):
    exe = ROOT / "bin" / f"{target}{EXT}"
    if not exe.is_file():
        pytest.skip(f"{exe.name} is not built")

    ws = tmp_path / "model"
    setup_model(ws)
    proc = subprocess.run([str(exe), NAM], cwd=ws, capture_output=True, text=True)
    assert "Normal termination" in proc.stdout, proc.stdout + proc.stderr

    heads = {key[1:3]: arr for key, arr in read_arrays(ws / CLN_HEADS)}
    ddn = [(key[1:3], arr) for key, arr in read_arrays(ws / CLN_DDN)]
    href = heads[(1, 1)]
    nchecked = 0
    for (kstp, kper), arr in ddn:
        if kper < 2:
            continue
        head = heads[(kstp, kper)]
        active = (np.abs(href) < HNOFLO) & (np.abs(head) < HNOFLO)
        diff = np.max(np.abs(arr[active] - (href[active] - head[active])))
        assert diff <= ATOL, (
            f"CLN drawdown in step {kstp} period {kper} differs by {diff}"
        )
        nchecked += 1
    assert nchecked > 0, "no CLN drawdown saved after stress period 1"
