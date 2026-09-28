"""The committed type stub, brille/_brille.pyi, matches the compiled module.

The stub ships in the wheel and the documentation's API reference is read from
it, so a binding change must come with a regenerated stub:
`python devtools/make_stub.py`.
"""
import difflib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.skipif(importlib.util.find_spec("pybind11_stubgen") is None, reason="needs pybind11-stubgen")
@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="C++ integer types differ between platforms; checked on Linux")
def test_committed_stub_is_current(tmp_path):
    committed = ROOT / "brille" / "_brille.pyi"
    if not committed.exists():
        pytest.skip(f"{committed} is not available")
    sys.path.insert(0, str(ROOT / "devtools"))
    try:
        from make_stub import generate
    finally:
        sys.path.pop(0)
    current = generate(tmp_path).read_text().splitlines()
    diff = list(difflib.unified_diff(committed.read_text().splitlines(), current, "committed", "generated", lineterm=""))
    assert not diff, "brille/_brille.pyi is stale; run `python devtools/make_stub.py`:\n" + "\n".join(diff[:40])
