"""Regenerate brille/_brille.pyi, the type stub for the compiled module.

Run it after changing the bindings, with the changed brille installed:

    python devtools/make_stub.py

The stub ships in the wheel, for type checkers and editors, and the
documentation's API reference is read from it. wrap/tests/test_24_type_stub.py
fails when it no longer matches the module.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

STUBGEN = "pybind11-stubgen==3.0.0"   # pinned: other versions format stubs differently
STUB = Path(__file__).resolve().parents[1] / "brille" / "_brille.pyi"


def generate(directory):
    """Write the stub for the installed brille._brille below directory; return its path"""
    # run outside the source tree, whose brille/ has no compiled module
    subprocess.run([sys.executable, "-m", "pybind11_stubgen", "brille._brille", "-o", str(directory)],
                   cwd=directory, check=True, capture_output=True)
    return Path(directory) / "brille" / "_brille.pyi"


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as directory:
        shutil.copyfile(generate(directory), STUB)
    print(f"wrote {STUB}")
