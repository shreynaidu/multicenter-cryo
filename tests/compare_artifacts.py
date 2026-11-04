"""List baseline figures whose regenerated copy differs.

PDFs are compared as 300-dpi pixels because their bytes include a creation date.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

FIGURE_SUFFIXES = {".png", ".tiff", ".tif", ".svg", ".pdf"}


def canonical_bytes(path):
    if path.suffix.lower() != ".pdf":
        return path.read_bytes()
    if not shutil.which("pdftoppm"):
        raise RuntimeError("PDF comparison requires pdftoppm (Poppler) on PATH")
    return subprocess.run(
        ["pdftoppm", "-r", "300", str(path)], check=True, capture_output=True
    ).stdout


def compare_figures(baseline, actual):
    frozen = [
        path
        for path in sorted(baseline.rglob("*"))
        if path.suffix.lower() in FIGURE_SUFFIXES
    ]
    if not frozen:
        raise AssertionError("No baseline figures found")
    differences = []
    for before in frozen:
        relative = before.relative_to(baseline)
        after = actual / relative
        if not after.is_file() or canonical_bytes(before) != canonical_bytes(after):
            differences.append(relative.as_posix())
    return {"figures_compared": len(frozen), "differences": differences}


if __name__ == "__main__":
    print(json.dumps(compare_figures(Path(sys.argv[1]), Path(sys.argv[2]))))
