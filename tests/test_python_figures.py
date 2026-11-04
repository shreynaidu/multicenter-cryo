#!/usr/bin/env python3
"""Regenerate the Python figures and compare them with the baseline.

Numbers have to match everywhere; figure files only on the baseline's platform.
"""

import argparse
import hashlib
import json
import math
import platform
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

PROJECT = Path(__file__).resolve().parents[1]


def default_baseline():
    private = PROJECT / "tests/baseline/python_figures"
    return (
        private
        if private.is_dir()
        else PROJECT / "tests/public_baseline/python_figures"
    )


def on_reference_platform():
    return platform.system() == "Darwin" and platform.machine() == "arm64"


def normalized_vector(path):
    content = path.read_bytes()
    if path.suffix == ".eps":
        return re.sub(
            rb"^%%(?:Title|CreationDate):[^\n]*(?:\n|$)",
            b"",
            content,
            flags=re.MULTILINE,
        )
    tree = ET.fromstring(content)
    dates = "{http://purl.org/dc/elements/1.1/}date"
    for parent in tree.iter():
        for child in list(parent):
            if child.tag == dates:
                parent.remove(child)
    identifiers = {
        element.attrib["id"]: f"id_{index}"
        for index, element in enumerate(tree.iter())
        if "id" in element.attrib
    }
    for element in tree.iter():
        for attribute, value in list(element.attrib.items()):
            if attribute == "id":
                element.attrib[attribute] = identifiers[value]
            else:
                for original, replacement in identifiers.items():
                    value = value.replace(f"url(#{original})", f"url(#{replacement})")
                    if value == "#" + original:
                        value = "#" + replacement
                element.attrib[attribute] = value
    return ET.tostring(tree)


def compare_csv(before, after):
    expected = pd.read_csv(before, float_precision="round_trip")
    actual = pd.read_csv(after, float_precision="round_trip")
    if list(expected.columns) != list(actual.columns) or expected.shape != actual.shape:
        raise AssertionError(f"Table structure differs: {before.name}")
    maximum = 0.0
    count = 0
    for column in expected:
        if pd.api.types.is_numeric_dtype(expected[column]):
            left = expected[column].to_numpy(dtype=float)
            right = actual[column].to_numpy(dtype=float)
            if not np.allclose(left, right, rtol=1e-12, atol=1e-12, equal_nan=True):
                raise AssertionError(f"Numerical inputs differ: {before.name}/{column}")
            finite = np.isfinite(left) & np.isfinite(right)
            maximum = max(
                maximum, np.max(np.abs(left[finite] - right[finite]), initial=0)
            )
            count += int(finite.sum())
        elif not expected[column].fillna("").equals(actual[column].fillna("")):
            raise AssertionError(f"Labels differ: {before.name}/{column}")
    return {
        "file": before.name,
        "numeric_cells": count,
        "maximum_absolute_difference": float(maximum),
    }


def compare_json_values(expected, actual, location="artists"):
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise AssertionError(f"JSON fields differ: {location}")
        for key in expected:
            compare_json_values(expected[key], actual[key], f"{location}/{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise AssertionError(f"JSON lengths differ: {location}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare_json_values(left, right, f"{location}/{index}")
    elif isinstance(expected, (int, float)):
        if not math.isclose(expected, actual, rel_tol=1e-12, abs_tol=1e-12):
            raise AssertionError(f"Plotted number differs: {location}")
    elif expected != actual:
        raise AssertionError(f"Figure annotation differs: {location}")


def compare_figure(before, after):
    if not after.is_file():
        raise AssertionError(f"Missing figure: {after.name}")
    if before.suffix in (".png", ".tiff"):
        with Image.open(before) as expected, Image.open(after) as actual:
            if (
                expected.size != actual.size
                or expected.mode != actual.mode
                or expected.tobytes() != actual.tobytes()
            ):
                raise AssertionError(f"Figure pixels differ: {before.name}")
            if expected.info.get("dpi") != actual.info.get("dpi"):
                raise AssertionError(f"Figure DPI differs: {before.name}")
    elif before.suffix == ".pdf":
        if not shutil.which("pdftoppm"):
            raise RuntimeError("PDF comparison requires pdftoppm (Poppler)")

        def render(path):
            return subprocess.run(
                ["pdftoppm", "-r", "300", str(path)], check=True, capture_output=True
            ).stdout

        if render(before) != render(after):
            raise AssertionError(f"PDF pixels differ: {before.name}")
    elif normalized_vector(before) != normalized_vector(after):
        raise AssertionError(f"Vector figure differs: {before.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT
        / "outputs/current/analysis_data/deidentified_corrected_analysis_dataset.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT / "outputs/python_figures_reproducibility",
    )
    parser.add_argument("--compare-only", type=Path)
    parser.add_argument("--analysis", type=Path, default=PROJECT / "outputs/current")
    parser.add_argument("--library", type=Path, default=PROJECT / "renv/library")
    parser.add_argument("--baseline", type=Path, default=default_baseline())
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--figures-report-only", action="store_true")
    args = parser.parse_args()
    baseline = args.baseline
    report_only = args.figures_report_only or not on_reference_platform()
    manifest = json.loads((baseline / "manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        path = baseline / entry["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise AssertionError(f"Python baseline changed: {entry['path']}")
    output = args.compare_only or args.output
    if args.compare_only is None:
        subprocess.run(
            [
                sys.executable,
                str(PROJECT / "analysis/scripts/generate_publication_figures.py"),
                "--input",
                str(args.input),
                "--output",
                str(output),
            ],
            check=True,
        )
    if not (output / "run_manifest.json").exists():
        raise AssertionError("Incomplete Python figure run")
    tables = [
        compare_csv(path, output / "numerical" / path.name)
        for path in sorted((baseline / "numerical").glob("*.csv"))
    ]
    compare_json_values(
        json.loads(
            (baseline / "numerical/figure_artists.json").read_text(encoding="utf-8")
        ),
        json.loads(
            (output / "numerical/figure_artists.json").read_text(encoding="utf-8")
        ),
    )
    figures = sorted((baseline / "results").iterdir())
    figure_differences = []
    for path in figures:
        try:
            compare_figure(path, output / path.name)
        except AssertionError as difference:
            figure_differences.append(str(difference))
    if figure_differences:
        message = "Figures differ from the baseline: " + "; ".join(figure_differences)
        if not report_only:
            raise AssertionError(message)
        print(message, file=sys.stderr)
    subprocess.run(
        [
            args.rscript,
            "--vanilla",
            str(PROJECT / "analysis/scripts/audit_python_figure_statistics.R"),
            "--figures",
            str(output),
            "--analysis",
            str(args.analysis),
            "--library",
            str(args.library),
        ],
        check=True,
    )
    consistency = json.loads(
        (output / "statistical_consistency.json").read_text(encoding="utf-8")
    )
    report = {
        "statistical_consistency": consistency,
        "status": "passed",
        "baseline_files_verified": len(manifest["files"]),
        "figure_files_compared": len(figures),
        "figure_differences_reported": figure_differences,
        "figures_required_to_match": not report_only,
        "tables_compared": len(tables),
        "numeric_cells_compared": sum(row["numeric_cells"] for row in tables),
        "maximum_absolute_difference": max(
            row["maximum_absolute_difference"] for row in tables
        ),
        "figure_rule": "All raster pixels/DPI and all PDF pixels at 300 dpi exact; EPS/SVG content exact after removing creation/title metadata and normalizing SVG random IDs",
        "artist_annotations_and_coordinates": "matched",
        "checks": tables,
    }
    (output / "reproducibility.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"PASS: {len(figures)} figure files, {len(tables)} tables, {report['numeric_cells_compared']} numerical cells; maximum difference {report['maximum_absolute_difference']}."
    )


if __name__ == "__main__":
    main()
