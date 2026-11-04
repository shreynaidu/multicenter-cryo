#!/usr/bin/env python3
"""Generate the supplied Python figures, separately from the R figures."""

import argparse
import importlib.metadata
import json
import os
import re
import shutil
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_DIR / ".cache/figure_matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from figure_inputs import SIZE_GROUPS, STAGE_GROUPS, load_figure_inputs, sha256
from figure_panels import (
    COMPOSITE_FIGURES,
    STANDALONE_FIGURES,
    build_composite,
    build_standalone,
    configure_style,
    figure_artist_snapshot,
)
from matplotlib import font_manager
from PIL import Image

OUTPUT_MARKER = ".python-figures"
SOURCE_MODULES = (
    "generate_publication_figures.py",
    "figure_estimators.py",
    "figure_inputs.py",
    "figure_panels.py",
)


def save_figure_formats(figure, output_dir, stem):
    paths = {}
    for extension in ("pdf", "eps", "svg", "png"):
        name = f"{stem}_600dpi.png" if extension == "png" else f"{stem}.{extension}"
        path = output_dir / name
        options = {"dpi": 600} if extension == "png" else {}
        if extension == "pdf":
            options["metadata"] = {"CreationDate": None, "ModDate": None}
        elif extension == "svg":
            options["metadata"] = {"Date": None}
        figure.savefig(path, format=extension, facecolor="white", **options)
        if extension == "eps":
            # Strip paths and timestamps so runs are byte-identical.
            content = re.sub(
                rb"^%%(?:Title|CreationDate):.*\n",
                b"",
                path.read_bytes(),
                flags=re.MULTILINE,
            )
            path.write_bytes(content)
        paths[extension] = path
    tiff = output_dir / f"{stem}_600dpi.tiff"
    with Image.open(paths["png"]) as image:
        image.convert("RGB").save(
            tiff, format="TIFF", compression="tiff_lzw", dpi=(600, 600)
        )
    paths["tiff"] = tiff
    return paths


def render_figures(inputs, output_dir):
    artists = {}
    for stem, alias, kind, height, bounds, letter in STANDALONE_FIGURES:
        figure = build_standalone(inputs, kind, height, bounds, letter)
        artists[stem] = figure_artist_snapshot(figure)
        artists[alias] = artists[stem]
        for path in save_figure_formats(figure, output_dir, stem).values():
            shutil.copyfile(path, path.with_name(path.name.replace(stem, alias, 1)))
        plt.close(figure)
    for supplementary, stem in COMPOSITE_FIGURES:
        figure = build_composite(inputs, supplementary)
        artists[stem] = figure_artist_snapshot(figure)
        save_figure_formats(figure, output_dir, stem)
        plt.close(figure)
    return artists


def export_numerical_inputs(inputs, output_dir):
    numerical = output_dir / "numerical"
    numerical.mkdir()

    def write_curve(table, name):
        table.to_csv(numerical / name, index=False, float_format="%.17g")

    write_curve(inputs.records, "figure_analysis_dataset.csv")
    for index, group in enumerate(["overall"] + SIZE_GROUPS):
        label = "overall" if index == 0 else f"size_{index}"
        for endpoint in ("os", "lpfs"):
            write_curve(inputs.curves[(group, endpoint)], f"{label}_{endpoint}_km.csv")
        for cause in (1, 2):
            write_curve(
                inputs.curves[(group, f"cif_{cause}")], f"{label}_cause_{cause}_cif.csv"
            )
    for stage in STAGE_GROUPS:
        write_curve(
            inputs.stage_curves[stage], f"{stage.replace(' ', '_')}_stage_steps.csv"
        )
    inputs.overrides.to_csv(numerical / "followup_override_audit.csv", index=False)
    summary = {
        "logrank_p": inputs.logrank_p,
        "median_os_months": inputs.median_os,
        "retained_reported_gray_p": inputs.reported_gray_p,
        "stage_curve_source": "Retained image-derived coordinates; not estimated from patient stage data",
    }
    (numerical / "computed_annotations.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


def prepare_output_directory(output_dir, input_path):
    if output_dir.is_symlink():
        raise ValueError("Output directory must not be a symbolic link")
    output_dir = output_dir.resolve()
    if output_dir == PROJECT_DIR or output_dir in input_path.resolve().parents:
        raise ValueError("Output directory overlaps project/source input")
    if output_dir.exists() and any(output_dir.iterdir()):
        if not (output_dir / OUTPUT_MARKER).is_file():
            raise ValueError(
                "Refusing a populated output directory without the figure marker"
            )
        previous = output_dir.with_name(output_dir.name + ".previous")
        if previous.exists():
            if not (previous / OUTPUT_MARKER).is_file():
                raise ValueError(
                    f"Refusing to replace an unmarked directory: {previous}"
                )
            shutil.rmtree(previous)
        output_dir.rename(previous)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / OUTPUT_MARKER).write_text(
        "Verified supplied Python publication figures.\n", encoding="utf-8"
    )
    return output_dir


def verify_environment():
    versions = {}
    lockfile = PROJECT_DIR / "requirements-figures.lock"
    for line in lockfile.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        package, expected = line.split("==", 1)
        versions[package] = importlib.metadata.version(package)
        if versions[package] != expected:
            raise RuntimeError(
                f"Package version differs: {package}; restore requirements-figures.lock"
            )
    return versions


def write_manifest(output_dir, input_hash, reference_dir, package_versions):
    font_path = Path(
        font_manager.findfont(font_manager.FontProperties(family=["sans-serif"]))
    )
    manifest = {
        "font": {"file": font_path.name, "sha256": sha256(font_path)},
        "input_sha256": input_hash,
        "python": sys.version.split()[0],
        "packages": package_versions,
        "source_sha256": {name: sha256(SCRIPT_DIR / name) for name in SOURCE_MODULES},
        "reference": json.loads(
            (reference_dir / "provenance.json").read_text(encoding="utf-8")
        ),
        "note": "Keeps values from the supplied figures; see statistical_consistency.csv.",
        "files": [
            {"path": path.relative_to(output_dir).as_posix(), "sha256": sha256(path)}
            for path in sorted(output_dir.rglob("*"))
            if path.is_file() and not path.name.startswith(".")
        ],
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_DIR
        / "outputs/current/analysis_data/deidentified_corrected_analysis_dataset.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_DIR / "outputs/python_publication_figures",
    )
    parser.add_argument(
        "--reference", type=Path, default=PROJECT_DIR / "analysis/figure_reference"
    )
    parser.add_argument(
        "--followup-overrides",
        type=Path,
        default=PROJECT_DIR / "data/reference/retained_followup_days.csv",
    )
    parser.add_argument("--check-environment", action="store_true")
    args = parser.parse_args()
    package_versions = verify_environment()
    if args.check_environment:
        print("Recorded Python figure packages are available.")
        return
    input_hash = sha256(args.input)
    inputs = load_figure_inputs(args.input, args.reference, args.followup_overrides)
    configure_style()
    output_dir = prepare_output_directory(args.output, args.input)
    artists = render_figures(inputs, output_dir)
    export_numerical_inputs(inputs, output_dir)
    (output_dir / "numerical/figure_artists.json").write_text(
        json.dumps(artists, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    if sha256(args.input) != input_hash:
        raise RuntimeError("Input CSV changed during generation")
    write_manifest(output_dir, input_hash, args.reference, package_versions)
    print(f"Figures written to {output_dir}")


if __name__ == "__main__":
    main()
