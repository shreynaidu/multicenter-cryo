"""Load and check the Python figures' inputs and fit their curves.

The follow-up overrides, stage curves and Gray P are values from the supplied figures,
kept in checksummed reference files. Overrides only touch an in-memory copy.
"""

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from figure_estimators import (
    cumulative_incidence_curve,
    kaplan_meier_curve,
    logrank_p_value,
    risk_counts,
    survival_median,
)

SIZE_GROUPS = ["<2.0 cm", "2.0 to <3.0 cm", "3.0 to <4.0 cm", "≥4.0 cm"]
STAGE_GROUPS = ["Stage I-II", "Stage III-IV", "Unknown stage"]
ENDPOINT_TIMES = ("os", "cr", "lpfs")
PLOTTED_COLUMNS = [
    "record_id",
    "tumor_size_group",
    "cr_status",
    "os_event",
    "lpfs_event",
] + [f"{endpoint}_time_days" for endpoint in ENDPOINT_TIMES]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plotted_columns_digest(records):
    # Other columns have decimals whose last digit can differ across platforms.
    rows = []
    for row in records[PLOTTED_COLUMNS].itertuples(index=False):
        record_id, group, cr_status, os_event, lpfs_event, *days = row
        whole_days = ",".join(f"{value:.6f}" for value in days)
        rows.append(
            f"{record_id},{group},{cr_status},{os_event},{lpfs_event},{whole_days}"
        )
    return hashlib.sha256("\n".join(rows).encode("utf-8")).hexdigest()


@dataclass
class FigureInputs:
    records: pd.DataFrame
    curves: dict
    stage_curves: dict
    stage_risks: dict
    logrank_p: dict
    median_os: tuple
    reported_gray_p: float
    overrides: pd.DataFrame


def verify_reference_files(provenance, reference_dir, overrides_path):
    if sha256(overrides_path) != provenance["followup_overrides_sha256"]:
        raise ValueError("Retained follow-up overrides changed")
    for name, expected in provenance["files"].items():
        if sha256(reference_dir / name) != expected:
            raise ValueError(f"Retained figure reference changed: {name}")


def validate_records(records, provenance):
    if set(PLOTTED_COLUMNS) - set(records):
        raise ValueError("Missing figure input columns")
    if plotted_columns_digest(records) != provenance["plotted_columns_sha256"]:
        raise ValueError(
            "Plotted columns differ from the verified figure inputs; review their provenance first."
        )
    if records.record_id.duplicated().any():
        raise ValueError("Expected unique source-row locators")
    if set(records.tumor_size_group) != set(SIZE_GROUPS):
        raise ValueError("Unexpected diameter groups")
    if (
        not records.cr_status.isin([0, 1, 2]).all()
        or not records[["os_event", "lpfs_event"]].isin([0, 1]).all().all()
    ):
        raise ValueError("Unexpected event codes")


def apply_followup_overrides(records, overrides):
    audit = []
    for row in overrides.itertuples(index=False):
        selected = records.record_id == row.record_id
        if (
            selected.sum() != 1
            or records.loc[selected, ["os_event", "lpfs_event", "cr_status"]]
            .to_numpy()
            .any()
        ):
            raise ValueError(
                "Follow-up override must identify exactly one censored record"
            )
        for endpoint in ENDPOINT_TIMES:
            field = f"{endpoint}_time_days"
            before = float(records.loc[selected, field].iloc[0])
            audit.append(
                {
                    "record_id": int(row.record_id),
                    "field": field,
                    "original_days": before,
                    "retained_figure_days": int(row.followup_days),
                    "changed": before != row.followup_days,
                }
            )
        records.loc[selected, ["os_time_days", "cr_time_days", "lpfs_time_days"]] = (
            row.followup_days
        )
    return pd.DataFrame(audit)


def add_endpoint_months(records, days_per_month):
    for endpoint in ENDPOINT_TIMES:
        months = records[f"{endpoint}_time_days"] / days_per_month
        if not np.isfinite(months).all() or (months <= 0).any():
            raise ValueError("Invalid endpoint duration")
        records[f"{endpoint}_time_months"] = months


def fit_curves(records):
    curves = {}
    subsets = [("overall", records)] + [
        (group, records[records.tumor_size_group == group]) for group in SIZE_GROUPS
    ]
    for group, subset in subsets:
        for endpoint in ("os", "lpfs"):
            curves[(group, endpoint)] = kaplan_meier_curve(
                subset[f"{endpoint}_time_months"], subset[f"{endpoint}_event"]
            )
        for cause in (1, 2):
            curves[(group, f"cif_{cause}")] = cumulative_incidence_curve(
                subset.cr_time_months, subset.cr_status, cause=cause
            )
    return curves


def load_figure_inputs(input_path, reference_dir, overrides_path):
    provenance = json.loads(
        (reference_dir / "provenance.json").read_text(encoding="utf-8")
    )
    verify_reference_files(provenance, reference_dir, overrides_path)
    records = pd.read_csv(input_path)
    validate_records(records, provenance)
    audit = apply_followup_overrides(records, pd.read_csv(overrides_path))
    add_endpoint_months(records, provenance["days_per_month"])
    expected_risks = provenance["expected_risk_counts"]
    if (
        risk_counts(records, "cr_time_months") != expected_risks["cr"]
        or risk_counts(records, "os_time_months") != expected_risks["os"]
    ):
        raise ValueError("Risk sets no longer reproduce the supplied figures")

    stage_table = pd.read_csv(
        reference_dir / "retained_stage_curves.csv", float_precision="round_trip"
    )
    stage_curves = {
        stage: stage_table.loc[
            stage_table.stage == stage, ["time", "estimate_pct"]
        ].reset_index(drop=True)
        for stage in STAGE_GROUPS
    }
    stage_counts = pd.read_csv(reference_dir / "retained_stage_risk_counts.csv")
    stage_risks = {
        stage: stage_counts.loc[stage_counts.stage == stage, "number_at_risk"].tolist()
        for stage in STAGE_GROUPS
    }
    curves = fit_curves(records)
    overall_os = curves[("overall", "os")]
    return FigureInputs(
        records,
        curves,
        stage_curves,
        stage_risks,
        {
            endpoint: logrank_p_value(records, endpoint, SIZE_GROUPS)
            for endpoint in ("os", "lpfs")
        },
        tuple(
            survival_median(overall_os, field)
            for field in ("estimate", "lower", "upper")
        ),
        provenance["reported_gray_p"],
        audit,
    )
