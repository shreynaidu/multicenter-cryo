"""Recalculate endpoints and estimates in Python and compare them with the R outputs.

Shares no code with the R analysis. Doesn't recompute CIF variances or model uncertainty.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
PRIOR_DATASET = PROJECT / "tests/baseline/source/analysis/reference/prior_dataset.csv"
DAYS_PER_MONTH = 365.25 / 12
MILESTONE_MONTHS = (12, 24, 36)
TOLERANCE = 1e-12
NORMAL_97_5 = 1.959963984540054

# Source field: prior dataset column.
SOURCE_TO_PRIOR = {
    "age_at_procedure": "age_years",
    "pack_years": "pack_years",
    "disease_stage": "stage_raw",
    "tumor_size_cm": "tumor_size_cm",
    "probe_count": "probe_count",
    "hospital_stay_days": "hospital_stay_days",
    "comp_pneumothorax": "pneumothorax",
    "comp_hemothorax": "hemothorax_hemorrhage",
    "comp_effusion": "pleural_effusion",
    "intervention_chest_tube": "chest_tube",
    "cirse_comp_grade": "sir_comp_grade",
    "local_recurrence": "local_progression_original",
    "repeat_ablation": "repeat_ablation",
    "death_status": "death_status",
    "time_to_recurrence_days": "supplied_recurrence_days",
    "time_to_death_days": "supplied_death_days",
    "time_to_lastfu_days": "supplied_followup_days",
    "pathology_result": "pathology_original",
}
# Source field: (code labels, prior dataset column).
SOURCE_CODES = {
    "sex": ({0: "Male", 1: "Female"}, "sex"),
    "smoking_history": ({0: "Never", 1: "Ever"}, "smoking_history"),
    "copd_status": ({0: "No", 1: "Yes"}, "copd"),
    "o2_requirement_baseline": (
        {0: "Room air", 1: "2 L/min", 2: "3 L/min", 3: ">3 L/min"},
        "oxygen_category",
    ),
    "pd_l1_status": ({0: "Negative", 1: "Positive", 99: "Unknown"}, "pd_l1"),
    "systemic_therapy": ({0: "No", 1: "Yes"}, "systemic_therapy"),
    "tumor_location_cp": ({0: "Peripheral", 1: "Central"}, "central_location"),
    "technical_success": ({0: "No", 1: "Yes"}, "technical_success"),
}
ENDPOINT_COLUMNS = {
    "crtime": "cr_time_days",
    "crstatus": "cr_status",
    "ostime": "os_time_days",
    "osevent": "os_event",
    "corrected_fu": "corrected_followup_days",
}


def number(value):
    return None if value is None or value == "" else float(value)


def same_value(left, right):
    if left in (None, "") and right in (None, ""):
        return True
    try:
        return float(left) == float(right)
    except (ValueError, TypeError):
        return str(left) == str(right)


def read_worksheet(worksheet):
    rows = list(worksheet.values)
    headers = rows[0]
    return [
        {key: value for key, value in zip(headers, row) if key is not None}
        for row in rows[1:]
        if any(value is not None for value in row)
    ]


def read_csv_rows(path):
    with path.open(encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file))


def read_source(path):
    if path.suffix.lower() != ".xlsx":
        return [
            {key: number(value) for key, value in row.items()}
            for row in read_csv_rows(path)
        ], None
    import openpyxl

    workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
    return read_worksheet(workbook["Sheet2"]), read_worksheet(workbook["Raw"])


def derive_endpoints(records):
    rows, extensions = [], []
    for record_id, record in enumerate(records, 1):
        progression = int(
            record["local_recurrence"] == 1
            or record["time_to_recurrence_days"] is not None
        )
        death = int(record["death_status"])
        recurrence_days = number(record["time_to_recurrence_days"])
        death_days = number(record["time_to_death_days"])
        followup_days = number(record["time_to_lastfu_days"])
        if (
            followup_days is None
            or followup_days <= 0
            or (progression and recurrence_days is None)
            or (death and death_days is None)
        ):
            raise ValueError("Invalid source endpoint coding or chronology.")
        corrected = max(
            [followup_days]
            + ([recurrence_days] if progression else [])
            + ([death_days] if death else [])
        )
        if corrected > followup_days:
            extensions.append(
                {
                    "record_id": record_id,
                    "supplied_followup_days": followup_days,
                    "extended_followup_days": corrected,
                    "recurrence_days": recurrence_days,
                    "death_days": death_days,
                }
            )
        status = (
            1
            if progression and (not death or recurrence_days <= death_days)
            else 2 if death else 0
        )
        rows.append(
            {
                "record_id": record_id,
                "lp": progression,
                "recode": int(progression and record["local_recurrence"] == 0),
                "death": death,
                "rec": recurrence_days,
                "dt": death_days,
                "fu": followup_days,
                "corrected_fu": corrected,
                "crtime": (
                    recurrence_days
                    if status == 1
                    else death_days if status == 2 else corrected
                ),
                "crstatus": status,
                "ostime": death_days if death else corrected,
                "osevent": death,
            }
        )
    return rows, extensions


def kaplan_meier(times, events):
    times = np.asarray(times, float)
    events = np.asarray(events, int)
    survival = 1.0
    greenwood = 0.0
    curve = []
    for time in sorted(set(times)):
        at_risk = int(np.sum(times >= time))
        deaths = int(np.sum((times == time) & (events == 1)))
        if deaths:
            survival *= 1 - deaths / at_risk
            greenwood = (
                greenwood + deaths / (at_risk * (at_risk - deaths))
                if at_risk > deaths
                else float("inf")
            )
        if 0 < survival < 1:
            error = math.sqrt(greenwood) / abs(math.log(survival))
            log_log = math.log(-math.log(survival))
            lower = math.exp(-math.exp(log_log + NORMAL_97_5 * error))
            upper = math.exp(-math.exp(log_log - NORMAL_97_5 * error))
        else:
            lower = upper = survival
        curve.append(
            {
                "days": float(time),
                "survival": survival,
                "lower": lower,
                "upper": upper,
                "at_risk": at_risk,
                "events": deaths,
            }
        )
    return curve


def aalen_johansen(times, status):
    times = np.asarray(times, float)
    status = np.asarray(status, int)
    event_free = 1.0
    ltp_incidence = death_incidence = 0.0
    curve = []
    for time in sorted(set(times)):
        at_risk = int(np.sum(times >= time))
        ltp_events = int(np.sum((times == time) & (status == 1)))
        death_events = int(np.sum((times == time) & (status == 2)))
        ltp_incidence += event_free * ltp_events / at_risk
        death_incidence += event_free * death_events / at_risk
        event_free *= 1 - (ltp_events + death_events) / at_risk
        curve.append(
            {
                "days": float(time),
                "ltp_cif": ltp_incidence,
                "death_before_ltp_cif": death_incidence,
                "event_free": event_free,
                "at_risk": at_risk,
                "ltp_events": ltp_events,
                "death_events": death_events,
            }
        )
    return curve


def step_value(curve, time, key, default):
    reached = [row for row in curve if row["days"] <= time]
    return reached[-1][key] if reached else default


def median_months(curve, key="survival"):
    return next(
        (row["days"] / DAYS_PER_MONTH for row in curve if row[key] <= 0.5), None
    )


def milestone_estimates(curves, cr_times, os_times):
    estimates = []
    for month in MILESTONE_MONTHS:
        time = month * DAYS_PER_MONTH
        for key in ("ltp_cif", "death_before_ltp_cif"):
            estimates.append(
                {
                    "endpoint": key,
                    "month": month,
                    "estimate": step_value(curves["aj"], time, key, 0),
                    "at_risk": sum(value >= time for value in cr_times),
                }
            )
        for key in ("os", "composite"):
            estimates.append(
                {
                    "endpoint": key,
                    "month": month,
                    "estimate": step_value(curves[key], time, "survival", 1),
                    "lower": step_value(curves[key], time, "lower", 1),
                    "upper": step_value(curves[key], time, "upper", 1),
                    "at_risk": sum(
                        value >= time
                        for value in (os_times if key == "os" else cr_times)
                    ),
                }
            )
    return estimates


def reconcile_with_r(estimates, analysis_dir):
    comparisons = []
    for file_name, endpoint_column, endpoints in [
        (
            "cif_all_cohorts.csv",
            "event",
            {
                "Local tumor progression": "ltp_cif",
                "Death before local tumor progression": "death_before_ltp_cif",
            },
        ),
        ("km_all_cohorts.csv", "endpoint", {"OS": "os", "Composite": "composite"}),
    ]:
        for row in read_csv_rows(analysis_dir / "tables" / file_name):
            if row["cohort"] != "Full107" or row["stratum"] != "Overall":
                continue
            endpoint = endpoints[row[endpoint_column]]
            own = next(
                estimate
                for estimate in estimates
                if estimate["endpoint"] == endpoint
                and estimate["month"] == int(float(row["time_months"]))
            )
            pairs = [
                (field, own[field], float(row[field]))
                for field in ("estimate", "lower", "upper")
                if field in own
            ] + [("at_risk", own["at_risk"], int(row["number_at_risk"]))]
            comparisons += [
                {
                    "endpoint": endpoint,
                    "month": own["month"],
                    "quantity": quantity,
                    "independent": independent,
                    "R": from_r,
                    "absolute_difference": abs(independent - from_r),
                }
                for quantity, independent, from_r in pairs
            ]
    return comparisons


def differences_from_prior(records, prior):
    changes = []
    for record_id, (record, old) in enumerate(zip(records, prior), 1):
        if int(old["record_id"]) != record_id:
            raise ValueError("Prior dataset is not in source order.")
        for field, column in SOURCE_TO_PRIOR.items():
            if not same_value(record[field], old[column]):
                changes.append(
                    {
                        "record_id": record_id,
                        "field": field,
                        "old": old[column],
                        "new": record[field],
                    }
                )
        for field, (labels, column) in SOURCE_CODES.items():
            if labels.get(record[field]) != old[column]:
                changes.append(
                    {
                        "record_id": record_id,
                        "field": field,
                        "old": old[column],
                        "new": record[field],
                    }
                )
    return changes


def differences_between_sheets(records, raw):
    differences = []
    for record_id, (record, raw_record) in enumerate(zip(records, raw), 1):
        for field in (set(record) & set(raw_record)) - {"pt_id"}:
            if not same_value(record[field], raw_record[field]):
                differences.append(
                    {
                        "record_id": record_id,
                        "field": field,
                        "Raw": str(raw_record[field]),
                        "Sheet2": str(record[field]),
                    }
                )
    return differences


def count_by_field(rows):
    fields = sorted({row["field"] for row in rows})
    return {field: sum(row["field"] == field for row in rows) for field in fields}


def write_rows(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--output", required=True, type=Path, help="R results directory"
    )
    args = parser.parse_args()
    analysis_dir = args.output.resolve()
    report_dir = analysis_dir / "validation"
    report_dir.mkdir(parents=True, exist_ok=True)
    source = args.input.resolve()

    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    manifest = json.loads(
        (analysis_dir / "source_manifest.json").read_text(encoding="utf-8")
    )
    if source_hash != manifest["source_sha256"]:
        raise ValueError("Validation input does not match the analysis source.")

    records, raw = read_source(source)
    r_records = read_csv_rows(analysis_dir / "analysis_data/reanalysis_dataset.csv")
    if len(records) != len(r_records):
        raise ValueError("R and independent endpoint row counts differ")

    rows, extensions = derive_endpoints(records)
    cr_times = [row["crtime"] for row in rows]
    cr_status = [row["crstatus"] for row in rows]
    os_times = [row["ostime"] for row in rows]
    os_events = [row["osevent"] for row in rows]
    curves = {
        "aj": aalen_johansen(cr_times, cr_status),
        "os": kaplan_meier(os_times, os_events),
        "composite": kaplan_meier(cr_times, [int(status != 0) for status in cr_status]),
        "reverse_km": kaplan_meier(os_times, [1 - event for event in os_events]),
    }
    estimates = milestone_estimates(curves, cr_times, os_times)
    comparisons = reconcile_with_r(estimates, analysis_dir)
    mismatches = [
        {
            "record_id": row["record_id"],
            "field": own,
            "independent": row[own],
            "R": r_row[column],
        }
        for row, r_row in zip(rows, r_records)
        for own, column in ENDPOINT_COLUMNS.items()
        if not same_value(row[own], r_row[column])
    ]
    largest_difference = max(row["absolute_difference"] for row in comparisons)
    summary = {
        "source_sha256": source_hash,
        "records": len(records),
        "ltp_events": sum(row["lp"] for row in rows),
        "status_recodes": sum(row["recode"] for row in rows),
        "deaths": sum(os_events),
        "first_ltp_events": cr_status.count(1),
        "competing_deaths": cr_status.count(2),
        "composite_events": sum(status != 0 for status in cr_status),
        "recurrence_after_death": sum(
            bool(row["lp"] and row["death"] and row["rec"] > row["dt"]) for row in rows
        ),
        "ltp_death_ties": sum(
            bool(row["lp"] and row["death"] and row["rec"] == row["dt"]) for row in rows
        ),
        "followup_extensions": len(extensions),
        "supplied_fu_before_recurrence": sum(
            bool(row["lp"] and row["fu"] < row["rec"]) for row in rows
        ),
        "supplied_fu_before_death": sum(
            bool(row["death"] and row["fu"] < row["dt"]) for row in rows
        ),
        "median_os_months": median_months(curves["os"]),
        "median_os_lower_months": median_months(curves["os"], "lower"),
        "median_os_upper_months": median_months(curves["os"], "upper"),
        "median_composite_months": median_months(curves["composite"]),
        "median_composite_lower_months": median_months(curves["composite"], "lower"),
        "median_composite_upper_months": median_months(curves["composite"], "upper"),
        "reverse_km_followup_months": median_months(curves["reverse_km"]),
        "reverse_km_lower_months": median_months(curves["reverse_km"], "lower"),
        "reverse_km_upper_months": median_months(curves["reverse_km"], "upper"),
        "estimates": estimates,
        "R_reconciliation_max_absolute_difference": largest_difference,
        "R_endpoint_record_mismatches": mismatches,
    }
    if raw is not None:
        sheet_differences = differences_between_sheets(records, raw)
        summary["unique_pt_ids"] = len({record["pt_id"] for record in records})
        summary["raw_sheet2_disagreement_fields"] = count_by_field(sheet_differences)
        write_rows(report_dir / "raw_sheet2_disagreements.csv", sheet_differences)
    if PRIOR_DATASET.is_file():
        changes = differences_from_prior(records, read_csv_rows(PRIOR_DATASET))
        summary["changed_field_counts"] = count_by_field(changes)
        write_rows(report_dir / "source_changes.csv", changes)

    write_rows(report_dir / "chronology_extensions.csv", extensions)
    write_rows(report_dir / "independent_endpoint_records.csv", rows)
    write_rows(report_dir / "independent_estimates.csv", estimates)
    write_rows(report_dir / "R_reconciliation.csv", comparisons)
    for name, curve in curves.items():
        write_rows(report_dir / f"curve_{name}.csv", curve)
    (report_dir / "independent_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    if mismatches:
        raise ValueError("Derived endpoint records disagree with R outputs.")
    if largest_difference > TOLERANCE:
        raise ValueError(
            "Independent estimates, KM limits or risk sets disagree with R outputs."
        )
    print(
        f"Independent QC passed: {len(rows)} endpoint records and "
        f"{len(comparisons)} estimate, interval and risk-set comparisons."
    )


if __name__ == "__main__":
    main()
