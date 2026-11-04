"""Format the R results into publication tables, footnotes and a figure archive."""

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path

TERM_LABELS = {
    "tumor_size_cm": "Tumor diameter, per 1 cm",
    "age_per_10y": "Age, per 10 years",
    "central_locationCentral": "Central versus peripheral location",
    "oxygen_anyAny home oxygen": "Any home oxygen versus room air",
    "stage_groupStage III-IV": "Stage III-IV versus stage I-II",
    "stage_groupUnknown stage": "Unknown stage versus stage I-II",
}
PRIMARY_TERMS = list(TERM_LABELS)[:4]
SAFETY_TERM_LABELS = {
    "tumor_size_group2.0 to <3.0 cm": "Diameter 2.0 to <3.0 cm",
    "tumor_size_group3.0 to <4.0 cm": "Diameter 3.0 to <4.0 cm",
    "tumor_size_group≥4.0 cm": "Diameter ≥4.0 cm",
    "central_locationCentral": "Central versus peripheral location",
    "probe_count": "Probe count, per probe",
}
SAFETY_OUTCOME_LABELS = {
    "pneumothorax": "Pneumothorax",
    "chest_tube": "Chest-tube intervention",
    "major_complication": "SIR grade ≥3 adverse event",
}
ENDPOINT_LABELS = {"OS": "Overall survival", "Composite": "LTP-or-death-free survival"}
COUNT_WORDS = {"4": "four"}

# From the workbook update; can't be calculated from the current data.
NEWLY_SUPPLIED_FOLLOWUP_DURATIONS = 14

# Release name, internal stem, title. Stems keep their old numbers to match the baseline.
FIGURES = [
    (
        "MainFigure1",
        "figure1_overall_cif",
        "Overall cumulative incidence of LTP and death before LTP",
    ),
    (
        "MainFigure2",
        "figure2_ltp_cif_stage",
        "LTP cumulative incidence by recorded stage",
    ),
    ("MainFigure3", "figure3_overall_os", "Overall survival"),
    (
        "SupplementaryFigure1",
        "supplement_figure2_os_size",
        "Overall survival by diameter; log-rank P={os}",
    ),
    (
        "SupplementaryFigure2",
        "supplement_figure3_composite_size",
        "LTP-or-death-free survival by diameter; log-rank P={composite}",
    ),
    (
        "SupplementaryFigure3",
        "supplement_figure1_ltp_size",
        "LTP cumulative incidence by diameter; Gray P={ltp}",
    ),
]
EXPECTED_ARCHIVE_CONTENTS = {".pdf": 6, ".tiff": 6, ".csv": 14}


def format_p(value):
    p = float(value)
    return "<.001" if p < 0.001 else f"{p:.3f}".removeprefix("0")


def format_ratio(estimate, lower, upper, digits=2):
    return f"{float(estimate):.{digits}f} ({float(lower):.{digits}f}–{float(upper):.{digits}f})"


def ratio_cell(row, digits=2):
    return format_ratio(row["estimate"], row["lower"], row["upper"], digits)


def percent_cell(row):
    estimate, lower, upper = (
        100 * float(row[key]) for key in ("estimate", "lower", "upper")
    )
    return f"{estimate:.1f} ({lower:.1f}–{upper:.1f})"


def first_where(rows, **criteria):
    return next(
        row
        for row in rows
        if all(row[column] == value for column, value in criteria.items())
    )


class Release:
    def __init__(self, results):
        self.results = results
        self.tables_dir = results / "publication_tables"
        self.root = results / "release_assets"
        self.tables = []
        self.tables_dir.mkdir(exist_ok=True)
        if self.root.exists():
            raise FileExistsError(
                "Release directory already exists; run the main entry point to replace it."
            )
        (self.root / "Tables").mkdir(parents=True)
        (self.root / "Figures").mkdir()

    def read(self, relative_path):
        with (self.results / relative_path).open(encoding="utf-8", newline="") as file:
            return list(csv.DictReader(file))

    def save_rows(self, name, rows):
        with (self.tables_dir / name).open("w", encoding="utf-8", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        shutil.copyfile(self.tables_dir / name, self.root / "Tables" / name)

    def add_table(self, table_id, name, title, rows, footnotes, extra=None):
        self.save_rows(name, rows)
        entry = {
            "id": table_id,
            "title": title,
            "path": "publication_tables/" + name,
            "archive_path": "Tables/" + name,
            "columns": list(rows[0]),
            "footnotes": footnotes,
            "rows": len(rows),
        }
        if extra:
            entry["additional_csv"] = extra
        self.tables.append(entry)

    def group_p(self, endpoint, test):
        row = first_where(
            self.read("model_results/group_tests.csv"),
            grouping="tumor_size_group",
            endpoint=endpoint,
            test=test,
        )
        return format_p(row["p_value"])


class StudyCounts:
    def __init__(self, release):
        cohorts = release.read("tables/cohort_counts.csv")
        self.primary = first_where(cohorts, cohort="Full107")
        self.success = first_where(cohorts, cohort="Success103")
        self.known_stage = first_where(
            release.read("tables/cohort_flow.csv"), analysis="Known-stage sensitivity"
        )["records_after"]
        self.stages = {
            row["stage_group"]: row["n"]
            for row in release.read("tables/stage_group_counts.csv")
        }
        sizes = {}
        for row in release.read("tables/size_cif_milestones.csv"):
            sizes[row["stratum"]] = row["n"]
        self.size_groups = list(sizes.values())
        readiness = release.read("analysis_data/data_readiness_report.csv")
        self.administrative_derivations = int(
            first_where(readiness, check="Administrative follow-up derivations")[
                "observed"
            ]
        )


class Footnotes:
    def __init__(self, counts):
        primary = counts.primary
        self.base = [
            f"Primary cohort N={primary['n']}; supplied follow-up durations from updated Sheet2; no administrative follow-up derivation. Prior investigator event and grading decisions retained.",
            "LTP, local tumor progression; composite, LTP or all-cause death; OS, overall survival; HR, hazard ratio; sHR, subdistribution hazard ratio; CI, confidence interval.",
        ]
        self.regression = [
            "Main continuous models adjust for diameter per 1 cm, age per 10 years, central versus peripheral location, and any baseline home oxygen versus room air.",
            f"Fine–Gray models LTP with death before LTP competing ({primary['local_failures']} LTP and {primary['competing_deaths']} competing deaths). Composite Cox includes {primary['composite_events']} events; OS Cox includes {primary['os_deaths']} deaths. Cox ties use Efron; intervals and individual P values are Wald-based.",
        ]
        self.absolute = [
            "Estimates and 95% CIs are percentages. CIF intervals use pointwise normal approximation with cmprsk variance, bounded 0–100%; Kaplan–Meier intervals use log-log transformation.",
            "CIF/composite risk counts represent alive, LTP-free patients under observation; OS risk counts represent living patients under observation. Late subgroup estimates have sparse support; CIs are pointwise.",
        ]


def add_main_table_3(release, notes):
    models = release.read("model_results/main_continuous_models.csv")
    rows = []
    for term in PRIMARY_TERMS:
        row = {"Characteristic": TERM_LABELS[term]}
        for model, label in [
            ("Full107 Fine-Gray", "LTP sHR"),
            ("Full107 Composite Cox", "Composite HR"),
            ("Full107 OS Cox", "OS HR"),
        ]:
            estimate = first_where(models, model=model, term=term)
            row[label + " (95% CI)"] = ratio_cell(estimate)
            row[label.split()[0] + " P"] = format_p(estimate["p_value"])
        rows.append(row)
    release.add_table(
        "Main Table3",
        "MainTable3.csv",
        "Adjusted continuous-diameter models for LTP, composite outcome, and overall survival",
        rows,
        notes.base + notes.regression,
    )


def add_size_milestones(release, notes, counts):
    rows = [
        {
            "Endpoint": row["event"] + " CIF",
            "Tumor diameter": row["stratum"],
            "N": row["n"],
            "Month": row["time_months"],
            "Estimate % (95% CI)": percent_cell(row),
            "At risk": row["number_at_risk"],
        }
        for row in release.read("tables/size_cif_milestones.csv")
    ] + [
        {
            "Endpoint": ENDPOINT_LABELS[row["endpoint"]],
            "Tumor diameter": row["stratum"],
            "N": row["n"],
            "Month": row["time_months"],
            "Estimate % (95% CI)": percent_cell(row),
            "At risk": row["number_at_risk"],
        }
        for row in release.read("tables/size_km_milestones.csv")
    ]
    a, b, c, d = counts.size_groups
    release.add_table(
        "Supplementary Table1",
        "SupplementaryTable1.csv",
        "Size-specific absolute outcome estimates",
        rows,
        notes.base
        + notes.absolute
        + [
            f"Size groups contain{a},{b},{c},and {d} patients. Gray LTP P={release.group_p('LTP', 'Gray')}; size-group OS log-rank P={release.group_p('OS', 'Log-rank')}; composite log-rank P={release.group_p('Composite', 'Log-rank')}. These are unadjusted exploratory group comparisons."
        ],
    )


def add_cause_specific_ltp(release, notes, counts):
    rows = [
        {
            "Covariate": TERM_LABELS[row["term"]],
            "HR (95% CI)": ratio_cell(row),
            "P": format_p(row["p_value"]),
        }
        for row in release.read("model_results/continuous_models_all_cohorts.csv")
        if row["model"] == "Full107 LocalFailure Cox"
    ]
    release.add_table(
        "Supplementary Table2",
        "SupplementaryTable2.csv",
        "Exploratory cause-specific Cox model for LTP",
        rows,
        notes.base
        + [
            f"N={counts.primary['n']}, {counts.primary['local_failures']} LTP events. Death before LTP is censored at death. The model includes all four continuous-model covariates. This cause-specific hazard differs from the Fine–Gray subdistribution hazard."
        ],
    )


def add_technical_success_panels(release, notes, counts):
    primary, success = counts.primary, counts.success
    population = {
        "Full107": f"Primary {primary['n']}",
        "Success103": f"Technical success {success['n']}",
    }
    panel_a = [
        {
            "Population": population[row["cohort"]],
            "Endpoint": row["event"] + " CIF",
            "N": row["n"],
            "Events": row["events"],
            "Month": row["time_months"],
            "Estimate % (95% CI)": percent_cell(row),
            "At risk": row["number_at_risk"],
        }
        for row in release.read("tables/cif_all_cohorts.csv")
    ] + [
        {
            "Population": population[row["cohort"]],
            "Endpoint": ENDPOINT_LABELS[row["endpoint"]],
            "N": row["n"],
            "Events": row["events"],
            "Month": row["time_months"],
            "Estimate % (95% CI)": percent_cell(row),
            "At risk": row["number_at_risk"],
        }
        for row in release.read("tables/km_all_cohorts.csv")
    ]
    # Three decimals because some lower limits are close to 1.
    panel_b = [
        {
            "Endpoint": row["model"].replace("Success103 ", ""),
            "Covariate": TERM_LABELS[row["term"]],
            "N": row["n"],
            "Events": row["events"],
            "HR or sHR (95% CI)": ratio_cell(row, 3),
            "P": format_p(row["p_value"]),
        }
        for row in release.read("model_results/continuous_models_all_cohorts.csv")
        if row["model"].startswith("Success103 ")
    ]
    panel_b_name = "SupplementaryTable3_PanelB.csv"
    release.save_rows(panel_b_name, panel_b)
    failures = COUNT_WORDS[primary["technical_failures"]]
    release.add_table(
        "Supplementary Table3",
        "SupplementaryTable3_PanelA.csv",
        "Panel A: Primary and technical-success sensitivity absolute outcome estimates; Panel B: Technical-success continuous regressions",
        panel_a,
        notes.base
        + notes.absolute
        + [
            f"Technical success {success['n']} excludes the {failures} recorded initial technical failures;{success['local_failures']} LTP, {success['competing_deaths']} competing deaths,{success['composite_events']} composite events, and {success['os_deaths']} OS deaths remain. This restriction is a sensitivity analysis.",
            f"The accompanying Supplementary Table 3 Panel B CSV gives four-covariate continuous regressions in {success['n']}. Hazard-ratio intervals retain three decimals to display lower limits close to 1 accurately.",
        ],
        extra="Tables/" + panel_b_name,
    )


def add_categorical_models(release, notes):
    models = release.read("model_results/categorical_models_with_oxygen.csv")
    global_tests = release.read("model_results/categorical_global_tests.csv")
    specifications = [
        ("Full107 categorical Fine-Gray", "LTP sHR", "Fine-Gray"),
        ("Full107 categorical Composite", "Composite HR", "Composite"),
        ("Full107 categorical OS", "OS HR", "OS"),
    ]
    rows_spec = [
        ("reference", "Tumor diameter <2.0 cm (reference)"),
        ("tumor_size_group2.0 to <3.0 cm", "2.0 to <3.0 cm"),
        ("tumor_size_group3.0 to <4.0 cm", "3.0 to <4.0 cm"),
        ("tumor_size_group≥4.0 cm", "≥4.0 cm"),
        ("global", "Tumor diameter, global test (3 df)"),
    ] + [(term, TERM_LABELS[term]) for term in PRIMARY_TERMS[1:]]
    rows = []
    for term, label in rows_spec:
        row = {"Characteristic": label}
        for model, column, endpoint in specifications:
            if term == "reference":
                estimate, p_value = "1.00 (reference)", "—"
            elif term == "global":
                test = first_where(global_tests, endpoint=endpoint)
                estimate, p_value = "—", format_p(test["global_size_p"])
            else:
                fitted = first_where(models, model=model, term=term)
                digits = (
                    3 if term == "central_locationCentral" and endpoint == "OS" else 2
                )
                estimate, p_value = ratio_cell(fitted, digits), format_p(
                    fitted["p_value"]
                )
            row[column + " (95% CI)"] = estimate
            row[column.split()[0] + " P"] = p_value
        rows.append(row)
    release.add_table(
        "Supplementary Table4",
        "SupplementaryTable4.csv",
        "Supplementary categorical-diameter models with consistent oxygen adjustment",
        rows,
        notes.base
        + [
            "Every model includes diameter category,age,central location,and baseline oxygen. Global size tests use 3-df Wald test for Fine–Gray and nested-model likelihood-ratio tests for Cox.",
            "The central-location OS interval retains three decimals because its lower limit is close to 1. Categorical findings are supplementary; the continuous model is the primary adjusted analysis.",
        ],
    )


def add_safety_models(release, counts):
    models = release.read("model_results/safety_models.csv")
    rows = [
        {
            "Recorded outcome": row["outcome"],
            "Covariate": SAFETY_TERM_LABELS[row["term"]],
            "OR (95% CI)": ratio_cell(row),
            "P": format_p(row["p_value"]),
        }
        for row in models
    ]
    global_p = {row["outcome"]: format_p(row["global_size_p"]) for row in models}
    release.add_table(
        "Supplementary Table5",
        "SupplementaryTable5.csv",
        "Firth penalized logistic models for recorded safety outcomes",
        rows,
        [
            f"OR, odds ratio; CI, confidence interval. Each model includes tumor-diameter category, location, and probe count; diameter <2.0 cm is the reference. N={counts.primary['n']}.",
            "Intervals and individual tests use profile penalized likelihood. Global size tests constrain its three coefficients in the full model: "
            f"pneumothorax P={global_p['Pneumothorax']}; chest tube P={global_p['Chest-tube intervention']}; recorded SIR grade ≥3 P={global_p['SIR grade ≥3 adverse event']}.",
            "The grade outcome uses the retained recorded/adjudicated grades; chest-tube intervention does not automatically change a grade.",
            "SIR, Society of Interventional Radiology.",
        ],
    )


def convergence_suffix(row):
    return (
        "; converged" if row["converged"] == "TRUE" else "; convergence not confirmed"
    )


def add_diagnostics(release, notes):
    rows = []
    for row in release.read("diagnostics/location_leave_one_out_summary.csv"):
        rows.append(
            {
                "Analysis": row["endpoint"] + " central-location leave-one-out",
                "Result": f"HR {float(row['min_hr']):.3f}–{float(row['max_hr']):.3f}; P {float(row['min_p']):.3f}–{float(row['max_p']):.3f}; {row['omissions_p_ge_05']}/{row['successful']} omissions P≥.05; {row['successful']} successful fits",
            }
        )
    primary_cox = [
        "Full107 Composite Cox",
        "Full107 OS Cox",
        "Full107 LocalFailure Cox",
    ]
    for row in release.read("diagnostics/cox_ph_tests.csv"):
        if row["term"] == "GLOBAL" and row["model"] in primary_cox:
            rows.append(
                {
                    "Analysis": row["model"] + " global Schoenfeld",
                    "Result": "P=" + format_p(row["p_value"]),
                }
            )
    for row in release.read("diagnostics/primary_fine_gray_time_interactions.csv"):
        rows.append(
            {
                "Analysis": "Fine–Gray log-time interaction: "
                + TERM_LABELS[row["term"]],
                "Result": "P=" + format_p(row["p_value"]) + convergence_suffix(row),
            }
        )
    for row in release.read("diagnostics/continuous_size_nonlinearity.csv"):
        rows.append(
            {
                "Analysis": row["endpoint"] + " diameter nonlinearity",
                "Result": "P="
                + format_p(row["nonlinearity_p"])
                + convergence_suffix(row),
            }
        )
    competing_death = first_where(
        release.read("diagnostics/competing_death_ph_tests.csv"), term="GLOBAL"
    )
    rows.append(
        {
            "Analysis": "Exploratory competing-death global Schoenfeld",
            "Result": "P=" + format_p(competing_death["p_value"]),
        }
    )
    fitted = [
        row
        for name in (
            "main_continuous_models",
            "continuous_models_all_cohorts",
            "categorical_models_with_oxygen",
            "stage_adjusted_sensitivity",
            "competing_death_exploratory",
        )
        for row in release.read(f"model_results/{name}.csv")
    ]
    all_converged = all(
        row["converged"] == "TRUE" and not row.get("warnings") for row in fitted
    )
    rows.append(
        {
            "Analysis": "Convergence of main and supporting fitted models",
            "Result": (
                "All recorded fits converged; no captured Cox fitting warnings"
                if all_converged
                else "Review model_results for convergence flags and warnings"
            ),
        }
    )
    release.add_table(
        "Supplementary Table6",
        "SupplementaryTable6.csv",
        "Updated model diagnostics and location influence checks",
        rows,
        notes.base
        + [
            "Leave-one-out checks omit each patient once and assess coefficient influence; they are not predictive cross-validation.",
            "Nonlinearity uses 3-df natural spline versus linear diameter (2-df test); Fine–Gray uses nonlinear-contrast Wald test and Cox uses likelihood-ratio test.",
            "Fine–Gray proportionality checks add one covariate-by-log-time term at a time. Cox checks use Schoenfeld residuals. Diagnostics have limited power.",
        ],
    )


def add_stage_outcomes(release, notes, counts):
    rows = [
        {
            "Recorded stage": row["stratum"],
            "Endpoint": row["event"] + " CIF",
            "N": row["n"],
            "Events": row["events"],
            "Month": row["time_months"],
            "Estimate % (95% CI)": percent_cell(row),
            "At risk": row["number_at_risk"],
        }
        for row in release.read("tables/stage_cif_24months.csv")
    ] + [
        {
            "Recorded stage": row["stratum"],
            "Endpoint": "Overall survival",
            "N": row["n"],
            "Events": row["events"],
            "Month": row["time_months"],
            "Estimate % (95% CI)": percent_cell(row),
            "At risk": row["number_at_risk"],
        }
        for row in release.read("tables/stage_os_24months.csv")
    ]
    stages = counts.stages
    release.add_table(
        "Supplementary Table7",
        "SupplementaryTable7.csv",
        "Two-year outcomes by recorded stage group",
        rows,
        notes.base
        + notes.absolute
        + [
            f"Recorded stage I–II N={stages['Stage I-II']}, stage III–IV N={stages['Stage III-IV']}, and unknown N={stages['Unknown stage']}. Stage labels do not establish isolated recurrence,oligoprogression,or treatment intent. Events are total observed events of that endpoint in each stratum,not events restricted to 24 months."
        ],
    )


def add_stage_adjusted_models(release, notes, counts):
    rows = [
        {
            "Model": row["model"],
            "Covariate": TERM_LABELS[row["term"]],
            "N": row["n"],
            "Events": row["events"],
            "HR (95% CI)": ratio_cell(row),
            "P": format_p(row["p_value"]),
        }
        for row in release.read("model_results/stage_adjusted_sensitivity.csv")
    ]
    release.add_table(
        "Supplementary Table8",
        "SupplementaryTable8.csv",
        "Stage-adjusted sensitivity analyses for composite outcome and overall survival",
        rows,
        notes.base
        + [
            f"Each model includes the four main covariates plus recorded stage group. Unknown-category models retain {counts.primary['n']} patients; known-stage complete-case models retain {counts.known_stage}. No stage was imputed.",
            f"Reference stage group is I–II. The known-stage results remain unchanged because all {NEWLY_SUPPLIED_FOLLOWUP_DURATIONS} newly supplied follow-up durations belong to the unknown-stage group.",
        ],
    )


def add_competing_death_model(release, notes, counts):
    rows = [
        {
            "Covariate": TERM_LABELS[row["term"]],
            "HR (95% CI)": ratio_cell(row),
            "P": format_p(row["p_value"]),
        }
        for row in release.read("model_results/competing_death_exploratory.csv")
    ]
    release.add_table(
        "Supplementary Table9",
        "SupplementaryTable9.csv",
        "Exploratory cause-specific model for death before LTP",
        rows,
        notes.base
        + [
            f"N={counts.primary['n']}, {counts.primary['competing_deaths']} deaths before LTP. The model includes age per decade,central location,and oxygen,with censoring at LTP. Sparse-event estimates are imprecise.",
            "These are all-cause deaths before LTP; they are not classified as noncancer or procedure-related deaths.",
        ],
    )


def add_descriptive_tables(release):
    for number, source, title, note in [
        (
            1,
            "table1_baseline_characteristics.csv",
            "Baseline characteristics",
            "Table 1 sex labels conflict with the maintained analysis coding; owner adjudication is required. See docs/FLAGGED_ISSUES.md.",
        ),
        (
            2,
            "table2_procedure_safety.csv",
            "Procedure and recorded safety outcomes",
            "Retains recorded grades; deaths are all-cause. Causal procedure-related mortality is unavailable.",
        ),
    ]:
        name = f"MainTable{number}.csv"
        shutil.copyfile(release.results / "tables" / source, release.tables_dir / name)
        shutil.copyfile(release.tables_dir / name, release.root / "Tables" / name)
        release.tables.append(
            {
                "id": f"Main Table{number}",
                "title": title,
                "path": "publication_tables/" + name,
                "archive_path": "Tables/" + name,
                "footnotes": [
                    "Retains the original descriptive table, including its type-2 quartiles.",
                    note,
                ],
            }
        )


def add_probe_omitted_panel(release):
    omitted = release.read("supplemental/refitted_probe_omitted_models.csv")
    rows = [
        {
            "Recorded outcome": SAFETY_OUTCOME_LABELS[row["outcome"]],
            "Covariate": SAFETY_TERM_LABELS[row["term"]],
            "OR (95% CI)": format_ratio(row["OR"], row["lower"], row["upper"]),
            "P": format_p(row["p"]),
        }
        for row in omitted
    ]
    name = "SupplementaryTable5_PanelB.csv"
    release.save_rows(name, rows)
    global_p = "; ".join(
        label + "=" + format_p(first_where(omitted, outcome=outcome)["global_p"])
        for outcome, label in SAFETY_OUTCOME_LABELS.items()
    )
    release.tables.append(
        {
            "id": "Supplementary Table5 Panel B",
            "title": "Probe-omitted Firth safety models",
            "path": "publication_tables/" + name,
            "archive_path": "Tables/" + name,
            "rows": len(rows),
            "footnotes": [
                "Same outcomes, diameter groups, central-location adjustment, and profile penalized-likelihood intervals; probe count omitted.",
                "Global size P values: " + global_p,
            ],
        }
    )


def copy_figures(release, counts):
    p_values = {
        "os": release.group_p("OS", "Log-rank"),
        "composite": release.group_p("Composite", "Log-rank"),
        "ltp": release.group_p("LTP", "Gray"),
    }
    figure_dir = release.results / "publication_figures"
    entries = []
    for label, stem, title in FIGURES:
        shutil.copyfile(
            figure_dir / f"{stem}.pdf", release.root / "Figures" / f"{label}.pdf"
        )
        shutil.copyfile(
            figure_dir / f"{stem}_600dpi.tiff",
            release.root / "Figures" / f"{label}_600dpi.tiff",
        )
        entries.append(
            {
                "id": label,
                "title": title.format(**p_values),
                "pdf": f"Figures/{label}.pdf",
                "tiff": f"Figures/{label}_600dpi.tiff",
                "note": f"Primary N={counts.primary['n']}; pointwise 95% CIs; display 0–36 months; numbers at risk shown at 0, 12, 24, and 36 months. Size and stage groups are exploratory/descriptive.",
            }
        )
    return entries


def write_manifests(release, counts, figures):
    source = json.loads(
        (release.results / "source_manifest.json").read_text(encoding="utf-8")
    )
    primary_n, success_n = counts.primary["n"], counts.success["n"]
    package = {
        "source_workbook": source["source_name"],
        "source_sheet": source["source_sheet"],
        "source_sha256": source["source_sha256"],
        "primary_n": int(primary_n),
        "technical_success_sensitivity_n": int(success_n),
        "administrative_followup_derivations": counts.administrative_derivations,
        "tables": release.tables,
        "figures": figures,
        "quality_control": "See run_manifest.json for provenance and analysis/validation for independent checks. No manuscript or DOCX included.",
    }
    (release.tables_dir / "tables_manifest.json").write_text(
        json.dumps({"tables": release.tables}, indent=2), encoding="utf-8"
    )
    (release.root / "manifest.json").write_text(
        json.dumps(package, indent=2), encoding="utf-8"
    )
    (release.root / "README_manifest.md").write_text(
        "# Statistical tables and figures\n\nSource: "
        + source["source_name"]
        + f", {source['source_sheet']}. SHA-256: `"
        + source["source_sha256"]
        + f"`.\n\nAll {primary_n} supplied follow-up durations are used. The primary cohort contains {primary_n} records; technical-success sensitivity analyses contain {success_n}. Prior event and grading decisions are retained.\n\nMain Tables 1–3 and Supplementary Tables 1–9 are CSV files. Supplementary Tables 3 and 5 have two panels. Table footnotes and figure descriptions appear in manifest.json. Each of the six figures has a vector PDF and a 600-dpi TIFF.\n\nCIF intervals use bounded pointwise normal approximations; Kaplan–Meier intervals use the log-log transformation. See analysis/README.md for methods and validation instructions.\n",
        encoding="utf-8",
    )


def write_archive(release):
    archive_path = release.results / "Tables_and_figures.zip"
    with zipfile.ZipFile(
        archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6
    ) as archive:
        for path in sorted(release.root.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(release.root))
    with zipfile.ZipFile(archive_path) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("Archive checksum check failed")
        for suffix, expected in EXPECTED_ARCHIVE_CONTENTS.items():
            found = sum(name.endswith(suffix) for name in archive.namelist())
            if found != expected:
                raise RuntimeError(f"Expected {expected} {suffix} files; found {found}")
    return archive_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", required=True, type=Path, help="R results directory"
    )
    args = parser.parse_args()
    release = Release(args.output.resolve())
    counts = StudyCounts(release)
    notes = Footnotes(counts)
    add_main_table_3(release, notes)
    add_size_milestones(release, notes, counts)
    add_cause_specific_ltp(release, notes, counts)
    add_technical_success_panels(release, notes, counts)
    add_categorical_models(release, notes)
    add_safety_models(release, counts)
    add_diagnostics(release, notes)
    add_stage_outcomes(release, notes, counts)
    add_stage_adjusted_models(release, notes, counts)
    add_competing_death_model(release, notes, counts)
    add_descriptive_tables(release)
    add_probe_omitted_panel(release)
    figures = copy_figures(release, counts)
    write_manifests(release, counts, figures)
    print(write_archive(release))
    print(f"Tables: {len(release.tables)}; figures: {len(figures)} PDF and TIFF sets.")


if __name__ == "__main__":
    main()
