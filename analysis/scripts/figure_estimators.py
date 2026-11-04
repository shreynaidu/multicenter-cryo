"""Estimators for the supplied Python figures.

Kept as supplied so the figures reproduce. The CIF standard errors differ from cmprsk;
audit_python_figure_statistics.R reports by how much.
"""

import math

import numpy as np
import pandas as pd
from scipy import stats

RISK_TABLE_MONTHS = (0, 12, 24, 36)


def kaplan_meier_curve(times, events):
    """Kaplan-Meier with Greenwood SEs and log-log limits."""
    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=int)
    normal_quantile = stats.norm.ppf(0.975)
    rows = [
        {
            "time": 0.0,
            "estimate": 1.0,
            "se": 0.0,
            "lower": 1.0,
            "upper": 1.0,
            "n_risk": len(times),
            "n_event": 0,
            "n_censor": 0,
        }
    ]
    survival = 1.0
    greenwood_sum = 0.0
    for time in np.sort(np.unique(times)):
        n_risk = int(np.sum(times >= time))
        n_events = int(np.sum((times == time) & (events == 1)))
        n_censored = int(np.sum((times == time) & (events == 0)))
        if n_events > 0:
            survival = survival * (1.0 - n_events / n_risk)
            greenwood_sum += (
                n_events / (n_risk * (n_risk - n_events))
                if n_risk > n_events
                else np.inf
            )
        standard_error = (
            survival * math.sqrt(greenwood_sum)
            if survival > 0 and np.isfinite(greenwood_sum)
            else 0.0
        )
        if survival <= 0 or survival >= 1 or greenwood_sum == 0:
            lower, upper = (np.nan, np.nan) if survival <= 0 else (1.0, 1.0)
        else:
            log_log_error = math.sqrt(greenwood_sum) / abs(math.log(survival))
            lower = survival ** math.exp(normal_quantile * log_log_error)
            upper = survival ** math.exp(-normal_quantile * log_log_error)
        rows.append(
            {
                "time": float(time),
                "estimate": float(survival),
                "se": float(standard_error),
                "lower": float(lower),
                "upper": float(upper),
                "n_risk": n_risk,
                "n_event": n_events,
                "n_censor": n_censored,
            }
        )
    return pd.DataFrame(rows)


def cumulative_incidence_curve(times, status, cause=1):
    """Aalen-Johansen CIF for one cause, with log-log limits."""
    times = np.asarray(times, dtype=float)
    status = np.asarray(status, dtype=int)
    normal_quantile = stats.norm.ppf(0.975)
    rows = [
        {
            "time": 0.0,
            "estimate": 0.0,
            "se": 0.0,
            "lower": 0.0,
            "upper": 0.0,
            "n_risk": len(times),
            "n_event": 0,
            "n_competing": 0,
            "n_censor": 0,
        }
    ]
    event_free = 1.0
    incidence = 0.0
    risk_history, cause_events_history, all_events_history = [], [], []
    survival_history, incidence_history = [], []
    for time in np.sort(np.unique(times)):
        n_risk = int(np.sum(times >= time))
        n_cause_events = int(np.sum((times == time) & (status == cause)))
        n_competing_events = int(
            np.sum((times == time) & (status != 0) & (status != cause))
        )
        n_all_events = n_cause_events + n_competing_events
        n_censored = int(np.sum((times == time) & (status == 0)))
        if n_all_events > 0:
            incidence += event_free * (n_cause_events / n_risk)
            risk_history.append(n_risk)
            cause_events_history.append(n_cause_events)
            all_events_history.append(n_all_events)
            survival_history.append(event_free)
            incidence_history.append(incidence)
            event_free *= 1.0 - n_all_events / n_risk
            risk_values = np.array(risk_history, dtype=float)
            cause_event_values = np.array(cause_events_history, dtype=float)
            all_event_values = np.array(all_events_history, dtype=float)
            survival_values = np.array(survival_history, dtype=float)
            incidence_values = np.array(incidence_history, dtype=float)
            denominator = np.where(
                risk_values > all_event_values,
                risk_values * (risk_values - all_event_values),
                np.inf,
            )
            variance_term_1 = (
                (incidence - incidence_values) ** 2 * all_event_values / denominator
            )
            variance_term_2 = (
                survival_values**2
                * cause_event_values
                * (risk_values - cause_event_values)
                / risk_values**3
            )
            variance_term_3 = (
                -2.0
                * (incidence - incidence_values)
                * survival_values
                * cause_event_values
                / risk_values**2
            )
            variance = float(
                np.sum(variance_term_1 + variance_term_2 + variance_term_3)
            )
            standard_error = math.sqrt(max(0.0, variance))
        else:
            standard_error = rows[-1]["se"]
        if incidence <= 0.0 or incidence >= 1.0 or standard_error <= 0.0:
            lower, upper = (0.0, 0.0) if incidence <= 0.0 else (np.nan, np.nan)
        else:
            log_log_error = standard_error / (incidence * abs(math.log(incidence)))
            lower = incidence ** math.exp(normal_quantile * log_log_error)
            upper = incidence ** math.exp(-normal_quantile * log_log_error)
        rows.append(
            {
                "time": float(time),
                "estimate": float(incidence),
                "se": float(standard_error),
                "lower": float(lower),
                "upper": float(upper),
                "n_risk": n_risk,
                "n_event": n_cause_events,
                "n_competing": n_competing_events,
                "n_censor": n_censored,
            }
        )
    return pd.DataFrame(rows)


def curve_at(curve, time):
    return curve[curve["time"] <= time].iloc[-1]


def extend_curve(curve, horizon=36.0):
    shown = curve[curve["time"] <= horizon].copy()
    last = shown.iloc[-1].to_dict()
    if last["time"] < horizon:
        last["time"] = horizon
        last["n_event"] = 0
        last["n_censor"] = 0
        shown = pd.concat([shown, pd.DataFrame([last])], ignore_index=True)
    return shown


def risk_counts(records, time_column):
    times = records[time_column].values
    return [int(np.sum(times >= month)) for month in RISK_TABLE_MONTHS]


def logrank_p_value(records, endpoint, group_labels):
    """Same statistic as survival::survdiff."""
    times = records[f"{endpoint}_time_months"].to_numpy()
    events = records[f"{endpoint}_event"].to_numpy()
    groups = records["tumor_size_group"].to_numpy()
    n_groups = len(group_labels)
    observed = np.zeros(n_groups)
    expected = np.zeros(n_groups)
    variance = np.zeros((n_groups, n_groups))
    for time in np.unique(times[events == 1]):
        at_risk = np.array(
            [np.sum((times >= time) & (groups == group)) for group in group_labels]
        )
        deaths = np.array(
            [
                np.sum((times == time) & (events == 1) & (groups == group))
                for group in group_labels
            ]
        )
        total_risk, total_deaths = at_risk.sum(), deaths.sum()
        probabilities = at_risk / total_risk
        observed += deaths
        expected += probabilities * total_deaths
        if total_risk > 1:
            variance += (
                total_deaths
                * (total_risk - total_deaths)
                / (total_risk - 1)
                * (np.diag(probabilities) - np.outer(probabilities, probabilities))
            )
    difference = (observed - expected)[:-1]
    statistic = difference @ np.linalg.solve(variance[:-1, :-1], difference)
    return float(stats.chi2.sf(statistic, n_groups - 1))


def survival_median(curve, column="estimate"):
    candidates = curve.loc[curve[column] <= 0.5, "time"]
    return float(candidates.iloc[0]) if len(candidates) else math.nan
