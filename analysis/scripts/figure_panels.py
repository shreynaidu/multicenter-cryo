"""Panels, risk tables and layouts for the supplied Python figures."""

import matplotlib.pyplot as plt
from figure_estimators import RISK_TABLE_MONTHS, curve_at, extend_curve, risk_counts
from figure_inputs import SIZE_GROUPS, STAGE_GROUPS
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

PRIMARY_COLOR = "#1D4E56"
SECONDARY_COLOR = "#E68A00"
STAGE_COLORS = [PRIMARY_COLOR, SECONDARY_COLOR, "#666666"]
SIZE_COLORS = [PRIMARY_COLOR, "#3B7EA1", SECONDARY_COLOR, "#A63A2B"]
CONFIDENCE_FILL = "#BDC9CE"
BANNER_COLOR = "#2E7D32"
INK = "#111111"
SPINE_COLOR = "#555555"
GRID_COLOR = "#DCE2E5"

SIZE_LABEL_POSITIONS = {
    "os": [
        (35.8, 49.0, "top"),
        (35.8, 62.0, "top"),
        (35.8, 37.2, "top"),
        (35.8, 20.4, "top"),
    ],
    "lpfs": [
        (33.8, 36.6, "top"),
        (35.8, 45.6, "bottom"),
        (35.8, 11.2, "top"),
        (35.8, 22.8, "top"),
    ],
    "cif_1": [
        (33.8, 24.0, "top"),
        (35.8, 38.2, "top"),
        (35.8, 54.0, "bottom"),
        (35.8, 44.8, "top"),
    ],
}
STAGE_LABEL_POSITIONS = [
    (35.8, 31.2, "top"),
    (35.8, 55.4, "bottom"),
    (35.8, 39.6, "bottom"),
]
PANEL_TITLES = {
    "overall_cif": "Local progression and death before progression",
    "stage": "Local progression stratified by stage",
    "overall_os": "Overall survival",
    "os": "Overall survival stratified by tumor size",
    "lpfs": "Progression-or-death-free survival stratified by tumor size",
    "cif_1": "Local progression stratified by tumor size",
}
Y_LABELS = {
    "overall_cif": "Cumulative incidence, %",
    "stage": "Cumulative incidence, %",
    "overall_os": "Overall survival, %",
    "os": "Overall survival, %",
    "lpfs": "Progression-or-death-free probability, %",
    "cif_1": "Cumulative incidence, %",
}
# Name, alias, panel, height (in), axes bounds (in), letter.
STANDALONE_FIGURES = [
    (
        "MainFigure1",
        "Figure1_Overall_Competing_Risks_CIF",
        "overall_cif",
        4.15,
        (0.88, 4.98, 0.86, 3.62),
        "A",
    ),
    (
        "MainFigure2",
        "Figure2_Local_Progression_CIF_by_Stage",
        "stage",
        4.55,
        (1.06, 4.98, 1.25, 4.02),
        "B",
    ),
    (
        "MainFigure3",
        "Figure3_Overall_Survival",
        "overall_os",
        4.15,
        (0.88, 4.98, 0.86, 3.62),
        "C",
    ),
    (
        "SupplementaryFigure1",
        "Supplement_Figure1_OS_by_Tumor_Size",
        "os",
        4.75,
        (0.96, 4.98, 1.44, 4.22),
        "A",
    ),
    (
        "SupplementaryFigure2",
        "Supplement_Figure2_LPFS_by_Tumor_Size",
        "lpfs",
        4.75,
        (0.96, 4.98, 1.44, 4.22),
        "B",
    ),
    (
        "SupplementaryFigure3",
        "Figure2_Local_Progression_CIF_by_Tumor_Size",
        "cif_1",
        4.75,
        (0.96, 4.98, 1.44, 4.22),
        "C",
    ),
]
COMPOSITE_FIGURES = [
    (False, "Composite_JAMA_Figure1_Survival_and_Competing_Risks"),
    (True, "Composite_JAMA_Supplementary_Figures_1_to_3"),
]


def configure_style():
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "multicenter-cryo-publication-figures-v1",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "font.family": "sans-serif",
            "axes.unicode_minus": False,
        }
    )


def add_text(target, x, y, text, **style):
    style.setdefault("fontweight", "normal")
    target.text(x, y, text, color=INK, **style)


def draw_panel_header(figure, width, height, left, top, letter, title):
    box = 0.175
    figure.add_artist(
        Rectangle(
            (left / width, (top - box / 2.0) / height),
            box / width,
            box / height,
            transform=figure.transFigure,
            facecolor="none",
            edgecolor="#222222",
            linewidth=0.75,
            zorder=2,
        )
    )
    add_text(
        figure,
        (left + box / 2.0) / width,
        top / height,
        letter,
        fontsize=8.2,
        fontweight="bold",
        ha="center",
        va="center",
        zorder=10,
    )
    add_text(
        figure,
        (left + box + 0.08) / width,
        top / height,
        title,
        fontsize=9.2,
        ha="left",
        va="center",
        zorder=10,
    )


def style_axes(ax, ymax, ystep, ylabel):
    ax.set_xlim(-1.1, 36.8)
    ax.set_ylim(-0.038 * ymax, ymax * 1.025)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(SPINE_COLOR)
        ax.spines[side].set_linewidth(0.65)
    yticks = list(range(0, int(ymax) + 1, ystep))
    ax.set_yticks(yticks)
    ax.set_yticklabels(
        [str(y) for y in yticks], fontsize=8.2, fontweight="normal", color=INK
    )
    ax.set_xticks(RISK_TABLE_MONTHS)
    ax.set_xticklabels(
        [str(month) for month in RISK_TABLE_MONTHS],
        fontsize=8.2,
        fontweight="normal",
        color=INK,
    )
    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=3.5,
        width=0.65,
        color=SPINE_COLOR,
        pad=3.0,
    )
    for y in yticks:
        ax.axhline(y, color=GRID_COLOR, lw=0.6, ls="-", zorder=0)
    ax.set_ylabel(ylabel, fontsize=8.8, fontweight="normal", color=INK, labelpad=6)
    ax.set_xlabel("Time, mo", fontsize=8.8, fontweight="normal", color=INK, labelpad=5)


def line_width(composite):
    return 1.4 if composite else 1.45


def draw_step(ax, times, estimates_pct, color, composite):
    ax.step(
        times,
        estimates_pct,
        where="post",
        color=color,
        lw=line_width(composite),
        ls="-",
        zorder=4,
    )


def draw_censor_ticks(ax, curve, censored_times, color, horizon=36.0):
    for time in censored_times:
        if time <= horizon:
            estimate = curve_at(curve, time)["estimate"] * 100.0
            ax.plot(
                [time, time],
                [estimate, estimate + 2.2],
                color=color,
                lw=0.75,
                solid_capstyle="butt",
                zorder=5,
            )


def plot_overall_incidence(ax, inputs, composite):
    for cause, color in [(1, PRIMARY_COLOR), (2, SECONDARY_COLOR)]:
        curve = extend_curve(inputs.curves[("overall", f"cif_{cause}")])
        draw_step(ax, curve.time, curve.estimate * 100, color, composite)
    for y, label, vertical in [
        (43.2, "Local progression", "bottom"),
        (19.8, "Death before progression", "top"),
    ]:
        add_text(
            ax,
            35.8,
            y,
            label,
            fontsize=7.8 if composite else 8.2,
            ha="right",
            va=vertical,
        )
    progression = curve_at(inputs.curves[("overall", "cif_1")], 24)
    death = curve_at(inputs.curves[("overall", "cif_2")], 24)
    add_text(
        ax,
        1.2,
        57.5,
        f"24-mo local progression, {progression['estimate']*100:.1f}%; 95% CI, {progression['lower']*100:.1f}%-{progression['upper']*100:.1f}%\n"
        f"24-mo death before progression, {death['estimate']*100:.1f}%; 95% CI, {death['lower']*100:.1f}%-{death['upper']*100:.1f}%",
        fontsize=7.3 if composite else 7.8,
        ha="left",
        va="top",
        linespacing=1.35 if composite else 1.38,
    )


def plot_overall_survival(ax, inputs, composite):
    curve = inputs.curves[("overall", "os")]
    extended = extend_curve(curve)
    ax.fill_between(
        extended.time,
        extended.lower * 100,
        extended.upper * 100,
        step="post",
        facecolor=CONFIDENCE_FILL,
        edgecolor="none",
        zorder=2,
    )
    draw_step(ax, extended.time, extended.estimate * 100, PRIMARY_COLOR, composite)
    censored = inputs.records.loc[inputs.records.os_event == 0, "os_time_months"].values
    draw_censor_ticks(ax, curve, censored, PRIMARY_COLOR)
    milestone = curve_at(curve, 24)
    median, lower, upper = inputs.median_os
    add_text(
        ax,
        16.5 if composite else 17.2,
        96.5,
        f"24-mo OS, {milestone['estimate']*100:.1f}%; 95% CI, {milestone['lower']*100:.1f}%-{milestone['upper']*100:.1f}%\n"
        f"Median OS, {median:.1f} mo; 95% CI, {lower:.1f}-{upper:.1f}",
        fontsize=7.5 if composite else 7.9,
        ha="left",
        va="top",
        linespacing=1.35 if composite else 1.38,
    )


def plot_stage_incidence(ax, inputs, composite):
    for stage, color, (x, y, vertical) in zip(
        STAGE_GROUPS, STAGE_COLORS, STAGE_LABEL_POSITIONS
    ):
        curve = inputs.stage_curves[stage]
        draw_step(ax, curve.time, curve.estimate_pct, color, composite)
        add_text(
            ax,
            x,
            y,
            stage,
            fontsize=7.6 if composite else 7.9,
            ha="right",
            va=vertical,
        )


def plot_size_outcome(ax, inputs, endpoint, composite):
    for group, color, (x, y, vertical) in zip(
        SIZE_GROUPS, SIZE_COLORS, SIZE_LABEL_POSITIONS[endpoint]
    ):
        curve = inputs.curves[(group, endpoint)]
        extended = extend_curve(curve)
        draw_step(ax, extended.time, extended.estimate * 100, color, composite)
        if endpoint != "cif_1":
            subset = inputs.records[inputs.records.tumor_size_group == group]
            censored = subset.loc[
                subset[f"{endpoint}_event"] == 0, f"{endpoint}_time_months"
            ].values
            draw_censor_ticks(ax, curve, censored, color)
        add_text(
            ax,
            x,
            y,
            group,
            fontsize=7.4 if composite else 7.8,
            ha="right",
            va=vertical,
        )
    # The Gray P is a fixed value from the supplied figure, not calculated.
    label = (
        f"Gray test P = {inputs.reported_gray_p:.3f}"
        if endpoint == "cif_1"
        else f"Log-rank P = {inputs.logrank_p[endpoint]:.3f}"
    ).replace("= 0.", "= .")
    add_text(
        ax,
        21.5 if composite else 22.0,
        76.0 if endpoint == "cif_1" else 95.0,
        label,
        fontsize=7.6 if composite else 8.0,
        ha="left",
        va="top",
    )


def figure_x(ax, figure, x):
    display = ax.transData.transform((x, 0.0))
    return figure.transFigure.inverted().transform(display)[0]


def draw_risk_row(figure, ax, left_fraction, y, label, counts, label_size, count_size):
    add_text(
        figure, left_fraction, y, label, fontsize=label_size, ha="left", va="center"
    )
    for month, count in zip(RISK_TABLE_MONTHS, counts):
        add_text(
            figure,
            figure_x(ax, figure, month),
            y,
            str(count),
            fontsize=count_size,
            ha="center",
            va="center",
        )


def draw_risk_table(figure, ax, inputs, kind, width, height, left, bottom, composite):
    if kind in ("overall_cif", "overall_os"):
        times = "os_time_months" if kind == "overall_os" else "cr_time_months"
        size = 7.8 if composite else 8.0
        draw_risk_row(
            figure,
            ax,
            left / width,
            (bottom - (0.52 if composite else 0.56)) / height,
            "No. at risk",
            risk_counts(inputs.records, times),
            size,
            size,
        )
        return
    stage = kind == "stage"
    if composite:
        heading_offset = 0.45 if stage else 0.44
        offsets = [0.64, 0.82, 1.00] if stage else [0.62, 0.79, 0.96, 1.13]
        heading_size, row_size = (7.8, 7.6) if stage else (7.7, 7.5)
    else:
        heading_offset = 0.48
        offsets = [0.68, 0.87, 1.06] if stage else [0.68, 0.87, 1.06, 1.25]
        heading_size, row_size = 8.0, 8.0 if stage else 7.9
    add_text(
        figure,
        left / width,
        (bottom - heading_offset) / height,
        "No. at risk",
        fontsize=heading_size,
        ha="left",
        va="center",
    )
    for group, offset in zip(STAGE_GROUPS if stage else SIZE_GROUPS, offsets):
        if stage:
            counts = inputs.stage_risks[group]
        else:
            subset = inputs.records[inputs.records.tumor_size_group == group]
            counts = risk_counts(
                subset, "os_time_months" if kind == "os" else "cr_time_months"
            )
        draw_risk_row(
            figure,
            ax,
            left / width,
            (bottom - offset) / height,
            group,
            counts,
            row_size,
            row_size if composite or stage else 8.0,
        )


def add_panel(figure, inputs, kind, bounds, header, composite=False):
    width, height = figure.get_size_inches()
    left, right, bottom, top = bounds
    ax = figure.add_axes(
        [left / width, bottom / height, (right - left) / width, (top - bottom) / height]
    )
    header_left, header_y, letter = header
    title = PANEL_TITLES[kind]
    if composite and kind == "lpfs":
        title = "Progression-or-death-free survival by tumor size"
    draw_panel_header(figure, width, height, header_left, header_y, letter, title)
    maximum = 60 if kind == "overall_cif" else 80 if kind in ("stage", "cif_1") else 100
    style_axes(ax, maximum, 10 if maximum == 60 else 20, Y_LABELS[kind])
    if kind == "overall_cif":
        plot_overall_incidence(ax, inputs, composite)
    elif kind == "overall_os":
        plot_overall_survival(ax, inputs, composite)
    elif kind == "stage":
        plot_stage_incidence(ax, inputs, composite)
    else:
        plot_size_outcome(ax, inputs, kind, composite)
    draw_risk_table(
        figure, ax, inputs, kind, width, height, header_left, bottom, composite
    )
    return figure


def build_standalone(inputs, kind, height, bounds, letter):
    figure = plt.figure(figsize=(5.2, height), dpi=300, facecolor="white")
    return add_panel(figure, inputs, kind, bounds, (0.14, height - 0.24, letter))


def build_composite(inputs, supplementary=False):
    width, height = 9.6, 8.6 if supplementary else 8.1
    figure = plt.figure(figsize=(width, height), dpi=300, facecolor="white")
    figure.add_artist(
        Line2D(
            [0.15 / width, (width - 0.15) / width],
            [(height - 0.12) / height] * 2,
            color=BANNER_COLOR,
            lw=2.2,
        )
    )
    title = (
        "eFigure. Survival and Competing-Risks Outcomes Stratified by Baseline Tumor Size"
        if supplementary
        else "Figure. Competing-Risks Local Tumor Progression and Overall Survival After Lung Tumor Cryoablation"
    )
    add_text(
        figure,
        0.15 / width,
        (height - 0.30) / height,
        title,
        fontsize=9.8,
        fontweight="bold",
        ha="left",
        va="center",
    )
    figure.add_artist(
        Line2D(
            [0.15 / width, (width - 0.15) / width],
            [(height - 0.44) / height] * 2,
            color="#444444",
            lw=0.65,
        )
    )
    panels = (
        [
            ("os", (0.90, 4.58, 5.38, 7.66), (0.15, 7.86, "A")),
            ("lpfs", (5.70, 9.38, 5.38, 7.66), (4.95, 7.86, "B")),
            ("cif_1", (0.90, 4.58, 1.38, 3.66), (0.15, 3.86, "C")),
        ]
        if supplementary
        else [
            ("overall_cif", (0.84, 4.58, 4.92, 7.26), (0.15, 7.46, "A")),
            ("stage", (5.90, 9.38, 4.92, 7.26), (4.95, 7.46, "B")),
            ("overall_os", (0.84, 4.58, 1.02, 3.36), (0.15, 3.56, "C")),
        ]
    )
    for kind, bounds, header in panels:
        add_panel(figure, inputs, kind, bounds, header, composite=True)
    return figure


def figure_artist_snapshot(figure):
    return {
        "texts": [
            {"text": text.get_text(), "position": list(map(float, text.get_position()))}
            for text in list(figure.texts)
            + [text for ax in figure.axes for text in ax.texts]
        ],
        "lines": [
            {
                "x": list(map(float, line.get_xdata())),
                "y": list(map(float, line.get_ydata())),
            }
            for ax in figure.axes
            for line in ax.lines
        ],
    }
