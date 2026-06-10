#!/usr/bin/env python3
"""Compare PBMC relative frequencies for melanoma miraclib responders."""

from __future__ import annotations

import csv
import math
import sqlite3
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
STATS_PATH = ROOT / "miraclib_melanoma_pbmc_stats.csv"
PLOT_PATH = ROOT / "miraclib_melanoma_pbmc_boxplot.svg"

COHORT_QUERY = """
WITH sample_totals AS (
    SELECT
        sample_id,
        SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT
    s.sample_id AS sample,
    p.response,
    cp.name AS population,
    cc.count * 100.0 / st.total_count AS percentage
FROM patients AS p
JOIN samples AS s
    ON s.subject_id = p.subject_id
JOIN cell_counts AS cc
    ON cc.sample_id = s.sample_id
JOIN cell_populations AS cp
    ON cp.population_id = cc.population_id
JOIN sample_totals AS st
    ON st.sample_id = s.sample_id
WHERE p.indication = 'melanoma'
    AND p.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND p.response IN ('yes', 'no')
ORDER BY cp.population_id, s.sample_id;
"""


@dataclass(frozen=True)
class TestResult:
    population: str
    responder_n: int
    non_responder_n: int
    responder_mean: float
    non_responder_mean: float
    responder_median: float
    non_responder_median: float
    mean_difference: float
    median_difference: float
    mann_whitney_u: float
    p_value: float
    fdr_p_value: float
    rank_biserial: float
    significant: bool


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    middle = n // 2
    if n % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[int(position)]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def average_ranks(values: list[tuple[float, str]]) -> tuple[dict[str, float], list[int]]:
    indexed = sorted((value, label, index) for index, (value, label) in enumerate(values))
    ranks = [0.0] * len(values)
    tie_sizes: list[int] = []
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][0] == indexed[i][0]:
            j += 1
        average_rank = (i + 1 + j) / 2
        for _, _, original_index in indexed[i:j]:
            ranks[original_index] = average_rank
        if j - i > 1:
            tie_sizes.append(j - i)
        i = j

    rank_sums = {"yes": 0.0, "no": 0.0}
    for rank, (_, label) in zip(ranks, values):
        rank_sums[label] += rank
    return rank_sums, tie_sizes


def mann_whitney_u(responders: list[float], non_responders: list[float]) -> tuple[float, float, float]:
    n_yes = len(responders)
    n_no = len(non_responders)
    values = [(value, "yes") for value in responders]
    values.extend((value, "no") for value in non_responders)
    rank_sums, tie_sizes = average_ranks(values)

    u_yes = rank_sums["yes"] - n_yes * (n_yes + 1) / 2
    n_total = n_yes + n_no
    mean_u = n_yes * n_no / 2

    tie_term = sum(size**3 - size for size in tie_sizes)
    variance = n_yes * n_no / 12 * (
        (n_total + 1) - tie_term / (n_total * (n_total - 1))
    )
    if variance <= 0:
        return u_yes, 1.0, 0.0

    sd = math.sqrt(variance)
    if u_yes > mean_u:
        z = (u_yes - mean_u - 0.5) / sd
    elif u_yes < mean_u:
        z = (u_yes - mean_u + 0.5) / sd
    else:
        z = 0.0

    p_value = math.erfc(abs(z) / math.sqrt(2))
    rank_biserial = (2 * u_yes) / (n_yes * n_no) - 1
    return u_yes, p_value, rank_biserial


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    previous = 1.0
    total = len(p_values)
    for rank, (index, p_value) in reversed(list(enumerate(indexed, start=1))):
        value = min(previous, p_value * total / rank)
        adjusted[index] = min(value, 1.0)
        previous = value
    return adjusted


def load_percentages() -> dict[str, dict[str, list[float]]]:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {DB_PATH.name}. Run `python load_data.py` first."
        )

    grouped: dict[str, dict[str, list[float]]] = {}
    with sqlite3.connect(DB_PATH) as connection:
        connection.row_factory = sqlite3.Row
        for row in connection.execute(COHORT_QUERY):
            population = row["population"]
            response = row["response"]
            grouped.setdefault(population, {"yes": [], "no": []})
            grouped[population][response].append(float(row["percentage"]))
    return grouped


def analyze(grouped: dict[str, dict[str, list[float]]]) -> list[TestResult]:
    partial_results = []
    for population, by_response in grouped.items():
        responders = by_response["yes"]
        non_responders = by_response["no"]
        u_statistic, p_value, rank_biserial = mann_whitney_u(responders, non_responders)
        responder_mean = sum(responders) / len(responders)
        non_responder_mean = sum(non_responders) / len(non_responders)
        responder_median = median(responders)
        non_responder_median = median(non_responders)
        partial_results.append(
            {
                "population": population,
                "responder_n": len(responders),
                "non_responder_n": len(non_responders),
                "responder_mean": responder_mean,
                "non_responder_mean": non_responder_mean,
                "responder_median": responder_median,
                "non_responder_median": non_responder_median,
                "mean_difference": responder_mean - non_responder_mean,
                "median_difference": responder_median - non_responder_median,
                "mann_whitney_u": u_statistic,
                "p_value": p_value,
                "rank_biserial": rank_biserial,
            }
        )

    adjusted = benjamini_hochberg([result["p_value"] for result in partial_results])
    results = []
    for result, fdr_p_value in zip(partial_results, adjusted):
        results.append(
            TestResult(
                population=result["population"],
                responder_n=result["responder_n"],
                non_responder_n=result["non_responder_n"],
                responder_mean=result["responder_mean"],
                non_responder_mean=result["non_responder_mean"],
                responder_median=result["responder_median"],
                non_responder_median=result["non_responder_median"],
                mean_difference=result["mean_difference"],
                median_difference=result["median_difference"],
                mann_whitney_u=result["mann_whitney_u"],
                p_value=result["p_value"],
                fdr_p_value=fdr_p_value,
                rank_biserial=result["rank_biserial"],
                significant=fdr_p_value < 0.05,
            )
        )
    return results


def write_stats(results: list[TestResult]) -> None:
    headers = [field for field in TestResult.__dataclass_fields__]
    with STATS_PATH.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "population": result.population,
                    "responder_n": result.responder_n,
                    "non_responder_n": result.non_responder_n,
                    "responder_mean": f"{result.responder_mean:.6f}",
                    "non_responder_mean": f"{result.non_responder_mean:.6f}",
                    "responder_median": f"{result.responder_median:.6f}",
                    "non_responder_median": f"{result.non_responder_median:.6f}",
                    "mean_difference": f"{result.mean_difference:.6f}",
                    "median_difference": f"{result.median_difference:.6f}",
                    "mann_whitney_u": f"{result.mann_whitney_u:.3f}",
                    "p_value": f"{result.p_value:.6g}",
                    "fdr_p_value": f"{result.fdr_p_value:.6g}",
                    "rank_biserial": f"{result.rank_biserial:.6f}",
                    "significant": result.significant,
                }
            )


def box_stats(values: list[float]) -> dict[str, float | list[float]]:
    q1 = quantile(values, 0.25)
    q2 = quantile(values, 0.50)
    q3 = quantile(values, 0.75)
    iqr = q3 - q1
    lower_limit = q1 - 1.5 * iqr
    upper_limit = q3 + 1.5 * iqr
    inlier_values = [value for value in values if lower_limit <= value <= upper_limit]
    outliers = [value for value in values if value < lower_limit or value > upper_limit]
    return {
        "q1": q1,
        "median": q2,
        "q3": q3,
        "lower": min(inlier_values),
        "upper": max(inlier_values),
        "outliers": outliers,
    }


def svg_text(x: float, y: float, text: str, size: int = 13, anchor: str = "middle") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, sans-serif" text-anchor="{anchor}">{text}</text>'
    )


def write_boxplot(grouped: dict[str, dict[str, list[float]]]) -> None:
    populations = list(grouped)
    width = 1120
    height = 680
    left = 82
    right = 38
    top = 62
    bottom = 112
    plot_width = width - left - right
    plot_height = height - top - bottom
    max_value = max(
        max(values)
        for by_response in grouped.values()
        for values in by_response.values()
    )
    y_max = max(10, math.ceil(max_value / 5) * 5)

    def x_position(population_index: int, response: str) -> float:
        group_width = plot_width / len(populations)
        center = left + group_width * (population_index + 0.5)
        offset = -24 if response == "yes" else 24
        return center + offset

    def y_position(value: float) -> float:
        return top + plot_height - (value / y_max) * plot_height

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(width / 2, 30, "Melanoma miraclib PBMC relative frequencies", 20),
        svg_text(width / 2, 52, "Responders vs non-responders", 14),
    ]

    for tick in range(0, y_max + 1, 5):
        y = y_position(tick)
        elements.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width - right}" y2="{y:.1f}" '
            'stroke="#e6e6e6" stroke-width="1"/>'
        )
        elements.append(svg_text(left - 10, y + 4, str(tick), 12, "end"))
    elements.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" '
        'stroke="#333" stroke-width="1.5"/>'
    )
    elements.append(
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - right}" '
        f'y2="{height - bottom}" stroke="#333" stroke-width="1.5"/>'
    )

    colors = {"yes": "#2f78b7", "no": "#d95f02"}
    labels = {"yes": "Responder", "no": "Non-responder"}
    box_width = 34
    for population_index, population in enumerate(populations):
        group_width = plot_width / len(populations)
        center = left + group_width * (population_index + 0.5)
        elements.append(svg_text(center, height - 70, population, 13))

        for response in ("yes", "no"):
            stats = box_stats(grouped[population][response])
            x = x_position(population_index, response)
            y_q1 = y_position(float(stats["q1"]))
            y_q3 = y_position(float(stats["q3"]))
            y_med = y_position(float(stats["median"]))
            y_low = y_position(float(stats["lower"]))
            y_high = y_position(float(stats["upper"]))
            color = colors[response]

            elements.append(
                f'<line x1="{x:.1f}" y1="{y_high:.1f}" x2="{x:.1f}" '
                f'y2="{y_low:.1f}" stroke="{color}" stroke-width="2"/>'
            )
            elements.append(
                f'<line x1="{x - box_width / 3:.1f}" y1="{y_high:.1f}" '
                f'x2="{x + box_width / 3:.1f}" y2="{y_high:.1f}" '
                f'stroke="{color}" stroke-width="2"/>'
            )
            elements.append(
                f'<line x1="{x - box_width / 3:.1f}" y1="{y_low:.1f}" '
                f'x2="{x + box_width / 3:.1f}" y2="{y_low:.1f}" '
                f'stroke="{color}" stroke-width="2"/>'
            )
            elements.append(
                f'<rect x="{x - box_width / 2:.1f}" y="{y_q3:.1f}" '
                f'width="{box_width}" height="{y_q1 - y_q3:.1f}" '
                f'fill="{color}" fill-opacity="0.28" stroke="{color}" stroke-width="2"/>'
            )
            elements.append(
                f'<line x1="{x - box_width / 2:.1f}" y1="{y_med:.1f}" '
                f'x2="{x + box_width / 2:.1f}" y2="{y_med:.1f}" '
                'stroke="#111" stroke-width="2"/>'
            )

            for outlier in stats["outliers"][:80]:
                y_outlier = y_position(float(outlier))
                elements.append(
                    f'<circle cx="{x:.1f}" cy="{y_outlier:.1f}" r="2" '
                    f'fill="{color}" fill-opacity="0.5"/>'
                )

    elements.append(svg_text(24, top + plot_height / 2, "Relative frequency (%)", 13, "middle"))
    elements.append(
        '<g transform="translate(840,24)">'
        '<rect x="0" y="0" width="14" height="14" fill="#2f78b7" fill-opacity="0.28" stroke="#2f78b7"/>'
        '<text x="22" y="12" font-size="13" font-family="Arial, sans-serif">Responder</text>'
        '<rect x="118" y="0" width="14" height="14" fill="#d95f02" fill-opacity="0.28" stroke="#d95f02"/>'
        '<text x="140" y="12" font-size="13" font-family="Arial, sans-serif">Non-responder</text>'
        '</g>'
    )
    elements.append("</svg>")
    PLOT_PATH.write_text("\n".join(elements), encoding="utf-8")


def print_summary(results: list[TestResult]) -> None:
    print("Cohort: melanoma patients receiving miraclib, PBMC samples only")
    print(f"Statistics: Mann-Whitney U test, Benjamini-Hochberg FDR across {len(results)} tests")
    print()
    print(
        "population     resp_n  nonresp_n  resp_median  nonresp_median  "
        "median_diff  p_value   fdr_p    significant"
    )
    print(
        "------------  ------  ---------  -----------  --------------  "
        "-----------  --------  --------  -----------"
    )
    for result in results:
        print(
            f"{result.population:<12}  {result.responder_n:>6}  "
            f"{result.non_responder_n:>9}  {result.responder_median:>11.3f}  "
            f"{result.non_responder_median:>14.3f}  "
            f"{result.median_difference:>11.3f}  {result.p_value:>8.3g}  "
            f"{result.fdr_p_value:>8.3g}  {str(result.significant):>11}"
        )
    print()
    print(f"Wrote statistics to {STATS_PATH.name}")
    print(f"Wrote boxplot to {PLOT_PATH.name}")


def main() -> None:
    grouped = load_percentages()
    results = analyze(grouped)
    write_stats(results)
    write_boxplot(grouped)
    print_summary(results)


if __name__ == "__main__":
    main()
