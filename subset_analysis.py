#!/usr/bin/env python3
"""Summarize baseline melanoma PBMC samples from miraclib-treated patients."""

from __future__ import annotations

import csv
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
BASELINE_SAMPLES_PATH = ROOT / "baseline_melanoma_miraclib_pbmc_samples.csv"
PROJECT_SUMMARY_PATH = ROOT / "baseline_project_summary.csv"
RESPONSE_SUMMARY_PATH = ROOT / "baseline_response_summary.csv"
GENDER_SUMMARY_PATH = ROOT / "baseline_gender_summary.csv"


BASELINE_QUERY = """
SELECT
    s.sample_id AS sample,
    p.subject_id,
    p.project,
    p.indication,
    p.treatment,
    s.sample_type,
    s.time_from_treatment_start,
    p.response,
    p.gender
FROM patients AS p
JOIN samples AS s
    ON s.subject_id = p.subject_id
WHERE p.indication = 'melanoma'
    AND p.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
ORDER BY p.project, p.subject_id, s.sample_id;
"""

PROJECT_SUMMARY_QUERY = """
SELECT
    p.project,
    COUNT(*) AS sample_count
FROM patients AS p
JOIN samples AS s
    ON s.subject_id = p.subject_id
WHERE p.indication = 'melanoma'
    AND p.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
GROUP BY p.project
ORDER BY p.project;
"""

RESPONSE_SUMMARY_QUERY = """
SELECT
    p.response,
    COUNT(DISTINCT p.subject_id) AS subject_count
FROM patients AS p
JOIN samples AS s
    ON s.subject_id = p.subject_id
WHERE p.indication = 'melanoma'
    AND p.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
GROUP BY p.response
ORDER BY p.response;
"""

GENDER_SUMMARY_QUERY = """
SELECT
    p.gender,
    COUNT(DISTINCT p.subject_id) AS subject_count
FROM patients AS p
JOIN samples AS s
    ON s.subject_id = p.subject_id
WHERE p.indication = 'melanoma'
    AND p.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
GROUP BY p.gender
ORDER BY p.gender;
"""


def print_table(title: str, rows: list[sqlite3.Row], headers: tuple[str, ...]) -> None:
    print(title)
    if not rows:
        print("No rows found.")
        print()
        return

    widths = {
        header: max(len(header), *(len(str(row[header])) for row in rows))
        for header in headers
    }
    print("  ".join(header.ljust(widths[header]) for header in headers))
    print("  ".join("-" * widths[header] for header in headers))
    for row in rows:
        print("  ".join(str(row[header]).ljust(widths[header]) for header in headers))
    print()


def write_baseline_samples(rows: list[sqlite3.Row]) -> None:
    headers = (
        "sample",
        "subject_id",
        "project",
        "indication",
        "treatment",
        "sample_type",
        "time_from_treatment_start",
        "response",
        "gender",
    )
    with BASELINE_SAMPLES_PATH.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({header: row[header] for header in headers})


def write_summary(path: Path, rows: list[sqlite3.Row], headers: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({header: row[header] for header in headers})


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {DB_PATH.name}. Run `python load_data.py` first."
        )

    with sqlite3.connect(DB_PATH) as connection:
        connection.row_factory = sqlite3.Row
        baseline_samples = connection.execute(BASELINE_QUERY).fetchall()
        project_summary = connection.execute(PROJECT_SUMMARY_QUERY).fetchall()
        response_summary = connection.execute(RESPONSE_SUMMARY_QUERY).fetchall()
        gender_summary = connection.execute(GENDER_SUMMARY_QUERY).fetchall()

    write_baseline_samples(baseline_samples)
    write_summary(PROJECT_SUMMARY_PATH, project_summary, ("project", "sample_count"))
    write_summary(RESPONSE_SUMMARY_PATH, response_summary, ("response", "subject_count"))
    write_summary(GENDER_SUMMARY_PATH, gender_summary, ("gender", "subject_count"))

    print("Subset: melanoma PBMC baseline samples from miraclib-treated patients")
    print(f"Baseline samples identified: {len(baseline_samples)}")
    print(f"Wrote sample list to {BASELINE_SAMPLES_PATH.name}")
    print(
        "Wrote summaries to "
        f"{PROJECT_SUMMARY_PATH.name}, {RESPONSE_SUMMARY_PATH.name}, "
        f"{GENDER_SUMMARY_PATH.name}"
    )
    print()

    print_table(
        "Samples from each project",
        project_summary,
        ("project", "sample_count"),
    )
    print_table(
        "Responder/non-responder subjects",
        response_summary,
        ("response", "subject_count"),
    )
    print_table(
        "Male/female subjects",
        gender_summary,
        ("gender", "subject_count"),
    )


if __name__ == "__main__":
    main()
