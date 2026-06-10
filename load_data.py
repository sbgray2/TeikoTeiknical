#!/usr/bin/env python3
"""Create and populate a SQLite database from cell-count.csv."""

from __future__ import annotations

import csv
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

CELL_POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")


SCHEMA = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS cell_populations;
DROP TABLE IF EXISTS patients;

CREATE TABLE patients (
    subject_id TEXT PRIMARY KEY,
    project TEXT NOT NULL,
    indication TEXT NOT NULL,
    age INTEGER,
    gender TEXT,
    treatment TEXT NOT NULL,
    response TEXT NOT NULL
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES patients(subject_id)
);

CREATE TABLE cell_populations (
    population_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL,
    population_id INTEGER NOT NULL,
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    FOREIGN KEY (population_id) REFERENCES cell_populations(population_id)
);

CREATE INDEX idx_samples_subject_id ON samples(subject_id);
CREATE INDEX idx_samples_time_from_treatment_start
    ON samples(time_from_treatment_start);
CREATE INDEX idx_cell_counts_population_id ON cell_counts(population_id);
"""


def require_columns(fieldnames: list[str] | None) -> None:
    required = {
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample",
        "sample_type",
        "time_from_treatment_start",
        *CELL_POPULATIONS,
    }
    missing = required - set(fieldnames or [])
    if missing:
        missing_columns = ", ".join(sorted(missing))
        raise ValueError(f"{CSV_PATH.name} is missing required columns: {missing_columns}")


def to_int(value: str, column_name: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"Invalid integer in column '{column_name}': {value!r}") from exc


def initialize_database(connection: sqlite3.Connection) -> dict[str, int]:
    connection.executescript(SCHEMA)
    connection.executemany(
        "INSERT INTO cell_populations (name) VALUES (?)",
        ((population,) for population in CELL_POPULATIONS),
    )
    rows = connection.execute("SELECT population_id, name FROM cell_populations").fetchall()
    return {name: population_id for population_id, name in rows}


def load_csv(connection: sqlite3.Connection, population_ids: dict[str, int]) -> int:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Could not find {CSV_PATH}")

    inserted_rows = 0
    with CSV_PATH.open(newline="", encoding="utf-8") as csv_file:
        reader = csv.DictReader(csv_file)
        require_columns(reader.fieldnames)

        for row in reader:
            connection.execute(
                """
                INSERT OR IGNORE INTO patients (
                    subject_id, project, indication, age, gender, treatment, response
                )
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    row["subject"],
                    row["project"],
                    row["condition"],
                    to_int(row["age"], "age"),
                    row["sex"],
                    row["treatment"],
                    row["response"],
                ),
            )
            connection.execute(
                """
                INSERT INTO samples (
                    sample_id, subject_id, sample_type, time_from_treatment_start
                )
                VALUES (?, ?, ?, ?)
                """,
                (
                    row["sample"],
                    row["subject"],
                    row["sample_type"],
                    to_int(
                        row["time_from_treatment_start"],
                        "time_from_treatment_start",
                    ),
                ),
            )
            connection.executemany(
                """
                INSERT INTO cell_counts (sample_id, population_id, count)
                VALUES (?, ?, ?)
                """,
                (
                    (
                        row["sample"],
                        population_ids[population],
                        to_int(row[population], population),
                    )
                    for population in CELL_POPULATIONS
                ),
            )
            inserted_rows += 1

    return inserted_rows


def main() -> None:
    with sqlite3.connect(DB_PATH) as connection:
        connection.execute("PRAGMA foreign_keys = ON")
        population_ids = initialize_database(connection)
        sample_count = load_csv(connection, population_ids)

    print(f"Created {DB_PATH.name} with {sample_count} samples.")


if __name__ == "__main__":
    main()
