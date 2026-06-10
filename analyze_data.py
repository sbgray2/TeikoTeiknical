#!/usr/bin/env python3
"""Display relative cell population frequencies for each sample."""

from __future__ import annotations

import sqlite3
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
FREQUENCY_OUTPUT_PATH = ROOT / "cell_population_frequencies.csv"


FREQUENCY_QUERY = """
WITH sample_totals AS (
    SELECT
        sample_id,
        SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT
    cc.sample_id AS sample,
    st.total_count,
    cp.name AS population,
    cc.count,
    ROUND(cc.count * 100.0 / st.total_count, 2) AS percentage
FROM cell_counts AS cc
JOIN sample_totals AS st
    ON st.sample_id = cc.sample_id
JOIN cell_populations AS cp
    ON cp.population_id = cc.population_id
ORDER BY cc.sample_id, cp.population_id;
"""


def print_table(rows: list[sqlite3.Row]) -> None:
    headers = ("sample", "total_count", "population", "count", "percentage")
    widths = {
        header: max(len(header), *(len(str(row[header])) for row in rows))
        for header in headers
    }

    print(
        "  ".join(
            header.ljust(widths[header])
            if header in {"sample", "population"}
            else header.rjust(widths[header])
            for header in headers
        )
    )
    print(
        "  ".join(
            "-" * widths[header]
            for header in headers
        )
    )

    for row in rows:
        print(
            "  ".join(
                (
                    str(row[header]).ljust(widths[header])
                    if header in {"sample", "population"}
                    else str(row[header]).rjust(widths[header])
                )
                for header in headers
            )
        )


def write_frequencies(rows: list[sqlite3.Row]) -> None:
    headers = ("sample", "total_count", "population", "count", "percentage")
    with FREQUENCY_OUTPUT_PATH.open("w", newline="", encoding="utf-8") as csv_file:
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
        rows = connection.execute(FREQUENCY_QUERY).fetchall()

    write_frequencies(rows)
    print(f"Wrote {len(rows)} rows to {FREQUENCY_OUTPUT_PATH.name}")
    print("Preview:")
    print_table(rows[:20])


if __name__ == "__main__":
    main()
