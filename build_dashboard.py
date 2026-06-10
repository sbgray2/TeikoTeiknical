#!/usr/bin/env python3
"""Build a lightweight HTML dashboard from pipeline outputs."""

from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
DASHBOARD_DIR = ROOT / "dashboard"
DASHBOARD_PATH = DASHBOARD_DIR / "index.html"
STATS_PATH = ROOT / "miraclib_melanoma_pbmc_stats.csv"
PLOT_PATH = ROOT / "miraclib_melanoma_pbmc_boxplot.svg"
BASELINE_SAMPLES_PATH = ROOT / "baseline_melanoma_miraclib_pbmc_samples.csv"
PROJECT_SUMMARY_PATH = ROOT / "baseline_project_summary.csv"
RESPONSE_SUMMARY_PATH = ROOT / "baseline_response_summary.csv"
GENDER_SUMMARY_PATH = ROOT / "baseline_gender_summary.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def database_metrics() -> dict[str, int]:
    with sqlite3.connect(DB_PATH) as connection:
        return {
            "patients": connection.execute("SELECT COUNT(*) FROM patients").fetchone()[0],
            "samples": connection.execute("SELECT COUNT(*) FROM samples").fetchone()[0],
            "cell_counts": connection.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0],
            "populations": connection.execute(
                "SELECT COUNT(*) FROM cell_populations"
            ).fetchone()[0],
        }


def json_script(name: str, data: object) -> str:
    payload = json.dumps(data).replace("<", "\\u003c").replace("&", "\\u0026")
    return (
        f'<script id="{name}" type="application/json">'
        f"{payload}</script>"
    )


def build_html() -> str:
    metrics = database_metrics()
    stats = read_csv(STATS_PATH)
    project_summary = read_csv(PROJECT_SUMMARY_PATH)
    response_summary = read_csv(RESPONSE_SUMMARY_PATH)
    gender_summary = read_csv(GENDER_SUMMARY_PATH)
    baseline_samples = read_csv(BASELINE_SAMPLES_PATH)
    boxplot_svg = PLOT_PATH.read_text(encoding="utf-8")

    significant_count = sum(row["significant"] == "True" for row in stats)
    baseline_count = len(baseline_samples)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Loblaw Bio Immune Cell Analysis</title>
  <style>
    :root {{
      --ink: #17202a;
      --muted: #5d6d7e;
      --line: #d9e2ec;
      --panel: #f7f9fb;
      --accent: #1f77b4;
      --orange: #d95f02;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      font-family: Arial, sans-serif;
      background: #ffffff;
    }}
    header {{
      padding: 28px clamp(18px, 5vw, 56px);
      border-bottom: 1px solid var(--line);
      background: #f9fbfd;
    }}
    h1 {{ margin: 0 0 8px; font-size: 30px; }}
    h2 {{ margin: 0 0 14px; font-size: 22px; }}
    p {{ line-height: 1.5; }}
    main {{
      width: min(1180px, calc(100% - 32px));
      margin: 24px auto 42px;
    }}
    section {{
      margin: 28px 0;
      padding-bottom: 26px;
      border-bottom: 1px solid var(--line);
    }}
    .muted {{ color: var(--muted); }}
    .kpis {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
      gap: 12px;
      margin-top: 20px;
    }}
    .kpi {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      background: white;
    }}
    .kpi strong {{ display: block; font-size: 26px; margin-bottom: 4px; }}
    .controls {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      align-items: center;
      margin: 12px 0;
    }}
    button {{
      border: 1px solid var(--line);
      border-radius: 8px;
      background: white;
      color: var(--ink);
      padding: 8px 12px;
      cursor: pointer;
    }}
    button.active {{
      border-color: var(--accent);
      background: #eaf4fc;
      color: #0f568a;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 12px;
      font-size: 14px;
    }}
    th, td {{
      border-bottom: 1px solid var(--line);
      padding: 9px 8px;
      text-align: left;
      white-space: nowrap;
    }}
    th {{ background: var(--panel); }}
    .plot {{
      overflow-x: auto;
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 8px;
      background: white;
    }}
    .plot svg {{
      width: 100%;
      min-width: 760px;
      height: auto;
      display: block;
    }}
    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 18px;
    }}
    .panel {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      background: white;
    }}
    .tag {{
      display: inline-block;
      border-radius: 999px;
      padding: 3px 8px;
      background: #edf2f7;
      color: var(--muted);
      font-size: 12px;
    }}
    .sig-false {{ color: #7b241c; }}
    .sig-true {{ color: #0b6b3a; font-weight: 700; }}
    footer {{
      width: min(1180px, calc(100% - 32px));
      margin: 0 auto 28px;
      color: var(--muted);
      font-size: 13px;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Loblaw Bio Immune Cell Analysis</h1>
    <p class="muted">SQLite-backed analysis of immune cell populations, response to miraclib, and melanoma PBMC baseline subsets.</p>
    <div class="kpis">
      <div class="kpi"><strong>{metrics["patients"]:,}</strong><span>patients</span></div>
      <div class="kpi"><strong>{metrics["samples"]:,}</strong><span>samples</span></div>
      <div class="kpi"><strong>{metrics["cell_counts"]:,}</strong><span>cell count records</span></div>
      <div class="kpi"><strong>{metrics["populations"]}</strong><span>immune populations</span></div>
      <div class="kpi"><strong>{baseline_count:,}</strong><span>baseline melanoma PBMC miraclib samples</span></div>
      <div class="kpi"><strong>{significant_count}</strong><span>FDR-significant populations</span></div>
    </div>
  </header>

  <main>
    <section>
      <h2>Responder vs Non-responder Frequencies</h2>
      <p class="muted">Melanoma patients receiving miraclib, PBMC samples only. Percentages are relative frequencies within each sample.</p>
      <div class="plot">{boxplot_svg}</div>
    </section>

    <section>
      <h2>Statistical Results</h2>
      <div class="controls">
        <button class="active" data-filter="all">All populations</button>
        <button data-filter="significant">FDR significant only</button>
        <span class="tag">Mann-Whitney U + Benjamini-Hochberg FDR</span>
      </div>
      <table id="stats-table"></table>
    </section>

    <section>
      <h2>Baseline Subset</h2>
      <p class="muted">Melanoma PBMC baseline samples at time 0 from miraclib-treated patients.</p>
      <div class="summary-grid">
        <div class="panel">
          <h2>Projects</h2>
          <table id="project-table"></table>
        </div>
        <div class="panel">
          <h2>Response</h2>
          <table id="response-table"></table>
        </div>
        <div class="panel">
          <h2>Gender</h2>
          <table id="gender-table"></table>
        </div>
      </div>
    </section>
  </main>

  <footer>
    Generated by <code>make pipeline</code>. Source data: <code>cell-count.csv</code>.
  </footer>

  {json_script("stats-data", stats)}
  {json_script("project-data", project_summary)}
  {json_script("response-data", response_summary)}
  {json_script("gender-data", gender_summary)}
  <script>
    const stats = JSON.parse(document.getElementById("stats-data").textContent);
    const projectSummary = JSON.parse(document.getElementById("project-data").textContent);
    const responseSummary = JSON.parse(document.getElementById("response-data").textContent);
    const genderSummary = JSON.parse(document.getElementById("gender-data").textContent);

    function renderTable(id, rows, columns) {{
      const table = document.getElementById(id);
      const header = "<thead><tr>" + columns.map(([key, label]) => `<th>${{label}}</th>`).join("") + "</tr></thead>";
      const body = rows.map(row => {{
        const cells = columns.map(([key]) => {{
          const value = row[key];
          const cls = key === "significant" ? ` class="sig-${{String(value).toLowerCase()}}"` : "";
          return `<td${{cls}}>${{value}}</td>`;
        }}).join("");
        return `<tr>${{cells}}</tr>`;
      }}).join("");
      table.innerHTML = header + `<tbody>${{body}}</tbody>`;
    }}

    function renderStats(filter) {{
      const rows = filter === "significant"
        ? stats.filter(row => row.significant === "True")
        : stats;
      renderTable("stats-table", rows, [
        ["population", "Population"],
        ["responder_n", "Responder n"],
        ["non_responder_n", "Non-responder n"],
        ["responder_median", "Responder median %"],
        ["non_responder_median", "Non-responder median %"],
        ["median_difference", "Median diff"],
        ["p_value", "p-value"],
        ["fdr_p_value", "FDR p-value"],
        ["significant", "Significant"]
      ]);
    }}

    document.querySelectorAll("button[data-filter]").forEach(button => {{
      button.addEventListener("click", () => {{
        document.querySelectorAll("button[data-filter]").forEach(item => item.classList.remove("active"));
        button.classList.add("active");
        renderStats(button.dataset.filter);
      }});
    }});

    renderStats("all");
    renderTable("project-table", projectSummary, [["project", "Project"], ["sample_count", "Samples"]]);
    renderTable("response-table", responseSummary, [["response", "Response"], ["subject_count", "Subjects"]]);
    renderTable("gender-table", genderSummary, [["gender", "Gender"], ["subject_count", "Subjects"]]);
  </script>
</body>
</html>
"""


def main() -> None:
    required_paths = [
        DB_PATH,
        STATS_PATH,
        PLOT_PATH,
        BASELINE_SAMPLES_PATH,
        PROJECT_SUMMARY_PATH,
        RESPONSE_SUMMARY_PATH,
        GENDER_SUMMARY_PATH,
    ]
    missing = [path.name for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing pipeline outputs: "
            + ", ".join(missing)
            + ". Run `python statistical_analysis.py` and `python subset_analysis.py` first."
        )

    DASHBOARD_DIR.mkdir(exist_ok=True)
    DASHBOARD_PATH.write_text(build_html(), encoding="utf-8")
    print(f"Wrote dashboard to {DASHBOARD_PATH}")


if __name__ == "__main__":
    main()
