PYTHON ?= python3
PORT ?= 8000

.PHONY: setup pipeline dashboard

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analyze_data.py
	$(PYTHON) statistical_analysis.py
	$(PYTHON) subset_analysis.py
	$(PYTHON) build_dashboard.py

dashboard: pipeline
	$(PYTHON) -m http.server $(PORT) --directory dashboard
