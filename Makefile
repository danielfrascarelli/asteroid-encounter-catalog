# Catálogo de Encuentros Cercanos — atajos para los comandos más comunes.
# Todo se ejecuta dentro de Docker, así que no necesitas Python local.

.PHONY: help build pipeline pipeline-rebound \
        download-mpcorb download-mpcorb-2015 download-fienga download-galad \
        validate validate-fienga validate-galad validate-jpl validate-all \
        compare test lint format shell dashboard clean-cache

# ──────────────────────────────────────────────────────────────────────────────
# Default: print this help
# ──────────────────────────────────────────────────────────────────────────────
help:
	@echo "Common targets:"
	@echo "  make build              Build the Docker image (only needed after pyproject.toml changes)"
	@echo "  make pipeline           Run the full encounter detection pipeline (config.yaml)"
	@echo "  make pipeline-rebound   Run the pipeline with N-body propagation (writes config.local.yaml)"
	@echo ""
	@echo "  make download-mpcorb           Fetch the current MPCORB.DAT from MPC"
	@echo "  make download-mpcorb-2015      Fetch a 2015 MPCORB snapshot via Wayback Machine"
	@echo "  make download-fienga           Fetch Fienga 2003 catalog from VizieR"
	@echo "  make download-galad            Scrape Galád 2002 tables from the A&A HTML"
	@echo ""
	@echo "  make validate-all              Run all four validators (Fienga, Galád, JPL, compare)"
	@echo "  make validate-fienga           Cross-match against Fienga 2003"
	@echo "  make validate-galad            Cross-match against Galád 2002"
	@echo "  make validate-jpl              Three-way cross-check against JPL Horizons"
	@echo "  make compare                   Kepler vs rebound side-by-side"
	@echo ""
	@echo "  make test               Run the pytest suite (194 tests)"
	@echo "  make lint               ruff + black --check + mypy"
	@echo "  make format             Apply ruff --fix + black"
	@echo "  make shell              Drop into a bash shell inside the pipeline container"
	@echo "  make dashboard          Bring up the Streamlit dashboard on localhost:8501"
	@echo "  make clean-cache        Delete data/cache/ contents"

# ──────────────────────────────────────────────────────────────────────────────
# Build
# ──────────────────────────────────────────────────────────────────────────────
build:
	docker compose build

# ──────────────────────────────────────────────────────────────────────────────
# Pipeline runs
# ──────────────────────────────────────────────────────────────────────────────
pipeline:
	docker compose run --rm pipeline python -m scripts.run_pipeline

pipeline-rebound:
	@printf 'propagation:\n  method: "rebound"\n  rebound:\n    integrator: "whfast"\n    include_planets:\n      - "sun"\n      - "jupiter"\n      - "saturn"\n    include_major_asteroids: false\noutput:\n  filename: "encounters_catalog_rebound"\n' > config.local.yaml
	docker compose run --rm pipeline python -m scripts.run_pipeline

# ──────────────────────────────────────────────────────────────────────────────
# Data downloads
# ──────────────────────────────────────────────────────────────────────────────
download-mpcorb:
	docker compose run --rm pipeline python -m scripts.download_mpcorb

download-mpcorb-2015:
	docker compose run --rm pipeline python -m scripts.download_mpcorb_historical --year 2015 --month 6

download-fienga:
	docker compose run --rm pipeline python -m scripts.download_fienga_2003

download-galad:
	docker compose run --rm pipeline python -m scripts.download_galad_2002

# ──────────────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────────────
validate: validate-fienga validate-galad validate-jpl

validate-all: validate compare

validate-fienga:
	docker compose run --rm pipeline python -m scripts.validate_fienga_2003

validate-galad:
	docker compose run --rm pipeline python -m scripts.validate_galad_2002

validate-jpl:
	docker compose run --rm pipeline python -m scripts.validate_jpl_horizons

compare:
	docker compose run --rm pipeline python -m scripts.compare_kepler_vs_rebound

# ──────────────────────────────────────────────────────────────────────────────
# Development
# ──────────────────────────────────────────────────────────────────────────────
test:
	docker compose run --rm test pytest tests/ -v

lint:
	docker compose run --rm pipeline sh -c "ruff check . && black --check . && mypy src scripts"

format:
	docker compose run --rm pipeline sh -c "ruff check . --fix && black ."

shell:
	docker compose run --rm pipeline bash

dashboard:
	docker compose up dashboard

# ──────────────────────────────────────────────────────────────────────────────
# Cleanup
# ──────────────────────────────────────────────────────────────────────────────
clean-cache:
	rm -rf data/cache/*.npy data/cache/*.json
	@echo "Trajectory cache cleared (data/cache/)"
