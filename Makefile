.PHONY: dev test lint fmt clean docs-api docs-figures docs-dev docs-build docs-preview

VENV   := .venv
PYTHON := $(VENV)/bin/python
PIP    := $(VENV)/bin/pip

$(VENV)/bin/activate:
	python3 -m venv $(VENV)

## Install package in editable mode with dev extras
dev: $(VENV)/bin/activate
	$(PIP) install -e ".[dev]"

## Run test suite
test: dev
	$(VENV)/bin/pytest

## Check linting (no fixes)
LINTED := src/shore/ tests/ examples/tutorials/ examples/cookbook/ docs/scripts/api_docs.py docs/scripts/figures/

lint: dev
	$(VENV)/bin/ruff check $(LINTED)
	$(VENV)/bin/ruff format --check $(LINTED)
	$(PYTHON) docs/scripts/api_docs.py --check

## Auto-fix lint and format
fmt: dev
	$(VENV)/bin/ruff check --fix $(LINTED)
	$(VENV)/bin/ruff format $(LINTED)

## Remove build artifacts
clean:
	rm -rf dist/ build/ *.egg-info .pytest_cache .ruff_cache .coverage coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null; true

## Regenerate the Python API reference from the docstrings
docs-api:
	$(PYTHON) docs/scripts/api_docs.py

## Redraw the docs figures (light and dark; needs pyvista, cad and gui extras)
docs-figures:
	QT_QPA_PLATFORM=offscreen $(PYTHON) docs/scripts/figures/make.py

## Start VitePress dev server (hot-reload)
docs-dev:
	cd docs && npm run docs:dev

## Build static VitePress site to docs/.vitepress/dist/
docs-build: docs-api
	cd docs && npm run docs:build

## Preview the built static site
docs-preview:
	cd docs && npm run docs:preview
