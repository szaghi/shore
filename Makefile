.PHONY: dev test lint fmt clean docs-dev docs-build docs-preview

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
	$(VENV)/bin/pytest --cov=shore --cov-report=term-missing

## Check linting (no fixes)
lint: dev
	$(VENV)/bin/ruff check src/shore/ tests/
	$(VENV)/bin/ruff format --check src/shore/ tests/

## Auto-fix lint and format
fmt: dev
	$(VENV)/bin/ruff check --fix src/shore/ tests/
	$(VENV)/bin/ruff format src/shore/ tests/

## Remove build artifacts
clean:
	rm -rf dist/ build/ *.egg-info .pytest_cache .ruff_cache .coverage coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null; true

## Start VitePress dev server (hot-reload)
docs-dev:
	cd docs && npm run docs:dev

## Build static VitePress site to docs/.vitepress/dist/
docs-build:
	cd docs && npm run docs:build

## Preview the built static site
docs-preview:
	cd docs && npm run docs:preview
