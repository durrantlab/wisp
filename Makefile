SHELL := /usr/bin/env bash
PYTHON_VERSION := 3.11
PYTHON_VERSION_CONDENSED := 311
PACKAGE_NAME := wisp
REPO_PATH := $(shell git rev-parse --show-toplevel)
CONDA_PATH := $(REPO_PATH)/.venv
CONDA := conda run -p $(CONDA_PATH)

###   ENVIRONMENT   ###

.PHONY: conda-setup
conda-setup:
	conda create -y -p $(CONDA_PATH) python=$(PYTHON_VERSION)
	conda install -y conda-lock -p $(CONDA_PATH)
	conda install -y -c conda-forge poetry pre-commit tomli tomli-w -p $(CONDA_PATH)
	$(CONDA) pip install conda_poetry_liaison

.PHONY: write-conda-lock
write-conda-lock:
	$(CONDA) conda env export --from-history | grep -v "^prefix" | grep -v "^name" > environment.yml
	$(CONDA) conda-lock -f environment.yml -p linux-64 -p osx-64 -p win-64
	$(CONDA) cpl-deps $(REPO_PATH)/pyproject.toml --env_path $(CONDA_PATH)
	$(CONDA) cpl-clean $(CONDA_PATH)

.PHONY: from-conda-lock
from-conda-lock:
	$(CONDA) conda-lock install --prefix $(REPO_PATH)/.venv $(REPO_PATH)/conda-lock.yml
	$(CONDA) cpl-clean $(CONDA_PATH)

# notify_poetry_of_conda.py
.PHONY: pre-commit-install
pre-commit-install:
	$(CONDA) pre-commit install

# Reads `pyproject.toml`, solves environment, then writes lock file.
.PHONY: poetry-lock
poetry-lock:
	$(CONDA) poetry lock --no-interaction
	$(CONDA) poetry export --without-hashes > requirements.txt

.PHONY: install
install:
	$(CONDA) poetry install --no-interaction

.PHONY: validate
validate:
	$(CONDA) pre-commit run --all-files

.PHONY: codestyle
codestyle:
	$(CONDA) pyupgrade --exit-zero-even-if-changed --py311-plus **/*.py
	$(CONDA) isort --settings-path pyproject.toml ./
	$(CONDA) black --config pyproject.toml ./

.PHONY: formatting
formatting: codestyle

.PHONY: test
test:
	$(CONDA) pytest -c pyproject.toml --cov=$(PACKAGE_NAME) --cov-report=xml tests/