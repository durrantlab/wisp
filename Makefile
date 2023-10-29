SHELL := /usr/bin/env bash
PYTHON_VERSION := 3.11
PYTHON_VERSION_CONDENSED := 311
PACKAGE_NAME := wisp
REPO_PATH := $(shell git rev-parse --show-toplevel)
CONDA_NAME := $(PACKAGE_NAME)-dev
CONDA_BASE_PATH = $(shell conda info --base)
CONDA_PATH := $(CONDA_BASE_PATH)/envs/$(CONDA_NAME)
CONDA := conda run -n $(CONDA_NAME)

###   ENVIRONMENT   ###

.PHONY: conda-setup
conda-setup:
	conda remove -y --name $(CONDA_NAME) --all
	conda create -y -n $(CONDA_NAME) python=$(PYTHON_VERSION)
	conda install -y conda-lock -n $(CONDA_NAME)
	conda install -y -c conda-forge poetry pre-commit tomli tomli-w -n $(CONDA_NAME)
	$(CONDA) pip install conda_poetry_liaison

.PHONY: write-conda-lock
write-conda-lock:
	$(CONDA) conda env export --from-history | grep -v "^prefix" > environment.yml
	$(CONDA) conda-lock -f environment.yml -p linux-64 -p osx-64 -p win-64
	$(CONDA) cpl-deps $(REPO_PATH)/pyproject.toml --env_path $(CONDA_PATH)
	$(CONDA) cpl-clean $(CONDA_PATH)

.PHONY: from-conda-lock
from-conda-lock:
	$(CONDA) conda-lock install --prefix $(REPO_PATH)/.venv $(REPO_PATH)/conda-lock.yml
	$(CONDA) pip install conda_poetry_liaison
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
	- $(CONDA) pre-commit run --all-files

.PHONY: formatting
formatting:
	- $(CONDA) pyupgrade --exit-zero-even-if-changed --py311-plus **/*.py
	- $(CONDA) isort --settings-path pyproject.toml ./
	- $(CONDA) black --config pyproject.toml ./

.PHONY: test
test:
	$(CONDA) pytest -c pyproject.toml --cov=$(PACKAGE_NAME) --cov-report=xml tests/

.PHONY: refresh
refresh: conda-setup write-conda-lock pre-commit-install poetry-lock install formatting validate
