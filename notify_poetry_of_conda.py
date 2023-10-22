#!/usr/bin/env python3

"""
Add all packages in the conda environment to pyproject.toml.
This will overwrite these files.
"""

import os
import subprocess

import tomli
import tomli_w

REPO_ROOT = os.path.dirname(os.path.realpath(__file__))
CONDA_PATH = os.path.join(REPO_ROOT, ".venv/")
PYPROJECT_TOML = os.path.join(REPO_ROOT, "pyproject.toml")
CONDA_GROUP_LABEL = "tool.poetry.group.conda.dependencies"

results = subprocess.run(
    ["conda", "run", "-p", CONDA_PATH, "pip", "list"],
    check=True,
    capture_output=True,
    text=True,
)
conda_packages = results.stdout.split("\n")[3:]


def check_if_include_package(conda_str):
    r"""Checks if we should include a package in pyproject.toml"""
    if "" == conda_str.strip():
        return False
    if "_" == conda_str[0]:
        return False
    return True


with open(PYPROJECT_TOML, mode="rb") as f:
    toml_file = tomli.load(f)

# Overwrite (or assign) conda group dependencies
toml_file[CONDA_GROUP_LABEL] = {}
for conda_package in conda_packages:
    if not check_if_include_package(conda_package):
        continue
    package_name, package_version = conda_package.strip().split()
    toml_file[CONDA_GROUP_LABEL][package_name] = f"^{package_version}"

with open(PYPROJECT_TOML, "wb") as f:
    tomli_w.dump(toml_file, f)
