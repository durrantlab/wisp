<h1 align="center">WISP</h1>

<h4 align="center">Calculate and visualize allosteric pathways from molecular trajectories.</h4>

<h4 align="center" style="padding-bottom: 0.5em;"><a href="https://durrantlab.github.io/wisp/">Documentation</a></h4>

<p align="center">
    <a href="https://github.com/durrantlab/wisp/actions/workflows/build.yml">
        <img src="https://github.com/durrantlab/wisp/actions/workflows/build.yml/badge.svg" alt="Build Status ">
    </a>
    <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/wisp">
    <a href="https://codecov.io/gh/durrantlab/wisp">
        <img src="https://codecov.io/gh/durrantlab/wisp/branch/main/graph/badge.svg?token=74wLrsOMTD" alt="codecov">
    </a>
    <a href="https://github.com/durrantlab/wisp/releases">
        <img src="https://img.shields.io/github/v/release/durrantlab/wisp" alt="GitHub release (latest by date)">
    </a>
    <a href="https://github.com/durrantlab/wisp/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/durrantlab/wisp" alt="License">
    </a>
    <a href="https://github.com/durrantlab/wisp/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/durrantlab/wisp" alt="GitHub repo size">
    </a>
    <a href="https://github.com/psf/black" target="_blank">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
    </a>
    <a href="https://github.com/PyCQA/pylint" target="_blank">
        <img src="https://img.shields.io/badge/linting-pylint-yellowgreen" alt="Black style">
    </a>
</p>

## Installation

First, you must obtain the WISP code by cloning the [GitHub repository](https://github.com/durrantlab/wisp)

```bash
git clone https://github.com/durrantlab/wisp && cd wisp
```

and installing with the following command.

```bash
pip install .
```

This will install the `wisp` Python package in addition to the `wisp` command-line tool.

TODO: check that the VMD plugin installation works.

1. To use the VMD plugin, add these two lines of code to your `.vmdrc` file
   (usually located in your `~/` directory)
   - `set auto_path "$auto_path <your directory here>"  ; #(NOTE: this may require the full pathname)`
   - `vmd_install_extension wisp wisp_tk_cb "Analysis/Wisp"`
2. Now open VMD: Click Extensions > Analysis > WISP

## Usage

You may use the `wisp` command-line interface as shown below.

```bash
wisp tests/files/trajectory_20_frames.pdb --source_residues C_LEU_10 --sink_residues C_ASP_11
```

Or, you may use wisp as a library in Python.

```python
from wisp.run import run_wisp
from wisp.contexts import ContextManager

# Update context
context_manager = ContextManager()
context_manager.pdb_path = "tests/files/trajectory_20_frames.pdb"
context_manager.source_residues = ["C_LEU_10"]
context_manager.sink_residues = ["C_ASP_11"]

# Run wisp
paths = run_wisp(context_manager)
```

## Development

We use [pixi](https://pixi.sh/latest/) to manage Python environments and simplify the developer workflow.
Once you have [pixi](https://pixi.sh/latest/) installed, move into `wisp` directory (e.g., `cd wisp`) and install the  environment using the command

```bash
pixi install
```

Now you can activate the new virtual environment using

```sh
pixi shell
```

## Program Output

The output directory contains multiple output files. The best way to visualize
the output is to use a free program called VMD, which can be downloaded from
[http://www.ks.uiuc.edu/Research/vmd/](http://www.ks.uiuc.edu/Research/vmd/).

The WISP output can be automatically loaded into VMD using the TCL script named `visualize.tcl`. Assuming `vmd` is the full path to your installed VMD
executable, just run the following from the command line:

`vmd -e visualize.tcl`

If you prefer not to use the command line, simply run the `vmd` executable and
load the `visualize.vmd` file using "File->Load Visualization State..." from the main menu.

The above methods are very slow. If your output is so large that a faster option is required, the Tk Console can be used. Use "Extensions->Tk Console" from the VMD main menu to pull up the Tk Console. Then run the following command, with
the full path to `visualize.tcl` included if necessary:

`source visualize.tcl`

Regardless of the method you use to load in the WISP output, the visualization
will be the same. Individual pathways are shown as tubes (i.e., "wisps"), the
protein is shown in ribbon representation, and protein residues that participate in any path are shown in licorice representation.

The WISP output directory contains a number of other files as well. Here are
descriptions of each:

- `log.txt`: Details describing WISP execution.
- `parameters_used.txt`: The WISP parameters used to generate the output.
- `average_structure.pdb`: The average structure of your PDB trajectory.
- `draw_frame.pdb`: If the user requests that a separate single-structure PDB
  file be used for calculating node and wisp positions, that file is saved as
  "draw_frame.pdb". Otherwise, the average structure is used.
- `functionalized_matrix_with_contact_map_applied.pickle`: A python pickle
  file that contains the matrix obtained by multiplying a functionalized
  correlation matrix and a contact map. This file is not human readable but
  can be loaded into WISP for use in subsequent runs with the
  `-load_wisp_saved_matrix` and `-wisp_saved_matrix_filename` parameters.
  Thus, the matrix needs only to be calculated once for each trajectory,
  rather than every time WISP is executed. Use `wisp -help` for more
  information.
- `contact_map_matrix.txt`: A human readable representation of the contact
  map. If the user wishes to generate their own contact map rather than
  letting WISP generate one automatically, a custom contact map formatted like
  this one can be loaded into WISP using the
  `-contact_map_path parameter`.
- `functionalized_correlation_matrix.txt`: A human readable representation of
  the functionalized correlation matrix, prior to multiplication by the
  contact map. If the user wishes to generate their own functionalized
  correlation matrix rather than letting WISP generate one automatically, a
  custom matrix formatted like this one can be loaded into WISP using the
  `-functionalized_matrix_filename parameter`.
- `simply_formatted_paths.txt`: A simple list of path lengths and nodes. The
  first column contains the lengths, and all following columns contain node
  indices. This file may be helpful for subsequent statistical analyses of the
  WISP output. Note that the `simply_formatted_paths.txt` output file reindexes the residues. See the `visualize.tcl` file instead for a more human-readable output.

## Deploying

We use [bump-my-version](https://github.com/callowayproject/bump-my-version) to release a new version.
This will create a git tag that is used by [poetry-dynamic-version](https://github.com/mtkennerly/poetry-dynamic-versioning) to generate version strings and update `CHANGELOG.md`.

For example, to bump the `minor` version you would run the following command.

```bash
poetry run bump-my-version bump minor
```

After releasing a new version, you need to push and include all tags.

```bash
git push --follow-tags
```

## Citation

If you use WISP in your work, please cite:

> A.T. Van Wart, J.D. Durrant, L. Votapka, R.E. Amaro. Weighted implementation of suboptimal paths (WISP): An optimized algorithm and tool for dynamical network analysis, J. Chem. Theory Comput. 10 (2014) 511-517

## License

It is licensed under the [Academic Free License 3.0](http://opensource.org/licenses/AFL-3.0).
