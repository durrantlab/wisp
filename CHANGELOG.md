Changes
=======

1.3
---

* Minor updates to deal with PDB files that have both ENDMDL and END tags.
* Formatting improvements to `README.md`.
* Changed development-branch name to `main`.

1.2
---

* Reformatted code with [Black Formatter](https://github.com/psf/black).
* Updated email in copyright header.
* Updated `contact_map_distance_limit` to 4.5 in example scripts, help
  descriptions, etc., to avoid confusion (given new default, see version 1.1
  changes).
* Added download and citation information to the standard output.
* Investigated Python2/3 differences. In brief, we recently noticed that WISP
  gives slightly different output when running under Python2 (numpy 1.16.0,
  networkx 2.2, scipy 1.2.0) vs. Python3 (numpy 1.17.4, networkx 2.2, scipy
  1.3.1). We believe [differences in the two numpy versions are
  responsible](https://numpy.org/doc/stable/release/1.17.0-notes.html#float16-subnormal-rounding),
  though it is possible that similar differences may arise when running on
  different operating systems, etc. Regardless, the top-ranked paths are the
  same, but lower-ranked paths sometimes aren't. After some investigation and
  test updates to the code, we now believe subtle rounding-error differences
  to be the cause. For consistency's sake, we recommend using the same
  Python/OS/dependencies when comparing multiple WISP analyses.

1.1
---

* Changed the default `contact_map_distance_limit` to 4.5 so `wisp.py` now
  matches `Wisp.tcl`.
* Updated the `README.md` file.
* Ensured Python2/3 compatibility.
* Modernized/cleaned up the code.

1.0
---

* The original version.
