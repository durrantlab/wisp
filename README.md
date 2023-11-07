# WISP

WISP is a trajectory analysis tool that calculates and visualizes allosteric pathways.
It is licensed under the Academic Free License 3.0.
For more information, please see http://opensource.org/licenses/AFL-3.0

WISP is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

Copyright 2012 Adam VanWart and Jacob D. Durrant.
If you have any questions, comments, or suggestions, please contact durrantj [at] pitt [dot] edu.

The latest version of WISP can be downloaded from [http://git.durrantlab.com/jdurrant/wisp](http://git.durrantlab.com/jdurrant/wisp).

If you use WISP in your work, please cite:

> A.T. Van Wart, J.D. Durrant, L. Votapka, R.E. Amaro. Weighted implementation of suboptimal paths (WISP): An optimized algorithm and tool for dynamical network analysis, J. Chem. Theory Comput. 10 (2014) 511-517

## How to Install the VMD WISP Plugin (Linux/Mac)

1. Unzip the file to a directory of your choice using this command: `tar -xzf wisp.tgz <your directory here>`
2. To use the VMD plugin, add these two lines of code to your `.vmdrc` file
   (usually located in your `~/` directory)
   - `set auto_path "$auto_path <your directory here>"  ; #(NOTE: this may require the full pathname)`
   - `vmd_install_extension wisp wisp_tk_cb "Analysis/Wisp"`
3. Now open VMD: Click Extensions > Analysis > WISP

## How to Use WISP from the Command Line

To learn how to use wisp from the command line, see the example in the
`./example_commandline/` directory. Here is a simple example:

```python
python wisp.py -pdb_trajectory_filename multi_frame_pdb.pdb -source_residues "X_SER_1 X_LEU_4" -sink_residues X_ARG_37`
```

## Program Output

The output directory contains multiple output files. The best way to visualize
the output is to use a free program called VMD, which can be downloaded from
[http://www.ks.uiuc.edu/Research/vmd/](http://www.ks.uiuc.edu/Research/vmd/).

The WISP output can be automatically loaded into VMD using the TCL script
named `visualize.tcl`. Assuming `vmd` is the full path to your installed VMD
executable, just run the following from the command line:

`vmd -e visualize.tcl`

If you prefer not to use the command line, simply run the `vmd` executable and
load the `visualize.vmd` file using "File->Load Visualization State..." from
the main menu.

The above methods are very slow. If your output is so large that a faster
option is required, the Tk Console can be used. Use "Extensions->Tk Console"
from the VMD main menu to pull up the Tk Console. Then run the following
command, with the full path to `visualize.tcl` included if necessary:

`source visualize.tcl`

Regardless of the method you use to load in the WISP output, the visualization
will be the same. Individual pathways are shown as tubes (i.e., "wisps"), the
protein is shown in ribbon representation, and protein residues that
participate in any path are shown in licorice representation.

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
  rather than every time WISP is executed. Use `python wisp.py -help` for more
  information.
- `contact_map_matrix.txt`: A human readable representation of the contact
  map. If the user wishes to generate their own contact map rather than
  letting WISP generate one automatically, a custom contact map formatted like
  this one can be loaded into WISP using the
  `-user_specified_contact_map_filename parameter`.
- `functionalized_correlation_matrix.txt`: A human readable representation of
  the functionalized correlation matrix, prior to multiplication by the
  contact map. If the user wishes to generate their own functionalized
  correlation matrix rather than letting WISP generate one automatically, a
  custom matrix formatted like this one can be loaded into WISP using the
  `-user_specified_functionalized_matrix_filename parameter`.
- `simply_formatted_paths.txt`: A simple list of path lengths and nodes. The
  first column contains the lengths, and all following columns contain node
  indices. This file may be helpful for subsequent statistical analyses of the
  WISP output. Note that the `simply_formatted_paths.txt` output file
  reindexes the residues. See the `visualize.tcl` file instead for a more
  human-readable output.

## Parameter Description

`python wisp.py -help` displays the following text:

```text
FILE-SYSTEM PARAMETERS
----------------------
pdb_trajectory_filename: The filename of the multi-frame PDB to
    analyze. Individual frames should be separated by "END" or
    "ENDMDL" lines.
output_directory: A new directory where the WISP output should be
    written. If this parameter is not specified, a default output
    directory is created whose name includes the current date for
    future reference. The default value is
    wisp_output__Sep_12_2019__03_51_AM.

COVARIANCE-MATRIX PARAMETERS
----------------------------
node_definition: WISP calculates the covariance matrix by defining
    nodes associated with each protein residue. If node_definition is
    set to "CA," the alpha carbon will be used. If set to
    "RESIDUE_COM,", "SIDECHAIN_COM,", or "BACKBONE_COM," the whole-
    residue, side-chain, or backbone center of mass will be used,
    respectively. The default value is RESIDUE_COM.
contact_map_distance_limit: If you use WISP's default contact-map
    generator, node pairs with average inter-node distances greater
    than this value will not be considered in calculating the
    covariance matrix. The default value is 4.5.
load_wisp_saved_matrix: If the covariance matrix (appropriately
    modifed by a contact map) has been previously saved to a file, set
    this parameter to "TRUE" to load the matrix instead of generating
    it from scratch. WISP automatically saves a copy of this matrix to
    the file "functionalized_matrix_with_contact_map_applied.pickle"
    in the output directory every time it is run. The default value is
    FALSE.
wisp_saved_matrix_filename: If load_wisp_saved_matrix is set to
    "TRUE," this parameter specifies the file to load. If it is set to
    "FALSE," this parameter specifies the file to which the matrix
    should be saved.

PATH-SEARCHING PARAMETERS
-------------------------
desired_number_of_paths: One of the advantages of WISP is that it can
    calculate not only the optimal path between residues, but multiple
    good paths. This parameter specifies the desired number of paths.
    The default value is 1.
source_residues: This parameter specifies the source residues for path
    generation. A list of residues should be constructed of the form
    "CHAIN_RESNAME_RESID," separated by spaces. For example: "X_SER_1
    X_LEU_4." For unix to treat a space-containing command-line
    parameter as a single parameter, it must be enclosed in quotes.
    If your PDB file does not have a chain, use `A`: `A_LEU_4`.
sink_residues: This parameter specifies the sink residues for path
    generation. The format is the same as for the source_residues
    parameter.

MULTI-PROCESSOR PARAMETERS
--------------------------
number_processors: On unix-like machines, WISP can use multiple
    processors to significantly increase speed. This parameter
    specifies the number of processors to use. The default value is 1.
num_frames_to_load_before_processing: When WISP is run with multiple
    processors, the frames from the PDB are loaded in chunks before
    being distributed to the many processors. This parameter specifies
    the number of frames to load before distribution. The default
    value is 96.

VISUALIZATION PARAMETERS
------------------------
shortest_path_radius: WISP outputs a VMD state file to facilitate
    visualization. The shortest path is represented by a strand with
    the largest radius. Longer paths have progressively smaller radii.
    This parameter specifies the radius of the shortest path, in
    Angstroms. The default value is 0.1.
longest_path_radius: This parameter specifies the radius of the
    longest path visualized, in Angstroms. The default value is 0.01.
spline_smoothness: The paths are represented by splines connecting the
    nodes. This parameter indicates the smoothness of the splines.
    Smaller values produce smoother splies, but take longer to render.
    The default value is 0.01.
vmd_resolution: When visualizing in VMD, a number of cylinders and
    spheres are drawn. This parameter specifies the resolution to use.
    The default value is 6.
node_sphere_radius: When visualizing in VMD, spheres are placed at the
    locations of the nodes. This parameter specifies the radius of
    these spheres. The default value is 1.0.
shortest_path_r: The color of the shortest path is given by an RGB
    color code. This parameter specifies the R value, ranging from 0.0
    to 1.0. The default value is 0.0.
shortest_path_g: The color of the shortest path is given by an RGB
    color code. This parameter specifies the G value, ranging from 0.0
    to 1.0. The default value is 0.0.
shortest_path_b: The color of the shortest path is given by an RGB
    color code. This parameter specifies the B value, ranging from 0.0
    to 1.0. The default value is 1.0.
longest_path_r: The color of the longest path is given by an RGB color
    code. This parameter specifies the R value, ranging from 0.0 to
    1.0. The default value is 1.0.
longest_path_g: The color of the longest path is given by an RGB color
    code. This parameter specifies the G value, ranging from 0.0 to
    1.0. The default value is 0.0.
longest_path_b: The color of the longest path is given by an RGB color
    code. This parameter specifies the B value, ranging from 0.0 to
    1.0. The default value is 0.0.
node_sphere_r: The color of the node spheres is given by an RGB color
    code. This parameter specifies the R value, ranging from 0.0 to
    1.0. The default value is 1.0.
node_sphere_g: The color of the node spheres is given by an RGB color
    code. This parameter specifies the G value, ranging from 0.0 to
    1.0. The default value is 1.0.
node_sphere_b: The color of the node spheres is given by an RGB color
    code. This parameter specifies the B value, ranging from 0.0 to
    1.0. The default value is 1.0.
shortest_path_opacity: The opacity of the shortest path, ranging from
    0.0 (transparent) to 1.0 (fully opaque). Note that if
    --shortest_path_opacity, --longest_path_opacity, and
    --node_sphere_opacity are not all identical, the output TCL file
    will contain many materials, which may be less-than-desirable for
    some users. The default value is 1.0.
longest_path_opacity: The opacity of the longest path, ranging from
    0.0 (transparent) to 1.0 (fully opaque). The default value is 1.0.
node_sphere_opacity: The opacity of the node spheres, ranging from 0.0
    (transparent) to 1.0 (fully opaque). The default value is 1.0.
pdb_single_frame_filename: By default, WISP uses the trajectory-
    average structure for positioning the nodes, visualizing the paths
    and protein, etc. However, if desired, a separate PDB structure
    with the same residue order and number can be specified for this
    purpose using the "pdb_single_frame_filename" parameter.

ADVANCED FEATURES
-----------------
seconds_to_wait_before_parallelizing_path_finding: WISP identifies
    paths from the source to the sink by recursively visiting node
    neighbors. The program begins the recursion algorithm on a single
    processor before distributing the search efforts to multiple
    processors. This parameter specifies how long WISP should search
    for source-sink paths using a single processor before distributing
    the search effort over multiple processors. By waiting longer
    before distribution, the search efforts are ultimately distributed
    more evenly over the multiple processors, potentially increasing
    speed in the long run. On the other hand, specifiying a lower
    value for this parameter means the program will spend more time
    running on multiple processors, also potentially increasing speed.
    A balance must be struck. The default value is 5.0.
user_specified_functionalized_matrix_filename: A text file containing
    a user-specified functionalized correlation matrix. If not given,
    WISP's default functionalized correlation matrix, as described in
    the WISP publication, will be automatically calculated. For
    convenience, WISP automatically saves a human-readable copy of the
    matrix used to the file "functionalized_correlation_matrix.txt" in
    the output directory every time it is run.
user_specified_contact_map_filename: A text file containing a user-
    specified contact map. If given, each element of the
    functionalized matrix will be multiplied by the corresponding
    value specified in the file. If not given, WISP's default contact
    map, based on the distances between average node locations, will
    be automatically applied. For convenience, WISP automatically
    saves a human-readable copy of the contact-map matrix to the file
    "contact_map_matrix.txt" in the output directory every time it is
    run.

Notes:
1) To visualize in VMD, first load the output TCL file, then load the PDB file.
2) WISP ignores PDB segnames. Every residue in your PDB trajectory must be
   uniquely identifiable by the combination of its chain, resname, and resid.

Example:
     python wisp.py -pdb_trajectory_filename multi_frame_pdb.pdb
         -node_definition CA -contact_map_distance_limit 4.5
         -load_wisp_saved_matrix false -wisp_saved_matrix_filename
         matrix.file -desired_number_of_paths 30 -source_residues
         "X_SER_1 X_LEU_4" -sink_residues X_ARG_37 -number_processors
         24 -num_frames_to_load_before_processing 96
         -seconds_to_wait_before_parallelizing_path_finding 10.0
         -shortest_path_radius 0.2 -longest_path_radius 0.05
         -spline_smoothness 0.05 -vmd_resolution 6 -node_sphere_radius
         1.0
```
