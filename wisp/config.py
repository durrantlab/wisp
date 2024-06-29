import time
from collections.abc import Iterable

from pydantic import BaseModel, Field


class WispConfig(BaseModel):
    source_residues: Iterable[str] = []
    """This parameter specifies the source residues for path generation.
    A list of residues should be constructed of the form 'CHAIN_RESNAME_RESID',
    separated by spaces. For example: 'X_SER_1 X_LEU_4'.
    """

    sink_residues: Iterable[str] = []
    """This parameter specifies the sink residues for path generation.
    The format is the same as for the source_residues parameter.
    """

    n_paths: int = Field(default=1, ge=1)
    """The desired number of paths."""

    output_dir: str = f'wisp_output__{time.strftime("%b_%d_%Y__%I_%M_%p")}'
    """A new directory where the WISP output should be written. If this parameter
    is not specified, a default output directory is created whose name includes
    the current date for future reference."""

    n_cores: int = 1
    """On unix-like machines, WISP can use multiple processors to significantly
    increase speed. This parameter specifies the number of processors to use."""

    n_paths_max: int = 100000
    """Terminate calculation if the number of estimated paths is greater than
    this."""

    check_max_paths: bool = True
    """Estimate the maximum number of paths and, if larger than `n_paths_max`,
    terminate the calculation.
    """

    frame_chunks: int = 96
    """When WISP is run with multiple processors, the frames from the PDB are
    loaded in chunks before being distributed to the many processors. This
    parameter specifies the number of frames to load before distribution."""

    node_definition: str = "RESIDUE_COM"
    """WISP calculates the covariance matrix by defining nodes associated with
    each protein residue. If node_definition is set to 'CA', the alpha carbon
    will be used. If set to 'RESIDUE_COM,', 'SIDECHAIN_COM,', or 'BACKBONE_COM,'
    the whole-residue, side-chain, or backbone center of mass will be used,
    respectively."""

    contact_map_distance_limit: float = 4.5
    """If you use WISP's default contact-map generator, node pairs with average
    inter-node distances greater than this value will not be considered in
    calculating the covariance matrix. Use a value of 999999.999 to deactivate."""

    seconds_to_wait_before_parallelizing_path_finding: float = 5.0
    """WISP identifies paths from the source to the sink by recursively visiting
    node neighbors. The program begins the recursion algorithm on a single processor
    before distributing the search efforts to multiple processors. This parameter
    specifies how long WISP should search for source-sink paths using a single
    processor before distributing the search effort over multiple processors.
    By waiting longer before distribution, the search efforts are ultimately
    distributed more evenly over the multiple processors, potentially increasing
    speed in the long run. On the other hand, specifiying a lower value for
    this parameter means the program will spend more time running on multiple
    processors, also potentially increasing speed. A balance must be struck."""

    path_saved_matrix: str | None = None
    """If the covariance matrix (appropriately modified by a contact map) has been
    previously saved to a file, set this parameter to 'TRUE' to load the matrix
    instead of generating it from scratch. WISP automatically saves a copy of this
    matrix to the file 'functionalized_matrix_with_contact_map_applied.pickle'
    in the output directory every time it is run."""

    contact_map_path: str | None = None
    """A text file containing a user-specified contact map. If given, each element
    of the functionalized matrix will be multiplied by the corresponding value
    specified in the file. If not given, WISP's default contact map, based on the
    distances between average node locations, will be automatically applied. For
    convenience, WISP automatically saves a human-readable copy of the contact-map
    matrix to the file contact_map_matrix.txt in the output directory every time
    it is run."""

    functionalized_matrix_path: str | None = None
    """A text file containing a user-specified functionalized correlation matrix.
    If not given, WISP's default functionalized correlation matrix, as described
    in the WISP publication, will be automatically calculated. For convenience,
    WISP automatically saves a human-readable copy of the matrix used to the file
    functionalized_correlation_matrix.txt in the output directory every time
    it is run."""

    pdb_single_frame_path: str | None = None
    """By default, WISP uses the trajectory-average structure for positioning the
    nodes, visualizing the paths and protein, etc. However, if desired, a separate
    PDB structure with the same residue order and number can be specified for this
    purpose using the pdb_single_frame_path parameter."""

    write_formatted_paths: bool = False
    """Write a text file containing a simply formatted list of paths.
    """

    shortest_path_radius: float = 0.1
    """WISP outputs a VMD state file to facilitate visualization. The shortest path
    is represented by a strand with the largest radius. Longer paths have
    progressively smaller radii. This parameter specifies the radius of the
    shortest path, in Angstroms."""

    longest_path_radius: float = 0.01
    """This parameter specifies the radius of the longest path visualized,
    in Angstroms."""

    spline_smoothness: float = 0.01
    """The paths are represented by splines connecting the nodes. This parameter
    indicates the smoothness of the splines. Smaller values produce smoother
    splines, but take longer to render."""

    vmd_resolution: int = 6
    """When visualizing in VMD, a number of cylinders and spheres are drawn.
    This parameter specifies the resolution to use."""

    node_sphere_radius: float = 1.0
    """When visualizing in VMD, a number of cylinders and spheres are drawn.
    This parameter specifies the resolution to use."""

    shortest_path_r: float = 0.0
    """The color of the shortest path is given by an RGB color code. This parameter
    specifies the R value, ranging from 0.0 to 1.0."""

    shortest_path_g: float = 0.0
    """The color of the shortest path is given by an RGB color code. This parameter
    specifies the G value, ranging from 0.0 to 1.0."""

    shortest_path_b: float = 0.0
    """The color of the shortest path is given by an RGB color code. This parameter
    specifies the B value, ranging from 0.0 to 1.0."""

    longest_path_r: float = 0.0
    """The color of the longest path is given by an RGB color code. This parameter
    specifies the R value, ranging from 0.0 to 1.0."""

    longest_path_g: float = 0.0
    """The color of the longest path is given by an RGB color code. This parameter
    specifies the G value, ranging from 0.0 to 1.0."""

    longest_path_b: float = 0.0
    """The color of the longest path is given by an RGB color code. This parameter
    specifies the B value, ranging from 0.0 to 1.0."""

    node_sphere_r: float = 0.0
    """The color of the node spheres is given by an RGB color code. This parameter
    specifies the R value, ranging from 0.0 to 1.0."""

    node_sphere_g: float = 0.0
    """The color of the node spheres is given by an RGB color code. This parameter
    specifies the G value, ranging from 0.0 to 1.0."""

    node_sphere_b: float = 0.0
    """The color of the node spheres is given by an RGB color code. This parameter
    specifies the B value, ranging from 0.0 to 1.0."""

    longest_path_opacity: float = 1.0
    """The opacity of the longest path, ranging from 0.0 (transparent) to 1.0
    (fully opaque)."""

    shortest_path_opacity: float = 1.0
    """The opacity of the shortest path, ranging from 0.0 (transparent) to 1.0
    (fully opaque)."""

    node_sphere_opacity: float = 1.0
    """The opacity of the node spheres, ranging from 0.0 (transparent) to
    1.0 (fully opaque)."""
