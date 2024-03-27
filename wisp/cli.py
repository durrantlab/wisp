import argparse
import os
import time

from . import __version__


class UserInput:

    def __init__(self) -> None:
        # Initialize the argument parser
        parser = argparse.ArgumentParser(
            description="WISP: A Tool for Protein Analysis"
        )

        parser.add_argument("pdb_path", type=str, help="Path to PDB file to analyze.")
        parser.add_argument(
            "--source_residues",
            nargs="+",
            help="This parameter specifies the source residues for path generation. A list of residues should be constructed of the form 'CHAIN_RESNAME_RESID', separated by spaces. For example: 'X_SER_1 X_LEU_4'. For unix to treat a space-containing command-line parameter as a single parameter, it must be enclosed in quotes.",
        )
        parser.add_argument(
            "--sink_residues",
            nargs="+",
            help="This parameter specifies the sink residues for path generation. The format is the same as for the source_residues parameter.",
        )
        parser.add_argument(
            "--n_paths", type=int, default=1, help="The desired number of paths."
        )
        parser.add_argument(
            "--output_dir",
            type=str,
            default=f'wisp_output__{time.strftime("%b_%d_%Y__%I_%M_%p")}',
            required=False,
            help="A new directory where the WISP output should be written. If this parameter is not specified, a default output directory is created whose name includes the current date for future reference.",
        )

        # Define arguments
        parser.add_argument(
            "--n_cores",
            type=int,
            default=1,
            help="On unix-like machines, WISP can use multiple processors to significantly increase speed. This parameter specifies the number of processors to use.",
        )
        parser.add_argument("--n_paths_max", type=int, default=100000, help="TODO")
        parser.add_argument(
            "--frame_chunks",
            type=int,
            default=96,
            help="When WISP is run with multiple processors, the frames from the PDB are loaded in chunks before being distributed to the many processors. This parameter specifies the number of frames to load before distribution.",
        )
        parser.add_argument(
            "--node_definition",
            type=str,
            default="RESIDUE_COM",
            choices=["CA", "SIDECHAIN_COM", "BACKBONE_COM", "RESIDUE_COM"],
            help="WISP calculates the covariance matrix by defining nodes associated with each protein residue. If node_definition is set to 'CA', the alpha carbon will be used. If set to 'RESIDUE_COM,', 'SIDECHAIN_COM,', or 'BACKBONE_COM,' the whole-residue, side-chain, or backbone center of mass will be used, respectively.",
        )
        parser.add_argument(
            "--contact_map_distance_limit",
            type=float,
            default=4.5,
            help="If you use WISP's default contact-map generator, node pairs with average inter-node distances greater than this value will not be considered in calculating the covariance matrix. Use a value of 999999.999 to deactivate.",
        )
        # TODO: Need to handle boolean logic when a wisp matrix path is provided
        parser.add_argument(
            "--wisp_saved_matrix_path",
            type=str,
            default=None,
            required=False,
            help="If the covariance matrix (appropriately modifed by a contact map) has been previously saved to a file, set this parameter to 'TRUE' to load the matrix instead of generating it from scratch. WISP automatically saves a copy of this matrix to the file 'functionalized_matrix_with_contact_map_applied.pickle' in the output directory every time it is run.",
        )
        parser.add_argument(
            "--shortest_path_radius",
            type=float,
            default=0.1,
            help="WISP outputs a VMD state file to facilitate visualization. The shortest path is represented by a strand with the largest radius. Longer paths have progressively smaller radii. This parameter specifies the radius of the shortest path, in Angstroms.",
        )
        parser.add_argument(
            "--longest_path_radius",
            type=float,
            default=0.01,
            help="This parameter specifies the radius of the longest path visualized, in Angstroms.",
        )
        parser.add_argument(
            "--spline_smoothness",
            type=float,
            default=0.01,
            help="The paths are represented by splines connecting the nodes. This parameter indicates the smoothness of the splines. Smaller values produce smoother splies, but take longer to render.",
        )
        # TODO: Add VMD flag
        parser.add_argument(
            "--vmd_resolution",
            type=float,
            default=6,
            help="When visualizing in VMD, a number of cylinders and spheres are drawn. This parameter specifies the resolution to use.",
        )
        parser.add_argument(
            "--node_sphere_radius",
            type=float,
            default=1.0,
            help="When visualizing in VMD, spheres are placed at the locations of the nodes. This parameter specifies the radius of these spheres.",
        )
        parser.add_argument(
            "--seconds_to_wait_before_parallelizing_path_finding",
            type=float,
            default=5.0,
            help="WISP identifies paths from the source to the sink by recursively visiting node neighbors. The program begins the recursion algorithm on a single processor before distributing the search efforts to multiple processors. This parameter specifies how long WISP should search for source-sink paths using a single processor before distributing the search effort over multiple processors. By waiting longer before distribution, the search efforts are ultimately distributed more evenly over the multiple processors, potentially increasing speed in the long run. On the other hand, specifiying a lower value for this parameter means the program will spend more time running on multiple processors, also potentially increasing speed. A balance must be struck.",
        )
        parser.add_argument(
            "--contact_map_path",
            type=str,
            default=None,
            required=False,
            help="A text file containing a user-specified contact map. If given, each element of the functionalized matrix will be multiplied by the corresponding value specified in the file. If not given, WISP's default contact map, based on the distances between average node locations, will be automatically applied. For convenience, WISP automatically saves a human-readable copy of the contact-map matrix to the file contact_map_matrix.txt in the output directory every time it is run.",
        )
        parser.add_argument(
            "--functionalized_matrix_path",
            type=str,
            default=None,
            required=False,
            help="A text file containing a user-specified functionalized correlation matrix. If not given, WISP's default functionalized correlation matrix, as described in the WISP publication, will be automatically calculated. For convenience, WISP automatically saves a human-readable copy of the matrix used to the file functionalized_correlation_matrix.txt in the output directory every time it is run.",
        )
        parser.add_argument(
            "--shortest_path_r",
            type=float,
            default=0.0,
            help="The color of the shortest path is given by an RGB color code. This parameter specifies the R value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--shortest_path_g",
            type=float,
            default=0.0,
            help="The color of the shortest path is given by an RGB color code. This parameter specifies the G value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--shortest_path_b",
            type=float,
            default=1.0,
            help="The color of the shortest path is given by an RGB color code. This parameter specifies the B value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--longest_path_r",
            type=float,
            default=1.0,
            help="The color of the longest path is given by an RGB color code. This parameter specifies the R value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--longest_path_g",
            type=float,
            default=0.0,
            help="The color of the longest path is given by an RGB color code. This parameter specifies the G value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--longest_path_b",
            type=float,
            default=0.0,
            help="The color of the longest path is given by an RGB color code. This parameter specifies the B value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--node_sphere_r",
            type=float,
            default=1.0,
            help="The color of the node spheres is given by an RGB color code. This parameter specifies the R value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--node_sphere_g",
            type=float,
            default=1.0,
            help="The color of the node spheres is given by an RGB color code. This parameter specifies the G value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--node_sphere_b",
            type=float,
            default=1.0,
            help="The color of the node spheres is given by an RGB color code. This parameter specifies the B value, ranging from 0.0 to 1.0.",
        )
        parser.add_argument(
            "--longest_path_opacity",
            type=float,
            default=1.0,
            help="The opacity of the longest path, ranging from 0.0 (transparent) to 1.0 (fully opaque).",
        )
        parser.add_argument(
            "--shortest_path_opacity",
            type=float,
            default=1.0,
            help="The opacity of the shortest path, ranging from 0.0 (transparent) to 1.0 (fully opaque). Note that if --shortest_path_opacity, --longest_path_opacity, and --node_sphere_opacity are not all identical, the output TCL file will contain many materials, which may be less-than-desirable for some users.",
        )
        parser.add_argument(
            "--node_sphere_opacity",
            type=float,
            default=1.0,
            help="The opacity of the node spheres, ranging from 0.0 (transparent) to 1.0 (fully opaque).",
        )
        parser.add_argument(
            "--pdb_single_frame_path",
            type=str,
            default=None,
            required=False,
            help="By default, WISP uses the trajectory-average structure for positioning the nodes, visualizing the paths and protein, etc. However, if desired, a separate PDB structure with the same residue order and number can be specified for this purpose using the pdb_single_frame_filename parameter.",
        )

        # Parse the arguments
        args = parser.parse_args()

        # Store arguments in a dictionary or directly access them using args.attribute
        self.parameters = vars(args)  # Convert Namespace to dict

        # Additional setup like checking for help, setting up directories, etc.
        # For example, setting up the output directory
        self.setup_output_directory()

        # TODO: Log parsed parameters or any other necessary information

    def setup_output_directory(self) -> None:
        # Assuming self.parameters is already populated
        output_dir = self.parameters.get(
            "output_directory", f'wisp_output__{time.strftime("%b_%d_%Y__%I_%M_%p")}'
        )
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.parameters["output_directory"] = output_dir


if __name__ == "__main__":
    user_input = UserInput()
