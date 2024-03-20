import os
import sys
import textwrap
import time

from loguru import logger

from . import __version__


class UserInput:
    """Process and store user-specified command-line parameters"""

    def __init__(self):
        """Receives, processes, and stores command-line parameters"""

        # Display program information.
        logger.info("WISP v{}", __version__)
        logger.info(
            "If you use WISP in your work, please cite: https://doi.org/10.1021/ct4008603"
        )

        # get the user input
        self.parameters = {}

        # set defaults
        self.parameters["number_processors"] = 1

        # only relevant if number_processors > 1
        self.parameters["num_frames_to_load_before_processing"] = 96

        # can be CA, SIDECHAIN_COM, BACKBONE_COM, or RESIDUE_COM
        self.parameters["node_definition"] = "RESIDUE_COM"
        self.parameters["pdb_trajectory_filename"] = ""  # The trajectory to be analyzed

        # If specified, correlations between residues whose average node
        # locations are greater than this value will be ignored
        self.parameters["contact_map_distance_limit"] = 4.5

        self.parameters["desired_number_of_paths"] = 1  # how many paths to consider
        self.parameters["source_residues"] = []
        self.parameters["sink_residues"] = []
        self.parameters["load_wisp_saved_matrix"] = "FALSE"
        self.parameters["wisp_saved_matrix_filename"] = ""
        self.parameters["shortest_path_radius"] = 0.1
        self.parameters["longest_path_radius"] = 0.01
        self.parameters["spline_smoothness"] = 0.01
        self.parameters["vmd_resolution"] = 6
        self.parameters["node_sphere_radius"] = 1.0
        self.parameters["seconds_to_wait_before_parallelizing_path_finding"] = 5.0
        self.parameters["user_specified_contact_map_filename"] = ""
        self.parameters["user_specified_functionalized_matrix_filename"] = ""
        self.parameters["shortest_path_r"] = 0.0
        self.parameters["shortest_path_g"] = 0.0
        self.parameters["shortest_path_b"] = 1.0
        self.parameters["longest_path_r"] = 1.0
        self.parameters["longest_path_g"] = 0.0
        self.parameters["longest_path_b"] = 0.0
        self.parameters["node_sphere_r"] = 1.0
        self.parameters["node_sphere_g"] = 1.0
        self.parameters["node_sphere_b"] = 1.0
        self.parameters["longest_path_opacity"] = 1.0
        self.parameters["shortest_path_opacity"] = 1.0
        self.parameters["node_sphere_opacity"] = 1.0
        self.parameters[
            "output_directory"
        ] = f'wisp_output__{time.strftime("%b_%d_%Y__%I_%M_%p")}'
        self.parameters["pdb_single_frame_filename"] = ""

        # first, check if the help file has been requested
        for t in sys.argv:
            if t.replace("-", "").lower() == "help":
                self.get_help()

        # load the parameters
        parameters_that_are_floats = [
            "node_sphere_opacity",
            "node_sphere_r",
            "node_sphere_g",
            "node_sphere_b",
            "longest_path_opacity",
            "shortest_path_opacity",
            "shortest_path_r",
            "shortest_path_g",
            "shortest_path_b",
            "longest_path_r",
            "longest_path_g",
            "longest_path_b",
            "contact_map_distance_limit",
            "shortest_path_radius",
            "longest_path_radius",
            "spline_smoothness",
            "seconds_to_wait_before_parallelizing_path_finding",
            "node_sphere_radius",
        ]
        parameters_that_are_ints = [
            "number_processors",
            "num_frames_to_load_before_processing",
            "desired_number_of_paths",
            "vmd_resolution",
        ]
        parameters_that_are_strings = [
            "pdb_single_frame_filename",
            "output_directory",
            "user_specified_contact_map_filename",
            "user_specified_functionalized_matrix_filename",
            "node_definition",
            "pdb_trajectory_filename",
            "load_wisp_saved_matrix",
            "wisp_saved_matrix_filename",
            "simply_formatted_paths_filename",
        ]
        parameters_that_are_lists = ["source_residues", "sink_residues"]

        for t in range(len(sys.argv)):
            key = sys.argv[t].lower().replace("-", "")
            if key in parameters_that_are_floats:
                self.parameters[key] = float(sys.argv[t + 1])
                sys.argv[t] = ""
                sys.argv[t + 1] = ""
            if key in parameters_that_are_ints:
                self.parameters[key] = int(sys.argv[t + 1])
                sys.argv[t] = ""
                sys.argv[t + 1] = ""
            if key in parameters_that_are_strings:
                self.parameters[key] = sys.argv[t + 1]
                sys.argv[t] = ""
                sys.argv[t + 1] = ""
            if (
                key in parameters_that_are_lists
            ):  # format: "CHAIN_RESNAME_RESID CHAIN_RESNAME_RESID".
                residues = sys.argv[t + 1].strip()
                residues = residues.replace("\t", " ")
                while "  " in residues:
                    residues = residues.replace("  ", " ")
                residues = residues.split(" ")
                self.parameters[key] = residues
                sys.argv[t] = ""
                sys.argv[t + 1] = ""

        # some parameters need to always be caps
        tocap = ["node_definition", "load_wisp_saved_matrix"]
        for param in tocap:
            self.parameters[param] = self.parameters[param].upper()

        # The paths from A to B are the same as the paths from B to A. If the
        # same residue is in both source_residues and sink_residues there will
        # be redundancies. Furthermore, the path from A to A will produce an
        # error. So these two lists need to be made mutually exclusive.
        # parameters_that_are_lists = ['source_residues', 'sink_residues']

        # what if not all required parameters have been specified?
        if (
            self.parameters["pdb_trajectory_filename"] == ""
            or self.parameters["source_residues"] == []
            or self.parameters["sink_residues"] == []
        ):
            logger.critical(
                "You have failed to provide all the required parameters. In its simplest form, WISP can be used like this:"
            )
            sys.exit(0)

        # make the output directory
        if self.parameters["output_directory"][-1:] != os.sep:
            self.parameters["output_directory"] = (
                self.parameters["output_directory"] + os.sep
            )
        if os.path.exists(self.parameters["output_directory"]):
            logger.critical(
                "The output directory, "
                + self.parameters["output_directory"]
                + ", already exists. Please delete this directory or select a different one for output before proceeding."
            )
            sys.exit()
        else:
            os.mkdir(self.parameters["output_directory"])

        # some parameters are auto generated
        autogenerated_parameters = ["logfile", "simply_formatted_paths_filename"]
        self.parameters["logfile"] = open(
            self.parameters["output_directory"] + "log.txt", "w"
        )
        self.parameters["simply_formatted_paths_filename"] = (
            self.parameters["output_directory"] + "simply_formatted_paths.txt"
        )

        # inform what parameters will be used
        with open(
            self.parameters["output_directory"] + "parameters_used.txt", "w"
        ) as parameters_file:
            logger.info("Wisp {}", __version__)

            logger.info(
                "# Command-line Parameters:",
                [self.parameters["logfile"], parameters_file],
            )
            somekeys = self.parameters.keys()
            somekeys = sorted(somekeys)
            for key in somekeys:
                if not key in autogenerated_parameters:
                    logger.info(
                        "#\t" + key + ": " + str(self.parameters[key]),
                        [self.parameters["logfile"], parameters_file],
                    )

            logger.info(
                "# A command like the following should regenerate this output:",
                [self.parameters["logfile"], parameters_file],
            )
            prog = "# " + sys.executable + " " + os.path.basename(sys.argv[0]) + " "
            for key in somekeys:
                if not key in autogenerated_parameters and self.parameters[key] != "":
                    if not key in ["sink_residues", "source_residues"]:
                        prog = prog + "-" + key + " " + str(self.parameters[key]) + " "
                    else:
                        prog = (
                            prog
                            + "-"
                            + key
                            + ' "'
                            + " ".join(self.parameters[key])
                            + '" '
                        )

            prog = prog.strip()
            logger.info(prog, [self.parameters["logfile"], parameters_file])

    def __getitem__(self, key):
        return self.parameters[key.lower()]

    def get_help(self):
        """Returns a help file describing program usage"""

        description = []

        # File-system related
        description.append(("title", "FILE-SYSTEM PARAMETERS"))
        description.append(
            (
                "pdb_trajectory_filename",
                'The filename of the multi-frame PDB to analyze. Individual frames should be separated by "END" or "ENDMDL" lines.',
            )
        )
        description.append(
            (
                "output_directory",
                "A new directory where the WISP output should be written. If this parameter is not specified, a default output directory is created whose name includes the current date for future reference.",
            )
        )

        # Covariance-matrix related
        description.append(("title", "COVARIANCE-MATRIX PARAMETERS"))
        description.append(
            (
                "node_definition",
                'WISP calculates the covariance matrix by defining nodes associated with each protein residue. If node_definition is set to "CA," the alpha carbon will be used. If set to "RESIDUE_COM,", "SIDECHAIN_COM,", or "BACKBONE_COM," the whole-residue, side-chain, or backbone center of mass will be used, respectively.',
            )
        )
        description.append(
            (
                "contact_map_distance_limit",
                "If you use WISP's default contact-map generator, node pairs with average inter-node distances greater than this value will not be considered in calculating the covariance matrix. Use a value of 999999.999 to deactivate.",
            )
        )
        description.append(
            (
                "load_wisp_saved_matrix",
                'If the covariance matrix (appropriately modifed by a contact map) has been previously saved to a file, set this parameter to "TRUE" to load the matrix instead of generating it from scratch. WISP automatically saves a copy of this matrix to the file "functionalized_matrix_with_contact_map_applied.pickle" in the output directory every time it is run.',
            )
        )
        description.append(
            (
                "wisp_saved_matrix_filename",
                'If load_wisp_saved_matrix is set to "TRUE," this parameter specifies the file to load. If it is set to "FALSE," this parameter specifies the file to which the matrix should be saved.',
            )
        )

        # Path related
        description.append(("title", "PATH-SEARCHING PARAMETERS"))
        description.append(
            (
                "desired_number_of_paths",
                "One of the advantages of WISP is that it can calculate not only the optimal path between residues, but multiple good paths. This parameter specifies the desired number of paths.",
            )
        )
        description.append(
            (
                "source_residues",
                'This parameter specifies the source residues for path generation. A list of residues should be constructed of the form "CHAIN_RESNAME_RESID," separated by spaces. For example: "X_SER_1 X_LEU_4." For unix to treat a space-containing command-line parameter as a single parameter, it must be enclosed in quotes.',
            )
        )
        description.append(
            (
                "sink_residues",
                "This parameter specifies the sink residues for path generation. The format is the same as for the source_residues parameter.",
            )
        )

        # Multi-processor related
        description.append(("title", "MULTI-PROCESSOR PARAMETERS"))
        description.append(
            (
                "number_processors",
                "On unix-like machines, WISP can use multiple processors to significantly increase speed. This parameter specifies the number of processors to use.",
            )
        )
        description.append(
            (
                "num_frames_to_load_before_processing",
                "When WISP is run with multiple processors, the frames from the PDB are loaded in chunks before being distributed to the many processors. This parameter specifies the number of frames to load before distribution.",
            )
        )

        # Visualization
        description.append(("title", "VISUALIZATION PARAMETERS"))
        description.append(
            (
                "shortest_path_radius",
                "WISP outputs a VMD state file to facilitate visualization. The shortest path is represented by a strand with the largest radius. Longer paths have progressively smaller radii. This parameter specifies the radius of the shortest path, in Angstroms.",
            )
        )
        description.append(
            (
                "longest_path_radius",
                "This parameter specifies the radius of the longest path visualized, in Angstroms.",
            )
        )
        description.append(
            (
                "spline_smoothness",
                "The paths are represented by splines connecting the nodes. This parameter indicates the smoothness of the splines. Smaller values produce smoother splies, but take longer to render.",
            )
        )
        description.append(
            (
                "vmd_resolution",
                "When visualizing in VMD, a number of cylinders and spheres are drawn. This parameter specifies the resolution to use.",
            )
        )
        description.append(
            (
                "node_sphere_radius",
                "When visualizing in VMD, spheres are placed at the locations of the nodes. This parameter specifies the radius of these spheres.",
            )
        )
        description.append(
            (
                "shortest_path_r",
                "The color of the shortest path is given by an RGB color code. This parameter specifies the R value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "shortest_path_g",
                "The color of the shortest path is given by an RGB color code. This parameter specifies the G value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "shortest_path_b",
                "The color of the shortest path is given by an RGB color code. This parameter specifies the B value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "longest_path_r",
                "The color of the longest path is given by an RGB color code. This parameter specifies the R value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "longest_path_g",
                "The color of the longest path is given by an RGB color code. This parameter specifies the G value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "longest_path_b",
                "The color of the longest path is given by an RGB color code. This parameter specifies the B value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "node_sphere_r",
                "The color of the node spheres is given by an RGB color code. This parameter specifies the R value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "node_sphere_g",
                "The color of the node spheres is given by an RGB color code. This parameter specifies the G value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "node_sphere_b",
                "The color of the node spheres is given by an RGB color code. This parameter specifies the B value, ranging from 0.0 to 1.0.",
            )
        )
        description.append(
            (
                "shortest_path_opacity",
                "The opacity of the shortest path, ranging from 0.0 (transparent) to 1.0 (fully opaque). Note that if --shortest_path_opacity, --longest_path_opacity, and --node_sphere_opacity are not all identical, the output TCL file will contain many materials, which may be less-than-desirable for some users.",
            )
        )
        description.append(
            (
                "longest_path_opacity",
                "The opacity of the longest path, ranging from 0.0 (transparent) to 1.0 (fully opaque).",
            )
        )
        description.append(
            (
                "node_sphere_opacity",
                "The opacity of the node spheres, ranging from 0.0 (transparent) to 1.0 (fully opaque).",
            )
        )
        description.append(
            (
                "pdb_single_frame_filename",
                'By default, WISP uses the trajectory-average structure for positioning the nodes, visualizing the paths and protein, etc. However, if desired, a separate PDB structure with the same residue order and number can be specified for this purpose using the "pdb_single_frame_filename" parameter.',
            )
        )

        # Advanced features
        description.append(("title", "ADVANCED FEATURES"))
        description.append(
            (
                "seconds_to_wait_before_parallelizing_path_finding",
                "WISP identifies paths from the source to the sink by recursively visiting node neighbors. The program begins the recursion algorithm on a single processor before distributing the search efforts to multiple processors. This parameter specifies how long WISP should search for source-sink paths using a single processor before distributing the search effort over multiple processors. By waiting longer before distribution, the search efforts are ultimately distributed more evenly over the multiple processors, potentially increasing speed in the long run. On the other hand, specifiying a lower value for this parameter means the program will spend more time running on multiple processors, also potentially increasing speed. A balance must be struck.",
            )
        )
        description.append(
            (
                "user_specified_functionalized_matrix_filename",
                'A text file containing a user-specified functionalized correlation matrix. If not given, WISP\'s default functionalized correlation matrix, as described in the WISP publication, will be automatically calculated. For convenience, WISP automatically saves a human-readable copy of the matrix used to the file "functionalized_correlation_matrix.txt" in the output directory every time it is run.',
            )
        )
        description.append(
            (
                "user_specified_contact_map_filename",
                'A text file containing a user-specified contact map. If given, each element of the functionalized matrix will be multiplied by the corresponding value specified in the file. If not given, WISP\'s default contact map, based on the distances between average node locations, will be automatically applied. For convenience, WISP automatically saves a human-readable copy of the contact-map matrix to the file "contact_map_matrix.txt" in the output directory every time it is run.',
            )
        )

        for item in description:
            if item[0] == "title":
                logger.info(item[1])
                logger.info("-" * len(item[1]))
            else:
                towrap = f"{item[0]}: {item[1]}"
                towrap = (
                    f"{towrap}"
                    if self.parameters[item[0]] in [[], ""]
                    else f"{towrap} The default value is {str(self.parameters[item[0]])}."
                )
                wrapper = textwrap.TextWrapper(
                    initial_indent="", subsequent_indent="    "
                )
                logger.info(wrapper.fill(towrap))
        logger.info("")
        logger.info("Notes:")
        logger.info(
            "1) To visualize in VMD, first load the output TCL file, then load the PDB file."
        )
        logger.info(
            "2) WISP ignores PDB segnames. Every residue in your PDB trajectory must be uniquely identifiable by the combination of its chain, resname, and resid."
        )
        logger.info("")
        logger.info("Example:")
        wrapper = textwrap.TextWrapper(
            initial_indent="     ", subsequent_indent="         "
        )
        logger.info(
            wrapper.fill(
                'python wisp.py -pdb_trajectory_filename multi_frame_pdb.pdb -node_definition CA -contact_map_distance_limit 4.5 -load_wisp_saved_matrix false -wisp_saved_matrix_filename matrix.file -desired_number_of_paths 30 -source_residues "X_SER_1 X_LEU_4" -sink_residues X_ARG_37 -number_processors 24 -num_frames_to_load_before_processing 96 -seconds_to_wait_before_parallelizing_path_finding 10.0 -shortest_path_radius 0.2 -longest_path_radius 0.05 -spline_smoothness 0.05 -vmd_resolution 6 -node_sphere_radius 1.0'
            )
        )
        logger.info("")
        sys.exit(0)
