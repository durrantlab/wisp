"""WISP is licensed under the Academic Free License 3.0. For more
information, please see http://opensource.org/licenses/AFL-3.0

WISP is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Copyright 2012 Adam VanWart and Jacob D. Durrant. If you have any questions,
comments, or suggestions, please don't hesitate to contact durrantj [at] pitt
[dot] edu.

The latest version of WISP can be downloaded from
http://git.durrantlab.com/jdurrant/wisp

If you use WISP in your work, please cite A.T. Van Wart, J.D. Durrant, L.
Votapka, R.E. Amaro. Weighted implementation of suboptimal paths (WISP): An
optimized algorithm and tool for dynamical network analysis, J. Chem. Theory
Comput. 10 (2014) 511-517."""

import copy
import gc
import multiprocessing
import os
import sys
import time
import textwrap
import networkx
import numpy
from scipy import interpolate
from scipy.spatial.distance import cdist

try:
    # Python 2
    import cPickle as pickle
except ImportError:
    # Python 3
    import pickle

opened_log_files = {}

######################## For logging ########################
def log(astring, fileobjects):  # prints to screen and to log file
    """Outputs WISP messages

    Arguments:
    astring -- a string containing the message
    fileobjects -- a list of python file objects (or filenames) specifying 
        where the messages should be saved
    """

    global opened_log_files

    if not isinstance(fileobjects, list):
        # it's not a list, so make it one
        fileobjects = [fileobjects]

    # If it's a string, open it as a file with append. Below is because file
    # objects cannot be pickled for distribution to multiple processors, but
    # file names (strings) can be.
    for t in range(len(fileobjects)):
        if isinstance(fileobjects[t], str):
            if fileobjects[t] in opened_log_files:
                fileobjects[t] = opened_log_files[fileobjects[t]]
            else:
                fileobjects[t] = open(fileobjects[t], "w")
                opened_log_files[fileobjects[t]] = fileobjects[t]

    print(astring)

    for fileobject in fileobjects:
        fileobject.write(astring + "\n")



######################## Data Structures to Store PDB Info ########################


class Atom:
    """A class containing atomic information."""

    def read_pdb_line(self, Line):
        """Reads atomic information from a string formatted according to the PDB standard.

        Arguments:
        Line -- A string formatted according to the PDB standard.
        """

        self.line = Line
        self.atomname = Line[11:16].strip()

        if len(self.atomname) == 1:
            self.atomname = self.atomname + "  "
        elif len(self.atomname) == 2:
            self.atomname = self.atomname + " "
        elif len(self.atomname) == 3:
            # This line is necessary for babel to work, though many PDBs in
            # the PDB would have this line commented out
            self.atomname = self.atomname + " "

        # now get the chain
        self.chain = Line[21:22]

        # now get the resid
        try:
            self.resid = int(Line[22:26])
        except Exception:
            self.resid = 0

        self.coordinates_numpy = numpy.array(
            [float(Line[30:38]), float(Line[38:46]), float(Line[46:54])], numpy.float64
        )

        self.element = ""
        if len(Line) >= 79:
            # element specified explicitly at end of life
            self.element = Line[76:79].strip().upper()
        if self.element == "":  # try to guess at element from name
            two_letters = self.atomname[:2].strip().upper()
            if two_letters == "BR":
                self.element = "BR"
            elif two_letters == "CL":
                self.element = "CL"
            elif two_letters == "BI":
                self.element = "BI"
            elif two_letters == "AS":
                self.element = "AS"
            elif two_letters == "AG":
                self.element = "AG"
            elif two_letters == "LI":
                self.element = "LI"
            elif two_letters == "MG":
                self.element = "MG"
            elif two_letters == "RH":
                self.element = "RH"
            elif two_letters == "ZN":
                self.element = "ZN"
            else:
                # So, just assume it's the first letter.
                self.element = self.atomname[:1].strip().upper()

        # Any number needs to be removed from the element name
        self.element = self.element.replace("0", "")
        self.element = self.element.replace("1", "")
        self.element = self.element.replace("2", "")
        self.element = self.element.replace("3", "")
        self.element = self.element.replace("4", "")
        self.element = self.element.replace("5", "")
        self.element = self.element.replace("6", "")
        self.element = self.element.replace("7", "")
        self.element = self.element.replace("8", "")
        self.element = self.element.replace("9", "")

        self.molecule_index = Line[6:12].strip()
        self.resname = Line[16:20]
        if self.resname.strip() == "":
            self.resname = " MOL"


class Molecule:
    """Loads, saves, and manupulates molecuar models."""

    def load_pdb_from_list(self, alist):
        """Loads a list of PDB ATOM/HETATM lines into the current Molecule object

        Arguments:
        alist -- the list of PDB lines
        """

        gc.disable()

        # have to use python lists initially because not sure of size
        self.atomnames = []
        self.chains = []
        self.resids = []
        self.elements = []
        self.resnames = []
        self.coordinates = []

        for line in alist:
            if len(line) >= 7 and (line[:4] == "ATOM" or line[:6] == "HETATM"):
                temp_atom = Atom()
                temp_atom.read_pdb_line(line)
                self.atomnames.append(temp_atom.atomname)
                self.chains.append(temp_atom.chain)
                self.resids.append(temp_atom.resid)
                self.elements.append(temp_atom.element)
                self.resnames.append(temp_atom.resname)
                self.coordinates.append(temp_atom.coordinates_numpy)

        # convert them into numpy arrays
        self.atomnames = numpy.array(self.atomnames)
        self.chains = numpy.array(self.chains)
        self.resids = numpy.array(self.resids)
        self.elements = numpy.array(self.elements)
        self.resnames = numpy.array(self.resnames)
        self.coordinates = numpy.array(self.coordinates, numpy.float64)

        gc.enable()

    def save_pdb(self, filename):
        """Saves a pdb file

        Arguments:
        filename -- a string specifying the file name
        """

        with open(filename, "w") as f:
            for index in range(len(self.atomnames)):
                line = (
                    "ATOM  "
                    + str(index + 1).rjust(5)
                    + self.atomnames[index].rjust(5)
                    + self.resnames[index].strip().rjust(4)
                    + self.chains[index].strip().rjust(2)
                    + str(self.resids[index]).rjust(4)
                    + "    "
                    + ("%.3f" % self.coordinates[index][0]).rjust(8)
                    + ("%.3f" % self.coordinates[index][1]).rjust(8)
                    + ("%.3f" % self.coordinates[index][2]).rjust(8)
                    + " " * 24
                )
                f.write(line + "\n")

    def map_atoms_to_residues(self):
        """Sets up self.residue_identifier_to_atom_indices, which matches
        chain_resname_resid to associated atom indices"""

        # each residue is uniquely identified by its chain, resname, resid
        # triplet
        residue_identifiers_for_all_atoms = [
            f"{self.chains[index].strip()}_{self.resnames[index].strip()}_{str(self.resids[index])}"
            for index in range(len(self.coordinates))
        ]
        self.residue_identifiers_in_order = residue_identifiers_for_all_atoms[:]
        for t in range(len(self.residue_identifiers_in_order) - 1, 0, -1):
            if (
                self.residue_identifiers_in_order[t]
                == self.residue_identifiers_in_order[t - 1]
            ):
                self.residue_identifiers_in_order.pop(t)

        residue_identifiers_for_all_atoms = numpy.array(
            residue_identifiers_for_all_atoms
        )
        self.residue_identifiers_in_order = numpy.array(
            self.residue_identifiers_in_order
        )

        self.residue_identifier_to_atom_indices = {}
        for the_id in self.residue_identifiers_in_order:
            self.residue_identifier_to_atom_indices[the_id] = numpy.nonzero(
                residue_identifiers_for_all_atoms == the_id
            )[0]

    def get_indices_of_atoms_in_a_residue_by_atom_name(
        self, residue_identifier, atom_names_list, not_selection=False
    ):
        """Gets the indices of atoms in a specified residue

        Arguments:
        residue_identifier -- a string (chain_resname_resid) specifying the residue
        atom_names_list -- a list of strings containing the names of the atoms to keep
        not_selection -- a optional boolean. if False, match the atom_names_list items.
                         if True, match the items not in atom_names_list

        Returns a numpy array containing the indices of the atoms to keep
        """

        # first, get the indices of the residue
        residue_indices = self.residue_identifier_to_atom_indices[residue_identifier]

        # now find which of these correspond to the desired atom names
        indices_to_keep = []
        for indx in residue_indices:
            if not not_selection:
                if self.atomnames[indx].strip() in atom_names_list:
                    indices_to_keep.append(indx)
            elif self.atomnames[indx].strip() not in atom_names_list:
                indices_to_keep.append(indx)
        indices_to_keep = numpy.array(indices_to_keep, numpy.float64)

        if not indices_to_keep:
            raise Exception(
                f"No atoms found in residue {residue_identifier} with atom names {str(atom_names_list)}"
            )

        return indices_to_keep

    def get_center_of_mass_from_selection_by_atom_indices(self, indices_selection):
        """Gets the center of mass of a set of atoms

        Arguments:
        indices_selection -- a numpy array containing the indices of the
                             atoms to consider

        Returns a numpy array containing the 3D coordinates of the center of
        mass of those atoms
        """

        coors = self.coordinates[indices_selection]
        eles = self.elements[indices_selection]
        masses = numpy.empty(coors.shape)

        for t in range(masses.shape[0]):
            themass = self.get_mass(eles[t])
            masses[t][0] = themass
            masses[t][1] = themass
            masses[t][2] = themass

        center_of_mass = ((coors * masses).sum(axis=0)) * (1.0 / masses.sum(axis=0))
        return numpy.array(
            [center_of_mass[0], center_of_mass[1], center_of_mass[2]], numpy.float64
        )

    def get_mass(self, element_name):
        """A library to provide the mass of a given element

        Arguments:
        element_name -- a string that specifies the element

        Returns a float, the mass of the specified element. If the element is
        not in the library, returns None.
        """

        element_name = element_name.upper()

        mass = {
            "H": 1.00794,
            "C": 12.0107,
            "CL": 35.453,
            "N": 14.0067,
            "O": 15.9994,
            "P": 30.973762,
            "S": 32.065,
            "BR": 79.904,
            "I": 126.90447,
            "F": 18.9984032,
            "B": 24.3051,
            "HG": 200.59,
            "BI": 208.9804,
            "AS": 74.9216,
            "AG": 107.8682,
            "K": 39.0983,
            "LI": 6.941,
            "MG": 24.305,
            "RH": 102.9055,
            "ZN": 65.38,
        }

        try:
            return mass[element_name]
        except Exception:
            return None

    def map_nodes_to_residues(self, node_definition):
        """For each residue in the molecule, define the node

        Arguments:
        node_definition -- a string describing the definition of the node: CA,
                           RESIDUE_COM, BACKBONE_COM, or SIDECHAIN_COM
        """

        self.nodes = numpy.empty((len(self.residue_identifiers_in_order), 3))
        for index, residue_iden in enumerate(self.residue_identifiers_in_order):
            if node_definition == "CA":  # the node is at the alpha carbon
                indices_to_consider = self.get_indices_of_atoms_in_a_residue_by_atom_name(
                    residue_iden, ["CA"]
                )
                node_loc = self.coordinates[int(indices_to_consider[0])]
            elif (
                node_definition == "RESIDUE_COM"
            ):  # the node is the residue center of mass
                node_loc = self.get_center_of_mass_from_selection_by_atom_indices(
                    self.residue_identifier_to_atom_indices[residue_iden]
                )
            elif (
                node_definition == "BACKBONE_COM"
            ):  # the node is the residue center of mass
                indices_to_consider = self.get_indices_of_atoms_in_a_residue_by_atom_name(
                    residue_iden,
                    [
                        "C",
                        "CA",
                        "H",
                        "H1",
                        "H2",
                        "H3",
                        "HA",
                        "HA2",
                        "HH1",
                        "HN",
                        "HT1",
                        "HT2",
                        "HT3",
                        "HW",
                        "N",
                        "O",
                        "O1",
                        "O2",
                        "OT1",
                        "OT2",
                        "OXT",
                    ],
                )
                node_loc = self.get_center_of_mass_from_selection_by_atom_indices(
                    indices_to_consider
                )
            elif (
                node_definition == "SIDECHAIN_COM"
            ):  # the node is the residue center of mass
                indices_to_consider = self.get_indices_of_atoms_in_a_residue_by_atom_name(
                    residue_iden,
                    [
                        "C",
                        "CA",
                        "H",
                        "H1",
                        "H2",
                        "H3",
                        "HA",
                        "HA2",
                        "HH1",
                        "HN",
                        "HT1",
                        "HT2",
                        "HT3",
                        "HW",
                        "N",
                        "O",
                        "O1",
                        "O2",
                        "OT1",
                        "OT2",
                        "OXT",
                    ],
                    True,
                )
                node_loc = self.get_center_of_mass_from_selection_by_atom_indices(
                    indices_to_consider
                )
            self.nodes[index][0] = node_loc[0]
            self.nodes[index][1] = node_loc[1]
            self.nodes[index][2] = node_loc[2]


######################## Get User Input ########################


class UserInput:
    """Process and store user-specified command-line parameters"""

    def __init__(self):
        """Receives, processes, and stores command-line parameters"""

        # Display program information.
        print("WISP 1.4\n")
        print(
            "The latest version of WISP can be downloaded from\nhttp://git.durrantlab.com/jdurrant/wisp\n"
        )

        print(
            "If you use WISP in your work, please cite:\nJ. Chem. Theory Comput. 10 (2014) 511-517.\n"
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
        self.parameters["output_directory"] = f'wisp_output__{time.strftime("%b_%d_%Y__%I_%M_%p")}'
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
            print(
                "\nYou have failed to provide all the required parameters. In its simplest form, WISP can be used like this:"
            )
            print(
                '     python wisp.py -pdb_trajectory_filename multi_frame_pdb.pdb -source_residues "X_SER_1 X_LEU_4" -sink_residues X_ARG_37'
            )
            print("")
            print(
                "For more detailed help, use the -help command-line parameter: python wisp.py -help\n"
            )
            sys.exit(0)

        # make the output directory
        if self.parameters["output_directory"][-1:] != os.sep:
            self.parameters["output_directory"] = (
                self.parameters["output_directory"] + os.sep
            )
        if os.path.exists(self.parameters["output_directory"]):
            print(
                "The output directory, "
                + self.parameters["output_directory"]
                + ", already exists. Please delete this directory or select a different one for output before proceeding."
            )
            sys.exit()
        else:
            os.mkdir(self.parameters["output_directory"])

        # some parameters are auto generated
        autogenerated_parameters = ["logfile", "simply_formatted_paths_filename"]

        self.parameters["logfile"] = self.parameters["output_directory"] + "log.txt"
        # open(
        #     self.parameters["output_directory"] + "log.txt", "w"
        # )
        self.parameters["simply_formatted_paths_filename"] = (
            self.parameters["output_directory"] + "simply_formatted_paths.txt"
        )

        # inform what parameters will be used
        with open(
            self.parameters["output_directory"] + "parameters_used.txt", "w"
        ) as parameters_file:
            log(
                "# Wisp 1.4\n# ========\n",
                [self.parameters["logfile"], parameters_file],
            )

            log(
                "# Command-line Parameters:",
                [self.parameters["logfile"], parameters_file],
            )
            somekeys = self.parameters.keys()
            somekeys = sorted(somekeys)
            for key in somekeys:
                if not key in autogenerated_parameters:
                    log(
                        "#\t" + key + ": " + str(self.parameters[key]),
                        [self.parameters["logfile"], parameters_file],
                    )

            log(
                "\n# A command like the following should regenerate this output:",
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
            log(prog, [self.parameters["logfile"], parameters_file])

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
                print("\n" + item[1])
                print("-" * len(item[1]))
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
                print(wrapper.fill(towrap))
        print("")
        print("Notes:")
        print(
            "1) To visualize in VMD, first load the output TCL file, then load the PDB file."
        )
        print(
            "2) WISP ignores PDB segnames. Every residue in your PDB trajectory must be uniquely identifiable by the combination of its chain, resname, and resid."
        )
        print("")
        print("Example:")
        wrapper = textwrap.TextWrapper(
            initial_indent="     ", subsequent_indent="         "
        )
        print(
            wrapper.fill(
                'python wisp.py -pdb_trajectory_filename multi_frame_pdb.pdb -node_definition CA -contact_map_distance_limit 4.5 -load_wisp_saved_matrix false -wisp_saved_matrix_filename matrix.file -desired_number_of_paths 30 -source_residues "X_SER_1 X_LEU_4" -sink_residues X_ARG_37 -number_processors 24 -num_frames_to_load_before_processing 96 -seconds_to_wait_before_parallelizing_path_finding 10.0 -shortest_path_radius 0.2 -longest_path_radius 0.05 -spline_smoothness 0.05 -vmd_resolution 6 -node_sphere_radius 1.0'
            )
        )
        print("")
        sys.exit(0)


######################## Objects to Handle Usign Multiple Processors ########################

# --- first, to generate a covariance matrix ---#
class multi_threading_to_collect_data_from_frames:
    """Launch PDB-frame processing on multiple processors"""

    combined_results = None

    def __init__(self, inputs, num_processors):
        """Launches PDB-frame data processing on multiple processors

        Arguments:
        inputs -- the data to be processed, in a list
        num_processors -- the number of processors to use to process this data, an integer

        """

        self.results = []

        # first, if num_processors <= 0, determine the number of processors to
        # use programatically
        if num_processors <= 0:
            num_processors = multiprocessing.cpu_count()

        # reduce the number of processors if too many have been specified
        if len(inputs) < num_processors:
            num_processors = len(inputs)

        # now, divide the inputs into the appropriate number of processors
        inputs_divided = {t: [] for t in range(num_processors)}

        for t in range(0, len(inputs), num_processors):
            for t2 in range(num_processors):
                index = t + t2
                if index < len(inputs):
                    inputs_divided[t2].append(inputs[index])

        # now, run each division on its own processor
        running = multiprocessing.Value("i", num_processors)
        mutex = multiprocessing.Lock()

        arrays = []
        threads = []
        for _ in range(num_processors):
            threads.append(collect_data_from_frames())
            arrays.append(multiprocessing.Array("i", [0, 1]))

        results_queue = multiprocessing.Queue()  # to keep track of the results

        processes = []
        for i in range(num_processors):
            p = multiprocessing.Process(
                target=threads[i].runit,
                args=(running, mutex, results_queue, inputs_divided[i]),
            )
            p.start()
            processes.append(p)

        while running.value > 0:
            continue  # wait for everything to finish

        # compile all results
        total_summed_coordinates = None
        dictionary_of_node_lists = {}
        for _ in threads:
            chunk = results_queue.get()

            if total_summed_coordinates is None:
                total_summed_coordinates = chunk[0]
            else:
                total_summed_coordinates = total_summed_coordinates + chunk[0]

            for key in chunk[1].keys():
                try:
                    dictionary_of_node_lists[key].extend(chunk[1][key])
                except Exception:
                    dictionary_of_node_lists[key] = chunk[1][key]

        self.combined_results = (total_summed_coordinates, dictionary_of_node_lists)


class collect_data_from_frames:
    """PDB-frame data processing on a single processor"""

    summed_coordinates = None
    nodes = {}

    def runit(self, running, mutex, results_queue, items):
        """PDB-frame data processing on a single processor

        Arguments:
        running -- a multiprocessing.Value() object
        mutex -- a multiprocessing.Lock() object
        results_queue -- where the results will be stored [multiprocessing.Queue()]
        items -- the data to be processed, in a list

        """
        for item in items:
            self.value_func(item)  # , results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put((self.summed_coordinates, self.nodes))

    def value_func(
        self, params_and_res_keys_and_pdb_lines_and_res_maps
    ):  # , results_queue): # so overwriting this function
        """Process a single PDB frame: identify the relevant nodes

        Arguments:
        params_and_res_keys_and_pdb_lines_and_res_maps -- a tuple containing required information.
             The first item contains user-defined parameters (a UserInput object)
             The second item is a list containing string representations of each residue ("CHAIN_RESNAME_RESID")
             The third item is a list of strings representing the PDB frame to be processed, where each string
                  contains a PDB ATOM or HETATM entry
             The fourth item is a dictionary that maps residue string identifiers ("CHAIN_RESNAME_RESID") to a list
                  of the indices of the atoms that correspond to that residue
        """

        # user-defined parameters
        params = params_and_res_keys_and_pdb_lines_and_res_maps[
            0
        ]

        # make sure this is not empty
        pdb_lines = params_and_res_keys_and_pdb_lines_and_res_maps[
            1
        ]

        # now load the frame into its own Molecule object
        pdb = Molecule()
        pdb.load_pdb_from_list(pdb_lines)

        if self.summed_coordinates is None:
            self.summed_coordinates = pdb.coordinates
        else:
            self.summed_coordinates = self.summed_coordinates + pdb.coordinates

        pdb.map_atoms_to_residues()
        pdb.map_nodes_to_residues(params["node_definition"])

        for index, residue_iden in enumerate(pdb.residue_identifiers_in_order):
            try:
                self.nodes[residue_iden].append(pdb.nodes[index])
            except:
                self.nodes[residue_iden] = [pdb.nodes[index]]


# --- second, to identify possible paths ---#
class multi_threading_find_paths:
    """Launches path finding on multiple processors"""

    results = []

    def __init__(self, inputs, num_processors):
        """Launches path finding on multiple processors

        Arguments:
        inputs -- the data to be processed, in a list
        num_processors -- the number of processors to use to process this data, an integer
        """

        self.results = []

        # first, if num_processors <= 0, determine the number of processors to
        # use programatically
        if num_processors <= 0:
            num_processors = multiprocessing.cpu_count()

        # reduce the number of processors if too many have been specified
        if len(inputs) < num_processors:
            num_processors = len(inputs)

        # now, divide the inputs into the appropriate number of processors
        inputs_divided = {t: [] for t in range(num_processors)}
        for t in range(0, len(inputs), num_processors):
            for t2 in range(num_processors):
                index = t + t2
                if index < len(inputs):
                    inputs_divided[t2].append(inputs[index])

        # now, run each division on its own processor
        running = multiprocessing.Value("i", num_processors)
        mutex = multiprocessing.Lock()

        arrays = []
        threads = []
        for _ in range(num_processors):
            threads.append(find_paths())
            arrays.append(multiprocessing.Array("i", [0, 1]))

        results_queue = multiprocessing.Queue()  # to keep track of the results

        processes = []
        for i in range(num_processors):
            p = multiprocessing.Process(
                target=threads[i].runit,
                args=(running, mutex, results_queue, inputs_divided[i]),
            )
            p.start()
            processes.append(p)

        while running.value > 0:
            continue  # wait for everything to finish

        # compile all results
        for _ in threads:
            chunk = results_queue.get()
            for chun in chunk:
                self.results.extend(chun)


class find_paths:  # other, more specific classes with inherit this one
    """Path-finding data processing on a single processor"""

    results = []

    def runit(self, running, mutex, results_queue, items):
        """Path-finding data processing on a single processor

        Arguments:
        running -- a multiprocessing.Value() object
        mutex -- a multiprocessing.Lock() object
        results_queue -- where the results will be stored [multiprocessing.Queue()]
        items -- the data to be processed, in a list
        """

        for item in items:
            self.value_func(item, results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    def value_func(
        self, item, results_queue
    ):  # this is the function that changes through inheritance
        """Process a single path-finding "branch"

        Arguments:
        item -- a tuple containing required information.
             The first is a numpy array containing a single float, the path-length cutoff
             The second is an index corresponding to the ultimate path sink
             The third is a networkx.Graph object describing the connectivity of the different nodes
             The fourth is a list corresponding to a path. The first item is the length of the path (float).
                  The remaining items are the indices of the nodes in the path (int).
        results_queue -- where the results will be stored [multiprocessing.Queue()]
        """

        cutoff = item[0]
        sink = item[1]
        G = item[2]

        pths_growing_frm_src = [item[3]]
        full_pth_strt_to_snk = []

        while pths_growing_frm_src:
            pths_growing_frm_src, full_pth_strt_to_snk = self.expand_pths_one_step(
                pths_growing_frm_src,
                full_pth_strt_to_snk,
                cutoff,
                sink,
                G,
            )

        # here save the results for later compilation
        self.results.append(full_pth_strt_to_snk)

    def expand_pths_one_step(
        self,
        pths_growing_frm_src,
        full_pth_strt_to_snk,
        cutoff,
        sink,
        G,
    ):
        """Expand the paths growing out from the source to the sink by one step
           (to the neighbors of the terminal node) of the expanding paths

        Arguments:
        pths_growing_frm_src -- a list of paths, where each path is represented by a list. The first item in each path
             is the length of the path (float). The remaining items are the indices of the nodes in the path (int).
        full_pth_strt_to_snk -- a growing list of identified paths that connect the source and the sink, where each
             path is formatted as above.
        cutoff -- a numpy array containing a single element (float), the length cutoff. Paths with lengths greater than the cutoff
             will be ignored.
        sink -- the index of the sink (int)
        G -- a networkx.Graph object describing the connectivity of the different nodes
        """

        for i, pth_growing_frm_src in enumerate(pths_growing_frm_src):
            if pth_growing_frm_src[0] > cutoff:
                # Because if the path is already greater than the cutoff, no
                # use continuing to branch out, since subsequent branhes will
                # be longer.
                pths_growing_frm_src.pop(i)
                break
            elif pth_growing_frm_src[-1] == sink:
                # so the sink has been reached
                full_pth_strt_to_snk.append(pth_growing_frm_src)
                pths_growing_frm_src.pop(i)
                break
            elif pth_growing_frm_src[-1] != sink:
                # sink not yet reached, but paths still short enough. So add
                # new paths, same as old, but with neighboring element
                # appended.
                # print("====================")
                # # print(pths_growing_frm_src)
                # print("\n".join([str(p) for p in pths_growing_frm_src]))
                node_neighbors = list(G.neighbors(pth_growing_frm_src[-1]))
                for j, node_neighbor in enumerate(node_neighbors):
                    if not node_neighbor in pth_growing_frm_src:
                        temp = pth_growing_frm_src[:]
                        temp.append(node_neighbor)
                        temp[0] = temp[0] + G.edges[temp[-2], temp[-1]]["weight"]
                        pths_growing_frm_src.insert((i + j + 1), temp)
                pths_growing_frm_src.pop(i)
                # print("")
                # # print(pths_growing_frm_src)
                # print("\n".join([str(p) for p in pths_growing_frm_src]))
                # print("")
                # print("")
                # import pdb; pdb.set_trace()
                break
            else:
                print("SOMETHING IS WRONG")

        # print(sorted(pths_growing_frm_src, reverse=True)[0])
        # import pdb; pdb.set_trace()
        return pths_growing_frm_src, full_pth_strt_to_snk

    # def expand_pths_one_step(
    #     self,
    #     pths_growing_frm_src,
    #     full_pth_strt_to_snk,
    #     cutoff,
    #     sink,
    #     G,
    # ):
    #     """Expand the paths growing out from the source to the sink by one step
    #        (to the neighbors of the terminal node) of the expanding paths

    #     Arguments:
    #     pths_growing_frm_src -- a list of paths, where each path is represented by a list. The first item in each path
    #          is the length of the path (float). The remaining items are the indices of the nodes in the path (int).
    #     full_pth_strt_to_snk -- a growing list of identified paths that connect the source and the sink, where each
    #          path is formatted as above.
    #     cutoff -- a numpy array containing a single element (float), the length cutoff. Paths with lengths greater than the cutoff
    #          will be ignored.
    #     sink -- the index of the sink (int)
    #     G -- a networkx.Graph object describing the connectivity of the different nodes
    #     """

    #     # ****

    #     # First, remove all paths that are already greater thanthe cutoff. No
    #     # use contrinuing to branch out, since subsequent branches will be
    #     # longer.
    #     pths_growing_frm_src = [
    #         path for path in pths_growing_frm_src if path[0] <= cutoff
    #     ]

    #     # Now, remove all paths that have already reached the sink. Add these
    #     # to the full paths list.
    #     full_pth_strt_to_snk.extend(
    #         [path for path in pths_growing_frm_src if path[-1] == sink]
    #     )
    #     pths_growing_frm_src = [
    #         path for path in pths_growing_frm_src if path[-1] != sink
    #     ]

    #     # Now, add new paths, same as old, but with neighboring element
    #     # appended.
    #     expanded_pths = []
    #     for pth_growing_frm_src in pths_growing_frm_src:
    #         # sink not yet reached, but paths still short enough. So add
    #         # new paths, same as old, but with neighboring element
    #         # appended.
    #         # print("====================")
    #         # # print(pths_growing_frm_src)
    #         # print("\n".join([str(p) for p in pths_growing_frm_src]))
    #         last_node_in_pth = pth_growing_frm_src[-1]
    #         node_neighbors = list(G.neighbors(last_node_in_pth))
    #         for node_neighbor in node_neighbors:
    #             if node_neighbor not in pth_growing_frm_src:
    #                 # It's a copy
    #                 updated_pth = pth_growing_frm_src[:]

    #                 # Add the neighbor
    #                 updated_pth.append(node_neighbor)

    #                 # Update the path length
    #                 updated_pth[0] = updated_pth[0] + G.edges[updated_pth[-2], updated_pth[-1]]["weight"]

    #                 # Insert this new path into the list of growing paths
    #                 # insert_idx = path_idx + neighbor_idx + 1

    #                 # print(insert_idx, len(pths_growing_frm_src))
    #                 # pths_growing_frm_src.insert(insert_idx, updated_pth)
    #                 expanded_pths.append(updated_pth)

    #         # print("")
    #         # # print(pths_growing_frm_src)
    #         # print("\n".join([str(p) for p in pths_growing_frm_src]))
    #         # print("")
    #         # print("")
    #         # import pdb; pdb.set_trace()

    #     # print(len(expanded_pths))
    #     return expanded_pths, full_pth_strt_to_snk


######################## To Generate Covariant Matrix ########################


class GetCovarianceMatrix:
    """Calculate and store the covariance matrix"""

    def __init__(self, params):
        """Calculates a covariance matrix

        Arguments:
        params -- user-specified command-line parameters (a UserInput object)
        """

        # first, split the file into frames. ^END matches both VMD and ENDMDL
        # formats.
        afile = open(params["pdb_trajectory_filename"], "r")
        this_frame = []
        first_frame = True
        number_of_frames = 0

        log(
            "\n# Loading frames from the PDB file and building the covariance matrix...",
            params["logfile"],
        )

        if params["number_processors"] == 1:

            load_frames_data = collect_data_from_frames()

            # a pdb object that will eventually contain the average structure
            self.average_pdb = None

            while 1:
                line = afile.readline()
                if not line:
                    break  # until eof

                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    this_frame.append(line)
                if line[:3] == "END":  # so reached end of frame
                    if len(this_frame) == 0:
                        # Happens sometimes with VMD-formatted files that end in
                        # ENDMDL, then END, so a frame has no coordinates in it.
                        continue

                    if first_frame:
                        self.average_pdb = Molecule()
                        self.average_pdb.load_pdb_from_list(this_frame)
                        first_frame = False

                    load_frames_data.value_func((params, this_frame))

                    this_frame = []  # so deleted for next time

                    log(
                        "#      Loading frame " + str(number_of_frames) + "...",
                        params["logfile"],
                    )
                    number_of_frames = number_of_frames + 1

            total_coordinate_sum = load_frames_data.summed_coordinates
            dictionary_of_node_lists = load_frames_data.nodes

            # because the previous coordinates belonged to the first frame
            self.average_pdb.coordinates = total_coordinate_sum / float(
                number_of_frames
            )

        else:
            # so more than one processor. Load in 100 frames, work on those.
            multiple_frames = []

            # this will keep a tallied sum of the coordinates of each frame
            # for subsequently calculating the average structure
            total_coordinate_sum = None

            dictionary_of_node_lists = {}

            # a pdb object that will eventually contain the average structure
            self.average_pdb = None

            while 1:
                line = afile.readline()
                if not line:
                    break  # until eof

                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    this_frame.append(line)
                if line[:3] == "END":  # so reached end of frame

                    if first_frame:
                        # self.get_residue_mappings(this_frame)
                        self.average_pdb = Molecule()
                        self.average_pdb.load_pdb_from_list(this_frame)
                        first_frame = False

                    if len(this_frame) > 0:
                        # Because sometimes ENDMDL followed by END => empty
                        # frame.
                        multiple_frames.append((params, this_frame))
                    this_frame = []  # so deleted for next time

                    if (
                        number_of_frames
                        % params["num_frames_to_load_before_processing"]
                        == 0
                    ):
                        # so you've collected 100 frames. Time to send them
                        # off to the multiple processes. note that the results
                        # are cumulative within the object.
                        tmp = multi_threading_to_collect_data_from_frames(
                            multiple_frames, params["number_processors"]
                        ).combined_results

                        if total_coordinate_sum is None:
                            total_coordinate_sum = tmp[0]
                            dictionary_of_node_lists = tmp[1]
                        else:
                            total_coordinate_sum = total_coordinate_sum + tmp[0]
                            for key in tmp[1]:
                                try:
                                    dictionary_of_node_lists[key].extend(tmp[1][key])
                                except:
                                    dictionary_of_node_lists[key] = tmp[1][key]

                        # so you're done processing the 100 frames, start over
                        # with the next 100
                        multiple_frames = []

                    log(
                        "#      Loading frame " + str(number_of_frames) + "...",
                        params["logfile"],
                    )
                    number_of_frames = number_of_frames + 1

            log("\n# Analyzing frames...", params["logfile"])

            # you need to get the last chunck
            tmp = multi_threading_to_collect_data_from_frames(
                multiple_frames, params["number_processors"]
            ).combined_results  # note that the results are cumulative within the object
            if total_coordinate_sum is None:
                total_coordinate_sum = tmp[0]
                dictionary_of_node_lists = tmp[1]
            else:
                total_coordinate_sum = total_coordinate_sum + tmp[0]
                for key in tmp[1]:
                    try:
                        dictionary_of_node_lists[key].extend(tmp[1][key])
                    except:
                        dictionary_of_node_lists[key] = tmp[1][key]

            self.average_pdb.coordinates = total_coordinate_sum / float(
                number_of_frames
            )  # because the previous coordinates belonged to the first frame

        afile.close()

        # numpyify dictionary_of_node_lists
        for res_iden in dictionary_of_node_lists:
            dictionary_of_node_lists[res_iden] = numpy.array(
                dictionary_of_node_lists[res_iden], numpy.float64
            )

        # now process the data that has been loaded
        # now get the average location of each node

        log("#      Saving the average PDB file...", params["logfile"])
        self.average_pdb.save_pdb(params["output_directory"] + "average_structure.pdb")

        log(
            "#      Calculating the average location of each node...", params["logfile"]
        )
        self.average_pdb.map_atoms_to_residues()
        self.average_pdb.map_nodes_to_residues(params["node_definition"])

        # now compute a set of deltas for each node, stored in a big array. delta = distance from node to average node location
        # so note that the nodes do need to be computed for each frame
        log(
            "#      Calculating the correlation for each node-node pair...",
            params["logfile"],
        )
        set_of_deltas = {}
        for index, residue_iden in enumerate(
            self.average_pdb.residue_identifiers_in_order
        ):
            set_of_deltas[residue_iden] = (
                dictionary_of_node_lists[residue_iden] - self.average_pdb.nodes[index]
            )

        # now compute the ensemble average of the deltas dot-producted with themselves
        ensmeble_average_deltas_self_dotproducted = {}
        for residue_iden in self.average_pdb.residue_identifiers_in_order:
            dot_products = (
                set_of_deltas[residue_iden] * set_of_deltas[residue_iden]
            ).sum(axis=1)
            ensmeble_average_deltas_self_dotproducted[residue_iden] = numpy.average(
                dot_products
            )

        # now build the correlation matrix
        if params["user_specified_functionalized_matrix_filename"] == "":
            log("#      Building the correlation matrix...", params["logfile"])
            self.correlations = numpy.empty(
                (
                    len(self.average_pdb.residue_identifiers_in_order),
                    len(self.average_pdb.residue_identifiers_in_order),
                )
            )

            for x in range(len(self.average_pdb.residue_identifiers_in_order)):
                residue1_key = self.average_pdb.residue_identifiers_in_order[x]
                for y in range(len(self.average_pdb.residue_identifiers_in_order)):
                    residue2_key = self.average_pdb.residue_identifiers_in_order[y]

                    # first, you need to calculate the dot products in the numerator
                    residue1_deltas = set_of_deltas[residue1_key]
                    residue2_deltas = set_of_deltas[residue2_key]

                    if len(residue1_deltas) != len(residue2_deltas):
                        log(
                            "ERROR: There were "
                            + str(len(residue1_deltas))
                            + ' residues that matched "'
                            + residue1_key
                            + '", but '
                            + str(len(residue2_deltas))
                            + ' residues that matched "'
                            + residue2_key
                            + '". Are each of your residues uniquely defined?',
                            params["logfile"],
                        )
                        sys.exit(0)

                    # generate a list of the dot products for all frames
                    dot_products = (residue1_deltas * residue2_deltas).sum(axis=1)

                    ensemble_average_dot_products = numpy.average(dot_products)

                    C = ensemble_average_dot_products / numpy.power(
                        ensmeble_average_deltas_self_dotproducted[residue1_key]
                        * ensmeble_average_deltas_self_dotproducted[residue2_key],
                        0.5,
                    )

                    self.correlations[x][y] = -numpy.log(
                        numpy.fabs(C)
                    )  # functionalizing the covariances
        else:  # so the user has specified a filename containing the covariance matrix
            log(
                "#      Loading the user-specified functionalized correlation matrix from the file "
                + params["user_specified_functionalized_matrix_filename"],
                params["logfile"],
            )
            self.correlations = numpy.loadtxt(
                params["user_specified_functionalized_matrix_filename"], dtype=float
            )

        # save the correlation matrix in a human-readable format
        numpy.savetxt(
            params["output_directory"] + "functionalized_correlation_matrix.txt",
            self.correlations,
        )

        # now modify the coorelation matrix, setting to 0 wherever the average distance between nodes is greater than a given cutoff
        contact_map = numpy.ones(self.correlations.shape)
        if params["user_specified_contact_map_filename"] == "":
            if params["contact_map_distance_limit"] != 999999.999:
                log(
                    "#      Applying the default WISP distance-based contact-map filter to the matrix so that distant residues will never be considered correlated...",
                    params["logfile"],
                )
                for index1 in range(
                    len(self.average_pdb.residue_identifiers_in_order) - 1
                ):
                    residue_iden1 = self.average_pdb.residue_identifiers_in_order[
                        index1
                    ]
                    residue1_pts = self.average_pdb.coordinates[
                        self.average_pdb.residue_identifier_to_atom_indices[
                            residue_iden1
                        ]
                    ]
                    for index2 in range(
                        index1 + 1, len(self.average_pdb.residue_identifiers_in_order)
                    ):
                        residue_iden2 = self.average_pdb.residue_identifiers_in_order[
                            index2
                        ]
                        residue2_pts = self.average_pdb.coordinates[
                            self.average_pdb.residue_identifier_to_atom_indices[
                                residue_iden2
                            ]
                        ]
                        min_dist_between_residue_atoms = numpy.min(
                            cdist(residue1_pts, residue2_pts)
                        )
                        if (
                            min_dist_between_residue_atoms
                            > params["contact_map_distance_limit"]
                        ):
                            # so they are far apart
                            self.correlations[index1][index2] = 0.0
                            self.correlations[index2][index1] = 0.0
                            contact_map[index1][index1] = 0.0
                            contact_map[index2][index1] = 0.0
        else:  # so the user has specified a contact map
            log(
                "#      Loading and applying the user-specified contact map from the file "
                + params["user_specified_contact_map_filename"],
                params["logfile"],
            )
            contact_map = numpy.loadtxt(
                params["user_specified_contact_map_filename"], dtype=float
            )
            self.correlations = self.correlations * contact_map

        # save the contact map in a human-readable format
        numpy.savetxt(
            params["output_directory"] + "contact_map_matrix.txt", contact_map
        )

        # now save the matrix if needed
        if (
            params["wisp_saved_matrix_filename"] != ""
        ):  # because it only would have gotten here if load_wisp_saved_matrix = FALSE
            pickle.dump(self, open(params["wisp_saved_matrix_filename"], "wb"))

    def convert_list_of_residue_keys_to_residue_indices(self, list_residue_keys):
        """Identify the indices in a networkx.Graph object corresponding to the identified residue string ids (CHAIN_RESNAME_RESID).

        Arguments:
        list_residue_keys -- a list of strings representing protein residues (format: CHAIN_RESNAME_RESID)

        Returns a list of ints, the networkx.Graph indices corresponding to the residue string ids in list_residue_keys
        """

        networkx_residue_indices = []
        for key in list_residue_keys:
            try:
                index_of_key = numpy.nonzero(
                    self.average_pdb.residue_identifiers_in_order == key
                )[0][0]
            except:
                raise Exception(
                    "ERROR: The residue "
                    + key
                    + " was not found in the average structure. Are you sure that the residue is present in the PDB file?"
                )
            networkx_residue_indices.append(index_of_key)
        return networkx_residue_indices

    def __getitem__(self, _):
        return self.correlations


######################### To Identify Paths ##############################


class GetPaths:
    """Get the paths from a list of sources to a list of sinks"""

    def __init__(self, corr_matrix, srcs, snks, params, residue_keys):
        """Identify paths that link the source and the sink and order them by their lengths

        Arguments:
        corr_matrix -- a numpy.array, the calculated correlation matrix
        srcs -- a list of ints, the indices of the sources for path finding
        snks -- a list of ints, the indices of the sinks for path finding
        params -- the user-specified command-line parameters, a UserInput object
        residue_keys -- a list containing string representations of each residue
        """

        # populate graph nodes and weighted edges
        G = networkx.Graph(incoming_graph_data=corr_matrix)

        # first calculate length of shortest path between any source and sink
        log("\n# Calculating paths...", params["logfile"])
        log(
            "#       Calculating the shortest path between any of the specified sources and any of the specified sinks...",
            params["logfile"],
        )
        shortest_length, shortest_path = self.get_shortest_path_length(
            corr_matrix, srcs, snks, G
        )
        log(
            f"#           The shortest path has length {str(shortest_length)}",
            params["logfile"],
        )

        path = [shortest_length]
        path.extend(shortest_path)
        
        # need to create this initial path in case only one path is requrested
        pths = [path]  

        cutoff = shortest_length

        cutoff_yields_max_num_paths_below_target = 0
        cutoff_yields_min_num_paths_above_target = 1000000.0

        # first step, keep incrementing a little until you have more than the
        # desired number of paths
        log(
            "#      Identifying the cutoff required to produce "
            + str(params["desired_number_of_paths"])
            + " paths...",
            params["logfile"],
        )
        num_paths = 1
        while num_paths < params["desired_number_of_paths"]:
            log(f"#          Testing the cutoff {str(cutoff)}...", params["logfile"])
            cutoff_in_array = numpy.array([cutoff], numpy.float64)
            pths_redundant = self.get_paths_between_multiple_endpoints(
                cutoff_in_array, corr_matrix, srcs, snks, G, params
            )
            pths = self.remove_redundant_paths(pths_redundant)
            num_paths = len(pths)

            log(
                f"#                The cutoff {str(cutoff)} produces {num_paths} paths...",
                params["logfile"],
            )

            if (
                num_paths < params["desired_number_of_paths"]
                and cutoff > cutoff_yields_max_num_paths_below_target
            ):
                cutoff_yields_max_num_paths_below_target = cutoff
            if (
                num_paths > params["desired_number_of_paths"]
                and cutoff < cutoff_yields_min_num_paths_above_target
            ):
                cutoff_yields_min_num_paths_above_target = cutoff

            # Original code adds .1 each time... but this may be to fast for
            # some systems and very slow for others... lets try increasing by
            # a percentage of the minimum path length instead... ideally this
            # could be an input parameter in the future.
            cutoff = cutoff + shortest_length * 0.1

        pths = self.remove_redundant_paths(pths)

        pths.sort()  # sort the paths by length

        if (
            num_paths != params["desired_number_of_paths"]
        ):  # so further refinement is needed
            pths = pths[: params["desired_number_of_paths"]]
            log(
                "#          Keeping the first "
                + str(params["desired_number_of_paths"])
                + " of these paths...",
                params["logfile"],
            )

        self.paths_description = ""

        self.paths_description = (
            self.paths_description + "\n# Output identified paths" + "\n"
        )
        index = 1

        if params["simply_formatted_paths_filename"] != "":
            simp = open(params["simply_formatted_paths_filename"], "w")
        for path in pths:
            self.paths_description = (
                f"{self.paths_description}#     Path {str(index)}:" + "\n"
            )
            self.paths_description = (
                f"{self.paths_description}#          Length: {str(path[0])}" + "\n"
            )
            self.paths_description = (
                f"{self.paths_description}#          Nodes: "
                + " - ".join([residue_keys[item] for item in path[1:]])
                + "\n"
            )
            if params["simply_formatted_paths_filename"] != "":
                simp.write(" ".join([str(item) for item in path]) + "\n")
            index = index + 1
        if params["simply_formatted_paths_filename"] != "":
            simp.close()

        self.paths = pths

    def remove_redundant_paths(self, pths):
        """Removes redundant paths

        Arguments:
        pths -- a list of paths

        Returns a list of paths with the redundant ones eliminated
        """

        if len(pths) == 1:
            # no reason to check if there's only one
            return pths

        for indx1 in range(len(pths) - 1):
            path1 = pths[indx1]
            if path1 is not None:
                for indx2 in range(indx1 + 1, len(pths)):
                    path2 = pths[indx2]
                    if path2 is not None and len(path1) == len(
                        path2
                    ):  # paths are the same length
                        pth1 = copy.deepcopy(path1[1:])
                        pth2 = copy.deepcopy(path2[1:])

                        if pth1[0] < pth1[-1]:
                            pth1.reverse()
                        if pth2[0] < pth2[-1]:
                            pth2.reverse()

                        if pth1 == pth2:
                            pths[indx2] = None

        while None in pths:
            pths.remove(None)

        return pths

    def get_shortest_path_length(
        self, corr_matrix, srcs, snks, G
    ):  # where sources and sinks are lists
        """Identify the length of the shortest path connecting any of the sources and any of the sinks

        Arguments:
        corr_matrix -- a numpy.array, the calculated correlation matrix
        srcs -- a list of ints, the indices of the sources for path finding
        snks -- a list of ints, the indices of the sinks for path finding
        G -- a networkx.Graph object describing the connectivity of the different nodes

        Returns a float, the length of the shortest path, and a list of ints corresponding to
             the nodes of the shortest path
        """

        shortest_length = 99999999.999
        shortest_path = []

        for source in srcs:
            for sink in snks:
                if source != sink:  # important to avoid this situation
                    short_path = networkx.dijkstra_path(
                        G, source, sink, weight="weight"
                    )
                    length = self.get_length_of_path(short_path, corr_matrix)
                    if length < shortest_length:
                        shortest_length = length
                        shortest_path = short_path
        return shortest_length, shortest_path

    def get_length_of_path(self, path, corr_matrix):
        """Calculate the length of a path

        Arguments:
        path -- a list of ints, the indices of the path
        corr_matrix -- a numpy.array, the calculated correlation matrix

        Returns a float, the length of the path
        """

        length = 0.0
        for t in range(len(path) - 1):
            length = length + corr_matrix[path[t], path[t + 1]]
        return length

    def get_paths_between_multiple_endpoints(
        self, cutoff, corr_matrix, srcs, snks, G, params
    ):  # where sources and sinks are lists
        """Get paths between sinks and sources

        Arguments:
        cutoff -- a numpy.array containing a single float, the cutoff specifying the maximum permissible path length
        corr_matrix -- a numpy.array, the calculated correlation matrix
        srcs -- a list of ints, the indices of the sources for path finding
        snks -- a list of ints, the indices of the sinks for path finding
        G -- a networkx.Graph object describing the connectivity of the different nodes
        params -- the user-specified command-line parameters, a UserInput object

        Returns a list of paths, where each path is represented by a list. The first item in each path is the length
             of the path (float). The remaining items are the indices of the nodes in the path (int).
        """

        pths = []
        for source in srcs:
            for sink in snks:
                if source != sink:  # avoid this situation
                    pths_to_add = self.get_paths_fixed_endpoints(
                        cutoff, corr_matrix, source, sink, G, params
                    )
                    pths.extend(pths_to_add)
        return pths

    def get_paths_fixed_endpoints(self, cutoff, corr_matrix, source, sink, G, params):
        """Get paths between a single sink and a single source

        Arguments:
        cutoff -- a numpy.array containing a single float, the cutoff specifying the maximum permissible path length
        corr_matrix -- a numpy.array, the calculated correlation matrix
        source -- the index of the source for path finding
        sink -- the index of the sink for path finding
        G -- a networkx.Graph object describing the connectivity of the different nodes
        params -- the user-specified command-line parameters, a UserInput object

        Returns a list of paths, where each path is represented by a list. The first item in each path is the length
             of the path (float). The remaining items are the indices of the nodes in the path (int).
        """

        if source == sink:
            return []

        source_lengths, source_paths = networkx.single_source_dijkstra(
            G, source, target=None, cutoff=None, weight="weight"
        )
        sink_lengths, sink_paths = networkx.single_source_dijkstra(
            G, sink, target=None, cutoff=None, weight="weight"
        )

        so_l = [source_lengths[key] for key in source_lengths.keys()]
        so_p = [source_paths[key] for key in source_paths.keys()]
        si_l = [sink_lengths[key] for key in sink_lengths.keys()]
        si_p = [sink_paths[key] for key in sink_paths.keys()]

        check_list_1 = []
        check_list_2 = []
        for i in range(len(so_l)):
            check_list_1.extend([so_p[i][-1]])
            check_list_2.extend([si_p[i][-1]])

        node_list = []
        dijkstra_list = []
        upper_minimum_length = 0
        if not set(check_list_1).difference(check_list_2):
            for i, _ in enumerate(so_l):
                if so_l[i] + si_l[i] <= cutoff:
                    node_list.extend(so_p[i][:])
                    node_list.extend(si_p[i][:])
                    si_pReversed = si_p[i][:]
                    si_pReversed.reverse()
                    temp_path = so_p[i][:] + si_pReversed[1:]
                    temp_length = so_l[i] + si_l[i]
                    dijkstra_list.append(temp_path)
                    if (so_l[i] + si_l[i]) > upper_minimum_length:
                        upper_minimum_length = temp_length
        else:
            print("paths do not match up")

        unique_nodes = list(set(node_list))
        unique_nodes.sort()

        node_length = len(unique_nodes)
        new_matrix = numpy.zeros((len(corr_matrix), len(corr_matrix)))

        for i in range(node_length):
            for j in range(node_length):
                new_matrix[unique_nodes[i]][unique_nodes[j]] = corr_matrix[
                    unique_nodes[i]
                ][unique_nodes[j]]

        corr_matrix = new_matrix
        G = networkx.Graph(incoming_graph_data=corr_matrix, labels=unique_nodes)

        length = 0.0
        pths_growing_frm_src = [[length, source]]
        full_pth_strt_to_snk = []

        # This is essentially this list-addition replacement for a recursive
        # algorithm you've envisioned. To parallelize, just get the first N
        # branches, and send them off to each node. Rest of branches filled out
        # in separate processes.

        find_paths_object = find_paths()
        if params["number_processors"] == 1:
            while pths_growing_frm_src:
                pths_growing_frm_src, full_pth_strt_to_snk = find_paths_object.expand_pths_one_step(
                    pths_growing_frm_src,
                    full_pth_strt_to_snk,
                    cutoff,
                    sink,
                    G,
                )
        else:
            # just get some of the initial paths on a single processor
            log(
                "#                Starting serial portion of path-finding algorithm (will run for "
                + str(params["seconds_to_wait_before_parallelizing_path_finding"])
                + " seconds)...",
                params["logfile"],
            )
            atime = time.time()
            while (
                pths_growing_frm_src
                and time.time() - atime
                < params["seconds_to_wait_before_parallelizing_path_finding"]
            ):
                pths_growing_frm_src, full_pth_strt_to_snk = find_paths_object.expand_pths_one_step(
                    pths_growing_frm_src,
                    full_pth_strt_to_snk,
                    cutoff,
                    sink,
                    G,
                )

            # ok, so having generated just a first few, divy up those among multiple processors
            if pths_growing_frm_src:  # in case you've already finished
                log(
                    "#                Starting parallel portion of path-finding algorithm running on "
                    + str(params["number_processors"])
                    + " processors...",
                    params["logfile"],
                )
                pths_growing_frm_src = [
                    (cutoff, sink, G, path) for path in pths_growing_frm_src
                ]
                additional_full_pth_strt_to_snk = multi_threading_find_paths(
                    pths_growing_frm_src, params["number_processors"]
                )
                full_pth_strt_to_snk.extend(
                    additional_full_pth_strt_to_snk.results
                )
            else:
                log(
                    "#                     (All paths found during serial path finding; parallelization not required)",
                    params["logfile"],
                )

        full_pth_strt_to_snk.sort()

        pths = []

        for full_path_from_start_to_sink in full_pth_strt_to_snk:
            pths.append(full_path_from_start_to_sink)

        return pths


class Visualize:
    """A class to facilitate the visualization of the identified paths in VMD"""

    def __init__(self, params, corr_matrix_object, pths):
        """Get paths between a single sink and a single source

        Arguments:
        params -- the user-specified command-line parameters, a UserInput object
        residue_keys -- a list containing string representations of each residue ("CHAIN_RESNAME_RESID")
        node_locs -- a dictionary, mapping string representations of each residue to a numpy.array representation
             of the node location
        pths -- a GetPaths object
        """

        log_files = [
            params["logfile"],
            open(params["output_directory"] + "visualize.tcl", "w"),
            open(params["output_directory"] + "visualize.vmd", "w"),
        ]

        # output a easy-to-read-representation of the paths
        log(pths.paths_description, log_files)

        if (
            params["longest_path_opacity"] == params["shortest_path_opacity"]
            and params["longest_path_opacity"] == params["node_sphere_opacity"]
        ):
            opacity_required = False
        else:
            opacity_required = True

        log("", params["logfile"])
        log("# Creating VMD state (TCL) file for visualization...", log_files)

        # get the color range
        lengths = [t[:1][0] for t in pths.paths]
        min_length = min(lengths)
        max_length = max(lengths)

        # Get the sizes of the various strands
        ratios = (
            []
        )  # so the shortest (best) path has a ratio of 0, and longest has a ratio of 1
        for length in lengths:
            if len(pths.paths) > 1:
                ratio = (length - min_length) / float(max_length - min_length)
            else:
                ratio = 0.5

            ratios.append(ratio)

        # define the colors
        a1 = numpy.array(
            [
                params["shortest_path_r"],
                params["shortest_path_g"],
                params["shortest_path_b"],
            ],
            numpy.float64,
        )
        a2 = numpy.array(
            [
                params["longest_path_r"],
                params["longest_path_g"],
                params["longest_path_b"],
            ],
            numpy.float64,
        )
        color_defs = {}
        for coloridd in range(23, 33):  # 33 because I want it to go to 32
            thiscolor = (
                ((a2 - a1) / 9.0) * coloridd + (32.0 / 9.0) * a1 - (23.0 / 9.0) * a2
            )
            color_defs[coloridd] = (
                "color change rgb "
                + str(coloridd)
                + " "
                + str(thiscolor[0])
                + " "
                + str(thiscolor[1])
                + " "
                + str(thiscolor[2])
            )
        color_defs[coloridd] = (
            "color change rgb 22 "
            + str(params["node_sphere_r"])
            + " "
            + str(params["node_sphere_g"])
            + " "
            + str(params["node_sphere_b"])
        )  # node sphere color

        # determine whether the average structure or a user-specified structure will be used for drawing
        if params["pdb_single_frame_filename"] != "":  # use a user-specified structure
            molecule_object_to_use = Molecule()
            molecule_object_to_use.load_pdb_from_list(
                open(params["pdb_single_frame_filename"], "r").readlines()
            )
            molecule_object_to_use.map_atoms_to_residues()
            molecule_object_to_use.map_nodes_to_residues(params["node_definition"])
            molecule_object_to_use.save_pdb(
                params["output_directory"] + "draw_frame.pdb"
            )
            molecule_filename = "draw_frame.pdb"
        else:  # so use the average structure
            molecule_object_to_use = corr_matrix_object.average_pdb
            molecule_filename = "average_structure.pdb"

        # get all the nodes
        nodes_used = []
        for path in pths.paths:
            for index in path[1:]:
                if not index in nodes_used:
                    nodes_used.append(index)

        node_ids = [
            corr_matrix_object.average_pdb.residue_identifiers_in_order[i]
            for i in nodes_used
        ]
        node_ids = [item.split("_") for item in node_ids]
        selection = "".join(
            [
                "(resid " + item[2] + " and chain " + item[0] + ") or "
                for item in node_ids
            ]
        )[:-4]

        log("\n# Load in the protein structure", log_files)
        log(
            "mol new "
            + molecule_filename
            + " type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            log_files,
        )
        log("mol delrep 0 top", log_files)
        log("mol representation NewCartoon 0.300000 10.000000 4.100000 0", log_files)
        log("mol color Name", log_files)
        log("mol selection {all}", log_files)
        log("mol material Opaque", log_files)
        log("mol addrep top", log_files)
        log("mol representation Licorice 0.300000 10.000000 10.000000", log_files)
        log("mol color Name", log_files)
        log("mol selection {" + selection + "}", log_files)
        log("mol material Opaque", log_files)
        log("mol addrep top", log_files)
        log("mol rename top " + molecule_filename, log_files)

        if (
            params["node_sphere_radius"] != 0.0
        ):  # if the radius is 0.0, don't even draw the spheres
            # draw spheres
            log("\n# Draw spheres at the nodes", log_files)
            if opacity_required:
                log(
                    "if {[lsearch [material list] node_spheres] == -1} {material add node_spheres}",
                    log_files,
                )
                log(
                    "material change opacity node_spheres "
                    + str(params["node_sphere_opacity"]),
                    log_files,
                )
                log("draw material node_spheres", log_files)
                log('mol rename top "Node Spheres"', log_files)

            nodes = [molecule_object_to_use.nodes[i] for i in nodes_used]
            log("graphics top color 22", log_files)

            for node in nodes:
                log(
                    "draw sphere {"
                    + str(node[0])
                    + " "
                    + str(node[1])
                    + " "
                    + str(node[2])
                    + "} resolution "
                    + str(params["vmd_resolution"])
                    + " radius "
                    + str(params["node_sphere_radius"]),
                    log_files,
                )

        color_index = 0
        log(
            "set wisp_num_paths %i" % len(pths.paths), log_files
        )  # tell the WISP plugin how many paths there are total
        for path in pths.paths:

            log("\n# Draw a new path", log_files)
            log("# \tLength: " + str(path[0]), log_files)
            log("# \tNodes: " + str(path[1:]) + "\n", log_files)

            # make the spline from the paths

            nodes = [corr_matrix_object.average_pdb.nodes[i] for i in path[1:]]

            x_vals = []
            y_vals = []
            z_vals = []

            for node in nodes:
                x_vals.append(node[0])
                y_vals.append(node[1])
                z_vals.append(node[2])

            ratio = ratios[color_index]
            color = (
                numpy.floor(ratio * 9) + 23
            )  # so 0.95 => 9.5 => 9 => 32 (red); 0.05 => 0.5 > 0 > 23 (pretty close to blue, but not perfect blue)

            radius = (
                ratio * (params["longest_path_radius"] - params["shortest_path_radius"])
                + params["shortest_path_radius"]
            )
            opacity = (
                ratio
                * (params["longest_path_opacity"] - params["shortest_path_opacity"])
                + params["shortest_path_opacity"]
            )

            if opacity_required:
                log("mol new", log_files)
            log(color_defs[color], log_files)
            log("graphics top color " + str(int(color)), log_files)

            mat_name = "wisp_material" + str(color_index)
            if opacity_required:
                log(
                    "if {[lsearch [material list] %s] == -1} {material add %s}"
                    % (mat_name, mat_name),
                    log_files,
                )
                log(
                    "material change opacity " + mat_name + " " + str(opacity),
                    log_files,
                )
                log('mol rename top "Path Length: ' + str(path[0]) + '"', log_files)
                log("draw material " + mat_name, log_files)

            try:

                degree = len(x_vals) - 1
                if degree > 3:
                    # so at most degree 3
                    degree = 3

                tck, _ = interpolate.splprep([x_vals, y_vals, z_vals], s=0, k=degree)

                # now interpolate
                unew = numpy.arange(0, 1.01, params["spline_smoothness"])

                out = interpolate.splev(unew, tck)

                for t in range(len(out[0]) - 1):
                    x1 = str(out[0][t])
                    y1 = str(out[1][t])
                    z1 = str(out[2][t])

                    x2 = str(out[0][t + 1])
                    y2 = str(out[1][t + 1])
                    z2 = str(out[2][t + 1])

                    log(
                        "draw cylinder {"
                        + x1
                        + " "
                        + y1
                        + " "
                        + z1
                        + "} {"
                        + x2
                        + " "
                        + y2
                        + " "
                        + z2
                        + "} radius "
                        + str(radius)
                        + " resolution "
                        + str(params["vmd_resolution"])
                        + " filled 0",
                        log_files,
                    )

            except:  # so just draw a single cylinder as a backup
                log(
                    "draw cylinder {"
                    + str(x_vals[0])
                    + " "
                    + str(y_vals[0])
                    + " "
                    + str(z_vals[0])
                    + "} {"
                    + str(x_vals[-1])
                    + " "
                    + str(y_vals[-1])
                    + " "
                    + str(z_vals[-1])
                    + "} radius "
                    + str(radius)
                    + " resolution "
                    + str(params["vmd_resolution"])
                    + " filled 0",
                    log_files,
                )

            color_index = color_index + 1


def output_directory_info(params):
    """Create a README.txt file in the output directory describing the directory contents

    Arguments:
    params -- the user-specified command-line parameters, a UserInput object
    """

    f = open(params["output_directory"] + "README.txt", "w")
    f.write(
        "This directory contains output from the program WISP. The best way to visualize the output is to use a free program called VMD, which can be downloaded from http://www.ks.uiuc.edu/Research/vmd/ ."
        + "\n\n"
    )
    f.write(
        'The WISP output can be automatically loaded into VMD using the TCL script named "visualize.tcl". Assuming "vmd" is the full path to your installed VMD executable, just run the following from the command line:'
        + "\n\n"
    )
    f.write("vmd -e visualize.tcl" + "\n\n")
    f.write(
        'If you prefer not to use the command line, simply run the vmd executable and load the "visualize.vmd" file using "File->Load Visualization State..." from the main menu.'
        + "\n\n"
    )
    f.write(
        'The above methods are very slow. If your output is so large that a faster option is required, the Tk Console can be used. Use "Extensions->Tk Console" from the VMD main menu to pull up the Tk Console. Then run the following command, with the full path to "visualize.tcl" included if necessary:'
        + "\n\n"
    )
    f.write("source visualize.tcl" + "\n\n")
    f.write(
        'Regardless of the method you use to load in the WISP output, the visualization will be the same. Individual pathways are shown as tubes (i.e., "wisps"), the protein is shown in ribbon representation, and protein residues that participate in any path are shown in licorice representation.'
        + "\n\n"
    )
    f.write(
        "The WISP output directory contains a number of other files as well. Here are descriptions of each:"
        + "\n\n"
    )
    f.write("log.txt: Details describing WISP execution." + "\n\n")
    f.write(
        "parameters_used.txt: The WISP parameters used to generate the output." + "\n\n"
    )
    f.write(
        "average_structure.pdb: The average structure of your PDB trajectory." + "\n\n"
    )
    f.write(
        'draw_frame.pdb: If the user requests that a separate single-structure PDB file be used for calculating node and wisp positions, that file is saved as "draw_frame.pdb". Otherwise, the average structure is used.'
        + "\n\n"
    )
    f.write(
        'functionalized_matrix_with_contact_map_applied.pickle: A python pickle file that contains the matrix obtained by multiplying a functionalized correlation matrix and a contact map. This file is not human readable but can be loaded into WISP for use in subsequent runs with the -load_wisp_saved_matrix and -wisp_saved_matrix_filename parameters. Thus, the matrix needs only to be calculated once for each trajectory, rather than every time WISP is executed. Use "python wisp.py -help" for more information.'
        + "\n\n"
    )
    f.write(
        "contact_map_matrix.txt: A human readable representation of the contact map. If the user wishes to generate their own contact map rather than letting WISP generate one automatically, a custom contact map formatted like this one can be loaded into WISP using the -user_specified_contact_map_filename parameter. "
        + "\n\n"
    )
    f.write(
        "functionalized_correlation_matrix.txt: A human readable representation of the functionalized correlation matrix, prior to multiplication by the contact map. If the user wishes to generate their own functionalized correlation matrix rather than letting WISP generate one automatically, a custom matrix formatted like this one can be loaded into WISP using the -user_specified_functionalized_matrix_filename parameter. "
        + "\n\n"
    )
    f.write(
        "simply_formatted_paths.txt: A simple list of path lengths and nodes. The first column contains the lengths, and all following columns contain node indices. This file may be helpful for subsequent statistical analyses of the WISP output."
        + "\n"
    )
    f.close()


if __name__ == "__main__":
    program_start_time = time.time()

    # get the commandline parameters
    parameters = UserInput()

    # compute the correlation matrix
    if parameters["load_wisp_saved_matrix"] == "TRUE":
        correlation_matrix_object = pickle.load(
            open(parameters["wisp_saved_matrix_filename"], "rb")
        )  # load the matrix instead of generating
    else:
        correlation_matrix_object = GetCovarianceMatrix(
            parameters
        )  # so generate the matrix instead of loading it
    correlation_matrix = correlation_matrix_object.correlations

    # always save a copy of the correlation matrix, regardless of how it was loaded/generated
    pickle.dump(
        correlation_matrix_object,
        open(
            parameters["output_directory"]
            + "functionalized_matrix_with_contact_map_applied.pickle",
            "wb",
        ),
    )

    # now get the source and sink locations from the parameters
    sources = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        parameters["source_residues"]
    )
    sinks = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        parameters["sink_residues"]
    )

    # compute the paths
    paths = GetPaths(
        correlation_matrix,
        sources,
        sinks,
        parameters,
        correlation_matrix_object.average_pdb.residue_identifiers_in_order,
    )

    # create the visualization
    vis = Visualize(parameters, correlation_matrix_object, paths)

    # provide the user information about the generated files
    output_directory_info(parameters)

    log(
        "\n# Program execution time: "
        + str(time.time() - program_start_time)
        + " seconds",
        parameters["logfile"],
    )

    # Only close if it isn't a string
    if type(parameters["logfile"]) != str:
        parameters["logfile"].close()
