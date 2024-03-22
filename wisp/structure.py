import gc

import numpy as np
from loguru import logger


class Atom:
    """A class containing atomic information."""

    def read_pdb_line(self, Line: str, default_chain_id: str = "A") -> None:
        """Reads atomic information from a string formatted according to the PDB
        standard.

        Args:
            Line: A string formatted according to the PDB standard.
            default_chain_id: If the chain ID is missing, replace it with this.
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

        # If chain is not filled out in PDB
        self.updated_chain = False
        if self.chain == " ":
            self.chain = default_chain_id
            self.updated_chain = True

        # now get the resid
        try:
            self.resid = int(Line[22:26])
        except Exception:
            self.resid = 0

        self.coordinates_numpy = np.array(
            [float(Line[30:38]), float(Line[38:46]), float(Line[46:54])], np.float64
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
    """Loads, saves, and manipulates molecular models."""

    def load_pdb_from_list(self, alist):
        """Loads a list of PDB ATOM/HETATM lines into the current Molecule object.

        Args:
            alist: the list of PDB lines
        """

        gc.disable()

        # have to use python lists initially because not sure of size
        self.atomnames = []
        self.chains = []
        self.resids = []
        self.elements = []
        self.resnames = []
        self.coordinates = []

        default_chain_id = "A"
        for line in alist:
            # If no chain IDs are specified and there are multiple chains, we need
            # to correctly account for multiple chains. We can do this by keeping track
            # of any `TER` lines and move it to the next letter.
            if line[:3] == "TER":
                default_chain_id = chr(ord(default_chain_id) + 1)
                # We also warn the user that we are updating chains.
                if temp_atom.updated_chain:
                    logger.warning(
                        "We found multiple chains in PDB frame, but no chain IDs are provided"
                    )
                    logger.warning("We assume first chain is A, second is B, etc.")

            if len(line) >= 7 and (line[:4] == "ATOM" or line[:6] == "HETATM"):
                temp_atom = Atom()
                temp_atom.read_pdb_line(line, default_chain_id=default_chain_id)
                self.atomnames.append(temp_atom.atomname)
                self.chains.append(temp_atom.chain)
                self.resids.append(temp_atom.resid)
                self.elements.append(temp_atom.element)
                self.resnames.append(temp_atom.resname)
                self.coordinates.append(temp_atom.coordinates_numpy)

        # convert them into numpy arrays
        self.atomnames = np.array(self.atomnames)
        self.chains = np.array(self.chains)
        self.resids = np.array(self.resids)
        self.elements = np.array(self.elements)
        self.resnames = np.array(self.resnames)
        self.coordinates = np.array(self.coordinates, np.float64)

        gc.enable()

    def save_pdb(self, filename):
        """Saves a pdb file

        Args:
            filename: a string specifying the file name
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

        residue_identifiers_for_all_atoms = np.array(residue_identifiers_for_all_atoms)
        self.residue_identifiers_in_order = np.array(self.residue_identifiers_in_order)

        self.residue_identifier_to_atom_indices = {}
        for the_id in self.residue_identifiers_in_order:
            self.residue_identifier_to_atom_indices[the_id] = np.nonzero(
                residue_identifiers_for_all_atoms == the_id
            )[0]

    def get_indices_of_atoms_in_a_residue_by_atom_name(
        self, residue_identifier, atom_names_list, not_selection=False
    ):
        """Gets the indices of atoms in a specified residue

        Args:
            residue_identifier: a string (chain_resname_resid) specifying the residue
            atom_names_list: a list of strings containing the names of the atoms to keep
            not_selection: a optional boolean. if False, match the atom_names_list items.
                if True, match the items not in atom_names_list

        Returns:
            a numpy array containing the indices of the atoms to keep
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
        indices_to_keep = np.array(indices_to_keep, np.float64)

        if not indices_to_keep:
            raise Exception(
                f"No atoms found in residue {residue_identifier} with atom names {str(atom_names_list)}"
            )

        return indices_to_keep

    def get_center_of_mass_from_selection_by_atom_indices(self, indices_selection):
        """Gets the center of mass of a set of atoms

        Args:
            indices_selection: a numpy array containing the indices of the
                atoms to consider

        Returns a numpy array containing the 3D coordinates of the center of
        mass of those atoms
        """

        coors = self.coordinates[indices_selection]
        eles = self.elements[indices_selection]
        masses = np.empty(coors.shape)

        for t in range(masses.shape[0]):
            themass = self.get_mass(eles[t])
            masses[t][0] = themass
            masses[t][1] = themass
            masses[t][2] = themass

        center_of_mass = ((coors * masses).sum(axis=0)) * (1.0 / masses.sum(axis=0))
        return np.array(
            [center_of_mass[0], center_of_mass[1], center_of_mass[2]], np.float64
        )

    def get_mass(self, element_name):
        """A library to provide the mass of a given element

        Args:
            element_name: a string that specifies the element

        Returns:
            a float, the mass of the specified element. If the element is not in the
            library, returns None.
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

        Args:
            node_definition: a string describing the definition of the node: `CA`,
                `RESIDUE_COM`, `BACKBONE_COM`, or `SIDECHAIN_COM`
        """

        self.nodes = np.empty((len(self.residue_identifiers_in_order), 3))
        for index, residue_iden in enumerate(self.residue_identifiers_in_order):
            if node_definition == "CA":  # the node is at the alpha carbon
                indices_to_consider = (
                    self.get_indices_of_atoms_in_a_residue_by_atom_name(
                        residue_iden, ["CA"]
                    )
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
                indices_to_consider = (
                    self.get_indices_of_atoms_in_a_residue_by_atom_name(
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
                )
                node_loc = self.get_center_of_mass_from_selection_by_atom_indices(
                    indices_to_consider
                )
            elif (
                node_definition == "SIDECHAIN_COM"
            ):  # the node is the residue center of mass
                indices_to_consider = (
                    self.get_indices_of_atoms_in_a_residue_by_atom_name(
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
                )
                node_loc = self.get_center_of_mass_from_selection_by_atom_indices(
                    indices_to_consider
                )
            self.nodes[index][0] = node_loc[0]
            self.nodes[index][1] = node_loc[1]
            self.nodes[index][2] = node_loc[2]
