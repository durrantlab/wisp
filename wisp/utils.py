import os
import pickle
import sys
from collections.abc import Collection

import numpy as np
from loguru import logger
from scipy.spatial.distance import cdist

from .structure import Molecule
from .traj import collect_data_from_frames, multi_threading_to_collect_data_from_frames


class GetCovarianceMatrix:
    """Calculate and store the covariance matrix"""

    def __init__(self, context):
        """
        Args:
            context: WISP context for computing the covariance matrix.
        """

        # first, split the file into frames. ^END matches both VMD and ENDMDL
        # formats.
        afile = open(context["pdb_path"], mode="r", encoding="utf-8")
        this_frame = []
        first_frame = True
        number_of_frames = 0

        logger.info(
            "Loading frames from the PDB file and building the covariance matrix"
        )

        if context["n_cores"] == 1:
            load_frames_data = collect_data_from_frames()

            # a pdb object that will eventually contain the average structure
            self.average_pdb = None

            for line in afile:
                if line[:4] == "ATOM" or line[:6] == "HETATM" or line[:3] == "TER":
                    this_frame.append(line)
                if line.startswith(("END", "ENDMDL")):  # so reached end of frame
                    if first_frame:
                        self.average_pdb = Molecule()
                        self.average_pdb.load_pdb_from_list(this_frame)
                        first_frame = False

                    # Happens when you have ENDMDL and then END.
                    if len(this_frame) == 0:
                        break

                    load_frames_data.value_func((context, this_frame))

                    this_frame = []  # so deleted for next time

                    logger.trace("Loading frame {}", str(number_of_frames))
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

            for line in afile:
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
                        multiple_frames.append((context, this_frame))
                    this_frame = []  # so deleted for next time

                    if number_of_frames % context["frame_chunks"] == 0:
                        # so you've collected 100 frames. Time to send them
                        # off to the multiple processes. note that the results
                        # are cumulative within the object.
                        tmp = multi_threading_to_collect_data_from_frames(
                            multiple_frames, context["n_cores"]
                        ).combined_results

                        if total_coordinate_sum is None:
                            total_coordinate_sum = tmp[0]
                            dictionary_of_node_lists = tmp[1]
                        else:
                            total_coordinate_sum = total_coordinate_sum + tmp[0]
                            for key in tmp[1]:
                                try:
                                    dictionary_of_node_lists[key].extend(tmp[1][key])
                                except Exception:
                                    dictionary_of_node_lists[key] = tmp[1][key]

                        # so you're done processing the 100 frames, start over
                        # with the next 100
                        multiple_frames = []

                    logger.trace("Loading frame " + str(number_of_frames))
                    number_of_frames = number_of_frames + 1

            logger.info("Analyzing frames...")

            # you need to get the last chunk
            tmp = multi_threading_to_collect_data_from_frames(
                multiple_frames, context["n_cores"]
            ).combined_results  # note that the results are cumulative within the object
            if total_coordinate_sum is None:
                total_coordinate_sum = tmp[0]
                dictionary_of_node_lists = tmp[1]
            else:
                total_coordinate_sum = total_coordinate_sum + tmp[0]
                for key in tmp[1]:
                    try:
                        dictionary_of_node_lists[key].extend(tmp[1][key])
                    except Exception:
                        dictionary_of_node_lists[key] = tmp[1][key]

            self.average_pdb.coordinates = total_coordinate_sum / float(
                number_of_frames
            )  # because the previous coordinates belonged to the first frame

        afile.close()

        # numpyify dictionary_of_node_lists
        for res_iden in dictionary_of_node_lists:
            dictionary_of_node_lists[res_iden] = np.array(
                dictionary_of_node_lists[res_iden], np.float64
            )

        # now process the data that has been loaded
        # now get the average location of each node

        logger.info("Saving the average PDB file...")
        self.average_pdb.save_pdb(
            os.path.join(context["output_dir"], "average_structure.pdb")
        )

        logger.info("Calculating the average location of each node...")
        self.average_pdb.map_atoms_to_residues()
        self.average_pdb.map_nodes_to_residues(context["node_definition"])

        # now compute a set of deltas for each node, stored in a big array. delta = distance from node to average node location
        # so note that the nodes do need to be computed for each frame
        logger.info(
            "Calculating the correlation for each node-node pair...",
        )
        set_of_deltas = {}
        for index, residue_iden in enumerate(
            self.average_pdb.residue_identifiers_in_order
        ):
            set_of_deltas[residue_iden] = (
                dictionary_of_node_lists[residue_iden] - self.average_pdb.nodes[index]
            )

        # now compute the ensemble average of the deltas dot-producted with themselves
        ensemble_average_deltas_self_dotproducted = {}
        for residue_iden in self.average_pdb.residue_identifiers_in_order:
            dot_products = (
                set_of_deltas[residue_iden] * set_of_deltas[residue_iden]
            ).sum(axis=1)
            ensemble_average_deltas_self_dotproducted[residue_iden] = np.average(
                dot_products
            )

        # now build the correlation matrix
        if context["functionalized_matrix_path"] is None:
            logger.info("Building the correlation matrix...")
            self.correlations = np.empty(
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
                        logger.info(
                            "ERROR: There were "
                            + str(len(residue1_deltas))
                            + ' residues that matched "'
                            + residue1_key
                            + '", but '
                            + str(len(residue2_deltas))
                            + ' residues that matched "'
                            + residue2_key
                            + '". Are each of your residues uniquely defined?',
                        )
                        sys.exit(0)

                    # generate a list of the dot products for all frames
                    dot_products = (residue1_deltas * residue2_deltas).sum(axis=1)

                    ensemble_average_dot_products = np.average(dot_products)

                    C = ensemble_average_dot_products / np.power(
                        ensemble_average_deltas_self_dotproducted[residue1_key]
                        * ensemble_average_deltas_self_dotproducted[residue2_key],
                        0.5,
                    )

                    self.correlations[x][y] = -np.log(
                        np.fabs(C)
                    )  # functionalizing the covariances
        else:  # so the user has specified a filename containing the covariance matrix
            logger.info(
                "Loading the user-specified functionalized correlation matrix from the file "
                + context["functionalized_matrix_path"],
            )
            self.correlations = np.loadtxt(
                context["functionalized_matrix_path"], dtype=float
            )

        # save the correlation matrix in a human-readable format
        np.savetxt(
            os.path.join(
                context["output_dir"], "functionalized_correlation_matrix.txt"
            ),
            self.correlations,
        )

        # now modify the correlation matrix, setting to 0 wherever the average distance between nodes is greater than a given cutoff
        contact_map = np.ones(self.correlations.shape)
        if context["contact_map_path"] is None:
            if context["contact_map_distance_limit"] != 999999.999:
                logger.info(
                    "Applying the default WISP distance-based contact-map filter to the matrix so that distant residues will never be considered correlated...",
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
                        min_dist_between_residue_atoms = np.min(
                            cdist(residue1_pts, residue2_pts)
                        )
                        if (
                            min_dist_between_residue_atoms
                            > context["contact_map_distance_limit"]
                        ):
                            # so they are far apart
                            self.correlations[index1][index2] = 0.0
                            self.correlations[index2][index1] = 0.0
                            contact_map[index1][index1] = 0.0
                            contact_map[index2][index1] = 0.0
        else:  # so the user has specified a contact map
            logger.info(
                "Loading and applying the user-specified contact map from the file "
                + context["contact_map_path"],
            )
            contact_map = np.loadtxt(context["contact_map_path"], dtype=float)
            self.correlations = self.correlations * contact_map

        # save the contact map in a human-readable format
        np.savetxt(
            os.path.join(context["output_dir"], "contact_map_matrix.txt"),
            contact_map,
        )

        # now save the matrix if needed
        if context["wisp_saved_matrix_path"] is not None:
            pickle.dump(self, open(context["wisp_saved_matrix_path"], "wb"))

    def convert_list_of_residue_keys_to_residue_indices(
        self, list_residue_keys: Collection[str]
    ):
        """Identify the indices in a nx.Graph object corresponding to the identified residue string ids (CHAIN_RESNAME_RESID).

        Args:
            list_residue_keys: a list of strings representing protein residues
                (format: CHAIN_RESNAME_RESID)

        Returns:
            a list of ints, the nx.Graph indices corresponding to the residue
            string ids in list_residue_keys
        """
        networkx_residue_indices = []
        for key in list_residue_keys:
            try:
                index_of_key = np.nonzero(
                    self.average_pdb.residue_identifiers_in_order == key
                )[0][0]
            except IndexError as e:
                logger.critical(f"Cannot find {key} in the structure.")
                logger.critical(f"Please check your source and sink residues.")
                raise RuntimeError(f"Cannot find {key} in the structure") from e
            networkx_residue_indices.append(index_of_key)
        return networkx_residue_indices

    def __getitem__(self, _):
        return self.correlations
