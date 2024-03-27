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

import os
import pickle
import time

import numpy as np
from loguru import logger

from .cli import UserInput
from .io import output_directory_info
from .paths import GetPaths
from .utils import GetCovarianceMatrix
from .viz import Visualize


def run_wisp(config):
    program_start_time = time.time()

    # compute the correlation matrix
    if config["load_wisp_saved_matrix"] == "TRUE":
        correlation_matrix_object = pickle.load(
            open(config["wisp_saved_matrix_path"], "rb")
        )  # load the matrix instead of generating
    else:
        correlation_matrix_object = GetCovarianceMatrix(
            config
        )  # so generate the matrix instead of loading it
    correlation_matrix = correlation_matrix_object.correlations

    # Set all diagonal negative values to zero
    # See https://github.com/durrantlab/wisp/pull/2#issuecomment-1774429502
    correlation_matrix[np.diag_indices(correlation_matrix.shape[0])] = 0.0

    # always save a copy of the correlation matrix, regardless of how it was loaded/generated
    pickle.dump(
        correlation_matrix_object,
        open(
            os.path.join(
                config["output_directory"],
                "functionalized_matrix_with_contact_map_applied.pickle",
            ),
            "wb",
        ),
    )

    # now get the source and sink locations from the parameters
    sources = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        config["source_residues"]
    )
    sinks = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        config["sink_residues"]
    )

    # compute the paths
    paths = GetPaths(
        correlation_matrix,
        sources,
        sinks,
        config,
        correlation_matrix_object.average_pdb.residue_identifiers_in_order,
    )

    # create the visualization
    Visualize(config, correlation_matrix_object, paths)

    # provide the user information about the generated files
    output_directory_info(config)

    logger.info(
        "Program execution time: " + str(time.time() - program_start_time) + " seconds",
        config["logfile"],
    )

    config["logfile"].close()

    return paths.paths


def main():
    config = UserInput()
    run_wisp(config)


if __name__ == "__main__":
    main()
