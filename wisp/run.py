import os
import pickle
import time

import numpy as np
from loguru import logger

from .cli import run_cli
from .config import WispConfig
from .io import output_dir_info
from .paths import GetPaths
from .utils import GetCovarianceMatrix
from .viz import Visualize


def run_wisp(wisp_config: WispConfig) -> None:
    program_start_time = time.time()

    # compute the correlation matrix
    if wisp_config.path_saved_matrix is not None:
        correlation_matrix_object = pickle.load(
            open(wisp_config.path_saved_matrix, "rb")
        )  # load the matrix instead of generating
    else:
        correlation_matrix_object = GetCovarianceMatrix(
            wisp_config
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
                wisp_config.output_dir,
                "functionalized_matrix_with_contact_map_applied.pickle",
            ),
            "wb",
        ),
    )

    # now get the source and sink locations from the parameters
    sources = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        wisp_config.source_residues
    )
    sinks = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        wisp_config.sink_residues
    )

    # compute the paths
    paths = GetPaths(
        correlation_matrix,
        sources,
        sinks,
        wisp_config,
        correlation_matrix_object.average_pdb.residue_identifiers_in_order,
    )

    # create the visualization
    Visualize(wisp_config, correlation_matrix_object, paths)

    # provide the user information about the generated files
    output_dir_info(wisp_config)

    logger.info(
        "Program execution time: " + str(time.time() - program_start_time) + " seconds"
    )

    return paths.paths


def main():
    wisp_config = run_cli()
    run_wisp(wisp_config)
