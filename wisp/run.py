import os
import pickle
import time
from collections.abc import MutableMapping
from typing import Any

import numpy as np
from loguru import logger

from .contexts import ContextManager
from .io import output_directory_info
from .paths import GetPaths
from .utils import GetCovarianceMatrix
from .viz import Visualize


def run_wisp(context: MutableMapping[str, Any]) -> None:
    # Update default context
    context_manager = ContextManager()
    context = context_manager.update(context)

    program_start_time = time.time()

    # compute the correlation matrix
    if context["wisp_saved_matrix_path"] is not None:
        correlation_matrix_object = pickle.load(
            open(context["wisp_saved_matrix_path"], "rb")
        )  # load the matrix instead of generating
    else:
        correlation_matrix_object = GetCovarianceMatrix(
            context
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
                context["output_directory"],
                "functionalized_matrix_with_contact_map_applied.pickle",
            ),
            "wb",
        ),
    )

    # now get the source and sink locations from the parameters
    sources = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        context["source_residues"]
    )
    sinks = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        context["sink_residues"]
    )

    # compute the paths
    paths = GetPaths(
        correlation_matrix,
        sources,
        sinks,
        context,
        correlation_matrix_object.average_pdb.residue_identifiers_in_order,
    )

    # create the visualization
    Visualize(context, correlation_matrix_object, paths)

    # provide the user information about the generated files
    output_directory_info(context)

    logger.info(
        "Program execution time: " + str(time.time() - program_start_time) + " seconds"
    )

    return paths.paths
