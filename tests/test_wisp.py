import os
import shutil

import numpy as np
import pytest

from wisp.contexts import ContextManager
from wisp.run import run_wisp

# Ensures we execute from file directory (for relative paths).
current_dir = os.path.dirname(os.path.abspath(__file__))

FILE_DIR = os.path.join(current_dir, "files/")
WRITING_DIR = os.path.join(current_dir, "tmp/")


def test_example():
    pdb_path = os.path.join(FILE_DIR, "trajectory_20_frames.pdb")
    test_dir = os.path.join(WRITING_DIR, "test_example")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)
    context_manager = ContextManager()
    context_manager.pdb_path = pdb_path
    context_manager.n_paths = 15
    context_manager.contact_map_distance_limit = 4.5
    context_manager.output_dir = test_dir
    context_manager.source_residues = ["C_LEU_10"]
    context_manager.sink_residues = ["C_ASP_11"]
    context_manager.n_cores = 4
    context_manager.frame_chunks = 20
    context_manager.write_formatted_paths = True

    paths = run_wisp(context_manager)
    test_data = np.array([i[0] for i in paths], dtype=np.float64)

    ref_data = np.array(
        [
            1.1363589537262389,
            2.027418997284591,
            2.205901870787125,
            2.307449850568712,
            2.442623103273746,
            2.448087831142992,
            2.501932824731498,
            2.5026984744745864,
            2.621657347053458,
            2.6370708377784413,
            2.6804156982340315,
            2.7056328296310643,
            2.8932868716011986,
            2.977212301921493,
            3.0311375788312436,
        ],
        dtype=np.float64,
    )
    assert np.allclose(test_data, ref_data)


def test_issue_1_serial():
    pdb_path = os.path.join(FILE_DIR, "sample-issue-1.pdb")
    test_dir = os.path.join(WRITING_DIR, "test_issue_1")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)

    context_manager = ContextManager()
    context_manager.pdb_path = pdb_path
    context_manager.n_paths = 15
    context_manager.contact_map_distance_limit = 4.5
    context_manager.output_dir = test_dir
    context_manager.sink_residues = ["C_HIE_840"]
    context_manager.source_residues = ["B_ASN_692"]
    context_manager.n_cores = 1
    context_manager.frame_chunks = 1

    with pytest.raises(SystemExit):
        run_wisp(context_manager)


def test_issue_5():
    pdb_path = os.path.join(FILE_DIR, "issue5.pdb")
    test_dir = os.path.join(WRITING_DIR, "test_issue_5")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)

    context_manager = ContextManager()
    context_manager.pdb_path = pdb_path
    context_manager.n_paths = 15
    context_manager.contact_map_distance_limit = 4.5
    context_manager.output_dir = test_dir
    context_manager.sink_residues = ["C_PLP_1449"]
    context_manager.source_residues = ["C_PLP_645"]
    context_manager.n_cores = 1

    with pytest.raises(RuntimeError):
        run_wisp(context_manager)
