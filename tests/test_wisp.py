import os
import shutil

import numpy as np
import pytest

from wisp.config import WispConfig
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
    wisp_config = WispConfig()
    wisp_config.pdb_path = pdb_path
    wisp_config.n_paths = 15
    wisp_config.contact_map_distance_limit = 4.5
    wisp_config.output_dir = test_dir
    wisp_config.source_residues = ["C_LEU_10"]
    wisp_config.sink_residues = ["C_ASP_11"]
    wisp_config.n_cores = 4
    wisp_config.frame_chunks = 20
    wisp_config.write_formatted_paths = True

    paths = run_wisp(wisp_config)
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

    wisp_config = WispConfig()
    wisp_config.pdb_path = pdb_path
    wisp_config.n_paths = 15
    wisp_config.contact_map_distance_limit = 4.5
    wisp_config.output_dir = test_dir
    wisp_config.sink_residues = ["C_HIE_840"]
    wisp_config.source_residues = ["B_ASN_692"]
    wisp_config.n_cores = 1
    wisp_config.frame_chunks = 1

    with pytest.raises(SystemExit):
        run_wisp(wisp_config)


def test_issue_5():
    pdb_path = os.path.join(FILE_DIR, "issue5.pdb")
    test_dir = os.path.join(WRITING_DIR, "test_issue_5")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)

    wisp_config = WispConfig()
    wisp_config.pdb_path = pdb_path
    wisp_config.n_paths = 15
    wisp_config.contact_map_distance_limit = 4.5
    wisp_config.output_dir = test_dir
    wisp_config.sink_residues = ["C_PLP_1449"]
    wisp_config.source_residues = ["C_PLP_645"]
    wisp_config.n_cores = 1

    with pytest.raises(RuntimeError):
        run_wisp(wisp_config)
