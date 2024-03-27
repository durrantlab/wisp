import os
import shutil

import numpy as np
import pytest

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
    config = {
        "contact_map_distance_limit": 4.5,
        "n_paths": 15,
        "n_paths_max": 100000,
        "longest_path_b": 0.0,
        "longest_path_g": 0.0,
        "longest_path_opacity": 1.0,
        "longest_path_r": 1.0,
        "longest_path_radius": 0.01,
        "node_definition": "RESIDUE_COM",
        "node_sphere_b": 1.0,
        "node_sphere_g": 1.0,
        "node_sphere_opacity": 1.0,
        "node_sphere_r": 1.0,
        "node_sphere_radius": 1.0,
        "frame_chunks": 20,
        "n_cores": 4,
        "pdb_path": pdb_path,
        "seconds_to_wait_before_parallelizing_path_finding": 5.0,
        "shortest_path_b": 1.0,
        "shortest_path_g": 0.0,
        "shortest_path_opacity": 1.0,
        "shortest_path_r": 0.0,
        "shortest_path_radius": 1.0,
        "sink_residues": ["C_ASP_11"],
        "source_residues": ["C_LEU_10"],
        "spline_smoothness": 0.01,
        "vmd_resolution": 6,
        "output_directory": test_dir,
        "logfile": open(os.path.join(test_dir, "log.txt"), "w", encoding="utf-8"),
        "simply_formatted_paths_path": os.path.join(
            test_dir, "simply_formatted_paths.txt"
        ),
    }
    paths = run_wisp(config)
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
    config = {
        "contact_map_distance_limit": 4.5,
        "n_paths": 15,
        "load_wisp_saved_matrix": "FALSE",
        "wisp_saved_matrix_path": "",
        "n_paths_max": 100000,
        "longest_path_b": 0.0,
        "longest_path_g": 0.0,
        "longest_path_opacity": 1.0,
        "longest_path_r": 1.0,
        "longest_path_radius": 0.01,
        "node_definition": "RESIDUE_COM",
        "node_sphere_b": 1.0,
        "node_sphere_g": 1.0,
        "node_sphere_opacity": 1.0,
        "node_sphere_r": 1.0,
        "node_sphere_radius": 1.0,
        "frame_chunks": 1,
        "n_cores": 1,
        "pdb_path": pdb_path,
        "seconds_to_wait_before_parallelizing_path_finding": 5.0,
        "shortest_path_b": 1.0,
        "shortest_path_g": 0.0,
        "shortest_path_opacity": 1.0,
        "shortest_path_r": 0.0,
        "shortest_path_radius": 1.0,
        "sink_residues": ["C_HIE_840"],
        "source_residues": ["B_ASN_692"],
        "spline_smoothness": 0.01,
        "vmd_resolution": 6,
        "output_directory": test_dir,
        "logfile": open(os.path.join(test_dir, "log.txt"), "w", encoding="utf-8"),
        "contact_map_path": "",
        "functionalized_matrix_path": "",
        "simply_formatted_paths_path": os.path.join(
            test_dir, "simply_formatted_paths.txt"
        ),
        "pdb_single_frame_path": "",
    }
    with pytest.raises(SystemExit) as e:
        run_wisp(config)
