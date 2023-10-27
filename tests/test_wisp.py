import os
import shutil
import subprocess

import numpy as np

from wisp.run import run_wisp

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

FILE_DIR = "./files/"
WRITING_DIR = "./tmp/"


def test_example():
    pdb_path = os.path.join(FILE_DIR, "trajectory_20_frames.pdb")
    test_dir = os.path.join(WRITING_DIR, "test_example")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)
    config = {
        "contact_map_distance_limit": 4.5,
        "desired_number_of_paths": 15,
        "load_wisp_saved_matrix": "FALSE",
        "wisp_saved_matrix_filename": "",
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
        "num_frames_to_load_before_processing": 20,
        "number_processors": 4,
        "pdb_trajectory_filename": pdb_path,
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
        "logfile": open(os.path.join(test_dir, "log.txt"), "w"),
        "user_specified_contact_map_filename": "",
        "user_specified_functionalized_matrix_filename": "",
        "simply_formatted_paths_filename": os.path.join(
            test_dir, "simply_formatted_paths.txt"
        ),
        "pdb_single_frame_filename": "",
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