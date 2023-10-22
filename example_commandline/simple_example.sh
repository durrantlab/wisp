#!/usr/bin/env bash
# Run wisp
python ../wisp.py -contact_map_distance_limit 4.5 -desired_number_of_paths 15 -load_wisp_saved_matrix FALSE -longest_path_b 0.0 -longest_path_g 0.0 -longest_path_opacity 1.0 -longest_path_r 1.0 -longest_path_radius 0.01 -node_definition RESIDUE_COM -node_sphere_b 1.0 -node_sphere_g 1.0 -node_sphere_opacity 1.0 -node_sphere_r 1.0 -node_sphere_radius 1.0 -num_frames_to_load_before_proceesing 20 -number_processors 4 -pdb_trajectory_filename trajectory_20_frames.pdb -seconds_to_wait_before_parallelizing_path_finding 5.0 -shortest_path_b 1.0 -shortest_path_g 0.0 -shortest_path_opacity 1.0 -shortest_path_r 0.0 -shortest_path_radius 1.0 -sink_residues "C_ASP_11" -source_residues "C_LEU_10" -spline_smoothness 0.01 -vmd_resolution 6 -output_directory ./example_output/

# Get expected contents of example_output/simply_formatted_paths.txt
echo
echo "Expected contents of example_output/simply_formatted_paths.txt (Python 3, numpy 1.17.4, networkx 2.2, scipy 1.3.1):"

cat <<EOF
1.1363589537262389 9 10
2.027418997284591 9 11 10
2.205901870787125 9 49 11 10
2.307449850568712 9 49 10
2.442623103273746 9 49 51 10
2.448087831142992 9 8 10
2.501932824731498 9 11 51 10
2.5026984744745864 9 48 11 10
2.621657347053458 9 49 48 11 10
2.6370708377784413 9 31 15 10
2.6804156982340315 9 49 11 51 10
2.7056328296310643 9 31 11 10
2.8932868716011986 9 49 51 11 10
2.977212301921493 9 48 11 51 10
3.0311375788312436 9 31 16 10
EOF
echo

# Get actual contents of example_output/simply_formatted_paths.txt
echo "Actual contents of example_output/simply_formatted_paths.txt:"
cat example_output/simply_formatted_paths.txt
