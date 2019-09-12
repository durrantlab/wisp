# Run wisp
python ../wisp.py -contact_map_distance_limit 999999.999 -desired_number_of_paths 15 -load_wisp_saved_matrix FALSE -longest_path_b 0.0 -longest_path_g 0.0 -longest_path_opacity 1.0 -longest_path_r 1.0 -longest_path_radius 0.01 -node_definition RESIDUE_COM -node_sphere_b 1.0 -node_sphere_g 1.0 -node_sphere_opacity 1.0 -node_sphere_r 1.0 -node_sphere_radius 1.0 -num_frames_to_load_before_proceesing 20 -number_processors 4 -pdb_trajectory_filename trajectory_20_frames.pdb -seconds_to_wait_before_parallelizing_path_finding 5.0 -shortest_path_b 1.0 -shortest_path_g 0.0 -shortest_path_opacity 1.0 -shortest_path_r 0.0 -shortest_path_radius 1.0 -sink_residues "C_ASP_11" -source_residues "C_LEU_10" -spline_smoothness 0.01 -vmd_resolution 6 -output_directory ./example_output/

# Get expected contents of example_output/simply_formatted_paths.txt
echo
echo "Expected contents of example_output/simply_formatted_paths.txt (if Python 3):"
cat <<EOF
1.1363589537262389 9 10
1.920042178661312 9 31 10
2.027418997284591 9 11 10
2.0580187961654244 9 271 10
2.072362499435586 9 104 10
2.1979011747323174 9 260 10
2.205901870787125 9 49 11 10
2.2274854140539255 9 30 31 10
2.2709631887057427 9 104 271 10
2.3074498505687115 9 49 10
2.31609691842617 9 36 31 10
2.316805404722313 9 105 104 10
2.317004756199552 9 32 31 10
2.344862237079425 9 107 104 10
2.3571728828116276 9 107 329 10
EOF
echo

# Get actual contents of example_output/simply_formatted_paths.txt
echo "Actual contents of example_output/simply_formatted_paths.txt:"
cat example_output/simply_formatted_paths.txt
