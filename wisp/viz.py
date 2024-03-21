import os

import numpy as np
from scipy import interpolate

from .structure import Molecule


def log(astring, fileobjects, print_viz=False):  # prints to screen and to log file
    """Outputs WISP messages

    Args:
        astring: a string containing the message
        fileobjects: a list of python file objects specifying where the messages
            hould be saved.
        print_viz: Print viz commands as well.
    """

    if not isinstance(fileobjects, list):
        # it's not a list, so make it one
        fileobjects = [fileobjects]

    if print_viz:
        print(astring)

    for fileobject in fileobjects:
        fileobject.write(astring + "\n")


class Visualize:
    """A class to facilitate the visualization of the identified paths in VMD"""

    def __init__(self, params, corr_matrix_object, pths):
        """Get paths between a single sink and a single source

        Args:
            params: the user-specified command-line parameters, a UserInput object
            residue_keys: a list containing string representations of each residue ("CHAIN_RESNAME_RESID")
            node_locs: a dictionary, mapping string representations of each residue to a np.array representation
                of the node location
            pths: a GetPaths object
        """

        log_files = [
            params["logfile"],
            open(os.path.join(params["output_directory"], "visualize.tcl"), "w"),
            open(os.path.join(params["output_directory"], "visualize.vmd"), "w"),
        ]

        # output a easy-to-read-representation of the paths
        log(pths.paths_description, log_files)

        if (
            params["longest_path_opacity"] == params["shortest_path_opacity"]
            and params["longest_path_opacity"] == params["node_sphere_opacity"]
        ):
            opacity_required = False
        else:
            opacity_required = True

        log("", params["logfile"])
        log("# Creating VMD state (TCL) file for visualization...", log_files)

        # get the color range
        lengths = [t[:1][0] for t in pths.paths]
        min_length = min(lengths)
        max_length = max(lengths)

        # Get the sizes of the various strands
        ratios = (
            []
        )  # so the shortest (best) path has a ratio of 0, and longest has a ratio of 1
        for length in lengths:
            if len(pths.paths) > 1:
                ratio = (length - min_length) / float(max_length - min_length)
            else:
                ratio = 0.5

            ratios.append(ratio)

        # define the colors
        a1 = np.array(
            [
                params["shortest_path_r"],
                params["shortest_path_g"],
                params["shortest_path_b"],
            ],
            np.float64,
        )
        a2 = np.array(
            [
                params["longest_path_r"],
                params["longest_path_g"],
                params["longest_path_b"],
            ],
            np.float64,
        )
        color_defs = {}
        for coloridd in range(23, 33):  # 33 because I want it to go to 32
            thiscolor = (
                ((a2 - a1) / 9.0) * coloridd + (32.0 / 9.0) * a1 - (23.0 / 9.0) * a2
            )
            color_defs[coloridd] = (
                "color change rgb "
                + str(coloridd)
                + " "
                + str(thiscolor[0])
                + " "
                + str(thiscolor[1])
                + " "
                + str(thiscolor[2])
            )
        color_defs[coloridd] = (
            "color change rgb 22 "
            + str(params["node_sphere_r"])
            + " "
            + str(params["node_sphere_g"])
            + " "
            + str(params["node_sphere_b"])
        )  # node sphere color

        # determine whether the average structure or a user-specified structure will be used for drawing
        if params["pdb_single_frame_filename"] != "":  # use a user-specified structure
            molecule_object_to_use = Molecule()
            molecule_object_to_use.load_pdb_from_list(
                open(params["pdb_single_frame_filename"]).readlines()
            )
            molecule_object_to_use.map_atoms_to_residues()
            molecule_object_to_use.map_nodes_to_residues(params["node_definition"])
            molecule_object_to_use.save_pdb(
                os.path.join(params["output_directory"], "draw_frame.pdb")
            )
            molecule_filename = "draw_frame.pdb"
        else:  # so use the average structure
            molecule_object_to_use = corr_matrix_object.average_pdb
            molecule_filename = "average_structure.pdb"

        # get all the nodes
        nodes_used = []
        for path in pths.paths:
            for index in path[1:]:
                if not index in nodes_used:
                    nodes_used.append(index)

        node_ids = [
            corr_matrix_object.average_pdb.residue_identifiers_in_order[i]
            for i in nodes_used
        ]
        node_ids = [item.split("_") for item in node_ids]
        selection = "".join(
            [
                "(resid " + item[2] + " and chain " + item[0] + ") or "
                for item in node_ids
            ]
        )[:-4]

        log("\n# Load in the protein structure", log_files)
        log(
            "mol new "
            + molecule_filename
            + " type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            log_files,
        )
        log("mol delrep 0 top", log_files)
        log("mol representation NewCartoon 0.300000 10.000000 4.100000 0", log_files)
        log("mol color Name", log_files)
        log("mol selection {all}", log_files)
        log("mol material Opaque", log_files)
        log("mol addrep top", log_files)
        log("mol representation Licorice 0.300000 10.000000 10.000000", log_files)
        log("mol color Name", log_files)
        log("mol selection {" + selection + "}", log_files)
        log("mol material Opaque", log_files)
        log("mol addrep top", log_files)
        log("mol rename top " + molecule_filename, log_files)

        if (
            params["node_sphere_radius"] != 0.0
        ):  # if the radius is 0.0, don't even draw the spheres
            # draw spheres
            log("\n# Draw spheres at the nodes", log_files)
            if opacity_required:
                log(
                    "if {[lsearch [material list] node_spheres] == -1} {material add node_spheres}",
                    log_files,
                )
                log(
                    "material change opacity node_spheres "
                    + str(params["node_sphere_opacity"]),
                    log_files,
                )
                log("draw material node_spheres", log_files)
                log('mol rename top "Node Spheres"', log_files)

            nodes = [molecule_object_to_use.nodes[i] for i in nodes_used]
            log("graphics top color 22", log_files)

            for node in nodes:
                log(
                    "draw sphere {"
                    + str(node[0])
                    + " "
                    + str(node[1])
                    + " "
                    + str(node[2])
                    + "} resolution "
                    + str(params["vmd_resolution"])
                    + " radius "
                    + str(params["node_sphere_radius"]),
                    log_files,
                )

        color_index = 0
        log(
            "set wisp_num_paths %i" % len(pths.paths), log_files
        )  # tell the WISP plugin how many paths there are total
        for path in pths.paths:
            log("\n# Draw a new path", log_files)
            log("# \tLength: " + str(path[0]), log_files)
            log("# \tNodes: " + str(path[1:]) + "\n", log_files)

            # make the spline from the paths

            nodes = [corr_matrix_object.average_pdb.nodes[i] for i in path[1:]]

            x_vals = []
            y_vals = []
            z_vals = []

            for node in nodes:
                x_vals.append(node[0])
                y_vals.append(node[1])
                z_vals.append(node[2])

            ratio = ratios[color_index]
            color = (
                np.floor(ratio * 9) + 23
            )  # so 0.95 => 9.5 => 9 => 32 (red); 0.05 => 0.5 > 0 > 23 (pretty close to blue, but not perfect blue)

            radius = (
                ratio * (params["longest_path_radius"] - params["shortest_path_radius"])
                + params["shortest_path_radius"]
            )
            opacity = (
                ratio
                * (params["longest_path_opacity"] - params["shortest_path_opacity"])
                + params["shortest_path_opacity"]
            )

            if opacity_required:
                log("mol new", log_files)
            log(color_defs[color], log_files)
            log("graphics top color " + str(int(color)), log_files)

            mat_name = "wisp_material" + str(color_index)
            if opacity_required:
                log(
                    "if {[lsearch [material list] %s] == -1} {material add %s}"
                    % (mat_name, mat_name),
                    log_files,
                )
                log(
                    "material change opacity " + mat_name + " " + str(opacity),
                    log_files,
                )
                log('mol rename top "Path Length: ' + str(path[0]) + '"', log_files)
                log("draw material " + mat_name, log_files)

            try:
                degree = len(x_vals) - 1
                if degree > 3:
                    # so at most degree 3
                    degree = 3

                tck, _ = interpolate.splprep([x_vals, y_vals, z_vals], s=0, k=degree)

                # now interpolate
                unew = np.arange(0, 1.01, params["spline_smoothness"])

                out = interpolate.splev(unew, tck)

                for t in range(len(out[0]) - 1):
                    x1 = str(out[0][t])
                    y1 = str(out[1][t])
                    z1 = str(out[2][t])

                    x2 = str(out[0][t + 1])
                    y2 = str(out[1][t + 1])
                    z2 = str(out[2][t + 1])

                    log(
                        "draw cylinder {"
                        + x1
                        + " "
                        + y1
                        + " "
                        + z1
                        + "} {"
                        + x2
                        + " "
                        + y2
                        + " "
                        + z2
                        + "} radius "
                        + str(radius)
                        + " resolution "
                        + str(params["vmd_resolution"])
                        + " filled 0",
                        log_files,
                    )

            except Exception:  # so just draw a single cylinder as a backup
                log(
                    "draw cylinder {"
                    + str(x_vals[0])
                    + " "
                    + str(y_vals[0])
                    + " "
                    + str(z_vals[0])
                    + "} {"
                    + str(x_vals[-1])
                    + " "
                    + str(y_vals[-1])
                    + " "
                    + str(z_vals[-1])
                    + "} radius "
                    + str(radius)
                    + " resolution "
                    + str(params["vmd_resolution"])
                    + " filled 0",
                    log_files,
                )

            color_index = color_index + 1
