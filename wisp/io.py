import os


def output_directory_info(params):
    """Create a README.txt file in the output directory describing the directory
    contents.

    Args:
        params: The user-specified command-line parameters, a UserInput object
    """

    f = open(
        os.path.join(params["output_directory"], "README.txt"), "w", encoding="utf-8"
    )
    f.write(
        "This directory contains output from the program WISP. The best way to visualize the output is to use a free program called VMD, which can be downloaded from http://www.ks.uiuc.edu/Research/vmd/ ."
        + "\n\n"
    )
    f.write(
        'The WISP output can be automatically loaded into VMD using the TCL script named "visualize.tcl". Assuming "vmd" is the full path to your installed VMD executable, just run the following from the command line:'
        + "\n\n"
    )
    f.write("vmd -e visualize.tcl" + "\n\n")
    f.write(
        'If you prefer not to use the command line, simply run the vmd executable and load the "visualize.vmd" file using "File->Load Visualization State..." from the main menu.'
        + "\n\n"
    )
    f.write(
        'The above methods are very slow. If your output is so large that a faster option is required, the Tk Console can be used. Use "Extensions->Tk Console" from the VMD main menu to pull up the Tk Console. Then run the following command, with the full path to "visualize.tcl" included if necessary:'
        + "\n\n"
    )
    f.write("source visualize.tcl" + "\n\n")
    f.write(
        'Regardless of the method you use to load in the WISP output, the visualization will be the same. Individual pathways are shown as tubes (i.e., "wisps"), the protein is shown in ribbon representation, and protein residues that participate in any path are shown in licorice representation.'
        + "\n\n"
    )
    f.write(
        "The WISP output directory contains a number of other files as well. Here are descriptions of each:"
        + "\n\n"
    )
    f.write("log.txt: Details describing WISP execution." + "\n\n")
    f.write(
        "parameters_used.txt: The WISP parameters used to generate the output." + "\n\n"
    )
    f.write(
        "average_structure.pdb: The average structure of your PDB trajectory." + "\n\n"
    )
    f.write(
        'draw_frame.pdb: If the user requests that a separate single-structure PDB file be used for calculating node and wisp positions, that file is saved as "draw_frame.pdb". Otherwise, the average structure is used.'
        + "\n\n"
    )
    f.write(
        'functionalized_matrix_with_contact_map_applied.npy: A NumPy array that contains the matrix obtained by multiplying a functionalized correlation matrix and a contact map. This file is not human readable but can be loaded into WISP for use in subsequent runs with the -load_wisp_saved_matrix and -wisp_saved_matrix_filename parameters. Thus, the matrix needs only to be calculated once for each trajectory, rather than every time WISP is executed. Use "wisp -help" for more information.'
        + "\n\n"
    )
    f.write(
        "contact_map_matrix.txt: A human readable representation of the contact map. If the user wishes to generate their own contact map rather than letting WISP generate one automatically, a custom contact map formatted like this one can be loaded into WISP using the -contact_map_path parameter. "
        + "\n\n"
    )
    f.write(
        "functionalized_correlation_matrix.txt: A human readable representation of the functionalized correlation matrix, prior to multiplication by the contact map. If the user wishes to generate their own functionalized correlation matrix rather than letting WISP generate one automatically, a custom matrix formatted like this one can be loaded into WISP using the -functionalized_matrix_filename parameter. "
        + "\n\n"
    )
    f.write(
        "simply_formatted_paths.txt: A simple list of path lengths and nodes. The first column contains the lengths, and all following columns contain node indices. This file may be helpful for subsequent statistical analyses of the WISP output."
        + "\n"
    )
    f.close()
