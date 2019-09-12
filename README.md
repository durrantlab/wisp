README
======

WISP is a trajectory analysis tool that calculates and visualizes allosteric
pathways.

HOW TO INSTALL WISP (UNIX BASED OS: LINUX/MAC)
----------------------------------------------

1. Unzip the file to a directory of your choice using this command: `tar -xzf
   wisp.tgz <your directory here>`
2. To use the VMD plugin, add these two lines of code to your `.vmdrc` file
   (usually located in your `~/` directory)
   * `set auto_path "$auto_path <your directory here>"  ; #(NOTE: this may
     require the full pathname)`
   * `vmd_install_extension wisp wisp_tk_cb "Analysis/Wisp"`
3. Now open VMD: Click Extensions > Analysis > WISP
4. To use the command line, see the example in the `./example_commandline/`
   directory.

HOW TO USE WISP
---------------

`python wisp.py -pdb_trajectory_filename multi_frame_pdb.pdb -source_residues
"X_SER_1 X_LEU_4" -sink_residues X_ARG_37`
