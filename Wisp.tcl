#!/usr/bin/tcl

#Name:
# WISP GUI
#Synopsis:
# A GUI for the allosteric analysis program WISP
#Version
# 1.0
#

###################################
#package require Tcl 8.2
#package require struct::matrix 2.0.1

package provide wisp 1.11

set PI 3.14159265

variable w
namespace eval ::wisp:: {
  variable w
  variable mainwin
  variable sep "/"
  variable help_url "http://amarolab.ucsd.edu/wisp/"

  variable src_selection ""

  variable sink_selection ""
  variable load_into_vmd 1
  variable arglist {}
  array set descarray {} ;# a description of all arguments
  array set min_array {}
  array set max_array {}
  # wisp options & temporary variables

  variable tempdir "${sep}tmp${sep}"
  variable tempdir_temp "${sep}tmp${sep}"
  variable workdir "./"
  variable workdir_temp "$workdir"

  variable pdb_path ""
  variable source_residues "";# these three variables are dummy variables that don't actually hold anything, but keep the "other" menu option from appearing when it shouldn't
  variable sink_residues ""

  #if {![info exists wisp_directory]} {variable wisp_directory "./wisp.py"} ;# depending on whether the user has set the value in their .vmdrc
  variable wisp_directory_temp ""
  variable wisp_directory_format "file"
  variable output_dir ""
  variable output_dir_temp ""
  variable output_dir_format "directory"
  variable node_definition "CA" ;# RESIDUE_COM, SIDECHAIN_COM, BACKBONE_COM
  variable node_definition_temp ""
  variable node_definition_format "list {CA RESIDUE_COM SIDECHAIN_COM BACKBONE_COM}"
  variable contact_map_distance_limit "4.5"
  variable contact_map_distance_limit_temp ""
  variable load_wisp_saved_matrix FALSE
  variable load_wisp_saved_matrix_temp ""
  variable load_wisp_saved_matrix_format "list {TRUE FALSE}"
  variable wisp_saved_matrix_path "" ;# NOTE: need default here
  variable wisp_saved_matrix_path_temp ""
  variable wisp_saved_matrix_path_format "file"
  variable n_paths 20
  variable n_paths_temp ""
  variable n_cores 1
  variable n_cores_temp ""
  variable frame_chunks 96
  variable frame_chunks_temp ""
  variable shortest_path_radius 0.1
  variable shortest_path_radius_temp ""
  variable longest_path_radius 0.01
  variable longest_path_radius_temp ""
  variable spline_smoothness 0.01
  variable spline_smoothness_temp ""
  variable vmd_resolution 6
  variable vmd_resolution_temp ""
  variable node_sphere_radius 1.0
  variable node_sphere_radius_temp ""
  variable shortest_path_r 0.0; variable shortest_path_g 0.0; variable shortest_path_b 1.0;
  variable shortest_path_r_temp ""; variable shortest_path_g_temp ""; variable shortest_path_b_temp "";
  variable shortest_path_r_format "little_slider_begin 0.0 1.0"; variable shortest_path_g_format "little_slider_middle 0.0 1.0"; variable shortest_path_b_format "little_slider_end 0.0 1.0";
  variable longest_path_r 1.0; variable longest_path_g 0.0; variable longest_path_b 0.0;
  variable longest_path_r_temp ""; variable longest_path_g_temp ""; variable longest_path_b_temp "";
  variable longest_path_r_format "little_slider_begin 0.0 1.0"; variable longest_path_g_format "little_slider_middle 0.0 1.0"; variable longest_path_b_format "little_slider_end 0.0 1.0";
  variable node_sphere_r 1.0; variable node_sphere_g 0.0; variable node_sphere_b 0.0;
  variable node_sphere_r_temp ""; variable node_sphere_g_temp ""; variable node_sphere_b_temp "";
  variable node_sphere_r_format "little_slider_begin 0.0 1.0"; variable node_sphere_g_format "little_slider_middle 0.0 1.0"; variable node_sphere_b_format "little_slider_end 0.0 1.0";
  variable shortest_path_opacity 1.0; variable longest_path_opacity 1.0; variable node_sphere_opacity 1.0;
  variable shortest_path_opacity_temp ""; variable longest_path_opacity_temp ""; variable node_sphere_opacity_temp "";
  variable shortest_path_opacity_format "slider 0.0 1.0"; variable longest_path_opacity_format "slider 0.0 1.0"; variable node_sphere_opacity_format "slider 0.0 1.0";
  variable pdb_single_frame_path ""
  variable pdb_single_frame_path_temp ""
  variable pdb_single_frame_path_format "file"
  #variable simply_formatted_paths_path ""
  #variable simply_formatted_paths_path_temp ""
  variable seconds_to_wait_before_parallelizing_path_finding 5.0
  variable seconds_to_wait_before_parallelizing_path_finding_temp ""
  variable functionalized_matrix_path ""
  variable functionalized_matrix_path_temp ""
  variable functionalized_matrix_path_format "file"
  variable contact_map_path ""
  variable contact_map_path_temp ""
  variable contact_map_path_format "file"
  variable temp_pdb "wisp_temp.pdb"
  variable temp_pdb_temp ""
  variable path_contrast 1.0
  variable path_contrast_temp ""
  variable path_contrast_format "slider 1.0 30.0"
  variable num_paths 0 ;# variable will be modified once a WISP visualization is loaded

  variable arglist_files "wisp_directory output_dir"
  variable arglist_covariance "node_definition contact_map_distance_limit load_wisp_saved_matrix wisp_saved_matrix_path"
  variable arglist_pathsearching "n_paths"
  variable arglist_multiprocessor "n_cores frame_chunks"
  variable arglist_graphics "shortest_path_radius  longest_path_radius spline_smoothness vmd_resolution node_sphere_radius shortest_path_r shortest_path_g shortest_path_b longest_path_r longest_path_g longest_path_b node_sphere_r node_sphere_g node_sphere_b shortest_path_opacity longest_path_opacity node_sphere_opacity pdb_single_frame_path"
  variable arglist_advanced "seconds_to_wait_before_parallelizing_path_finding functionalized_matrix_path contact_map_path"

  variable arglist_other ""
  variable moltxt "(none)"
  variable which_mol 0
}
set ::wisp::auto_path $auto_path


proc ::wisp::UpdateMolecule {args} {
	# args is meaningless, but necessary
	set mainwin $::wisp::mainwin
	variable w
    	#variable moltxt
    	#variable molid
    	global vmd_molecule

    	# Update the molecule browser
    	set mollist [molinfo list]
    	$w.mollist.m configure -state disabled
        $w.mollist.m.menu delete 0 end
        set ::wisp::moltxt "(none)"
        #puts "Mainwin: $mainwin.mollist.m"
    	#if { [llength $mollist] > 0 } {
       # 	#$w.foot configure -state normal
       # 	$mainwin.mollist.m configure -state normal
       # 	#puts $mollist
       # 	foreach id $mollist {
       #     		$mainwin.mollist.m insert end "$id - [molinfo $id get name]"
       # 	}
#	}
    if { [llength $mollist] > 0 } {
        #$w.foot configure -state normal
        $w.mollist.m configure -state normal
        #puts $mollist
        foreach id $mollist {
            $w.mollist.m.menu add radiobutton -value $id \
                -command {global vmd_molecule ; if {[info exists vmd_molecule($::wisp::which_mol)]} {set ::wisp::moltxt "$::wisp::which_mol [molinfo $::wisp::which_mol get name]"} else {set ::wisp::moltxt "(none)" ; set ::wisp::which_mol -1} } \
                -label "$id [molinfo $id get name]" \
                -variable ::wisp::which_mol
            if {$id == $::wisp::which_mol} {
              puts "id: $id"
              puts "which_mol: $::wisp::which_mol"
              puts "exists: [info exists vmd_molecule($::wisp::which_mol)]"
                if {[info exists vmd_molecule($::wisp::which_mol)]} {
                    set ::wisp::moltxt "$::wisp::which_mol:[molinfo $::wisp::which_mol get name]"
                } else {
                    set ::wisp::moltxt "(none)"
                    set ::wisp::which_mol -1
                }
            }
        }
    }
}

proc ::wisp::wisp_mainwin {} {
  variable w
  if { [winfo exists .wisp] } {
    #destroy $::wisp::w
    wm deiconify $w
    #::wisp::UpdateMolecule $w
    return
  }
  set w [toplevel ".wisp"]
  wm title $w "WISP"
  wm resizable $w 0 0
  set ::wisp::mainwin $w

  wm protocol $w WM_DELETE_WINDOW {
    grab release $::wisp::mainwin
    after idle destroy $::wisp::mainwin
  }
	##
	## make menu bar
	##
  set ::wisp::w $w
  frame $w.menubar -relief raised -bd 2 ;# frame for menubar
  pack $w.menubar -padx 1 -fill x

  menubutton $w.menubar.help -text Help -underline 0 -menu $w.menubar.help.menu
  menubutton $w.menubar.settings -text Settings -underline 0 -menu $w.menubar.settings.menu
  menu $w.menubar.settings.menu -tearoff no
  $w.menubar.settings.menu add command -label "Files..." -command {::wisp::wisp_settings files}
  $w.menubar.settings.menu add command -label "Covariance Matrix..." -command {::wisp::wisp_settings covariance}
  #$w.menubar.settings.menu add command -label "Path Search..." -command {::wisp::wisp_settings pathsearching} ;# not currently needed
  $w.menubar.settings.menu add command -label "Multiprocessor..." -command {::wisp::wisp_settings multiprocessor}
  $w.menubar.settings.menu add command -label "Graphics..." -command {::wisp::wisp_settings graphics}
  $w.menubar.settings.menu add command -label "Advanced..." -command {::wisp::wisp_settings advanced}
  #puts "arglist_other: $::wisp::arglist_other"
  if {$::wisp::arglist_other != ""} { ;# if there are any subsequent settings, put them here
    $w.menubar.settings.menu add command -label "Other..." -command {::wisp::wisp_settings other}
  }

  menu $w.menubar.help.menu -tearoff no
  $w.menubar.help.menu add command -label "About.." -command {tk_messageBox -type ok -title "About WISP" -message "VMD GUI plugin for running allosteric network analyses using WISP.\n\nDeveloped in the:\n\tAmaro Laboratory\n\tUniversity of California, San Diego\n\thttp://amarolab.ucsd.edu\n\nDevelopers:\n\tAdam VanWart\n\tJacob Durrant\n\tLane Votapka\n\nPlease see README in the installation directory for additional information"}
  $w.menubar.help.menu add command -label "Help..." -command "vmd_open_url $::wisp::help_url"
  # XXX - set menubutton width to avoid truncation in OS X
  $w.menubar.help config -width 5
  $w.menubar.settings config -width 5

  pack $w.menubar.help -side right
  pack $w.menubar.settings -side left
  # Mol List
  frame $w.mollist
  pack [label $w.mollist.lbl -text "Select Molecule:"] -side top -anchor w
  #listbox $w.mollist.m -relief raised -selectmode single ;# molecule listbox
  #$w.mollist.m
  #scrollbar $w.mollist.scrollbar -command [list $w.mollist.m yview] ;# molecule listbox scrollbar
  #$w.mollist.m configure -yscrollcommand [list $w.mollist.scrollbar set]
  #pack $w.mollist.m -side left -fill x -expand 1  ;# pack the listbox
  #pack $w.mollist.scrollbar -side left -fill y  ;# pack the scrollbar
  #
  menubutton $w.mollist.m -relief raised -bd 2 -direction flush \
  	-textvariable ::wisp::moltxt \
	-menu $w.mollist.m.menu
  menu $w.mollist.m.menu -tearoff no
  pack $w.mollist.m -side left -fill x -expand yes
  pack $w.mollist -side top -fill x -padx 10 -pady 10
  #source/sink
  frame $w.src_sink
  pack [checkbutton $w.src_sink.load_into_vmd -text "Load into VMD when finished" -variable ::wisp::load_into_vmd] -side top
  pack [label $w.src_sink.lbl0 -text "Source selection:"] -side top
  pack [entry $w.src_sink.src_selection -textvariable ::wisp::src_selection] -side top
  pack [label $w.src_sink.lbl1 -text "Sink selection:"] -side top
  pack [entry $w.src_sink.sink_selection -textvariable ::wisp::sink_selection] -side top
  pack $w.src_sink -side top -anchor n -padx 10 -pady 10
  # num_paths
  frame $w.num_paths -relief groove -bd 2
  pack [label $w.num_paths.lbl1 -text "Desired Number of Paths:"] -side top
  pack [entry $w.num_paths.entry -textvariable ::wisp::n_paths] -side top
  pack $w.num_paths -side top -anchor n -padx 10 -pady 10
  # buttons
  frame $w.buttons
  button $w.buttons.run -text "Run WISP" -width 10 -command ::wisp::run_wisp
  pack $w.buttons.run -side left
  pack $w.buttons -side top -anchor n -padx 10 -pady 10
  ::wisp::UpdateMolecule $w
  global vmd_molecule
  trace variable vmd_molecule w ::wisp::UpdateMolecule;
}


proc ::wisp::run_wisp {} {
  #set selected_mol_index [$::wisp::mainwin.mollist.m curselection]
  if {$::wisp::which_mol == -1} {
    tk_messageBox -type ok -title "Error" -message "Error: Please select a molecule."; return
  }
  if {[file exists $::wisp::wisp_directory] == 0} {	;# wisp is not found
	#give an error message
	tk_messageBox -type ok -title "Wisp not found" -message "Error: The Wisp program location does not exist at the specified location: $::wisp::wisp_directory. Go to Settings > Files to change location"
	return
  }
  if {[file readable $::wisp::wisp_directory] == 0} {	;# wisp is not readable
	#give an error message
	tk_messageBox -type ok -title "Wisp not readable" -message "Error: The specified Wisp program does not have read permissions."
	return
  }
  set selected_mol $::wisp::which_mol ;#[lindex [molinfo list] $::wisp::which_mol]
  set pdb_struct [atomselect $selected_mol all]
  set src_response [catch {set src_sel [atomselect $selected_mol $::wisp::src_selection]}]
  if {$src_response || [$src_sel get index]==""} {
    tk_messageBox -type ok -title "Error" -message "Error: the Source selection string is not properly formatted: it must be formatted like a VDM atom selection.\n Example: \"resid 105 and chain B\""
    return
  }
  set sink_response [catch {set sink_sel [atomselect $selected_mol $::wisp::sink_selection]}]
  if { $sink_response || [$sink_sel get index]==""} {
    tk_messageBox -type ok -title "Error" -message "Error: the Sink selection string is not properly formatted: it must be formatted like a VDM atom selection.\n Example: \"resid 105 and chain B\""
    return
  }
  set src_str "[lindex [$src_sel get chain] 0]_[lindex [$src_sel get resname] 0]_[lindex [$src_sel get resid] 0]" ;# NOTE: incorrect! need to be a list
  set sink_str "[lindex [$sink_sel get chain] 0]_[lindex [$sink_sel get resname] 0]_[lindex [$sink_sel get resid] 0]"
  set other_args [::wisp::get_other_args]
  puts "Now writing pdb trajectory..."
  set dialog [toplevel .wait]
  pack [frame $dialog.frame -borderwidth 2 -relief raised] -fill both
  pack [message $dialog.frame.msg -text "Please wait\nWISP is running..." -width 100] -padx 100 -pady 100
  update
  grab $dialog
  animate write pdb $::wisp::temp_pdb sel $pdb_struct $selected_mol  ;# this writes the pdb trajectory explicitly

  set pdb_path $::wisp::temp_pdb ;#[molinfo $selected_mol get filename]
  set timename [clock format [clock seconds] -format "wisp_output__%h_%d_%Y__%H_%M_%p"]
  if {$::wisp::output_dir == ""} {
    set outdir $timename
  } else {
    set outdir [file join $::wisp::output_dir $timename]
  }
  set curdir [pwd]
  cd $::wisp::workdir
  set command "python $::wisp::wisp_directory -pdb_path $pdb_path -source_residues \"$src_str\" -sink_residues \"$sink_str\" -output_dir $outdir $other_args"
  puts "running command: $command"
  update

  set error [catch {exec {*}$command} results options]

  destroy $dialog
  puts "WISP std. output: $results"
  if {$error} {tk_messageBox -type ok -title "Failure" -message "Alert: WISP failed to run properly. See VMD standard output for more detailed information."; return}
  set error [catch {exec rm $::wisp::temp_pdb}] ;# delete the temporary trajectory
  if {$::wisp::load_into_vmd == 1} {
    cd $outdir ;# cd to the directory where everything is
    source "visualize.tcl"
    if {[info exists wisp_num_paths]} {set ::wisp::num_paths $wisp_num_paths}
  } ;# load the wisp tcl file for visualization
  cd $curdir
  ::wisp::UpdateMolecule
}

proc ::wisp::parse_helpfile {{help_command "-help"}} {
  set full_help_command "python $::wisp::wisp_directory $help_command"
  set helpstring [exec {*}$full_help_command]
  set helplist [split $helpstring "\n"]
  set ::wisp::desclist {}
  set append_desc False ;# whether we are appending the description contents to a particular argument
  set ::wisp::arglist {}
  set descvar ""
  set arg ""
  foreach line $helplist {
    set firstchar [lindex [split $line ""] 0] ;# get the first character in a line
    if {($firstchar == " ") || ($firstchar == "\t")} { ;# then the first character is a whitespace
      if {$append_desc == True} { ;# then we are appending this to the description list
        set descvar "${descvar}[string range $line 3 end]"
      }
      continue
    } elseif {([scan $firstchar %c] >= 97) && ([scan $firstchar %c] <= 122)} { ;# if the line begins with a lowercase letter
      if {$arg != ""} {set ::wisp::descarray($arg) $descvar}
      set arg [lindex [split $line ":"] 0]
      lappend ::wisp::arglist $arg
      #puts "descvar: $descvar, arg: $arg"
      set append_desc True
      set descvar [lindex [split $line ":"] 1]
    } else {
      set append_desc False
    }
  }
  #puts $::wisp::arglist
  #puts $::wisp::desclist
  foreach arg $::wisp::arglist {
    if {!([info exists "::wisp::$arg"])} { ;# if the variable doesn't already exist, create it
      set ::wisp::$arg ""
      lappend ::wisp::arglist_other $arg
    }
  }

}

proc ::wisp::find_py {} {
  # attempts to find the wisp.py python script
  if {[info exists ::wisp::wisp_directory] && [file exists $::wisp::wisp_directory]} {return True} ;# search the defined location
  foreach path $::wisp::auto_path {
    if {[file exists [file join $path wisp wisp.py]]} { ;# then we've found it
      set ::wisp::wisp_directory [file join $path wisp wisp.py]
      return True
    } elseif {[file exists [file join $path wisp.py]]} { ;# then we've found it
      set ::wisp::wisp_directory [file join $path wisp.py]
      return True
    }
  }
  set ::wisp::wisp_directory "./wisp.py" ;# then we couldn't find it
  return False

}

proc ::wisp::default_args {} { ;# in case there is some sort of problem with parsing the helpfile, this will assign the default arguments
  set ::wisp::arglist {}
  #lappend ::wisp::arglist output_dir
  lappend ::wisp::arglist node_definition
  lappend ::wisp::arglist contact_map_distance_limit
  lappend ::wisp::arglist load_wisp_saved_matrix
  lappend ::wisp::arglist wisp_saved_matrix_path
  lappend ::wisp::arglist n_paths
  lappend ::wisp::arglist n_cores
  lappend ::wisp::arglist frame_chunks
  lappend ::wisp::arglist shortest_path_radius
  lappend ::wisp::arglist longest_path_radius
  lappend ::wisp::arglist spline_smoothness
  lappend ::wisp::arglist vmd_resolution
  lappend ::wisp::arglist node_sphere_radius
  lappend ::wisp::arglist shortest_path_r
  lappend ::wisp::arglist shortest_path_g
  lappend ::wisp::arglist shortest_path_b
  lappend ::wisp::arglist longest_path_r
  lappend ::wisp::arglist longest_path_g
  lappend ::wisp::arglist longest_path_b
  lappend ::wisp::arglist node_sphere_r
  lappend ::wisp::arglist node_sphere_g
  lappend ::wisp::arglist node_sphere_b
  lappend ::wisp::arglist shortest_path_opacity
  lappend ::wisp::arglist longest_path_opacity
  lappend ::wisp::arglist node_sphere_opacity
  lappend ::wisp::arglist pdb_single_frame_path
  lappend ::wisp::arglist seconds_to_wait_before_parallelizing_path_finding
  lappend ::wisp::arglist functionalized_matrix_path
  lappend ::wisp::arglist contact_map_path

}

proc ::wisp::get_other_args {} {
  set argstr ""
  foreach arg $::wisp::arglist {
    if {$arg == "pdb_path" || $arg == "output_dir"} {continue}
    eval "set tmpvar \$::wisp::$arg"
    if {$tmpvar != {}} {set argstr [eval "concat $argstr \"-$arg \$::wisp::$arg\""]}
  }
  #puts $argstr
  return $argstr
}

proc ::wisp::wisp_settings {{mode graphics}} { ;# mode can be graphics, advanced
	if {$mode == "files"} {
	  set arglist $::wisp::arglist_files
	  set window_title "WISP Files Settings"
	} elseif { $mode == "covariance" } {
	  set arglist $::wisp::arglist_covariance
	  set window_title "WISP Covariance Matrix Settings"
	} elseif { $mode == "pathsearching" } {
	  set arglist $::wisp::arglist_pathsearching
	  set window_title "WISP Path-Searching Settings"
	} elseif { $mode == "multiprocessor" } {
	  set arglist $::wisp::arglist_multiprocessor
	  set window_title "WISP Multiprocessor Settings"
	} elseif { $mode == "graphics" } {
	  set arglist $::wisp::arglist_graphics
	  set window_title "WISP Graphics Settings"
	} elseif { $mode == "advanced" } {
	  set arglist $::wisp::arglist_advanced
	  set window_title "WISP Advanced Settings"
	} elseif { $mode == "other"} {
	  set arglist $::wisp::arglist_other
	  set window_title "WISP Other Settings"
	}

	variable w
	variable settings_win
	set win_name $w.settings_$mode
	if { [winfo exists $win_name] } {
		wm deiconify $settings_win
		return
	}
	set settings_win [toplevel "$w.settings_$mode"]
	wm transient $settings_win $w
	wm protocol $settings_win WM_DELETE_WINDOW {
		grab release $::wisp::settings_win
		after idle destroy $::wisp::settings_win
	}
	wm title $settings_win "$window_title"
	wm resizable $settings_win no no

	# temp variables
	#set ::wisp::workdir_temp $::wisp::workdir
	#set ::wisp::wispdir_temp $::wisp::wispdir

	# menu widgets


	foreach arg $arglist {
          eval "set ::wisp::${arg}_temp \$::wisp::$arg" ;# set all temporary argument variables
        }
	#puts "shortest_path_opacity: $::wisp::shortest_path_opacity_temp"
	set len_arglist [llength $arglist] ;# get the length of the argument list
	set half_len_arglist [expr "([llength $arglist] / 2) + 1"] ;# get half the length of the argument list for the two columns in the settings window

	frame $settings_win.leftcol ;# first do the left column
	frame $settings_win.rightcol


	#frame $settings_win.leftcol.temp_dir
	#pack [label $settings_win.leftcol.temp_dir.caption -text "Working Directory"] -side top -anchor w ; balloon $settings_win.leftcol.temp_dir.caption "description"
	#pack [entry $settings_win.leftcol.temp_dir.textbox -width 20 -textvariable ::wisp::workdir_temp] -side top -anchor w
	#pack $settings_win.leftcol.temp_dir -side top -anchor w -pady 10
#
	#frame $settings_win.leftcol.wisp_dir
	#pack [label $settings_win.leftcol.wisp_dir.caption -text "WISP location"] -side top -anchor w
	#pack [entry $settings_win.leftcol.wisp_dir.textbox -width 20 -textvariable ::wisp::wispdir_temp] -side top -anchor w
	#pack $settings_win.leftcol.wisp_dir -side top -anchor w -pady 10

	#puts "len_arglist: $len_arglist"

	for {set i 0} {$i < $len_arglist} {incr i} {
	  set arg [lindex $arglist $i]
	  if {([info exists "::wisp::${arg}_format"])} {
	    eval "set arg_format_list \$::wisp::${arg}_format"
	  } else {
	    set arg_format_list "normal"
	  }
	  set arg_format [lindex $arg_format_list 0] ;# get the first element of the format list, because the rest of the stuff are parameters
	  if {[llength $arg_format_list] > 2} { ;# then we have more parameters for this argument
	    set arg_min [lindex $arg_format_list 1] ;# the min value this argument can take
	    set ::wisp::min_array($arg) $arg_min ;# set the array that holds minimum allowed values
	    set arg_max [lindex $arg_format_list 2] ;# the max value this argument can take
	    set ::wisp::max_array($arg) $arg_max ;# set the array that holds maximum allowed values
	  } elseif {[llength $arg_format_list] == 2} { ;# then its a list
	    set arg_list_options [lindex $arg_format_list 1]
	  }
	  if {$arg == "pdb_path" || $arg == "source_residues" || $arg == "sink_residues" } { ;# or many other things...
	    continue ;# skip it, we already have it covered in the main menu
	  }
	  #if {($i < $half_len_arglist || $len_arglist < 6) } { ;# then we are on the left side
	    set column $settings_win.leftcol
	  #}
	  if {$i >= $half_len_arglist && $len_arglist > 6} { ;# then we are on the right side
	    if {$arg_format != "little_slider_middle"} {
	      if {$i > $half_len_arglist} {set column $settings_win.rightcol} ;# all this is very messy, but I'm trying to keep the little_sliders (which require the same frame) together on the same column
	    } elseif {$arg_format != "little_slider_end"} {
	      if {$i > [expr "$half_len_arglist + 1"]} {set column $settings_win.rightcol}
	    } else {
	      set column $settings_win.rightcol ;# if I'm not dealing with little_sliders, then no need to worry
	    }
	  }
	  set our_widget "$column.$arg"
	  #puts "now creating widget: $our_widget"
	  if {$arg_format == "little_slider_begin"}  {
	    set last_str_index [expr "[string length $arg] - 3"] ;# get the index of the third to last character in the string for the title
	    set widget_name [string range $arg 0 $last_str_index] ;# need a generic name for three subsequent widgets
	    set our_widget $column.$widget_name
	    frame $our_widget ;# need to create the frame for the title and entries
	    pack [label $our_widget.caption -text "$widget_name (RGB)"] -side top -anchor w ; balloon $our_widget.caption ::wisp::descarray($arg) ;# title
	    frame $our_widget.subframe ;# a subframe for the three entries
	    #pack [entry $our_widget.subframe.textbox_begin -width 6 -textvariable ::wisp::${arg}_temp] -side left -anchor w; balloon $our_widget.subframe.textbox_begin ::wisp::descarray($arg)
	    #::wisp::make_scale "$our_widget.subframe.slider_begin" $arg $arg_max $arg_min
	    eval "scale $our_widget.subframe.slider_begin -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 50 -sliderlength 20 -showvalue true -resolution 0.01 -command {set ::wisp::${arg}_temp }"
	    eval "$our_widget.subframe.slider_begin set \$::wisp::${arg}_temp"
	    pack $our_widget.subframe.slider_begin -side left -anchor w ; balloon $our_widget.subframe.slider_begin ::wisp::descarray($arg)
	  } elseif {$arg_format == "little_slider_middle"} { ;# for the middle field in a set of three small entries
	    set last_str_index [expr "[string length $arg] - 3"] ;# get the index of the third to last character in the string for the title
	    set widget_name [string range $arg 0 $last_str_index]
	    set our_widget $column.$widget_name
	    #pack [entry $our_widget.subframe.textbox_middle -width 6 -textvariable ::wisp::${arg}_temp] -side left -anchor w; balloon $our_widget.subframe.textbox_middle ::wisp::descarray($arg)
	    eval "scale $our_widget.subframe.slider_middle -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 50 -sliderlength 20 -showvalue true -resolution 0.01 -command {set ::wisp::${arg}_temp }"
	    eval "$our_widget.subframe.slider_middle set \$::wisp::${arg}_temp"
	    pack $our_widget.subframe.slider_middle -side left -anchor w ; balloon $our_widget.subframe.slider_middle ::wisp::descarray($arg)
	  } elseif {$arg_format == "little_slider_end"} { ;# for the last field in a set of three small entries
	    set last_str_index [expr "[string length $arg] - 3"] ;# get the index of the third to last character in the string for the title
	    set widget_name [string range $arg 0 $last_str_index]
	    set our_widget $column.$widget_name
	    #pack [entry $our_widget.subframe.textbox_end -width 6 -textvariable ::wisp::${arg}_temp] -side left -anchor w; balloon $our_widget.subframe.textbox_end ::wisp::descarray($arg)
	    eval "scale $our_widget.subframe.slider_end -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 50 -sliderlength 20 -showvalue true -resolution 0.01 -command {set ::wisp::${arg}_temp }"
	    eval "$our_widget.subframe.slider_end set \$::wisp::${arg}_temp"
	    pack $our_widget.subframe.slider_end -side left -anchor w ; balloon $our_widget.subframe.slider_end ::wisp::descarray($arg)
	    pack $our_widget.subframe -side top -anchor w -fill x -expand 1 ;# close the entry subframe
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 ;# close the entire frame
	  } elseif {$arg_format == "slider"} { ;# a slider widget
	    frame $our_widget ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg"] -side top -anchor w ; balloon $our_widget.caption ::wisp::descarray($arg)
	    frame $our_widget.cutoffslider
	    #pack [entry $our_widget.cutoffslider.textbox -width 6 -textvariable ::wisp::${arg}_temp] -side left -fill x

	    if {${arg} == "node_sphere_opacity"} { ;# this means that this argument is node_sphere opacity, and that a structure is loaded
	      eval "scale $our_widget.cutoffslider.slider -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -showvalue true -resolution 0.01 -command {if {[lsearch [material list] node_spheres] != -1} {material change opacity node_spheres \$::wisp::${arg}_temp}; set ::wisp::${arg}_temp }"
	    } elseif {${arg} == "shortest_path_opacity" || $arg == "longest_path_opacity"} {
	      ::wisp::make_scale $our_widget.cutoffslider.slider $arg $arg_max $arg_min

	    } else { ;# otherwise, just treat it normally
	      eval "scale $our_widget.cutoffslider.slider -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -showvalue true -resolution 0.01 -command {set ::wisp::${arg}_temp }"
	    }
	    eval "$our_widget.cutoffslider.slider set \$::wisp::${arg}_temp"
	    pack $our_widget.cutoffslider.slider -side left -fill x -expand 1
	    pack $our_widget.cutoffslider -side left -fill x -expand 1
	    pack $our_widget -side top -anchor w -pady 10 -padx 10
	  } elseif {$arg_format == "file" || $arg_format == "directory"} {
	    if {$arg_format == "file"} { set which_tk_getopen "tk_getOpenFile"
	    } elseif {$arg_format == "directory"} { set which_tk_getopen "tk_chooseDirectory" } ;# two different things the user might choose: files or directories
	    frame $our_widget -relief groove -bd 2  ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg"] -side top -anchor w ; balloon $our_widget.caption ::wisp::descarray($arg)
	    eval "pack \[entry $our_widget.textbox -width 25 -textvariable ::wisp::${arg}_temp -validate all -validatecommand {::wisp::isFile %P $our_widget.textbox}\] -side top -anchor w; balloon $our_widget.textbox ::wisp::descarray($arg)" ;# NOTE: will need to make an array in place of desclist
	    eval "pack \[button $our_widget.browse -width 10 -text \"Browse...\" -command {set ::wisp::${arg}_temp \[$which_tk_getopen\]}\] -side top -anchor w; balloon $our_widget.browse ::wisp::descarray($arg)" ;# NOTE: will need to make an array in place of desclist
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 -fill x
	  } elseif {$arg_format == "list"} {
	    frame $our_widget -relief groove -bd 2 ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg"] -side top -anchor w ; balloon $our_widget.caption ::wisp::descarray($arg)
	    #pack [entry $our_widget.textbox -width 20 -textvariable ::wisp::${arg}_temp] -side top -anchor w; balloon $our_widget.textbox ::wisp::descarray($arg);# NOTE: will need to make an array in place of desclist
	    foreach option $arg_list_options {
	      pack [radiobutton $our_widget.option_$option -text "$option" -value $option -variable ::wisp::${arg}_temp] -side top -anchor w ;
	    }
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 -fill x
	  } else {
	    frame $our_widget ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg"] -side top -anchor w ; balloon $our_widget.caption ::wisp::descarray($arg)
	    pack [entry $our_widget.textbox -width 20 -textvariable ::wisp::${arg}_temp] -side top -anchor w; balloon $our_widget.textbox ::wisp::descarray($arg);# NOTE: will need to make an array in place of desclist
	    pack $our_widget -side top -anchor w -pady 10 -padx 10
	  }
	}

	if { $mode == "graphics" } { ;# then place the contrast bar
	  set our_widget "$column.path_contrast"
	  set arg_max [lindex $::wisp::path_contrast_format 2]
	  set arg_min [lindex $::wisp::path_contrast_format 1]
	  frame $our_widget ;# balloon $our_widget "Test Help"
	  pack [label $our_widget.caption -text "path contrast"] -side top -anchor w ; balloon $our_widget.caption "Adjusts the contrast between paths of varying lengths"
	  frame $our_widget.cutoffslider
	  eval "scale $our_widget.cutoffslider.slider -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -showvalue true -resolution 0.1 -command {::wisp::adjust_path_contrast }"
	  # balloon?
	  $our_widget.cutoffslider.slider set $::wisp::path_contrast
	  pack $our_widget.cutoffslider.slider -side left -fill x -expand 1
	  pack $our_widget.cutoffslider -side left -fill x -expand 1
	  pack $our_widget -side top -anchor w -pady 10 -padx 10
	}

	grid $settings_win.leftcol -column 0 -row 0
	grid $settings_win.rightcol -column 1 -row 0

#frame $settings_win.leftcol.checkboxes
	#pack [checkbutton $settings_win.leftcol.checkboxes.sovereignty -text "No intrastructural clustering" -variable ::wisp::sovereignty_temp -onvalue 1 -offvalue 0] -side top -anchor w
	#pack [checkbutton $settings_win.leftcol.checkboxes.hydrophobic -text "Hydrophobic calculations" -variable ::wisp::find_hydrophobic_temp -onvalue 1 -offvalue 0] -side top -anchor w

	set arg shortest_path_b

	frame $settings_win.okaycancel
	# need to use an eval so that the $arglist doesn't have to be global
	eval "button \$settings_win.okaycancel.okay -text OK -width 6 \
	  -command {;\
	  	foreach arg [list $arglist] {\
	  	  eval \"set value \\\$::wisp::\${arg}_temp\";\
	  	  if {[info exists ::wisp::min_array(\$arg)] && [info exists ::wisp::max_array(\$arg)]} {\
	  	    if {\$value < \$::wisp::min_array(\$arg)} { \
	  	      tk_messageBox -type ok -title \"Error\" -message \"Error: the value you entered for \$arg is too low. It must be greater than \$::wisp::min_array(\$arg).\";\
	  	      return;\
	  	    } ;\
	  	    if {\$value > \$::wisp::max_array(\$arg)} {\
	  	      tk_messageBox -type ok -title \"Error\" -message \"Error: the value you entered for \$arg is too high. It must be less than \$::wisp::max_array(\$arg).\";\
	  	      return;\
	  	    };\
	  	  };\
     		  set ::wisp::\$arg \$value;\
     		} ;\
	  	grab release $::wisp::settings_win;\
	  	after idle destroy $::wisp::settings_win;\
	}"
	button $settings_win.okaycancel.cancel -text Cancel -width 6 \
	  -command {
	  	grab release $::wisp::settings_win
	  	after idle destroy $::wisp::settings_win
	}
	#button $settings_win.okaycancel.default -text Default -width 6 \
	#  -command {
	#  	set ::wisp::workdir_temp $::wisp::workdir
	#  	#...
	#
	#}
	pack $settings_win.okaycancel.okay $settings_win.okaycancel.cancel -side left ;# $settings_win.okaycancel.default
	grid $settings_win.okaycancel -column 0 -row 1

}

proc ::wisp::make_scale {name arg arg_max arg_min} {
  eval "scale $name -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -sliderlength 40 -showvalue true -resolution 0.01 -command {::wisp::adjust_wispmaterial ::wisp::shortest_path_opacity_temp ::wisp::longest_path_opacity_temp ::wisp::shortest_path_r ::wisp::longest_path_r ::wisp::num_paths; set ::wisp::${arg}_temp }"
  #eval "$name set \$::wisp::${arg}_temp"

}

proc ::wisp::adjust_path_contrast {value} {
 set ::wisp::path_contrast $value
 ::wisp::adjust_wispmaterial ::wisp::shortest_path_opacity_temp ::wisp::longest_path_opacity_temp ::wisp::shortest_path_r ::wisp::longest_path_r ::wisp::num_paths
 }

proc ::wisp::adjust_wispmaterial {shortest_path_opacity_var longest_path_opacity_var shortest_path_color_var longest_path_color_var num_paths_var} {
  set shortest_path_opacity [expr "\$$shortest_path_opacity_var"]
  set longest_path_opacity [expr "\$$longest_path_opacity_var"]
  set shortest_path_color [expr "\$$shortest_path_color_var"]
  set longest_path_color [expr "\$$longest_path_color_var"]
  set num_paths [expr "\$$num_paths_var"] ;# get the current value of the num_paths variable
  #puts "$shortest_path_opacity $longest_path_opacity $shortest_path_color $longest_path_color $num_paths"
  if {$num_paths == 0} {return} ;# if there are no paths, then get outta here
  set increment [expr "($longest_path_opacity - $shortest_path_opacity)/$num_paths"]


  foreach mat [material list] { ;# for all the materials in vmd
    #puts "mat: $mat"
    if {[string first wisp_material $mat] == 0} { ;# then this material is one we need to modify
      #puts "modifying material: $mat"
      set index [string range $mat [string length wisp_material] end] ;# this is the index at the end of the material name
      set opacity [expr "($shortest_path_opacity + $increment*$index)**($::wisp::path_contrast**2)"]
      material change opacity $mat $opacity
      # then do something about the color (???)
    }
  }
}

proc ::wisp::isFile {f w} {
    if { [file exists $f] && ([file isfile $f] || [file isdirectory $f]) } {
        $w configure -fg black
    } else {
        $w configure -fg red
    }
    return 1;
}

proc ::wisp::wisp_init {} {
  ;# the correct OS directory separator must be defined
  switch [vmdinfo arch] {
    WIN64 -
    WIN32 {
      set ::wisp::sep "\\"
      set ::wisp:tempdir ""
    } default {
      set ::wisp::sep "/"
    }
  }
  #set ::wisp::wisp_directory "[file dirname $argv0]${sep}wisp.py"
  if { [::wisp::find_py] == True } { ::wisp::parse_helpfile } ;# search for the python script
  ::wisp::wisp_mainwin ;# call the main window
}

### Balloon help ###
proc balloon {w help} {
    if {![info exists $help]} {return}
    eval "bind $w <Any-Enter> \"after 1000 [list balloon:show %W [list \$$help]]\""
    bind $w <Any-Leave> "destroy %W.balloon"
}
proc balloon:show {w arg} {
    if {[eval winfo containing  [winfo pointerxy .]]!=$w} {return}
    set top $w.balloon
    catch {destroy $top}
    toplevel $top -bd 1 -bg black
    wm overrideredirect $top 1
    if {[string equal [tk windowingsystem] aqua]}  {
        ::tk::unsupported::MacWindowStyle style $top help none
    }
    pack [message $top.txt -aspect 200 -bg lightyellow \
            -font fixed -text $arg]
    set wmx [winfo rootx $w]
    set wmy [expr [winfo rooty $w]+[winfo height $w]]
    wm geometry $top \
      [winfo reqwidth $top.txt]x[winfo reqheight $top.txt]+$wmx+$wmy
    raise $top
}


#this gets called by VMD the first time the menu is opened
proc wisp_tk_cb {} {
  ::wisp::wisp_init	;# start the tool
  return $::wisp::w
}
#::wisp::parse_helpfile
