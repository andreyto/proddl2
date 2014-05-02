
################################################################################
#
# batch thermodynamics calculation routines 
#
################################################################################

################################################################################
#
# display the batch thermodynamics screen
#
proc show_batch_screen {} {
	global Frame panvar stcVar

	set_val panvar(batch_file) . batch_file $stcVar(OUT_DIR)/stc.batch
	set_val stcVar(RUN_FILE) . run_file "$stcVar(OUT_DIR)/$stcVar(TMP).thermo.run"
	set f $Frame(batch)
	toplevel $f 
	init_system_type

	set labelLen 24
	set padx 10
	set pady 3

	# define title 
	label $f.title -text "\nCalculating Thermodynamics on Multiple PDB Files\n"
	pack $f.title

	# define file entry fields
	set panvar(pdbroot) $stcVar(OUT_DIR) 
	u_file_entry $f.input      "PDB directory:   " \
		$labelLen panvar(pdbroot) $panvar(file_width) $padx {} 2 0
	pack $f.input -pady 2

	u_file_entry $f.batchfile "Summary Results File:" \
		$labelLen panvar(batch_file) $panvar(file_width) $padx {} 0 0
	pack $f.batchfile -pady 3

	u_file_entry $f.batchvars "Display Variables:" \
		$labelLen panvar(batch_vars) 30 $padx {} 0 0
	pack $f.batchvars -pady 3

	u_file_entry $f.temp "Temperature:" \
		$labelLen panvar(temperature) 10 $padx {} 0 0
	pack $f.temp -pady 3

	pack [frame $f.nop1] -pady 10

	set panvar(sbb) 0
	pack [frame $f.sbb] -fill x -expand true -padx 1
	pack [frame $f.sbb.a] -side left
	pack [frame $f.sbb.b] -side left
	label $f.sbb.a.label -text "Sbb:   " -anchor e -width $labelLen 
	pack $f.sbb.a.label -side right
	radiobutton $f.sbb.b.zero -text "Set zero" \
		-variable panvar(sbb) -value 0
	radiobutton $f.sbb.b.data -text "Calculate from data files" \
		-variable panvar(sbb) -value 1
	pack $f.sbb.b.zero $f.sbb.b.data -padx 10 -pady 0 -side left

	set panvar(sexu) 0
	pack [frame $f.sexu] -fill x -expand true -padx 1
	pack [frame $f.sexu.a] -side left
	pack [frame $f.sexu.b] -side left
	label $f.sexu.a.label -text "Sex->u:" -anchor e -width $labelLen 
	pack $f.sexu.a.label -side right
	radiobutton $f.sexu.b.zero -text "Set zero" \
		-variable panvar(sexu) -value 0
	radiobutton $f.sexu.b.data -text "Calculate from data files" \
		-variable panvar(sexu) -value 1
	radiobutton $f.sexu.b.calc -text "Calculate from tables" \
		-variable panvar(sexu) -value 2
	pack $f.sexu.b.zero $f.sexu.b.data $f.sexu.b.calc \
		-padx 10 -pady 0 -side left

	pack [frame $f.nop2] -pady 10

	message $f.msg -text " " -width 500
	pack $f.msg
	show_pdb $f


	# calculate/exit buttons
	frame $f.buts
	pack $f.buts -pady 20 
	button $f.buts.calc -text "Calculate" -command "do_batch $f" 
	button $f.buts.exit -text "Exit" -command "destroy $f"
	pack $f.buts.calc $f.buts.exit -side left
}

proc show_pdb {f} {
	global panvar stcVar

	if [info exists panvar(pdbroot)] {
		if [file isdirectory $panvar(pdbroot)] {
			set stcVar(pdbfiles) [glob -nocomplain $panvar(pdbroot)/*.pdb]
			set filenum [llength $stcVar(pdbfiles)]
			$f.msg configure -text "Number of pdb files = $filenum"
			return
		}
	}
	$f.msg configure -text "Specify directory of pdb files"
}

proc do_batch { f } {
	global stcVar panvar

	set stcVar(pdbfiles) [glob -nocomplain $panvar(pdbroot)/*.pdb]
	set filenum [llength $stcVar(pdbfiles)]
	$f.msg configure -text "Number of pdb files = $filenum"

	if {$panvar(sys_mode) == "P_UNFOLD"} {
		error "Batch mode not implemented for unfolding"
	}

	if ![winfo exists $f.fr0] {
		frame $f.fr0
		pack $f.fr0 -pady 20 -padx 20
		set stcVar(resList) [listbox $f.fr0.listing \
			-setgrid false -selectmode single \
			-xscrollcommand "$f.fr0.xscroll set" \
			-yscrollcommand "$f.fr0.yscroll set" \
			-relief raised -bd 2 -width 70 -height 8]
		scrollbar $f.fr0.xscroll -orient horizontal -command "$f.fr0.listing xview"
		scrollbar $f.fr0.yscroll -command "$f.fr0.listing yview"
		pack $f.fr0.xscroll -side bottom -fill x
		pack $f.fr0.yscroll -side right -fill y
		pack $f.fr0.listing -side left -fill both -expand true
	}
	$stcVar(resList) delete 0 end
	$stcVar(resList) insert end "Results from '$panvar(batch_file)'"
	bind $stcVar(resList)  <ButtonRelease-1> {
		set root [lindex [selection get] 0]
		if { 0 != [string compare $root "Results"] } {
			show_thermo_results $root.basic $root.detail "a"
		}
	}
	update

	if [catch {open $panvar(batch_file) w} stcVar(batchfp)] {
		puts stdout "Cannot open batch results file '$panvar(batch_file)'."
		exit
	}
	foreach pdb $stcVar(pdbfiles) {
		set root [file tail [file rootname $pdb]]
		calc_asa_batch $f $pdb $root
	}
	check_batch_errors $f.msg

}


proc calc_asa_batch { f pdb root } {
	global stcVar panvar

	fix_pdb 2 $pdb $root.complex $root.ligand $root.enzyme
	calc_asa $root $root.complex $root.complex.acc
	calc_asa $root $root.ligand $root.ligand.acc
	calc_asa $root $root.enzyme $root.enzyme.acc
	calc_thermo_batch $f $root
}

################################################################################
#
# user clicks on the calculate button
#
proc calc_thermo_batch { f root } {
	global stcVar panvar

	if [catch {set fp [open $stcVar(RUN_FILE) w]} err] {
		set stcVar(errFlag) 1
		error "Cannot open file: $stcVar(RUN_FILE)"
	}

	if {$panvar(sys_mode) == "P_BIND"} {
		puts $fp "1"
	} elseif {$panvar(sys_mode) == "P_UNFOLD"} {
		puts $fp "2"
	} else { puts $fp "3" }

	puts $fp "$root.basic"
	puts $fp "$root.detail"
	puts $fp "$panvar(amino_table)"
	puts $fp "$root.complex.acc"
	puts $fp "$root.ligand.acc"
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		puts $fp "$root.enzyme.acc"
	}

	if {$panvar(sconf_free1) == ""} {
		puts $fp "n"
	} else {
		puts $fp "y"
		puts $fp "$panvar(sconf_free1)"
	}

	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		if {$panvar(sconf_free2) == ""} {
			puts $fp "n"
		} else {
			puts $fp "y"
			puts $fp "$panvar(sconf_free2)"
		}
	}

	if {$panvar(res_unf) == ""} {
		puts $fp "n"
	} else {
		puts $fp "y"
		puts $fp "$panvar(res_unf)"
	}

	if {$panvar(sbb) != 0} {
		puts $fp "3"
		puts $fp $panvar(sbb)
	}

	if {$panvar(sexu) != 0} {
		puts $fp "4"
		puts $fp $panvar(sexu)
	}

	puts $fp "5"
	puts $fp "0"
	close $fp

	set panvar(thermo_sconf) [string trim $panvar(thermo_sconf)]
	calc_thermo $panvar(thermo_sconf) $stcVar(RUN_FILE) $root.thermo.log

	set line "$root"
	foreach item $panvar(batch_vars) {
		set unixcmd "grep $item $root.basic"
		set ans [eval exec $unixcmd]
		set fans [format " : %s %4g %s " [lindex $ans 0] [lindex $ans 1] [lindex $ans 2]]
		set z [concat $line $fans ] 
		set line $z
	}
	$stcVar(resList) insert 1 $line
	puts $stcVar(batchfp) $line
	update
}

########################################################################
#
#  check log file for errors 
#
proc check_batch_errors { f } {
	global stcVar panvar

	# Are there errors in the log file? 
	set panvar(error_mode) 0
	set unixcmd "grep Error $stcVar(fn,log)"
	if ![catch { eval exec $unixcmd } err] {
		set panvar(error_mode) 2
	}
	set unixcmd "grep Unknown $stcVar(fn,log)"
	if ![catch { eval exec $unixcmd } err] {
		set panvar(error_mode) 1
	}
	if {$panvar(error_mode) == 0} {
		$f configure -text "Done calculations, no errors."
	}
	if {$panvar(error_mode) == 1} {
		$f configure -text "Done calculations, unknown atoms recorded in $stcVar(fn,log)."
	}
	if {$panvar(error_mode) == 2} {
		$f configure -text "Done calculations, errors recorded in $stcVar(fn,log)." }
}
