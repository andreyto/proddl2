
################################################################################
#
#  Accessible surface area routines 
#
################################################################################

################################################################################
#
# setup the calc asa screen 
#
proc show_asa_screen {} {
	global panvar stcVar Frame

	if ![file writable $stcVar(OUT_DIR)] {
		if ![make_dir $stcVar(OUT_DIR)] {
			ustc_err "Cannot write to output directory '$stcVar(OUT_DIR)'."
		}
	}

	set f $Frame(asa)
	toplevel $f 
	init_system_type

	# decide label names for file entry fields
	if {$panvar(sys_mode) == "P_OLIG"} {
		set n1 "Bound oligomer pdb file:"
		set n2 "(Optional) Free monomer pdb:"
	}
	if {$panvar(sys_mode) == "P_UNFOLD"} {
		set n1 "Folded pdb file:"
		set n2 "Unfolded pdb file:"
	}
	if {$panvar(sys_mode) == "P_BIND"} {
		set n1 "Bound complex pdb file (ligand first):"
		set n2 "(Optional) Free ligand pdb:"
		set n3 "(Optional) Free enzyme pdb:"
	}
	set labelLen 38

	# define asa fields 
	label $f.title -text "\nCalculating Accessible Surface Area\n"
	pack $f.title
	file_entry $f.pdb1file $n1 $labelLen panvar(pdb1) $panvar(file_width) \
		"focus $f.pdb2file.fname.entry" 
	pack $f.pdb1file -padx 10
	file_entry $f.pdb2file $n2 $labelLen panvar(pdb2) $panvar(file_width) "bind_pdb2_file"
	pack $f.pdb2file -padx 10
	if {$panvar(sys_mode) == "P_BIND"} {
		file_entry $f.pdb3file $n3 $labelLen panvar(pdb3) $panvar(file_width) \
			"focus $f.pdb1file.fname.entry" 
		pack $f.pdb3file -padx 10
	}
	label $f.status -text "\n"
	pack $f.status -pady 10

	# define asa buttons
	frame $f.buts
	pack $f.buts -pady 5
	button $f.buts.calc -text "Calculate" -command "calc_all_asa" 
	button $f.buts.disp -text "Display Results" -command "show_asa_results"
	button $f.buts.exit -text "Dismiss" -command "destroy $f"
	pack $f.buts.calc $f.buts.disp $f.buts.exit -side left

	focus $f.pdb1file.fname.entry
}

################################################################################
#
# user clicks on the button to calculate asa
#
proc calc_all_asa {} {
	global stcVar panvar Frame 

	# check input fields first
	set f $Frame(asa)
	set panvar(pdb1) [string trim $panvar(pdb1)]
	set panvar(pdb2) [string trim $panvar(pdb2)]
	set panvar(pdb3) [string trim $panvar(pdb3)]

	# reformat pdb input data
	$f.status configure -text "Reformatting PDB Data...\n"
	update
	if [catch {open $panvar(pdb1) r} fileId] {
		set stcVar(errFlag) 1
		ustc_show_log
		error "Cannot open 1st pdb data file: '$panvar(pdb1)'"
	}
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		fix_pdb 2 $panvar(pdb1) $panvar(pdb1Tmp) $panvar(pdb2Tmp) $panvar(pdb3Tmp)
		if {$panvar(pdb2) != ""} {
			if [catch {open $panvar(pdb2) r} fileId] {
				set stcVar(errFlag) 1
				ustc_show_log
				error "Cannot open 2nd pdb data file: '$panvar(pdb2)'"
			}
			fix_pdb 1 $panvar(pdb2) $panvar(pdb2Tmp) $stcVar(OUT_DIR)/tmp.1 $stcVar(OUT_DIR)/tmp.2 
		}
		if {$panvar(pdb3) != ""} {
			if [catch {open $panvar(pdb3) r} fileId] {
				set stcVar(errFlag) 1
				ustc_show_log
				error "Cannot open 3nd pdb data file: '$panvar(pdb3)'"
			}
			fix_pdb 1 $panvar(pdb3) $panvar(pdb3Tmp) $stcVar(OUT_DIR)/tmp.1 $stcVar(OUT_DIR)/tmp.2 
		}
	} elseif {$panvar(sys_mode) == "P_UNFOLD"} {
		if [catch {open $panvar(pdb2) r} fileId] {
			set stcVar(errFlag) 1
			ustc_show_log
			error "Cannot open 2nd pdb data file: '$panvar(pdb2)'"
		}
		fix_pdb 1 $panvar(pdb1) $panvar(pdb1Tmp) $stcVar(OUT_DIR)/tmp.1 $stcVar(OUT_DIR)/tmp.2 
		fix_pdb 1 $panvar(pdb2) $panvar(pdb2Tmp) $stcVar(OUT_DIR)/tmp.1 $stcVar(OUT_DIR)/tmp.2
	}

	# calculate ASA
	$f.status configure -text "Making surface accessible file $panvar(acc1)\n"
	update
	calc_asa $panvar(bound) $panvar(pdb1Tmp) $panvar(acc1)
	$f.status configure -text "Making surface accessible file $panvar(acc2)\n"
	update
	calc_asa $panvar(free1) $panvar(pdb2Tmp) $panvar(acc2)
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		$f.status configure -text "Making surface accessible file $panvar(acc3)...\n"
		update
		calc_asa $panvar(free2) $panvar(pdb3Tmp) $panvar(acc3)
	}

	check_asa_errors $f.status
	update
}

################################################################################
#
# display results for calculating accesible surface area
#
proc show_asa_results {} {
    global panvar Frame stcVar

	set f $Frame(showAsa)
	toplevel $f 
	frame $f.fr1
	
	# define screen widgets
	radiobutton $f.fr1.accb -text "$panvar(acc1)" -variable panvar(ofn) -value "a"
	radiobutton $f.fr1.acc1 -text "$panvar(acc2)" -variable panvar(ofn) -value "b"
	radiobutton $f.fr1.bound -text "$panvar(pdb1Tmp)" -variable panvar(ofn) -value "c"
	radiobutton $f.fr1.free1 -text "$panvar(pdb2Tmp)" -variable panvar(ofn) -value "d"
	radiobutton $f.fr1.log -text "Calculations log" -variable panvar(ofn) -value "e"

	bind $f.fr1.accb <Button-1> "u_load_file $panvar(acc1) $f.results.text"
	bind $f.fr1.acc1 <Button-1> "u_load_file $panvar(acc2) $f.results.text"
	bind $f.fr1.bound <Button-1> "u_load_file $panvar(pdb1Tmp) $f.results.text"
	bind $f.fr1.free1 <Button-1> "u_load_file $panvar(pdb2Tmp) $f.results.text"
	bind $f.fr1.log <Button-1> "u_load_file $stcVar(fn,log) $f.results.text"

	if {$panvar(sys_mode) != "P_UNFOLD"} {
		radiobutton $f.fr1.free2 -text "$panvar(pdb3Tmp)" -variable panvar(ofn) -value "y"
		radiobutton $f.fr1.acc2 -text "$panvar(acc3)" -variable panvar(ofn) -value "z"
		bind $f.fr1.free2 <Button-1> "u_load_file $panvar(pdb3Tmp) $f.results.text"
		bind $f.fr1.acc2 <Button-1> "u_load_file $panvar(acc3) $f.results.text"
	}

	create_file_panel $f.results 80 30
	button $f.quit -text "Dismiss" -command "destroy $f"

	# configure and pack widgets accordingly
	pack $f.fr1 -pady 10
	pack $f.fr1.accb $f.fr1.acc1 -anchor w -padx 15
	if {$panvar(sys_mode) != "P_UNFOLD"} {
		pack $f.fr1.acc2 -anchor w -padx 15
	}
	pack $f.fr1.bound $f.fr1.free1 -anchor w -padx 15 
	if {$panvar(sys_mode) != "P_UNFOLD"} {
		pack $f.fr1.free2 -anchor w -padx 15
	}
	pack $f.fr1.log -anchor w -padx 15 
	pack $f.results -pady 10 -side top -fill both -expand true
	pack $f.quit -pady 5

	if {$panvar(error_mode) > 0} {
		set panvar(ofn) "e"
		u_load_file $stcVar(fn,log) $f.results.text
	} else {
		# load first file if it exists 
		set panvar(ofn) "a"
		if [catch {open $panvar(acc1) r} fileId] { return }
		u_load_file $panvar(acc1) $f.results.text
	}
}

################################################################################
#
# user enters a pdb2file, decide where the focus is
#
proc bind_pdb2_file {} {
	global panvar Frame

	if {$panvar(sys_mode) == "P_BIND"} {
		focus $Frame(asa).pdb3file.fname.entry
	} else {
		focus $Frame(asa).pdb1file.fname.entry
	}
}

########################################################################
#
# run fix_pdb script
#
proc fix_pdb { type pdb1 bound1 free1 free2 } {
	global stcVar panvar

	set unixcmd "$stcVar(FIX_PDB) $type $pdb1 $bound1 $free1 $free2"
	if [catch {eval exec $unixcmd >> $stcVar(fn,log)} err] {
		set stcVar(errFlag) 1
		ustc_show_log
		error "Problems re-formatting $pdb1."
	}
}

########################################################################
#
# run calc asa script 
#
proc calc_asa { root pdb acc } {
	global stcVar panvar

	# standard parameters for calc_asa call
	set a "$stcVar(CALC_ASA) $stcVar(BIN_DIR) $panvar(amino_table) $stcVar(OUT_DIR)"

	# add variable parameters
	set cmd "$a $root $pdb $acc"

	if [catch {eval exec $cmd >> $stcVar(fn,log)} err] {
		set stcVar(errFlag) 1
		ustc_show_log
		error "Problems calculating asa for $pdb."
	}

	# we are successful, therefore we can remove temporary files 
	remove_file $stcVar(OUT_DIR)/$root.addradii
	remove_file $stcVar(OUT_DIR)/$root.edpdb
	remove_file $stcVar(OUT_DIR)/$root.dia
	remove_file $stcVar(OUT_DIR)/$root.seq
}

########################################################################
#
#  check log file for errors 
#
proc check_asa_errors { f } {
	global stcVar panvar

	# Are there errors in the log file? 
	set panvar(error_mode) 0
	set unixcmd "grep Unknown $stcVar(fn,log)"
	if ![catch { eval exec $unixcmd } err] {
		set panvar(error_mode) 1
	}
	set unixcmd "grep Error $stcVar(fn,log)"
	if ![catch { eval exec $unixcmd } err] {
		set panvar(error_mode) 2
	}
	if {$panvar(error_mode) == 0} {
		$f configure -text "Done calculations, no errors.\n"
	}
	if {$panvar(error_mode) == 1} {
		$f configure -text "Done calculations, display results, log file has unknown atoms.\n"
	}
	if {$panvar(error_mode) == 2} {
		$f configure -text "Done calculations, display results, log file has errors.\n"
	}
}
