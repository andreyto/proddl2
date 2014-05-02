
################################################################################
#
# thermodynamics calculation routines 
#
################################################################################

################################################################################
#
# display the thermodynamics screen
#
proc show_thermo_screen {} {
	global Frame panvar

	set f $Frame(thermo)
	toplevel $f 
	init_system_type

	set labelLen 30
	set padx 10

	# define title 
	label $f.title -text "\nCalculating Thermodynamics from ASA Files\n"
	pack $f.title

	# define file entry fields
	file_entry $f.t1file "ASA file - $panvar(bdesc) $panvar(bound):" \
		$labelLen panvar(acc1) $panvar(file_width) ""
	pack $f.t1file -padx $padx
	file_entry $f.t2file "ASA file - $panvar(fdesc) $panvar(free1):" \
		$labelLen panvar(acc2) $panvar(file_width) ""
	pack $f.t2file -padx $padx
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		file_entry $f.t3file "ASA file - $panvar(fdesc) $panvar(free2):" \
			$labelLen panvar(acc3) $panvar(file_width) ""
		pack $f.t3file -padx $padx
	}
	file_entry $f.sconf_free1 "(Optional) $panvar(free1) sConf File:" \
		$labelLen panvar(sconf_free1) $panvar(file_width) ""
	pack $f.sconf_free1 -padx $padx
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		file_entry $f.sconf_free2 "(Optional) $panvar(free2) sConf File:" \
			$labelLen panvar(sconf_free2) $panvar(file_width) ""
		pack $f.sconf_free2 -padx $padx
	}

	# define parameter fields
	set pady 3
	label $f.nop -text "\n"
	pack $f.nop
	parm_entry $f.temp "Temperature:" $labelLen panvar(temperature) $panvar(file_width) ""
	pack $f.temp -padx $padx -pady $pady -anchor ne
	parm_entry $f.resunf "(Optional) Residues Unfolded:" $labelLen panvar(res_unf) $panvar(file_width) ""
	pack $f.resunf -padx $padx -pady $pady -anchor ne
	parm_entry $f.sconf "(Optional) Total sConf:" $labelLen panvar(thermo_sconf) \
		$panvar(file_width) ""
	pack $f.sconf -padx $padx -pady $pady -anchor ne

	# radio buttons section
	set panvar(sbb) 0
	set panvar(sexu) 0
	label $f.nop2 -text "\n"
	pack $f.nop2

    frame $f.sbb
	pack $f.sbb -anchor w
	label $f.sbb.label -text "Sbb:   " -width 8
	radiobutton $f.sbb.zero -text "set zero" \
		-variable panvar(sbb) -value 0
	radiobutton $f.sbb.data -text "calculate from datafile" \
		-variable panvar(sbb) -value 1
	pack $f.sbb.label $f.sbb.zero $f.sbb.data -padx 10 -pady 0 -side left

    frame $f.sexu
	pack $f.sexu -anchor w
	label $f.sexu.label -text "Sex->u:" -width 8
	radiobutton $f.sexu.zero -text "set zero" \
		-variable panvar(sexu) -value 0
	radiobutton $f.sexu.data -text "calculate from datafile" \
		-variable panvar(sexu) -value 1
	radiobutton $f.sexu.calc -text "calculate from tables" \
		-variable panvar(sexu) -value 2
	pack $f.sexu.label $f.sexu.zero $f.sexu.data $f.sexu.calc \
		-padx 10 -pady 0 -side left

	# calculate/exit buttons
	frame $f.buts
	pack $f.buts -pady 20 
	button $f.buts.calc -text "Calculate" -command "do_thermo" 
	button $f.buts.exit -text "Exit" -command "destroy $f"
	pack $f.buts.calc $f.buts.exit -side left
}

################################################################################
#
# user clicks on the calculate button
#
proc do_thermo {} {
	global stcVar panvar

	set_val stcVar(RUN_FILE) . run_file "$stcVar(OUT_DIR)/$stcVar(TMP).thermo.run"
	set panvar(acc1) [string trim $panvar(acc1)]
	set panvar(acc2) [string trim $panvar(acc2)]
	set panvar(acc3) [string trim $panvar(acc3)]

	set panvar(thermo_sconf) [string trim $panvar(thermo_sconf)]
	if [catch {open $stcVar(RUN_FILE) w} err] {
		set stcVar(errFlag) 1
		error "Cannot open file: $stcVar(RUN_FILE)"
	}
	set fp [open $stcVar(RUN_FILE) w]

	if {$panvar(sys_mode) == "P_BIND"} {
		puts $fp "1"
	} elseif {$panvar(sys_mode) == "P_UNFOLD"} {
		puts $fp "2"
	} else { puts $fp "3" }

	puts $fp "$stcVar(OUT_DIR)/out.basic"
	puts $fp "$stcVar(OUT_DIR)/out.detail"
	puts $fp "$panvar(amino_table)"
	puts $fp "$panvar(acc1)"
	puts $fp "$panvar(acc2)"
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		puts $fp "$panvar(acc3)"
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

	calc_thermo $panvar(thermo_sconf) $stcVar(RUN_FILE) $stcVar(OUT_DIR)/out.thermo.log
	show_thermo_results $stcVar(OUT_DIR)/out.basic $stcVar(OUT_DIR)/out.detail "a"
}

################################################################################
#
# display thermodynamic results in a window
#
proc show_thermo_results {basic verbose showFirst} {
    global Frame panvar

	set f $Frame(showThermo)
	toplevel $f 

	# define widgets for the screen
    frame $f.fr0
	pack $f.fr0 -pady 10
	radiobutton $f.fr0.basic -text $basic -variable panvar(oFile) \
		-value "a" -command "load_results $basic $f.results.text"
	radiobutton $f.fr0.detail -text $verbose -variable panvar(oFile) \
		-value "b" -command "load_results $verbose $f.results.text"
	pack $f.fr0.basic $f.fr0.detail -padx 15 -side top

	frame $f.fr5
	pack $f.fr5
	button $f.fr5.graph -text "Histograms" -command "show_histo_panel $verbose"
	button $f.fr5.quit -text "Dismiss" -command "destroy $f"
	pack $f.fr5.graph $f.fr5.quit -padx 15 -pady 10 -side left

	create_file_panel $f.results 80 39
	pack $f.results -side top -fill both -expand true

	# load log file if it exists 
	set panvar(oFile) $showFirst
	if [catch {open $basic r} err] { return }
	load_results $basic $f.results.text
}

################################################################################
#
# run thermo script 
#
proc calc_thermo { sconf rundeck logfile } {
	global stcVar panvar

	set cmd "$stcVar(THERMO_PROC) -temp $panvar(temperature) -log $logfile"
	if {$sconf != ""} {
		set cmd "$cmd -sconf $sconf"
	}
	if [catch {eval exec cat $rundeck | $cmd } err] {
		set stcVar(errFlag) 1
		show_thermo_results $logfile $logfile "a"
		error "Problems with thermodynamics calculations. Press Ok to view log file."
	}
	remove_file $rundeck
}
