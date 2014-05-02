
######################################################################
#
# misc TCL utilities specific for stc
#
######################################################################


######################################################################
#
# write msg to the log file
#
proc ustc_log { myCode msg } {
	global stcVar

	# write message to the screen
	if { $myCode == "stdout" || $myCode == "all" } { 
		puts stdout "stc: $msg"
	}

	# put message in log file
	if { $myCode == "file" || $myCode == "all" } { 
		if [catch {open $stcVar(fn,log) a} fileId] {
			puts stdout "Cannot open log file '$stcVar(fn,log)'."
			return
		}
		puts $fileId $msg
		close $fileId
	}
	# put output line in batch log file
	if { $myCode == "batch" } { 
		if [catch {open $stcVar(fn,batch) a} fileId] {
			puts stdout "Cannot open batch file '$stcVar(fn,batch)'."
			return
		}
		puts $fileId $msg
		close $fileId
	}
}

######################################################################
#
# open log file
#
proc ustc_open_log { } {
	global stcVar

	if [catch {open $stcVar(fn,log) w} fileId] {
		puts stdout "Cannot open log file '$stcVar(fn,log)'."
		exit
	}
	set date [get_date] 
	ustc_log file ""
	ustc_log file "Installation directory for $stcVar(NAME) is: '$stcVar(INSTALL)'"
	ustc_log file "Date: '$date'"
	ustc_log file ""

}

######################################################################
#
# write err to the log file, show in tcl window
#
proc ustc_err { msg } {
	global stcVar 

	bell

	# This is because Macosx tcl error doesn't do what one would think here 
	if [info exists stcVar(stcmain,winName)] { 
		if [winfo exists $stcVar(stcmain,winName)] { 
			ustc_log all "Error: $msg"
		} else {	
			ustc_log file "Error: $msg"
		}
	}
	error "\n$msg"
}

######################################################################
#
# write err to the log file, show in tcl window
#
proc ustc_show_log { } {
	global stcVar Frame

	set f .showlog
	toplevel $f
	create_file_panel $f.fr 80 50
	pack $f.fr -pady 10 -side top -fill both -expand true
	u_load_file $stcVar(fn,log) $f.fr.text
	update
}
