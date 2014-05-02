
######################################################################
#
# misc TCL utilities
#
######################################################################

######################################################################
#
# set a variable from the defaults file
#
proc set_val {var struc entry val} {
	upvar $var a
	set val [string trim $val]
	if ![info exists a] {
		set a [string trim [option get $struc $entry {}]]
		if {$a == ""} { set a $val }
	}
}

######################################################################
#
# create a generic widget for file text display
#
proc create_file_panel { f width height } {
	global stcVar 
	frame $f
	text $f.text -relief raised -bd 2 \
		-width $width -height $height \
		-setgrid true -wrap none \
		-xscrollcommand [list $f.xscroll set] \
		-yscrollcommand [list $f.yscroll set]
	scrollbar $f.xscroll -orient horizontal \
		-command [list $f.text xview]
	scrollbar $f.yscroll -orient vertical \
		-command [list $f.text yview]
	pack $f.xscroll -side bottom -fill x
	pack $f.yscroll -side right -fill y
	pack $f.text -side left -fill both -expand true
}

######################################################################
#
# generic file field entry
#
proc file_entry {f labName labLen entVar entLen command} {
	frame $f
	button $f.but -text "Browse" -width 6 -command "openDataFile $entVar"
	frame $f.fname
	label $f.fname.label -text $labName -width $labLen -anchor w -relief flat
	entry $f.fname.entry -relief sunken -textvar $entVar -width $entLen
	bind $f.fname.entry <Return> $command
	pack $f.fname.label -side left
	pack $f.fname.entry -fill x -expand true 
	pack $f.but $f.fname -side left
}
proc openDataFile {aVar} {
	global panvar 
	set $aVar [tk_getOpenFile -title "Please select an input file"]
}

######################################################################
#
# parameter field entry
#
proc parm_entry {f labName labLen entVar entLen command} {
	frame $f
	label $f.label -text $labName -width $labLen -anchor w -relief flat
	entry $f.entry -relief sunken -textvar $entVar -width $entLen
	bind $f.entry <Return> $command
	pack $f.label -side left
	pack $f.entry -fill x -expand true
}

######################################################################
#
# generic command field entry
#
proc command_entry {name label width command args} {
	frame $name
	label $name.label -text $label -width $width -anchor w
	eval {entry $name.entry -relief sunken } $args
	pack $name.label -side left
	pack $name.entry -side right -fill x -expand true
	bind $name.entry <Return> $command
	return $name.entry
}

######################################################################
#
# save graphics to postscript file
#
proc pop_save_ps { f } {
	global Frame stcVar 

	toplevel $f 
	label $f.nop1
	label $f.nop2

	set savecmd {$Frame(fplotc) postscript -file $stcVar(OUT_DIR)/$psFname; exec chmod 600 $stcVar(OUT_DIR)/$psFname}

	command_entry $f.ps "Save to postscript file: " 26 \
		"$savecmd; destroy $f" -textvar psFname -width 20 
	pack $f.nop1 $f.ps $f.nop2 -side top

	button $f.but -text "Cancel" -command "destroy $f"
	pack $f.but -side top
}

######################################################################
#
# loads a file onto a canvas - code taken form "TCL and the TK toolkit"
#
proc u_load_file { file wind } {
	$wind delete 1.0 end
	if [catch {open $file r} f] {
   		$wind insert end "Not able to read selected file." 
		return -1
	}
	while {![eof $f]} {
   		$wind insert end [read $f 10000]
	}
	close $f
}

######################################################################
#
# loads result files with symbolic fonts added
#
proc load_results {file wt} {
	global stcVar panvar

	$wt tag configure symbol -font $panvar(symbol_font) 
	$wt tag configure subscript -font $panvar(subscript_font) 
	$wt tag configure superscript -font $panvar(subscript_font) \
		-offset $panvar(superscript_offset) 

	$wt delete 1.0 end
	set f [open $file ]
	while {![eof $f]} {
		$wt insert end [read $f 10000]
	}
	close $f

	set delta {D[GCSHA]}
	forAllMatches $wt $delta {
		$wt tag add symbol first "first + 1 chars"
	}

	set unitLabel {bind  |conf  |fold  |sol   |sc    |tot   |rt    |bu->ex|ex->u |bb    }
	forAllMatches $wt $unitLabel {
		$wt tag add subscript first last 
	}

	set angs [format "%c2" 197]
	forAllMatches $wt $angs {
		$wt tag add superscript "first + 1 chars" last
	}
}
proc forAllMatches {w pat script} {
	scan [$w index end] %d numLines
	for {set i 1} {$i < $numLines} {incr i} {
		$w mark set last $i.0
		while {[regexp -indices $pat [$w get last "last lineend"] indices]} {
			$w mark set first "last + [lindex $indices 0] chars"
			$w mark set last "last + 1 chars + [lindex $indices 1] chars"
			uplevel $script
		}
	}
}


######################################################################
#
# set a variable if it does not already exist
#
proc u_set_val {var val} {
	upvar $var a
	set val [string trim $val]
	if ![info exists a] {
		set a $val
	}
}

######################################################################
#
# purpose: to call a tk routine to search for a directory
#
proc u_dir_browse { dir } {
	global Frame stcVar panvar

	set chosen [tk_chooseDirectory]
	if { $chosen != "" } {
		set currentpwd [pwd]
		set cmd "regsub {$currentpwd\/} $chosen \"\" shortchosen"
		eval $cmd
		set $dir $shortchosen
	}
	show_pdb $Frame(batch)
}

######################################################################
#
# generic file entry widget
#
proc u_file_entry { f title tWidth tVar eWidth padx cmd browse view } {
    pack [frame $f] -fill x -expand true -padx 1
    pack [frame $f.a] -side left
    pack [frame $f.b] -side left
    pack [frame $f.c] -side left 
	label $f.a.label -text $title -anchor e -width $tWidth
	pack $f.a.label -side right
	entry $f.b.entry -bd 2 -textvariable $tVar -width $eWidth -relief sunken
#	$f.b.entry config -cursor {xterm white}
	bind $f.b.entry <Return> $cmd
	pack $f.b.entry -side left -anchor nw
	if { $browse > 0 } {
		if { $browse == 1 } {
			button $f.c.browse -text " Browse " \
				-command "u_file_browse $tVar" 
		} else {
			button $f.c.browse -text " Browse " \
				-command "u_dir_browse $tVar" 
		}
		pack $f.c.browse -side left -padx 8 -anchor nw
	}
	if { $view > 0 } {
		button $f.c.view -text " View " \
			-padx $padx \
			-command "u_file_view $tVar" 
		pack $f.c.view -side left -padx 8
	}
}

######################################################################
#
# show file view window
#
proc u_file_view { fVar } {
	global stcVar

	# get file name from variable
	upvar $fVar fn 

	# check state of popup window
	set f [ustc_show_win showfile "File: '$fn'"]
	if { $f == "FALSE" } { return }

	frame $f.win
	frame $f.buttons
	pack $f.win -fill both -expand true
	pack $f.buttons

	# show text window and load file
	u_textwin $f.win.panel 80 30
	u_load_file $fn $f.win.panel.textwin

	# buttons for this window
	button $f.buttons.quitButton \
		-text " Dismiss " \
		-command "destroy $f"
	pack $f.buttons.quitButton -side left -pady 10
}

######################################################################
#
# raise existing window on desktop or create it if it does not exist
#
proc ustc_show_win { winStr title } {

	set f .$winStr 

	# if window already exists, turf it 
	if [winfo exists $f] {
		destroy $f
	}

	# window does not exist, we should create it and set geometry
	toplevel $f
	wm title $f $title
	wm geometry $f +50+50
	return $f 
}


######################################################################
#   
# create a generic widget for file text display
#       
proc u_textwin { f width height } {
	frame $f
	pack $f -expand true -fill both -padx 10 -pady 10
    text $f.textwin -relief raised -bd 2 \
		-width $width -height $height \
		-setgrid false -wrap none \
		-xscrollcommand [list $f.xscroll set] \
		-yscrollcommand [list $f.yscroll set]
	scrollbar $f.xscroll -orient horizontal \
		-command [list $f.textwin xview]
	scrollbar $f.yscroll -orient vertical \
		-command [list $f.textwin yview]
	pack $f.xscroll -side bottom -fill x
	pack $f.yscroll -side right -fill y
	pack $f.textwin -side left -fill both -expand true
}
