 
###############################################################################
#
# Initialize and show main stc window
#
###############################################################################

###############################################################################
#
# Here starts the Main Program
#
proc stc_main {} {
	global stcVar panvar Frame

	init_vars

	remove_file $stcVar(OUT_DIR)/out.basic
	remove_file $stcVar(OUT_DIR)/out.detail 
	remove_file $stcVar(OUT_DIR)/out.asa.log
	remove_file $stcVar(OUT_DIR)/out.thermo.log
	set stcVar(errFlag)	0

	stc_show_win
}

###############################################################################
#
# init variables
#
proc init_vars { } {
	global NAME VERSION
	global stcVar gVar panvar Frame

	set stcVar(PROG)			stc
	set stcVar(VERSION)			$VERSION	
	set stcVar(VERSION_DATE) 	"Feb 3, 2006"
	set stcVar(NAME)			"$stcVar(PROG)-$stcVar(VERSION)"
	set stcVar(WEB_HOME)		"http://www.bionmr.ualberta.ca/software/stc"
	set stcVar(TMP)				".tmp"

	set outdir 					"/tmp"
	set USER_DEFAULTS			"stc.defaults"
	set SYS_DEFAULTS    		"$stcVar(INSTALL)/lib/defaults/standard"

	if ![catch {set stcVar(PWD) [exec printenv PWD]} err] {
		# output files written to current directory
		 if [file writable $stcVar(PWD)] {
			 set outdir $stcVar(PWD)
		 }
		# get copy of current settings
		set USER_DEFAULTS $stcVar(PWD)/$USER_DEFAULTS
		if ![file readable $USER_DEFAULTS ] {
			set unixcmd "cp -r $SYS_DEFAULTS $USER_DEFAULTS"
	    	if [ catch { eval exec $unixcmd } err] {
		        puts stderr  "Problems with copy command: $unixcmd \n Check your file permissions."
			}
		}
	}
	get_defaults $USER_DEFAULTS $SYS_DEFAULTS

	# find and open log files
	set_val stcVar(fn,log) 			. log	"stc.log"
	ustc_open_log

	# location of library files
	set_val stcVar(HEAD_IMAGE) . head_image	"$stcVar(INSTALL)/lib/header.gif"
	set_val stcVar(AMINO_DEF) . amino_def	"$stcVar(INSTALL)/lib/amino.def"
	set_val stcVar(HELP_FILE) . help_file	"$stcVar(INSTALL)/lib/docs/index.html"
	
	# location of binary files 
	set stcVar(BIN_DIR)			"$stcVar(INSTALL)/bin"
	set stcVar(FIX_PDB)			"$stcVar(BIN_DIR)/fix_pdb"
	set stcVar(CALC_ASA)		"$stcVar(BIN_DIR)/calc_asa"
	set stcVar(THERMO_PROC)		"$stcVar(BIN_DIR)/thermo"

	# system experiment type
	set_val panvar(sys_mode) 			. sys_mode				"P_BIND"

	# input pdb files
	set_val panvar(pdb1)	 			. pdb1					""
	set_val panvar(pdb2) 				. pdb2					""
	set_val panvar(pdb3) 				. pdb3					""
	set_val panvar(acc1) 				. acc1					""
	set_val panvar(acc2) 				. acc2					""
	set_val panvar(acc3) 				. acc3					""

	set_val stcVar(OUT_DIR)			. out_dir			$outdir	
	if ![file writable $stcVar(OUT_DIR)] {
		ustc_err "Cannot write to output directory '$stcVar(OUT_DIR)'"
	}
	ustc_log file "Output files are found in the directory: '$stcVar(OUT_DIR)'"
	
	# check if specified amino definitions file exists
	set_val panvar(amino_table_test)	. amino_def 	$stcVar(AMINO_DEF)
	set panvar(amino_table)	$stcVar(AMINO_DEF)
	if [file readable $panvar(amino_table_test)] {
		set panvar(amino_table) $panvar(amino_table_test)
	}

	# thermo calc variables
	set_val panvar(temperature)		. temperature		"25.0"
	set_val panvar(batch_vars)		. batch_vars		"DGbind Kd"
	set_val panvar(sconf_free1)		. sconf_free1		""	
	set_val panvar(sconf_free2)		. sconf_free2		""	
	set_val panvar(res_unf)			. res_unf			""
	set_val panvar(thermo_sconf)	. thermo_sconf		""


	# histogram settings
	set_val gVar(line1_color)		. line1_color		red
	set_val gVar(line2_color)		. line2_color		blue
	set_val gVar(line1_shape)		. line1_shape		oval
	set_val gVar(line2_shape)		. line2_shape		rect
	set_val gVar(line1_size)		. line1_size		3
	set_val gVar(line2_size)		. line2_size		3
	set_val gVar(canvas_geometry)	. canvas_geometry	800x400
	set_val gVar(axis_color)		. axis_color		black
	set_val gVar(title_color)		. title_color		black
	set_val gVar(bg)				. bg				white
	set_val gVar(xStart)			. xStart			90	
	set_val gVar(yStart)			. yStart			70	
	set_val gVar(yLabOffset)		. yLabOffset		20	
	set_val gVar(xScaleOffset)		. xScaleOffset		20
	set_val gVar(xTics)				. xTics				9	
	set_val gVar(yTics)				. yTics				7

	set_val panvar(main_geometry)		. main_geometry			+50+50
	set_val panvar(subscript_font)		. subscript_font		"8x13"
	set_val panvar(subscript_offset)	. subscript_offset		"-2"
	set_val panvar(superscript_offset)	. superscript_offset	"3"
	set_val panvar(symbol_font)			. symbol_font			"-adobe-symbol-medium-r-normal--17-*-*-*-*-*-adobe-fontspecific"

	set_val panvar(file_width)			. file_width 			50


	# other variables particular to this software 
	set_val stcVar(FIXED_FONT)	. fixed_font	"-bitstream-courier-medium-r-normal--0-0-0-0-m-0-iso8859-1"
	set_val stcVar(SUBSCRIPT_FONT)	. subscript_font	"8x13"
	set_val stcVar(FILE_WIDTH)		. file_width		40
	set_val stcVar(MAX_RES_ID)		. max_res_id		3000
	
}


###############################################################################
#
# draw main window for stc 
#
proc stc_show_win {} {
	global stcVar panvar Frame

	set f ".stcVar(PROG)"
	toplevel $f
	wm title $f "stc - $stcVar(VERSION)"
	wm geometry $f $panvar(main_geometry)

	# Standard tk window names for this software
	set_val Frame(stc)			. frame_stc		"$f.main"
	set_val Frame(asa)			. frame_asa		$Frame(stc)-calc-asa
	set_val Frame(showAsa)		. frame_showasa	$Frame(stc)-calc-asa-output
	set_val Frame(thermo)		. frame_thermo	$Frame(stc)-calc-thermo
	set_val Frame(batch)		. frame_batch	$Frame(stc)-multiple-pdfs
	set_val Frame(showThermo)	. frame_showthermo	$Frame(stc)-calc-thermo-output
	set_val Frame(histo)		. frame_histo	$Frame(stc)-histogram
	set_val Frame(fplotc)		. frame_fplotc	$Frame(stc)-histogram.c

	frame $Frame(stc)
	pack $Frame(stc) -anchor nw 
	show_menubar $Frame(stc)
	show_logo $Frame(stc)
	show_system_type $Frame(stc).type
	show_output_dir $Frame(stc)
	show_buts $Frame(stc).buts

	tkwait window $Frame(stc)
}

###############################################################################
#
# setup variables for specified system selected 
#
proc init_system_type {} {
	global panvar Frame stcVar

	set f $Frame(stc).type

	if {$panvar(sys_mode) == "P_OLIG"} {
		set panvar(bound) oligomer
		set panvar(free1) monomer1
		set panvar(free2) monomer2
		set panvar(bdesc) bound
		set panvar(fdesc) free 
	}

	if {$panvar(sys_mode) == "P_UNFOLD"} {
		set panvar(bound) folded
		set panvar(free1) unfolded
		set panvar(free2) "" 
		set panvar(bdesc) ""
		set panvar(fdesc) ""
	}

	if {$panvar(sys_mode) == "P_BIND"} {
		set panvar(bound) complex 
		set panvar(free1) ligand 
		set panvar(free2) enzyme 
		set panvar(bdesc) bound
		set panvar(fdesc) free 
	}

	# setup file names
	set panvar(pdb1Tmp) $stcVar(OUT_DIR)/$stcVar(TMP).$panvar(bound).pdb
	set panvar(pdb2Tmp) $stcVar(OUT_DIR)/$stcVar(TMP).$panvar(free1).pdb
	if {$panvar(acc1) == ""} {
		set panvar(acc1) $stcVar(OUT_DIR)/out.$panvar(bound).acc
	}
	if {$panvar(acc2) == ""} {
		set panvar(acc2) $stcVar(OUT_DIR)/out.$panvar(free1).acc 
	}
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		set panvar(pdb3Tmp) $stcVar(OUT_DIR)/$stcVar(TMP).$panvar(free2).pdb
		if {$panvar(acc3) == ""} {
			set panvar(acc3) $stcVar(OUT_DIR)/out.$panvar(free2).acc 
		}
	}
}

###############################################################################
#
# create menubar
#
proc show_menubar { f } {
	set mb $f.mbar
	frame $mb -relief raised -bd 2
	frame $f.dummy
	pack $mb $f -side top -fill x
	create_menus $mb
}

###############################################################################
#
# create menu items
#
proc create_menus {f} {
	global stcVar

	menubutton $f.file -padx 10 -text File -underline 0 -menu $f.file.menu
	menubutton $f.help -padx 10 -text Help -underline 0 -menu $f.help.menu
	
	menu $f.file.menu
	$f.file.menu add command -label "About the program..." \
		-command "pop_authors $f.authors"
	$f.file.menu add separator
	$f.file.menu add command -label "Exit" -command "exit_stc" 
	
	menu $f.help.menu
	$f.help.menu add command -label "see $stcVar(HELP_FILE)"
	
	pack $f.file $f.help -side left
	tk_menuBar $f $f.file $f.help
	focus $f
}

###############################################################################
#
# show main stc logo image in starting window 
#
proc show_logo { f } {
	global scr stcVar 

    canvas $f.header -height 150 -bg white -width 600
	pack $f.header -fill both -expand true
	set img [image create photo -file $stcVar(HEAD_IMAGE)]
	$f.header create image 0 0 -image $img -anchor nw
	set myTitle "Structural Thermodynamics Calculations"
	set text "$myTitle - $stcVar(VERSION) $stcVar(VERSION_DATE)"
	$f.header create text 290 120 -text $text \
		-justify center -anchor n -fill red
}

###############################################################################
#
# setup system type 
#
proc show_system_type {f} {
	global panvar

    frame $f
	pack $f -fill both -side top -pady 15 -padx 22
	label $f.syslab -text " System Type"
	pack $f.syslab

	radiobutton $f.pbind -variable panvar(sys_mode) -value P_BIND \
		-text "binding  (E + L <-> EL)"
	radiobutton $f.pfold -variable panvar(sys_mode) -value P_UNFOLD  \
		-text "folding (unfolded -> folded)"
	radiobutton $f.polig -variable panvar(sys_mode) -value P_OLIG \
		-text "oligomerization"

	foreach t "bind fold olig" {
		foreach j "1 2 3" {
			frame $f.$t$j
			text $f.$t$j.text -height 0 -padx 60 -width 40
			pack $f.$t$j.text -ipady 2
		}
	}
	pack $f.pbind $f.bind1 $f.bind2 $f.bind3 -side top -anchor w 
	pack $f.pfold $f.fold1 $f.fold2 $f.fold3 -side top -anchor w
	pack $f.polig $f.olig1 $f.olig2 $f.olig3 -side top -anchor w

	# formulate mathematical text - binding
	$f.bind1.text insert end {Kd <-> [E][L] / [EL]}
	$f.bind1.text tag add sub 1.1 1.2
	$f.bind1.text tag configure sub -font $panvar(subscript_font) \
		-offset $panvar(subscript_offset)
	$f.bind2.text insert end "DG = RT ln Kd"
	$f.bind2.text tag add delta 1.0
	$f.bind2.text tag add sub 1.12 1.13
	$f.bind2.text tag configure delta -font $panvar(symbol_font)
	$f.bind2.text tag configure sub -font $panvar(subscript_font) \
		-offset $panvar(subscript_offset)

	# formulate mathematical text - unfolding 
	$f.fold1.text insert end {Keq = [ ] / [ ]}
	$f.fold1.text tag add sub 1.1 1.3
	$f.fold1.text tag configure sub -font $panvar(subscript_font) \
		-offset $panvar(subscript_offset)
	$f.fold2.text insert end {DG = - RT ln Keq}
	$f.fold2.text tag add delta 1.0
	$f.fold2.text tag add sub 1.14 1.16
	$f.fold2.text tag configure delta -font $panvar(symbol_font)
	$f.fold2.text tag configure sub -font $panvar(subscript_font) \
		-offset $panvar(subscript_offset)

	# formulate mathematical text - oligomerization 
	$f.olig1.text insert end {E + E <-> E2}
	$f.olig1.text tag add sub 1.11 1.12
	$f.olig1.text tag configure sub -font $panvar(subscript_font) \
		-offset $panvar(subscript_offset)
}

###############################################################################
#
# show main output directory and amino acid table
#
proc show_output_dir {f} {
	global stcVar panvar

	parm_entry $f.outdir "Output directory:" 18 stcVar(OUT_DIR) 50 ""
	pack $f.outdir
	parm_entry $f.aminodef "Residue Library:" 18 panvar(amino_table) 50 ""
	pack $f.aminodef
}

###############################################################################
#
# show main window buttons 
#
proc show_buts {f} {

	frame $f
	pack $f -pady 20

	button $f.asacalcBut -text "Calc ASA" -command "show_asa_screen"
	button $f.thermocalcBut -text "Thermodynamics" -command "show_thermo_screen"
	button $f.batchBut -text "Batch" -command "show_batch_screen"
	button $f.quitBut -text "Exit" -command "exit_stc"

	pack $f.asacalcBut $f.thermocalcBut $f.batchBut $f.quitBut -side left
}

###############################################################################
#
# clean up before exiting windows 
#
proc exit_stc {} {
	global panvar stcVar stcVar

	if {$stcVar(errFlag) == 0} {
		set files [ glob -nocomplain $stcVar(OUT_DIR)/$stcVar(TMP).* ]
		foreach item $files {
			remove_file $item
		}
	}
	exit
}

###############################################################################
#
# About the program window
#
proc pop_authors { f } {
	global panvar stcVar 

	toplevel $f 

	label $f.vers -text "Structure Thermodynamic Calculations"
	label $f.date -text "$stcVar(VERSION) $stcVar(VERSION_DATE)"
	label $f.nop1
	label $f.aut \
		-text "Authors: Robert Boyko / Pierre Lavigne / John Bagu / Brian Sykes"
	label $f.nop2
	label $f.msg1 \
		-text " Note: This project is funded by the \n Protein Engineering Network Centres of Excellence and the \n Medical Research Council of Canada"

	label $f.nop3
	label $f.msg3 -text " Website - \
	 $stcVar(WEB_HOME) \n \
	 Email - webmaster@pence.ca"
	button $f.buts -text "Ok" -command "destroy $f"
	label $f.nop4
	pack $f.vers $f.date $f.nop1 $f.aut $f.nop2 $f.msg1 
	pack $f.nop3 $f.msg3 $f.nop4
	pack $f.buts -side top
}

#####################################################################
#
# read the defaults file (look in user's directory and system location)
#
proc get_defaults { f_user f_sys } {

#
# first read system defaults
#
    if ![file readable $f_sys ] {
		puts stderr "Unable to open installed defaults file: $f_sys"
		exit
	}

	if [catch {option readfile $f_sys interactive} defaults_err] {
		puts stderr "Problem reading option in $f_sys: $defaults_err"
		exit
	}

#
# then read user defaults if any exist
#
	if [file readable $f_user ] {
		if [catch {option readfile $f_user interactive} err] {
			puts stderr "Warning, problem reading option in $f_user: $err"
		}
	}
}
