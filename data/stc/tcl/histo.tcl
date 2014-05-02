
###############################################################################
#
# Histogram code for stc
#
###############################################################################

###############################################################################
#
# show main histogram panel
#
proc show_histo_panel { outFile } {
	global gVar panvar Frame

	set f $Frame(histo)
	toplevel $f
	wm title $f "Stc histograms" 
	wm geometry $f $gVar(canvas_geometry)

	frame $f.fr1
	pack $f.fr1 -anchor w -pady 10
	label $f.fr1.label -text "Structure:" -width 10
	radiobutton $f.fr1.ligand -text $panvar(free1) -variable gVar(type) \
		-value "a" -command "new_histo $f $outFile"
	if {$panvar(sys_mode) == "P_BIND" || $panvar(sys_mode) == "P_OLIG"} {
		radiobutton $f.fr1.enzyme -text $panvar(free2) -variable gVar(type) \
			-value "b" -command "new_histo $f $outFile"
		pack $f.fr1.label $f.fr1.ligand $f.fr1.enzyme -padx 5 -side left
	} else {
		pack $f.fr1.label $f.fr1.ligand -padx 5 -pady 1 -side left
	}

    frame $f.fr2
	pack $f.fr2 -anchor w
	label $f.fr2.label -text "Sconf:" -width 10
	radiobutton $f.fr2.asa -text "DASA" -variable gVar(func) \
		-value 3C3C4 -command "new_histo $f $outFile" 
	radiobutton $f.fr2.dasatot -text "DASATot" -variable gVar(func) \
		-value 4C3C3 -command "new_histo $f $outFile" 
	radiobutton $f.fr2.dasasc -text "DASASidechain" -variable gVar(func) \
		-value 4C4C4 -command "new_histo $f $outFile" 
	radiobutton $f.fr2.dsbuex -text "DSbu->ex" -variable gVar(func) \
		-value 4C5C5 -command "new_histo $f $outFile" 
	radiobutton $f.fr2.sConf -text "DsConf" -variable gVar(func)  \
		-value 4C8C8 -command "new_histo $f $outFile" 
	pack $f.fr2.label $f.fr2.asa $f.fr2.dasatot $f.fr2.dasasc $f.fr2.dsbuex \
		$f.fr2.sConf -padx 5 -side left

    frame $f.fr4
	pack $f.fr4 -anchor w
	label $f.fr4.label -text "Values:" -width 10
	radiobutton $f.fr4.cpbind -text "DCpBind" -variable gVar(func) \
		-value 5C3C3 -command "new_histo $f $outFile" 
	radiobutton $f.fr4.hbind -text "DhBind" -variable gVar(func) \
		-value 5C4C4 -command "new_histo $f $outFile"
	radiobutton $f.fr4.tdsol -text "TDSsol" -variable gVar(func) \
		-value 5C5C5 -command "new_histo $f $outFile"
	radiobutton $f.fr4.tdsconf -text "TDSconf" -variable gVar(func) \
		-value 5C6C6 -command "new_histo $f $outFile"
	radiobutton $f.fr4.tdsbind -text "TDSbind" -variable gVar(func) \
		-value 5C7C7 -command "new_histo $f $outFile"
	radiobutton $f.fr4.gbind -text "DGbind" -variable gVar(func) \
		-value 5C8C8 -command "new_histo $f $outFile"
	pack $f.fr4.label $f.fr4.cpbind $f.fr4.hbind $f.fr4.tdsol $f.fr4.tdsconf \
		$f.fr4.tdsbind $f.fr4.gbind -padx 5 -pady 1 -side left
	
	frame $f.buttons
	pack $f.buttons
	button $f.buttons.1 -text "Save Postscript" -command "pop_save_ps .ps"
	button $f.buttons.2 -text "Fit to Window" -command "draw_histo $f"
	button $f.buttons.3 -text "Dismiss" -command "destroy $f"
	pack $f.buttons.1 $f.buttons.2 $f.buttons.3 -side left -pady 10

	# create histogram canvas with movable elements
	canvas $Frame(fplotc) -bg $gVar(bg)
	pack $Frame(fplotc)
	$Frame(fplotc) bind movable <Button-1> {Mark %x %y %W}
	$Frame(fplotc) bind movable <B1-Motion> {Drag %x %y %W}
	update
	
	# default histogram shown
	set gVar(type) "a"
	set gVar(func) "3C3C4"
	new_histo $f $outFile
}

###############################################################################
#
# user has selected a different histogram to draw 
#
proc new_histo { f outFile } {
	global panvar gVar

	# read the data, set parameters and draw it.
	read_verbose $outFile $gVar(type) $gVar(func)
	set_histo_parms
	draw_histo $f
}

###############################################################################
#
# find specific parameters for the histogram we are to draw
#
proc set_histo_parms { } {
	global gVar gData

	set xLabText "Residue"
	set titleText "Ligand"
	if {$gVar(type) == "b"} { set titleText "Enzyme" }

	switch $gVar(func) {
		3C3C4 {
			set yLabText "dASA"
			set titleText " $titleText - Amino Acid Differences in dASA $gVar(line1_color)=polar $gVar(line2_color)=nonpolar"
		}
		4C3C3 {
			set yLabText "DASATot"
			set titleText " $titleText - Calculating sConf per Residue"
		}
		4C4C4 {
			set yLabText "DASA\nSideChain"
			set titleText " $titleText - Calculating sConf per Residue"
		}
		4C5C5 {
			set yLabText "DSbu->ex"
			set titleText " $titleText - Calculating sConf per Residue"
		}
		4C8C8 {
			set yLabText "DSconf"
			set titleText " $titleText - Calculating sConf per Residue"
		}
		5C3C3 {
			set yLabText "DcpBind"
			set titleText " $titleText - Thermodynamic Values per Residue"
		}
		5C4C4 {
			set yLabText "DhBind"
			set titleText " $titleText - Thermodynamic Values per Residue"
		}
		5C5C5 {
			set yLabText "TDSsol"
			set titleText " $titleText - Thermodynamic Values per Residue"
		}
		5C6C6 {
			set yLabText "TDSconf"
			set titleText " $titleText - Thermodynamic Values per Residue"
		}
		5C7C7 {
			set yLabText "TDSbind"
			set titleText " $titleText - Thermodynamic Values per Residue"
		}
		5C8C8 {
			set yLabText "DGbind"
			set titleText " $titleText - Thermodynamic Values per Residue"
		}
	}
	set xaxis "middle"
	if {$gData(yMax) == 0} { set xaxis "top" }
	if {$gData(yMin) == 0} { set xaxis "bottom" }
	if {$xaxis == "middle"} { 
		set yTmp [expr 0 - $gData(yMin)]
		if {$yTmp > $gData(yMax)} { set gData(yMax) $yTmp }
		if {$yTmp < $gData(yMax)} { set gData(yMin) [expr 0 - $gData(yMax)] }
	}

	set gVar(xLabText) $xLabText
	set gVar(yLabText) $yLabText
	set gVar(titleText) $titleText
	set gVar(xaxisLoc) $xaxis
	set gVar(xMin) $gData(xMin)
	set gVar(xMax) $gData(xMax)
	set gVar(yMin) $gData(yMin)
	set gVar(yMax) $gData(yMax)
}

###############################################################################
#
# how to read detailed output file for histogram data to show
#
proc read_verbose { file aStruct func } {
	global stcVar gData
	global x1 y1 y2

	set section [lindex [split $func C] 0]
	set y1Col [lindex [split $func C] 1]
	set y2Col [lindex [split $func C] 2]

	set pat "*$section$aStruct:*"
	set gData(xMax) 0 
	set gData(xMin) $stcVar(MAX_RES_ID)
	set gData(yMax) 0.0 
	set gData(yMin) 0.0 
	set gData(numData) 0
	set numLines 0
	for {set i 0} {$i <= $stcVar(MAX_RES_ID)} {incr i} {
		set y1($i) 0.0
		set y2($i) 0.0
	}
	set finp [open $file]
	while {[gets $finp line] >= 0} {
		if { 1 == [string match $pat $line] } {
			gets $finp line
			gets $finp line
			gets $finp line
			gets $finp line
			gets $finp line

			set nd [scan $line "%d %s %f %f %f %f %f %f "  id res a(3) a(4) a(5) a(6) a(7) a(8)]
			while { $nd > 2 } {
				set y1($id) $a($y1Col) 
				set y2($id) $a($y2Col)
				if {$id > $gData(xMax)} { set gData(xMax) $id }
				if {$id < $gData(xMin)} { set gData(xMin) $id }
				if {$y1($id) > $gData(yMax)} { set gData(yMax) $y1($id) }
				if {$y1($id) < $gData(yMin)} { set gData(yMin) $y1($id) }
				if {$y2($id) > $gData(yMax)} { set gData(yMax) $y2($id) }
				if {$y2($id) < $gData(yMin)} { set gData(yMin) $y2($id) }
				incr gData(numData)
				gets $finp line
				set nd [scan $line "%d %s %f %f %f %f %f %f "  id res a(3) a(4) a(5) a(6) a(7) a(8)]
			}
			break
		}
	}
	set gData(numData) 0
	for {set i $gData(xMin)} {$i <= $gData(xMax)} {incr i} {
		set x1($gData(numData)) $i
		set y1($gData(numData)) $y1($i)
		set y2($gData(numData)) $y2($i)
		incr gData(numData)
	}
}

###############################################################################
###############################################################################
#
# update histogram canvas, erase the old one and draw the elements of the
# new one. These are the low level routines
#
proc draw_histo { wind } {
	global gVar Frame gData 
	global x1 y1 y2 

	# update histogram variables
	set gVar(canvas_geometry) [wm geometry $wind] 
	set newHeight [lindex [split $gVar(canvas_geometry) x+-] 1]
	set gVar(width) [lindex [split $gVar(canvas_geometry) x+-] 0]
	set gVar(height) [expr $newHeight - 120]

	# update canvas 
	$Frame(fplotc) configure -width $gVar(width)
	$Frame(fplotc) configure -height $gVar(height)

	# erase old histogram
	set items [$Frame(fplotc) find all]
	foreach item $items {
		$Frame(fplotc) delete $item
	}

	# re-draw histogram for this canvas 
	draw_axis $Frame(fplotc) 
	draw_title $Frame(fplotc)
	draw_xScale $Frame(fplotc)
	draw_yScale $Frame(fplotc)
	draw_data $Frame(fplotc) x1 y1 $gData(numData) \
		$gVar(line1_color) $gVar(line1_shape) $gVar(line1_size)
	draw_data $Frame(fplotc) x1 y2 $gData(numData) \
		$gVar(line2_color) $gVar(line2_shape) $gVar(line2_size)
}

###############################################################################
#
# draw axis on screen 
# 
proc draw_axis { fc } {
	global gVar

	# find axis lengths
	set gVar(yLen) [expr $gVar(height) - $gVar(yStart) - $gVar(yStart)]
	set gVar(xLen) [expr $gVar(width) - $gVar(xStart) - $gVar(xStart)]

	# find yOrigin for xaxis
	set gVar(yEnd) [expr $gVar(yStart) + $gVar(yLen)]
	set gVar(xEnd) [expr $gVar(xStart) + $gVar(xLen)]
	switch $gVar(xaxisLoc) {
		top { set gVar(yOrigin) $gVar(yStart) }
		bottom { set gVar(yOrigin) $gVar(yEnd) }
		middle { set gVar(yOrigin) [expr ($gVar(yStart)+$gVar(yEnd)) / 2] }
	}

	# y axis
	$fc create line $gVar(xStart) $gVar(yStart) \
		$gVar(xStart) $gVar(yEnd) -width 1 -fill $gVar(axis_color)

	# x axis
	$fc create line $gVar(xStart) $gVar(yOrigin) \
		$gVar(xEnd) $gVar(yOrigin) -width 1 -fill $gVar(axis_color)
}

###############################################################################
#
# draw plot title 
# 
proc draw_title { fc } {
	global gVar

	set xTitle [expr $gVar(width) / 2]
	set yTitle [expr $gVar(yStart) / 2]

	$fc create text $xTitle $yTitle \
		-text $gVar(titleText) -fill $gVar(title_color) -tag movable
}

###############################################################################
#
# draw x axis scale 
# 
proc draw_xScale { fc } {
	global gVar 

	set xLabInc [expr 1 + ($gVar(xMax) - $gVar(xMin)) / $gVar(xTics)]
	set xLabStart [expr $gVar(xMin) + $xLabInc]
	set xTicTop [expr $gVar(yOrigin) - 5]
	set xTicBot [expr $gVar(yOrigin) + 5]

	if {[expr $xLabStart + $xLabInc] > $gVar(xMax)} { return }

	set xPixelRatio [expr $gVar(xLen) / double($gVar(xMax) - $gVar(xMin))]
	set xInc [expr $xLabStart - $gVar(xMin)]
	set px [expr $gVar(xStart) + $xInc * $xPixelRatio]
	set py [expr $gVar(yOrigin) + $gVar(xScaleOffset)]
	set xi $xLabStart 

	$fc create line $px $xTicTop $px $xTicBot -fill $gVar(axis_color)
	set sVal [format "%.4g" $xi]
	$fc create text $px $py -text $sVal -fill $gVar(axis_color) -tag movable

	set xEnd [expr $gVar(xStart) + $gVar(xLen)]
	set xi [expr $xi + $xLabInc]
	set px [expr $px + $xLabInc * $xPixelRatio]

	while {$px < $gVar(xEnd)} {
		$fc create line $px $xTicTop $px $xTicBot -fill $gVar(axis_color)
		set sVal [format "%.4g" $xi]
		$fc create text $px $py -text $sVal -fill $gVar(axis_color) -tag movable
		set xi [expr $xi + $xLabInc]
		set px [expr $px + $xLabInc * $xPixelRatio]
	}

	set px [expr $gVar(xLen) + $gVar(xStart)]
	set py [expr $py + 20]
	$fc create text $px $py -text $gVar(xLabText) \
		-fill $gVar(axis_color) -tag movable
}

###############################################################################
#
# draw y axis scale on graph 
# 
proc draw_yScale { fc } {
	global gVar 

	set yTicLeft [expr $gVar(xStart) - 5]
	set yTicRight [expr $gVar(xStart) + 5]

	set px [expr $gVar(xStart) - 40]

	if {$gVar(yTics) <= 0} { return }
	if {$gVar(yTics) == 1} { 
		$fc create line $yTicLeft $gVar(yOrigin) $yTicRight $yOrigin \
			-fill $aColor
		set sVal [format "%5.3f" 0.000]
		$fc create text $px $gVar(yOrigin) -text $sVal -fill $gVar(axis_color) \
			-tag movable
		return
	}

	set nTics [expr $gVar(yTics) - 1]
	set yInc [expr ($gVar(yMax) - $gVar(yMin)) / $nTics]
	set yi $gVar(yMin)
	set pInc [expr double($gVar(yLen)) / $nTics]
	for {set i 0} {$i <= $nTics} {incr i} {
		set py [expr $gVar(yStart) + $gVar(yLen) - round(0.50 + $pInc * $i)]
		$fc create line $yTicLeft $py $yTicRight $py -fill $gVar(axis_color)
		set sVal [format "%5.3f" $yi]
		$fc create text $px $py -text $sVal -fill $gVar(axis_color) -tag movable
		set yi [expr $yi + $yInc]
	}

	set py [expr $gVar(yStart) - $gVar(yLabOffset)]
	$fc create text $px $py -text $gVar(yLabText) \
		-fill $gVar(axis_color) -tag movable
}

###############################################################################
#
# draw data points on the screen 
# 
proc draw_data { fc x y nData aColor shape size } {
	global gVar

	upvar $x xd
	upvar $y yd

	if {$nData == 0} {return}

	set xPixelRatio [expr $gVar(xLen) / double($gVar(xMax) - $gVar(xMin))]
	set yPixelRatio [expr $gVar(yLen) / double($gVar(yMax) - $gVar(yMin))]

	for {set j 0} {$j < $nData} {incr j} {
		set px [expr $gVar(xStart) + ($xd($j) - $gVar(xMin)) * $xPixelRatio]
		if {$px >= $gVar(xStart) && $px <= [expr $gVar(xStart) + $gVar(xLen)]} {
			set py [expr $gVar(yOrigin) - $yd($j) * $yPixelRatio]
			$fc create $shape [expr $px-$size] [expr $py-$size] \
				[expr $px+$size] [expr $py+$size] -outline $aColor -width 1
			if [info exists oldpx] {
				set k [expr $j - 1]
				set q1 $xd($k)
				set q2 [expr $xd($j) - 1]
				if {$q1 == $q2} {
					$fc create line $oldpx $oldpy $px $py \
						-width 1 -fill $aColor
				}
			}
			set oldpx $px
			set oldpy $py
		}
	}
}

###############################################################################
#
# Allows user to drag elements of the histogram (eg, titles) to more
# appropriate places
#
proc Mark {x y w } {
	global state
	set state($w,obj) [$w find closest $x $y]
	set state($w,x) $x
	set state($w,y) $y
}
proc Drag {x y w} {
	global state
	set dx [expr $x - $state($w,x)]
	set dy [expr $y - $state($w,y)]
	$w move $state($w,obj) $dx $dy
	set state($w,x) $x
	set state($w,y) $y
}
