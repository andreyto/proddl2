########################################################################
#
# STC - structure thermodynamics calculations 
#
########################################################################
 
package require Tk
namespace eval alpha {}
tk appname Alphatk

# MacTk seems to have a bug if we exit using exit.
# see http://www.tcl.tk/software/mac/features.html
tk appname stc 
rename exit __exit
proc exit {} {
    destroy .
    # force the exit as well.
    __exit
}

##########################################################################

set VERSION		"v5.3.0"
set NAME		"stc-$VERSION"
set AQUA 		"$NAME.app/Contents/Resources/Scripts"

#
# look for installation directory if not set 
#
if [ catch {set stcVar(INSTALL) [exec printenv STC_HOME]} err] {
	set path "/Applications/$AQUA /Applications/nmr/$AQUA ~/Applications/$AQUA"
	if ![ catch {set MY_HOME [exec printenv HOME]} err] {
		set path "$MY_HOME/$AQUA $path"
		set path "$MY_HOME/Applications/nmr/$AQUA $path"
		set path "$MY_HOME/Applications/$AQUA $path"
		set path "$MY_HOME/Desktop/$AQUA $path"
	}

	# look for the software
	set stcVar(INSTALL) 0
	foreach loc $path {
		if [file isdirectory $loc] {
			set stcVar(INSTALL) $loc
			break
		}
	}

	# didn't find it
	if { $stcVar(INSTALL) == 0 } {
		error "Cannot find stc. Please place in Applications or on the desktop."
	}
}

# add tcl scripts
set stcVar(TCL_FILES) "util.tcl unix.tcl ustc.tcl init.tcl asa.tcl batch.tcl thermo.tcl histo.tcl"
foreach file $stcVar(TCL_FILES) {
	source $stcVar(INSTALL)/tcl/$file
}

wm withdraw .
stc_main
destroy .
