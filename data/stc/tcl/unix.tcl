
######################################################################
#
# Unix dependencies, tcl miscellaneous 
#
######################################################################

######################################################################
#
# remove file
#
proc remove_file { file } {
	set unixcmd "rm -f $file"
	if [ catch { eval exec $unixcmd } err] {
		return 0
	}
	return 1
}

######################################################################
#
# copy file
#
proc copy_file { from to } {
	set unixcmd "cp $from $to"
	if [ catch { eval exec $unixcmd } err] {
		return 0
	}
	return 1
}

######################################################################
#
# make dir 
#
proc make_dir { file } {
	set unixcmd "mkdir $file"
	if [ catch { eval exec $unixcmd } err] {
		return 0
	}
	return 1
}

######################################################################
# 
# get today's date
#
proc get_date { } {
	set unixcmd "date +%Y%m%d"
	if [ catch { set DATE [exec date] } err] {
		error "Problems executing unix command: $unixcmd"
	}
	return $DATE 
}

######################################################################
#
# grep file
#
proc grep_file { file string } {
	set unixcmd "grep $string $file"
	if [ catch { set ans [exec grep $string $file] } err] {
		error "Problems executing unix command: $unixcmd"
	}
	return $ans
}

######################################################################
#
# startup web browser 
#
proc start_browser { file } {
	global panvar 
	set unixcmd "$panvar(browser) $file &"
	if [ catch { eval exec $unixcmd } err] {
		return 0
	}
	return 1
}
