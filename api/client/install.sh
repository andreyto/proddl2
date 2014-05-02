#!/bin/bash

## This script will allow installing into INSTALL_PREFIX that is not
## currently in PYTHONPATH and is not configured to process .pth files.

install_prefix=$1
shift
if [ -z "$install_prefix" ]; then
	echo "Provide install prefix as the first argument" >&2
	exit 1
fi

python_exe=$1
shift
if [ -z "$python_exe" ]; then
	echo "Provide path to Python executable as the second argument" >&2
	exit 1
fi

mkdir -p "$install_prefix"
cat > "$install_prefix/sitecustomize.py"  <<End-of-message

#This file causes .pth files to be processed if the directory
##of this file is on PYTHONPATH. Placing this file in such a directory
##is the way to create an additional site-packages like directory in which
##eggs can be installed without having to edit any system or user wide Python
##config files.

import site
import os,inspect
site.addsitedir(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))

End-of-message

PYTHONPATH="$install_prefix":$PYTHONPATH INSTALL_PREFIX="$install_prefix" \
	"$python_exe" setup.py easy_install \
	--install-dir "$install_prefix" \
	--site-dirs "$install_prefix" \
	$@ \
	.

