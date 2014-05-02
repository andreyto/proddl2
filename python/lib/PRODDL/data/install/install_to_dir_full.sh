#!/bin/bash
set -e
target_dir=$1
[ -n "$target_dir" ] || exit 1
mkdir -p "$target_dir"
rc="PRODDL.rc"
python lib/PRODDL/data/install/helpers.py gen-rc install-dir "$target_dir" $rc
. "$rc"
python setup.py easy_install --install-dir "$target_dir" --always-copy .
#python setup.py easy_install --install-dir "$target_dir" .
mv "$rc" "$target_dir/"

