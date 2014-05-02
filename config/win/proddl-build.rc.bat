set PRODDL_BUILD_TOP=C:\proddl
set PRODDL_DEPS=%PRODDL_BUILD_TOP%\deps
set PRODDL_DEPS_BUILD=%PRODDL_BUILD_TOP%\deps_build
set PDODDL_SWIG_ROOT=%PRODDL_DEPS_BUILD%\swigwin-2.0.9
set PDODDL_DOXYGEN_ROOT=%PRODDL_DEPS_BUILD%\doxygen

set BOOST_ROOT=%PRODDL_DEPS%\boost_1_44

set PRODDL_SRC=%PRODDL_BUILD_TOP%\proddl2

set IMP_HOME=%PRODDL_DEPS%\IMP-2.1.1

REM adding cygwin only to have nice linux utils for development
REM adding VS path to have devenv in the path for starting debugger from command line
PATH=%PATH%;%PDODDL_SWIG_ROOT%;%PRODDL_DOXYGEN_ROOT%\bin;%PRODDL_DEPS%\bin;%PRODDL_DEPS%\HDF5_1.8.10\bin;%IMP_HOME%\bin;C:\cygwin32\bin;"C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE"

set PYTHONPATH=%PYTHONPATH%;%IMP_HOME%\python;%IMP_HOME%\bin
