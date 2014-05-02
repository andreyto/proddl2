./configure --prefix $PRODDL_PREFIX -with-cxx=gcc
make
## this take a while to build; some tests fail
#make check-testsuite
make install
