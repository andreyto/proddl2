#
#	On BSD machines, RANLIB should be 'ranlib'
#
#	On System V machines, RANLIB should be ':'
#
SHELL		= /bin/sh
RANLIB		= :

INSTALL         = install --owner=root --mode=644

PREFIX          = /usr/local

CXX		= g++ -I/home/andrey/science/src/python/include/
#CXX		= CC

.SUFFIXES:	.cc

.cc.o:
	$(CXX) -g -c $< -o $@

OPT		= -g
DEBUG		=
CFLAGS		= $(OPT) $(DEBUG)
LIBRARY		= pdb++

HDRS            = pdb++.hpp

LIBARCH		= lib$(LIBRARY).a
OBJS		= pdb_read.o pdb_sprntf.o pdb_sscanf.o pdb_chars.o \
		pdb_type.o pdb++.o pdbinput.o
SRCS		= pdb_read.cc pdb_sprntf.cc pdb_sscanf.cc pdb_chars.cc \
		pdb_type.cc pdb++.cc pdbinput.cc

all:		$(LIBARCH) test

install:	$(LIBARCH)
		$(INSTALL) $(LIBARCH) $(PREFIX)/lib
		$(INSTALL) $(HDRS) $(PREFIX)/include 

$(LIBARCH):     $(OBJS)
		@echo "Loading $(LIBARCH) ... "
		@ar cru $(LIBARCH) $(OBJS)
		@$(RANLIB) $(LIBARCH)
		@echo "done"

test:           $(LIBARCH)
		$(CXX) $(OPT) -o test_pdb++ test_pdb++.cpp $(LIBARCH) -lgzstream -lz

clean:;		@rm -f $(OBJS)

spotless:;	@rm -f $(OBJS) $(LIBARCH)

tags:           $(HDRS) $(SRCS); @ctags -w $(HDRS) $(SRCS)
