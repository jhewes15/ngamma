#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = ngamma_nnbarCint
SHLIB = libngamma_nnbar.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS  = -I../base
INCFLAGS += $(shell seltool-config --includes)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
LDFLAGS += $(shell seltool-config --libs)
LDFLAGS += -L$(LARLITE_LIBDIR) -lngamma_base
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE


