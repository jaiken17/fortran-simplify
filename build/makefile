# Directories break if there area any spaces in paths
.SUFFIXES:
PROJDIR := $(realpath $(CURDIR)/..)
SOURCEDIR = $(PROJDIR)/src
BUILDDIR = $(PROJDIR)/build
VPATH = $(SOURCEDIR)
TARGET = polytest

# Compiler defaults to gfortran
FORT ?= gfortran

# Flags are compiler dependent
ifeq ($(FORT), gfortran)
FLAGS = -Wall -Warray-bounds -std=f2008
endif
ifeq ($(FORT), ifort)
FLAGS = -warn -WB -stand 
endif

# Compilation commands for modules and main program
COMPOBJ = $(FORT) $(FLAGS) -c -o
COMPMAIN = $(FORT) $(FLAGS) -o

# Source files -> order of compilation is implicit from left to right
MAINSOURCE = polytest.f90
SOURCES = precision.f90 io.f90 simplify.f90


$(BUILDDIR)/$(TARGET): $(SOURCEDIR)/$(MAINSOURCE) $(subst .f90,.out, $(SOURCES))
	$(COMPMAIN) $(TARGET) $< $(subst .f90,.out, $(SOURCES))

%.out: $(SOURCEDIR)/%.f90
	$(COMPOBJ) $*.out $< 


.PHONY: clean
clean:
	rm -f *.out polytest *.mod