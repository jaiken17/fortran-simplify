.SUFFIXES:
PROJDIR := $(realpath $(CURDIR)/..)
SOURCEDIR = $(PROJDIR)/src
BUILDDIR = $(PROJDIR)/build
VPATH = $(SOURCEDIR)
TARGET = polytest

FORT ?= gfortran
COMPOBJ = $(FORT) -c -o
COMPMAIN = $(FORT) -o
MAINSOURCE = polytest.f90
SOURCES = $(BSPLINE) precision.f90 polysimplify.f90


$(BUILDDIR)/$(TARGET): $(SOURCEDIR)/$(MAINSOURCE) $(subst .f90,.out, $(SOURCES))
	$(COMPMAIN) $(TARGET) $< $(subst .f90,.out, $(SOURCES))

%.out: $(SOURCEDIR)/%.f90
	$(COMPOBJ) $*.out $< 


.PHONY: clean
clean:
	rm -f *.out polytest *.mod