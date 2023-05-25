# Derek Anderson
# 03.12.2015
#
# Makefile for a variety of experiments.  If an executable needs to be linked
# HepMC or likewise, uncomment the relevant block of code and add program to
# relevant list.
#
# Last updated: 01.08.2016
#
# Courtesy of: 
#    M. Kirsanov 07.04.2006

SHELL = /bin/sh

-include config.mk
ifeq (x$(PYTHIA8LOCATION),x)
 PYTHIA8LOCATION=../..
endif
-include $(PYTHIA8LOCATION)/config.mk

# Location of directories.
TOPDIR=$(shell \pwd)
INCDIR=include
SRCDIR=src
LIBDIR=lib
LIBDIRARCH=lib/archive
BINDIR=bin
ROOTFLAGS=-I`root-config --incdir` `root-config --libs`

# Libraries to include if GZIP support is enabled
ifeq (x$(ENABLEGZIP),xyes)
LIBGZIP=-L$(BOOSTLIBLOCATION) -lboost_iostreams -L$(ZLIBLOCATION) -lz
endif

# There is no default behaviour, so remind user.
all:
	@echo "Usage: for <file>, make <file>"

# Create an executable for the following programs
mymain01  GammaPi0study VectorExperiment TriggerCount FileExperiment GammaPi0study_hQ GammaPi0study_pG \
	  MomExperiment TreeExperiment GammaJetExperiment GammaJetExecute: \
	$(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc -o $(BINDIR)/$@.exe \
	-L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 -llhapdfdummy $(LIBGZIP) $(ROOTFLAGS)
	@ln -fs $(BINDIR)/$@.exe $@.exe

# Create an executable linked to HepMC (in principle).
# Owing to excessive warning output -Wshadow is not used for HepMC.
#ifneq (x$(HEPMCLOCATION),x)
#   <file>: \
#	$(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a $(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8tohepmc.a
#	@mkdir -p $(BINDIR)
#	$(CXX) $(CXXFLAGS) -Wno-shadow -I$(PYTHIA8LOCATION)/$(INCDIR) -I$(HEPMCLOCATION)/include \
#	$@.cc -o $(BINDIR)/$@.exe \
#	-L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 -llhapdfdummy $(LIBGZIP) \
#	-lpythia8tohepmc \
#	-L$(HEPMCLOCATION)/lib -lHepMC
#	@ln -fs $(BINDIR)/$@.exe $@.exe
#else
#   <file>:
#	@echo ERROR, this target needs HepMC, variable HEPMCLOCATION
#endif 

# Create an executable that links to LHAPDF
#<file>: $(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a
#	@mkdir -p $(BINDIR)
#	$(CXX) $(CXXFLAGS) -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc -o $(BINDIR)/$@.exe \
#	-L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 $(LIBGZIP) \
#	-L$(LHAPDFLOCATION) $(LHAPDFLIBNAME) \
#	$(FLIBS)
#	@ln -fs $(BINDIR)/$@.exe $@.exe

# Create an executable that links to LHAPDF and HepMC
#<file>: \
#	$(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a $(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8tohepmc.a
#	@mkdir -p $(BINDIR)
#	$(CXX) $(CXXFLAGS) -Wno-shadow -I$(PYTHIA8LOCATION)/$(INCDIR) -I$(HEPMCLOCATION)/include \
#	$@.cc -o $(BINDIR)/$@.exe \
#	-L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 -lpythia8tohepmc $(LIBGZIP) \
#	-L$(LHAPDFLOCATION) $(LHAPDFLIBNAME) \
#	-L$(HEPMCLOCATION)/lib -lHepMC \
#	$(FLIBS)
#	@ln -fs $(BINDIR)/$@.exe $@.exe

# Create an executable that links to Fastjet
# Owing to excessive warning output -Wshadow is not used for Fastjet.
# (Fixed as of Fastjet 3.0.1, so will be modified eventually.)
#ifneq (x$(FASTJETLOCATION),x)
#   <file>: $(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a
#	@mkdir -p $(BINDIR)
#       # Note: $(CXXFLAGS) is after Fastjet flags as Fastjet includes
#       # optimisation/debug flags which may be unwanted (e.g. -g -O2)
#	$(CXX) -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc \
#	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
#	$(CXXFLAGS) -Wno-shadow \
#	-o $(BINDIR)/$@.exe \
#	-L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 -llhapdfdummy $(LIBGZIP) \
#	-L$(FASTJETLOCATION)/lib \
#	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`
#	@ln -fs $(BINDIR)/$@.exe $@.exe
#	@rm -f $@.o
#else
#   <file>:
#	@echo ERROR, this target needs Fastjet, variable FASTJETLOCATION
#endif

# Create an executable that links to Fastjet, HepMC and LHApdf
# Owing to excessive warning output -Wshadow is not used for Fastjet.
# (Fixed as of Fastjet 3.0.1, so will be modified eventually.)
#ifneq (x$(FASTJETLOCATION),x)
#   <file>: \
#	$(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a  $(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8tohepmc.a
#	@mkdir -p $(BINDIR)
#       # Note: $(CXXFLAGS) is after Fastjet flags as Fastjet includes
#       # optimisation/debug flags which may be unwanted (e.g. -g -O2)
#	$(CXX) -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc \
#	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
#	$(CXXFLAGS) -Wno-shadow \
#	-I$(PYTHIA8LOCATION)/$(INCDIR) -I$(HEPMCLOCATION)/include \
#	-o $(BINDIR)/$@.exe \
#	-L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 \
#	-L$(LHAPDFLOCATION) $(LHAPDFLIBNAME) \
#	-lpythia8tohepmc \
#	-L$(HEPMCLOCATION)/lib -lHepMC \
#	-L$(FASTJETLOCATION)/lib \
#	-L$(LHAPDFLOCATION)/lib \
#	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`
#	@ln -fs $(BINDIR)/$@.exe $@.exe
#	@rm -f $@.o
#else
#   <file>:
#	@echo ERROR, this target needs Fastjet, variable FASTJETLOCATION
#endif

# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -rf $(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core*

# End
