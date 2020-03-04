# Micromorphic elasto-plasticity library
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   March 3, 2020
#
# This is the makefile for test_micromorphic_elasto_plasticity.cpp it tests the 
# implementation of the micromorphic elasto-plastic constitutive model.

include ../../../../config.mk

# List of executables
EXECUTABLE=test_micromorphic_elasto_plasticity
TESTOBJECT=test_micromorphic_elasto_plasticity.o
TESTSOURCE=test_micromorphic_elasto_plasticity.cpp

OBJECT=micromorphic_elasto_plasticity.o
SOURCE=../../micromorphic_elasto_plasticity.cpp
HEADER=../../micromorphic_elasto_plasticity.h

MICROMORPHICLINEARELASTICITYOBJECT=micromorphic_linear_elasticity.o
MICROMORPHICTOOLSOBJECT=micromorphic_tools.o
CONSTITUTIVEOBJECT=constitutive_tools.o
ERROROBJECT=error_tools.o

OBJECTS = $(TESTOBJECT)
OBJECTS += $(OBJECT) $(ERROROBJECT) $(CONSTITUTIVEOBJECT) $(MICROMORPHICTOOLSOBJECT) $(MICROMORPHICLINEARELASTICITYOBJECT)

INC+=-I../..
LOCAL_CFLAGS=-ggdb

#Set the mode due to the method variable
DBG =
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

# Makefile rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=c++11

$(TESTOBJECT): $(TESTSOURCE) $(HEADER) $(TESTHEADER) $(MICROMORPHICLINERELASTICITYHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $< -c -o $@ $(CFLAGS) -lstdc++ -std=c++11

$(OBJECT): $(SOURCE) $(HEADER) $(MICROMORPHICTOOLSHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(MICROMORPHICLINEARELASTICITYOBJECT): $(MICROMORPHICLINEARELASTICITYSOURCE) $(MICROMORPHICLINEARELASTICITYHEADER) $(MICROMORPHICTOOLSHEADER) $(CONSTITUTIVEHEADER) $(VECTORHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(MICROMORPHICTOOLSOBJECT): $(MICROMORPHICTOOLSSOURCE) $(MICROMORPHICTOOLSHEADER) $(CONSITUTIVEHEADER) $(VECTORHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(CONSTITUTIVEOBJECT): $(CONSTITUTIVESOURCE) $(CONSTITUTIVEHEADER) $(VECTORHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11


clean:
	rm -f $(EXECUTABLE)
	rm $(OBJECTS)

.PHONY: all clean
