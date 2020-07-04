# Micromorphic Elasto Plasticity, a implementation of a finite deformation Elasto-Plastic
# Material model in micromorphic theory.
#
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the common configuration file for all of the included makefiles

# Check for icc compiler or default to g++
ICC_EXIST:=$(shell which icc)
ifdef ICC_EXIST
    CXX=icc
    CFLAGS=-std=c++11 -Wall -ansi -pedantic -O3 -I. -fmax-errors=5
else
    CXX=g++
    CFLAGS=-std=gnu++11 -Wall -ansi -pedantic -O3 -I. -fmax-errors=5
endif

# Location of the Eigen library
EIGEN=-I$(abspath $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))/../eigen)

# Set the root directory
ROOTDIR:=$(abspath $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))/../)

# Add the location of the error_tools to the include and library
ERRORSOURCE = $(ROOTDIR)/error_tools/src/cpp/error_tools.cpp
ERRORHEADER = $(ROOTDIR)/error_tools/src/cpp/error_tools.h
INC=-I$(ROOTDIR)/error_tools/src/cpp
LIB=-L$(ROOTDIR)/error_tools/src/cpp

# Add the location of the vector_tools to the include and library
VECTORSOURCE = $(ROOTDIR)/vector_tools/src/cpp/vector_tools.cpp
VECTORHEADER = $(ROOTDIR)/vector_tools/src/cpp/vector_tools.h
INC+=-I$(ROOTDIR)/vector_tools/src/cpp
LIB+=-L$(ROOTDIR)/vector_tools/src/cpp

# Add the location of constitutive_tools to the include and library
CONSTITUTIVESOURCE = $(ROOTDIR)/constitutive_tools/src/cpp/constitutive_tools.cpp
CONSTITUTIVEHEADER = $(ROOTDIR)/constitutive_tools/src/cpp/constitutive_tools.h
INC+=-I$(ROOTDIR)/constitutive_tools/src/cpp
LIB+=-L$(ROOTDIR)/constitutive_tools/src/cpp

# Add the location of micromorphic_tools to the include and library
MICROMORPHICTOOLSSOURCE = $(ROOTDIR)/micromorphic_tools/src/cpp/micromorphic_tools.cpp
MICROMORPHICTOOLSHEADER = $(ROOTDIR)/micromorphic_tools/src/cpp/micromorphic_tools.h
INC+=-I$(ROOTDIR)/micromorphic_tools/src/cpp
LIB+=-L$(ROOTDIR)/micromorphic_tools/src/cpp

# Add the location of micromorphic_linear_elasticity to the include and library
MICROMORPHICLINEARELASTICITYSOURCE = $(ROOTDIR)/micromorphic_linear_elasticity/src/cpp/micromorphic_linear_elasticity.cpp
MICROMORPHICLINEARELASTICITYHEADER = $(ROOTDIR)/micromorphic_linear_elasticity/src/cpp/micromorphic_linear_elasticity.h
INC+=-I$(ROOTDIR)/micromorphic_linear_elasticity/src/cpp
LIB+=-L$(ROOTDIR)/micromorphic_linear_elasticity/src/cpp

# Add the location of solver_tools to the include and library
SOLVERSOURCE = $(ROOTDIR)/solver_tools/src/cpp/solver_tools.cpp
SOLVERHEADER = $(ROOTDIR)/solver_tools/src/cpp/solver_tools.h
INC+=-I$(ROOTDIR)/solver_tools/src/cpp
LIB+=-L$(ROOTDIR)/solver_tools/src/cpp

# Add the location of the micromorphic material library
MICROMORPHICMATERIALLIBSOURCE = $(ROOTDIR)/micromorphic_element/src/cpp/micromorphic_material_library.cpp
MICROMORPHICMATERIALLIBHEADER = $(ROOTDIR)/micromorphic_element/src/cpp/micromorphic_material_library.h
INC+=-I$(ROOTDIR)/micromorphic_element/src/cpp
LIB+=-L$(ROOTDIR)/micromorphic_element/src/cpp

# The python command
PYTHON=/apps/anaconda3/bin/python
