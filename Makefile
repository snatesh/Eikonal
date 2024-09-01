########################## BEGIN USER EDIT ##############################
# specify location of DoublyPeriodicStokes directory
export EIKONAL_ROOT    = $(PWD)
# specify desired location of shared libraries and test bin
export EIKONAL_INSTALL = $(EIKONAL_ROOT)/lib
export EIKONAL_TEST    = $(EIKONAL_ROOT)/bin
# if True, compilation will use debug mode for cpu
export DEBUG           ?= False

################################ END USER EDIT ##################################

CURDIR          = $(shell pwd)
CXX             = g++ 

ifneq ($(DEBUG), True)
  CXXFLAGS           = -I$(CURDIR)/include -w -O3 -march=native -shared -L$(EIKONAL_INSTALL)
  CXXFLAGS_test      = -I$(CURDIR)/include -w -O3 -march=native -L$(EIKONAL_INSTALL)
else
  CXXFLAGS      = -I$(CURDIR)/include -w -shared -g -O0 -DDEBUG -L$(EIKONAL_INSTALL)
endif

structureFactorSRC         = src/test.cpp 
structureFactorINC         = include/structure_factors.h include/common.h include/exceptions.h
LIBS_           = libstructureFactor.so 
LIBS            = $(patsubst %,$(EIKONAL_INSTALL)/%,$(LIBS_))

all: $(LIBS)

$(EIKONAL_INSTALL)/libstructureFactor.so: $(structureFactorSRC) $(structureFactorINC)
	@mkdir -p $(EIKONAL_INSTALL)
	@mkdir -p $(EIKONAL_TEST)
	$(CXX) -o $@ $(structureFactorSRC) $(CXXFLAGS) -fPIC
	$(CXX) -o $(EIKONAL_TEST)/test $(structureFactorSRC) $(CXXFLAGS_test)


