########################## BEGIN USER EDIT ##############################
# specify location of DoublyPeriodicStokes directory
export EIKONAL_ROOT    = $(PWD)
# specify desired location of shared libraries and test bin
export EIKONAL_LIB = $(EIKONAL_ROOT)/lib
export EIKONAL_BIN    = $(EIKONAL_ROOT)/bin
# if True, compilation will use debug mode for cpu
export DEBUG           ?= False

################################ END USER EDIT ##################################

CURDIR                = $(shell pwd)
CXX                   = g++ 

ifneq ($(DEBUG), True)
  CXXFLAGS_lib            = -I$(CURDIR)/include -w -O3 -march=native -shared -L$(EIKONAL_LIB)
  CXXFLAGS_bin       = -I$(CURDIR)/include -w -O3 -march=native -L$(EIKONAL_LIB)
else
  CXXFLAGS_lib            = -I$(CURDIR)/include -w -shared -g -O0 -DDEBUG -L$(EIKONAL_LIB)
  CXXFLAGS_bin            = -I$(CURDIR)/include -w -g -O0 -DDEBUG -L$(EIKONAL_LIB)
endif

structureFactorINC    = include/structure_factors.h
jPolyINC              = include/jPoly.h
testSRC               = src/test.cpp
testINC               = include/structure_factors.h include/jPoly.h
LIBS_                 = libstructureFactor.so libjPoly.so
EXEC_                 = test 
LIBS                  = $(patsubst %,$(EIKONAL_LIB)/%,$(LIBS_))
EXEC                  = $(patsubst %,$(EIKONAL_BIN)/%,$(EXEC_))

all: $(LIBS) $(EXEC)

$(EIKONAL_LIB)/libstructureFactor.so: $(structureFactorINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(structureFactorINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libjPoly.so: $(jPolyINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(jPolyINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_BIN)/test: $(testSRC) $(testINC) $(EIKONAL_LIB)/libstructureFactor.so $(EIKONAL_LIB)/libjPoly.so
	@mkdir -p $(EIKONAL_BIN)
	$(CXX) -o $(EIKONAL_BIN)/test $(testSRC) $(CXXFLAGS_bin)

clean: 
	rm -rf $(EIKONAL_LIB)/lib*.so $(EIKONAL_BIN)/*

