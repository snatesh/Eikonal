########################## BEGIN USER EDIT ##############################
# specify location of Eikonal directory
export EIKONAL_ROOT     = $(PWD)
# specify desired location of shared libraries, exec bin and test bin
export EIKONAL_LIB      = $(EIKONAL_ROOT)/lib
export EIKONAL_BIN      = $(EIKONAL_ROOT)/bin
export EIKONAL_TESTBIN     = $(EIKONAL_ROOT)/testing/bin
# if True, compilation will use debug mode for cpu
export DEBUG           ?= False

################################ END USER EDIT ##################################

CURDIR                = $(shell pwd)
CXX                   = g++ 

ifneq ($(DEBUG), True)
  CXXFLAGS_lib            = -I$(CURDIR)/include -w -O3 -march=native -shared -L$(EIKONAL_LIB)
  CXXFLAGS_bin            = -I$(CURDIR)/include -w -O3 -march=native -L$(EIKONAL_LIB)
  CXXFLAGS_test           = -I$(CURDIR)/include -I /usr/include/gtest/ -w -O3 -march=native -L$(EIKONAL_LIB)
else
  CXXFLAGS_lib            = -I$(CURDIR)/include -w -shared -g -O0 -DDEBUG -L$(EIKONAL_LIB)
  CXXFLAGS_bin            = -I$(CURDIR)/include -w -g -O0 -DDEBUG -L$(EIKONAL_LIB)
  CXXFLAGS_test           = -I$(CURDIR)/include -I /usr/include/gtest/ -w -g -O0 -DDEBUG -L$(EIKONAL_LIB)
endif

structureFactorINC        = include/structure_factors.h
jPolyINC                  = include/jPoly.h
testSRC                   = src/test.cpp
testINC                   = include/structure_factors.h include/jPoly.h
gtestSRC                  = testing/gtest.cpp
gtestINC                  = include/structure_factors.h include/jPoly.h /usr/include/gtest/gtest.h
LIBS_                     = libstructureFactor.so libjPoly.so 
gtestLIB                  = /usr/lib/x86_64-linux-gnu/libgtest.a
EXEC_                     = test 
GTEST_                    = gtest
LIBS                      = $(patsubst %,$(EIKONAL_LIB)/%,$(LIBS_))
EXEC                      = $(patsubst %,$(EIKONAL_BIN)/%,$(EXEC_))
GTEST                     = $(patsubst %,$(EIKONAL_TESTBIN)/%,$(GTEST_))

all: $(LIBS) $(EXEC) 
test: $(GTEST) 

$(EIKONAL_LIB)/libstructureFactor.so: $(structureFactorINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(structureFactorINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libjPoly.so: $(jPolyINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(jPolyINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_BIN)/test: $(testSRC) $(testINC) $(LIBS)
	@mkdir -p $(EIKONAL_BIN)
	$(CXX) -o $(EIKONAL_BIN)/test $(testSRC) $(CXXFLAGS_bin)

$(EIKONAL_TESTBIN)/gtest: $(gtestSRC) $(gtestINC) $(LIBS)
	@mkdir -p $(EIKONAL_TESTBIN)
	$(CXX) -o $(EIKONAL_TESTBIN)/gtest $(gtestSRC) $(gtestINC) $(gtestLIB) $(CXXFLAGS_test) -I/usr/include/python3.12 -l:libpython3.12.so -DWITHOUT_NUMPY
	@cd $(EIKONAL_TESTBIN) && ./gtest && open ../testdata/intconv.png

clean: 
	rm -rf $(EIKONAL_LIB) $(EIKONAL_BIN) $(EIKONAL_TESTBIN) $(EIKONAL_TESTBIN)/*.png


#g++ modern.cpp -I/usr/include/python3.12 -l:libpython3.12.so -w -DWITHOUT_NUMPY

