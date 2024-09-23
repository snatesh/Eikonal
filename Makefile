########################## BEGIN USER EDIT ##############################
# specify location of Eikonal directory
export EIKONAL_ROOT     = $(PWD)
# specify desired location of shared libraries, exec bin and test bin
export EIKONAL_LIB      = $(EIKONAL_ROOT)/lib
export EIKONAL_BIN      = $(EIKONAL_ROOT)/bin
export EIKONAL_TESTBIN     = $(EIKONAL_ROOT)/testing/bin
# if True, compilation will use debug mode for cpu
export DEBUG           ?= True
export PLOT            ?= False
################################ END USER EDIT ##################################

CURDIR                = $(shell pwd)
CXX                   = g++ 

ifneq ($(DEBUG), True)
  CXXFLAGS_lib            = -I$(CURDIR)/include -Ofast -w -march=native -shared -L$(EIKONAL_LIB) -ftree-vectorize -fopt-info-vec-all  -fopenmp  
  CXXFLAGS_bin            = -I$(CURDIR)/include -w -O3 -march=native -L$(EIKONAL_LIB) -fopenmp 
  CXXFLAGS_test           = -I$(CURDIR)/include -I /usr/include/gtest/ -w -O3 -march=native -L$(EIKONAL_LIB)
else
  CXXFLAGS_lib            = -I$(CURDIR)/include -w -shared -g -O0 -DDEBUG -L$(EIKONAL_LIB)
  CXXFLAGS_bin            = -I$(CURDIR)/include -w -g -O0 -DDEBUG -L$(EIKONAL_LIB)
  CXXFLAGS_test           = -I$(CURDIR)/include -I /usr/include/gtest/ -w -g -O0 -DDEBUG -L$(EIKONAL_LIB)
endif

cblasINC                = /usr/include/x86_64-linux-gnu/openblas-openmp/cblas.h
cblasLIBDIR             = /usr/lib/x86_64-linux-gnu/openblas-openmp/
cblasLIB                = libopenblasp-r0.3.26.so
lapackINC               = /usr/include/lapacke.h
lapackLIBDIR            = /usr/lib/x86_64-linux-gnu/
lapackLIB               = liblapacke.so.3.12.0

nloptLIBDIR             = /usr/local/lib/
nloptLIB                = libnlopt.so.0.12.0
nloptINC                = /usr/local/include/nlopt.h

sFactorsINC               = include/sFactors.h
kMatINC                   = include/kMat.h
dMatINC                   = include/dMat.h
jPolyINC                  = include/jPoly.h
jMatINC                   = include/jMat.h
jevdINC                   = include/jevd.h include/timer.h
ngjquadINC                = include/ngjquad.h

tetquadSRC                = src/tetquad.cpp
tetquadINC                = $(jevdINC) 

ngjquadSRC                = src/ngjquad.cpp
ngjquadINC                = $(sFactorsINC) $(jpolyINC) $(jMatINC) $(nloptINC) $(lapackINC) include/ngjquad.h
  
gtestSRC                  = testing/gtest.cpp
gtestINC                  = include/sFactors.h include/jPoly.h /usr/include/gtest/gtest.h
LIBS_                     = libsFactors.so libjPoly.so libkMat.so libdMat.so libjMat.so libjevd.so libngjquad.so
gtestLIB                  = /usr/lib/x86_64-linux-gnu/libgtest.a
EXEC_                     = ngjquad
GTEST_                    = gtest
LIBS                      = $(patsubst %,$(EIKONAL_LIB)/%,$(LIBS_))
EXEC                      = $(patsubst %,$(EIKONAL_BIN)/%,$(EXEC_))
GTEST                     = $(patsubst %,$(EIKONAL_TESTBIN)/%,$(GTEST_))

all: $(LIBS) $(EXEC) 
test: $(GTEST) 

$(EIKONAL_LIB)/libsFactors.so: $(sFactorsINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(sFactorsINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libjPoly.so: $(jPolyINC) $(sFactorsINC) $(EIKONAL_LIB)/libsFactors.so 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(jPolyINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libkMat.so: $(kMatINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(kMatINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libdMat.so: $(dMatINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(dMatINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libjMat.so: $(jMatINC) $(sFactorsINC) $(EIKONAL_LIB)/libsFactors.so
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(jMatINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libjevd.so: $(jevdINC)
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(jevdINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libngjquad.so: $(ngjquadINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(ngjquadINC) $(CXXFLAGS_lib) -fPIC -lm 


$(EXEC): $(ngjquadSRC) $(ngjquadINC) $(LIBS) $(lapackLIBDIR)/$(lapackLIB) $(cblasLIBDIR)/$(cblasLIB) $(nloptLIBDIR)/$(nloptLIB)
	@mkdir -p $(EIKONAL_BIN)
	$(CXX) -o $(EXEC) $(ngjquadSRC) $(ngjquadINC) $(CXXFLAGS_bin) -I$(lapackINC) -L$(lapackLIBDIR) $(lapackINC) -L$(cblasLIBDIR) $(cblasINC) -L$(nloptLIBDIR) $(nloptINC) -l:$(lapackLIB) -l:$(cblasLIB) -l:$(nloptLIB) -lm -I/usr/include/python3.12 -l:libpython3.12.so -DWITHOUT_NUMPY  

#$(EXEC): $(tetquadSRC) $(tetquadINC) $(EIKONAL_LIB)/libjevd.so $(cblasLIBDIR)/$(cblasLIB)
#	@mkdir -p $(EIKONAL_BIN)
#	$(CXX) -o $(EXEC) $(tetquadSRC) $(tetquadINC) $(CXXFLAGS_bin) -L$(cblasLIBDIR) $(cblasINC) -l:$(cblasLIB) 
#
$(EIKONAL_TESTBIN)/gtest: $(gtestSRC) $(gtestINC) $(LIBS) $(cblasLIBDIR)/$(cblasLIB) $(nloptLIBDIR)/$(nloptLIB)
	@mkdir -p $(EIKONAL_TESTBIN)
	$(CXX) -o $(EIKONAL_TESTBIN)/gtest $(gtestSRC) $(gtestINC) $(gtestLIB) $(CXXFLAGS_test) -L$(cblasLIBDIR) $(cblasINC) -L$(nloptLIBDIR) $(nloptINC) -l:$(cblasLIB) -l:$(nloptLIB) -lm 
	@cd $(EIKONAL_TESTBIN) && ./gtest 

clean: 
	rm -rf $(EIKONAL_LIB) $(EIKONAL_BIN) $(EIKONAL_TESTBIN) $(EIKONAL_TESTBIN)/*.png

snopt:
	cd $(SNOPT_INSTALL) && ./autogen.sh && ./configure --prefix=$(SNOPT_INSTALL) --with-blas="-L$(cblasLIBDIR) -l:$(cblasLIB)" --with-pic && make && make install


#g++ modern.cpp -I/usr/include/python3.12 -l:libpython3.12.so -w -DWITHOUT_NUMPY

