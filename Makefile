########################## BEGIN USER EDIT ##############################
# specify location of Eikonal directory
export EIKONAL_ROOT     = $(PWD)
# specify desired location of shared libraries, exec bin and test bin
export EIKONAL_LIB      = $(EIKONAL_ROOT)/lib
export EIKONAL_BIN      = $(EIKONAL_ROOT)/bin
export EIKONAL_TESTBIN     = $(EIKONAL_ROOT)/testing/bin
# if True, compilation will use debug mode for cpu
export DEBUG           ?= False
export PLOT            ?= False
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

cblasINC                = /usr/include/x86_64-linux-gnu/openblas-openmp/cblas.h
cblasLIBDIR             = /usr/lib/x86_64-linux-gnu/openblas-openmp/
cblasLIB                = libopenblasp-r0.3.26.so
lapackINC               = /usr/include/lapacke.h
lapackLIBDIR            = /usr/lib/x86_64-linux-gnu/
lapackLIB               = liblapacke.so.3.12.0

SNOPT_INSTALL           = $(EIKONAL_ROOT)/snopt-interface/
nloptLIBDIR             = /usr/local/lib/
nloptLIB                = libnlopt.so.0.12.0
nloptINC                = /usr/local/include/nlopt.h

structureFactorINC        = include/structure_factors.h
promotionMatINC           = include/promotion_mat_tri.h
dxdyINC                   = include/dxdy_mat_tri.h
jPolyINC                  = include/jPoly.h
jacobiMatINC              = include/jacobi_mat_ON_tri.h
testSRC                   = src/test.cpp
testINC                   = include/structure_factors.h include/jPoly.h
ngjquadoptSRC              = src/ngjquad_opt.cpp
ngjquadoptINC              = $(structureFactorINC) $(jpolyINC) $(jacobiMatINC) $(nloptINC) $(lapackINC) 

gtestSRC                  = testing/gtest.cpp
gtestINC                  = include/structure_factors.h include/jPoly.h /usr/include/gtest/gtest.h
LIBS_                     = libstructureFactor.so libjPoly.so libpromotionMat.so libdxdy.so libjacobiMat.so
gtestLIB                  = /usr/lib/x86_64-linux-gnu/libgtest.a
EXEC_                     = test ngjquad_opt 
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

$(EIKONAL_LIB)/libpromotionMat.so: $(promotionMatINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(promotionMatINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libdxdy.so: $(dxdyINC) 
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(dxdyINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_LIB)/libjacobiMat.so: $(jacobiMatINC) $(structureFactorINC) $(EIKONAL_LIB)/libstructureFactor.so
	@mkdir -p $(EIKONAL_LIB)
	$(CXX) -o $@ $(jacobiMatINC) $(CXXFLAGS_lib) -fPIC

$(EIKONAL_BIN)/test: $(testSRC) $(testINC) $(LIBS)
	@mkdir -p $(EIKONAL_BIN)
	$(CXX) -o $(EIKONAL_BIN)/test $(testSRC) $(CXXFLAGS_bin)

$(EIKONAL_BIN)/ngjquad_opt: $(ngjquadoptSRC) $(ngjquadoptINC) $(lapackLIBDIR)/$(lapackLIB) $(cblasLIBDIR)/$(cblasLIB) $(nloptLIBDIR)/$(nloptLIB)
	@mkdir -p $(EIKONAL_BIN)
	$(CXX) -o $(EIKONAL_BIN)/ngjquad_opt $(ngjquadoptSRC) $(ngjquadoptINC) $(CXXFLAGS_bin) -I$(lapackINC) -L$(lapackLIBDIR) $(lapackINC) -L$(cblasLIBDIR) $(cblasINC) -L$(nloptLIBDIR) $(nloptINC) -l:$(lapackLIB) -l:$(cblasLIB) -l:$(nloptLIB) -lm -I/usr/include/python3.12 -l:libpython3.12.so -DWITHOUT_NUMPY 


#ifeq ($(PLOT), True)
#$(EIKONAL_TESTBIN)/gtest: $(gtestSRC) $(gtestINC) $(LIBS) $(cblasLIBDIR)/$(cblasLIB)
#	@mkdir -p $(EIKONAL_TESTBIN)
#	$(CXX) -o $(EIKONAL_TESTBIN)/gtest $(gtestSRC) $(gtestINC) $(gtestLIB) $(CXXFLAGS_test) -I/usr/include/python3.12 -l:libpython3.12.so -DWITHOUT_NUMPY -L$(cblasLIBDIR) $(cblasINC) -l:$(cblasLIB) 
#	@cd $(EIKONAL_TESTBIN) && ./gtest && open ../testdata/intconv.png
#else
$(EIKONAL_TESTBIN)/gtest: $(gtestSRC) $(gtestINC) $(LIBS) $(cblasLIBDIR)/$(cblasLIB) $(nloptLIBDIR)/$(nloptLIB)
	@mkdir -p $(EIKONAL_TESTBIN)
	$(CXX) -o $(EIKONAL_TESTBIN)/gtest $(gtestSRC) $(gtestINC) $(gtestLIB) $(CXXFLAGS_test) -L$(cblasLIBDIR) $(cblasINC) -L$(nloptLIBDIR) $(nloptINC) -l:$(cblasLIB) -l:$(nloptLIB) -lm 
	@cd $(EIKONAL_TESTBIN) && ./gtest 

#endif

clean: 
	rm -rf $(EIKONAL_LIB) $(EIKONAL_BIN) $(EIKONAL_TESTBIN) $(EIKONAL_TESTBIN)/*.png

snopt:
	cd $(SNOPT_INSTALL) && ./autogen.sh && ./configure --prefix=$(SNOPT_INSTALL) --with-blas="-L$(cblasLIBDIR) -l:$(cblasLIB)" --with-pic && make && make install


#g++ modern.cpp -I/usr/include/python3.12 -l:libpython3.12.so -w -DWITHOUT_NUMPY

