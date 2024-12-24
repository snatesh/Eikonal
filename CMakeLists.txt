cmake_minimum_required(VERSION 3.6)
project(Eikonal CXX)
option(debug "Build with debugger" ON)
option(build_testing "compile and run tests" ON)
enable_language(C)

# lib and bin will be in top source
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
# set runtime path for shared libraries
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")


if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.5.0")
    message(FATAL_ERROR "Insufficient gcc version. Must be 7.5.0 or greater")
  endif()
else()
  message(FATAL_ERROR "Compiler ID must be GNU (Intel support coming soon!)")
endif()

if(NOT debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -Iinclude -O3 -march=native")
elseif(debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Iinclude -O3 -g -DDEBUG -march=native")
endif()

if (build_testing)
  enable_testing()
  add_subdirectory(testing)
endif()

set(BLA_VENDOR OpenBLAS)
find_package(OpenMP REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(VTK REQUIRED)
#find_package(VTK COMPONENTS 
#  CommonColor
#  CommonCore
#  IOImage
#  InteractionImage
#  InteractionStyle
#  RenderingContextOpenGL2
#  RenderingCore
#  RenderingFreeType
#  RenderingGL2PSOpenGL2
#  RenderingOpenGL2
#)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "JPEGReader: Unable to find the VTK build folder.")
endif()


include_directories(include)
include_directories(${Python_INCLUDE_DIRS})

# setting lib src files
set(sfactorSRC include/sFactors.hh wrapper/sFactorsWrapper.cpp)
set(kmatSRC include/kMat.hh wrapper/kMatWrapper.cpp)
set(dmatSRC include/dMat.hh wrapper/dMatWrapper.cpp)
set(lmatSRC include/lMat.hh wrapper/lMatWrapper.cpp)
set(jpolySRC include/jPoly.hh wrapper/jPolyWrapper.cpp)
set(tensorquadSRC include/mapTensorQuad.hh wrapper/tensorQuadWrapper.cpp)
set(jweightSRC include/jWeight.hh wrapper/jWeightWrapper.cpp)
set(jmatSRC include/jMat.hh wrapper/jMatWrapper.cpp)
set(jevdSRC include/jevd.hh wrapper/jevdWrapper.cpp)
#set(ngjquadSRC include/ngjquad.hh wrapper/ngjquadWrapper.cpp)
set(ngjquadSRC include/ngjquad_redobj.hh)
set(legquadSRC include/legQuad.hh wrapper/legQuadWrapper.cpp)
set(eikonalSRC include/eikonal.hh wrapper/eikonalWrapper.cpp)
set(tetquadSRC src/tetquad.cpp)
#set(gjquadoptSRC src/ngjquad.cpp)
set(gjquadoptSRC src/ngjquad_redobj.cpp)
set(EikonalSRC src/eikonal.cpp)
set(EikonalWeightedSRC src/eikonal_weighted.cpp)
set(nullspaceSRC  src/nullspace.cpp)
#set(eikonalSRC src/eikonal.cpp)
set(readjpgsrc src/readjpg.cpp)
set(readnefsrc src/readnef.cpp)

# set build output paths
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#building lib
add_library(sfactors SHARED ${sfactorSRC})
set_source_files_properties(${sfactorSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
target_link_libraries(sfactors gomp) 

add_library(kmat SHARED ${kmatSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${kmatSRC} PROPERTIES COMPILE_FLAGS "-fPIC")
else()
  set_source_files_properties(${kmatSRC} PROPERTIES COMPILE_FLAGS "-xHost -fpic")
endif()

add_library(dmat SHARED ${dmatSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${dmatSRC} PROPERTIES COMPILE_FLAGS "-fPIC")
else()
  set_source_files_properties(${dmatSRC} PROPERTIES COMPILE_FLAGS "-xHost -fpic ")
endif()

add_library(lmat SHARED ${lmatSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${lmatSRC} PROPERTIES COMPILE_FLAGS "-fPIC")
else()
  set_source_files_properties(${lmatSRC} PROPERTIES COMPILE_FLAGS "-xHost -fpic ")
endif()

add_library(tensorquad SHARED ${tensorquadSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${tensorquadSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
  target_link_libraries(tensorquad gomp)
else()
  set_source_files_properties(${tensorquadSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic")
  target_link_libraries(tensorquad ${OpenMP_iomp5_LIBRARY} ${OpenMP_pthread_LIBRARY})
endif()

add_library(jpoly SHARED ${jpolySRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${jpolySRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${jpolySRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(jpoly sfactors openblas)

add_library(jweight SHARED ${jweightSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${jweightSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${jweightSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(jweight tensorquad) 

add_library(jmat SHARED ${jmatSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${jmatSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${jmatSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(jmat jpoly jweight openblas) 

add_library(jevd SHARED ${jevdSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${jevdSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${jevdSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(jevd openblas gomp) 

add_library(ngjquad SHARED ${ngjquadSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${ngjquadSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${ngjquadSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(ngjquad jmat nlopt lapacke) 

add_library(legquad SHARED ${legquadSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${legquadSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${legquadSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(legquad ngjquad) 

add_library(eikonal SHARED ${eikonalSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${eikonalSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${eikonalSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(eikonal nlopt legquad) 

set_target_properties(sfactors kmat dmat jmat jweight eikonal
                      jpoly ngjquad legquad jevd tensorquad 
                      PROPERTIES LINKER_LANGUAGE CXX)

add_executable(gjquadopt ${gjquadoptSRC})
set_target_properties(gjquadopt PROPERTIES COMPILE_FLAGS "-DWITHOUT_NUMPY")
#target_link_libraries(gjquadopt dmat legquad ${Python_LIBRARIES})
target_link_libraries(gjquadopt dmat ngjquad ${Python_LIBRARIES})

add_executable(tetquad ${tetquadSRC})
target_link_libraries(tetquad jevd legquad)   

add_executable(Eikonal ${EikonalSRC})
set_target_properties(Eikonal PROPERTIES COMPILE_FLAGS "-DWITHOUT_NUMPY")
target_link_libraries(Eikonal eikonal)  

add_executable(nullspace ${nullspaceSRC})
target_link_libraries(nullspace openblas lapacke)  

add_executable(EikonalW ${EikonalWeightedSRC})
set_target_properties(EikonalW PROPERTIES COMPILE_FLAGS "-DWITHOUT_NUMPY")
target_link_libraries(EikonalW eikonal)  

add_executable(readjpg ${readjpgsrc})
target_link_libraries(readjpg ${VTK_LIBRARIES})

add_executable(readnef ${readnefsrc})
target_link_libraries(readnef ${VTK_LIBRARIES})

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS readjpg readnef
  MODULES ${VTK_LIBRARIES}
)

install(TARGETS sfactors ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS kmat ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS dmat ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS lmat ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS jmat ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS jweight ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS jpoly ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS ngjquad ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS legquad ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS jevd ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS tensorquad ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS eikonal ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS gjquadopt ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
install(TARGETS tetquad ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
install(TARGETS Eikonal ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
install(TARGETS nullspace ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
install(TARGETS EikonalW ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
install(TARGETS readjpg ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
install(TARGETS readnef ARCHIVE DESTINATION bin LIBRARY DESTINATION bin)
