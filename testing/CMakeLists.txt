find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(../include)
include_directories(${Python_INCLUDE_DIRS})
add_executable(runGtest gtest.cpp)
target_link_libraries(runGtest gtest gtest_main sfactors 
                      kmat dmat jmat jweight jpoly ngjquad legquad 
                      jevd tensorquad)
add_test(runGtest runGtest)
