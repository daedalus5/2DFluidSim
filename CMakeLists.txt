cmake_minimum_required(VERSION 3.5)
# Find eigen in system directory and treat it as an interface library
find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)# Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_cispba_executable(2DFluidSim KaminoCoreSolver/main.cpp
                        KaminoCoreSolver/include/KaminoQuantity.h
                        KaminoCoreSolver/src/KaminoQuantity.cpp
                        KaminoCoreSolver/src/KaminoSolver.cpp)
        
target_include_directories(2DFluidSim SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
set(CMAKE_CXX_FLAGS "-O3")
target_link_libraries(2DFluidSim partio)
