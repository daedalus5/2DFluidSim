cmake_minimum_required(VERSION 3.5)
# Find eigen in system directory and treat it as an interface library
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
### Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_cispba_executable(2DFluidSim
						KaminoCoreSolver/main.cpp
                        KaminoCoreSolver/include/Kamino.h
                        KaminoCoreSolver/include/KaminoQuantity.h
			            KaminoCoreSolver/include/CubicSolver.h
                        KaminoCoreSolver/src/Kamino.cpp
			            KaminoCoreSolver/src/KaminoSolver.cpp
			            KaminoCoreSolver/src/CubicSolver.cpp
			            KaminoCoreSolver/src/KaminoAdvection.cpp
			            KaminoCoreSolver/src/KaminoFFT.cpp
			            KaminoCoreSolver/src/KaminoGeometry.cpp
			            KaminoCoreSolver/src/KaminoInit.cpp
			            KaminoCoreSolver/src/KaminoProjection.cpp
                        KaminoCoreSolver/src/KaminoQuantity.cpp
                        KaminoCoreSolver/src/KaminoParticles.cpp
                        )
        
target_include_directories(2DFluidSim SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(2DFluidSim ${OpenCV_LIBS})
set(CMAKE_CXX_FLAGS "-O3")
target_link_libraries(2DFluidSim partio)
