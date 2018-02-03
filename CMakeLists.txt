cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)

project(phenotyping-genmodels3d)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(PCL 1.7 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/include")
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
set(PCL_BUILD_TYPE Release)

# the code was built using libvtk6 libraries
# there might be incompatibility issues when using libvtk4 libraries
# especially with the setInputDataObject()/setInputData()/setInput() functions

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# CMAES
include_directories("/usr/local/include/libcmaes")

set(srcs
	"src/main.cpp"
	"src/pcl-utils.cpp"
	"src/PhytomerGeneration.cpp"
	"src/PointsFromSurfaceSampler.cpp"
	"src/PairwiseRegistration.cpp"
	"src/CloudErrorMetrics.cpp"
	"src/CrossEntropyOptimization.cpp")

set(incs
	"include/pcl-utils.h"
	"include/EigenMultivariateNormal.h"
	"include/PhytomerGeneration.h"
	"include/PointsFromSurfaceSampler.h"
	"include/PairwiseRegistration.h"
	"include/CloudErrorMetrics.h"
	"include/CrossEntropyOptimization.h")

add_executable(main ${srcs})
target_link_libraries(main ${PCL_LIBRARIES} "/usr/local/lib/libcmaes.so")