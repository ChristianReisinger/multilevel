cmake_minimum_required(VERSION 3.10)
project(multilevel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -std=c++17")

set(NC 2 CACHE STRING "number of colors for SU(NC)")

add_definitions(-D__SUN_N__=${NC})

if(NC EQUAL 2)
	find_library(SUN_INTERFACE_LIB MCSU2Gaugefield "${LatticeTools_ROOT}/lib")
	set(SUN_ROOT "${SU2_ROOT}")
	
elseif(NC EQUAL 3)
	find_library(SUN_INTERFACE_LIB CL2QCDGaugefield "${LatticeTools_ROOT}/lib")
	
	find_package(Boost 1.59.0 REQUIRED COMPONENTS regex filesystem system program_options unit_test_framework)
	find_package(OpenCL REQUIRED)
	find_package(LibXml2 REQUIRED)
	list(APPEND CMAKE_MODULE_PATH "${CL2QCD_ROOT}/cmake")
	find_package(Nettle REQUIRED)
	find_package(LIME REQUIRED)
	find_package(LAPACK REQUIRED)

	include_directories(SYSTEM
		${Boost_INCLUDE_DIRS}
		${OPENCL_INCLUDE_DIRS}
		${LIBXML2_INCLUDE_DIR}
		${Nettle_INCLUDE_DIR}
		${LIME_INCLUDE_DIR}
		${LAPACK_INCLUDE_DIR}
	)

	set(SUN_LIB
		"${SU3_ROOT}/fields.o"
		"${Boost_LIBRARIES}"
		"${OPENCL_LIBRARIES}"
		"${LIBXML2_LIBRARIES}"
		"${Nettle_LIBRARIES}"
		"${LIME_LIBRARIES}"
		"${LAPACK_LIBRARIES}"
	)

	set(SUN_ROOT "${SU3_ROOT}")
	
else()
	message(FATAL_ERROR "NC must be 2 or 3")
	
endif()

add_library(SMEARING_LIB "${SUN_ROOT}/smearing_techniques.o")
set_target_properties(SMEARING_LIB PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	"${Tools_ROOT}/include"
	"${LatticeTools_ROOT}/include"
	"${SUN_ROOT}"
	include
)

find_library(HELPER_LIB helper_functions "${Tools_ROOT}/lib")
find_library(IOTOOLS_LIB io_tools "${Tools_ROOT}/lib")
find_library(TOOLS_LIB tools "${Tools_ROOT}/lib")
find_library(LINKPATH_LIB LinkPath "${LatticeTools_ROOT}/lib")

find_package(OpenMP REQUIRED)

add_library(mlvl_omp
	src/algorithm/MultilevelAnalyzer.cc
)

target_link_libraries(mlvl_omp
	OpenMP::OpenMP_CXX
)

add_library(mlvl
	src/algorithm/MultilevelConfig.cc
	src/algorithm/LevelDef.cc
	src/algorithm/TwolinkOperator.cc
	src/algorithm/TwolinkComputer.cc
	src/algorithm/T_field.cc
	src/physics/twolink_operator_functions.cc
	src/log/logger.cc
)

add_library(sublattice_algebra SHARED
	src/physics/sublattice_algebra.cc
)

target_compile_definitions(mlvl PUBLIC
	F_
)

add_executable(multilevel
	src/bin/multilevel.cc
	src/parameters/parse_parameters.cc
)
target_link_libraries(multilevel
	mlvl_omp
	mlvl
	sublattice_algebra
	"${SUN_INTERFACE_LIB}"
	"${LINKPATH_LIB}"
	"${SUN_LIB}"
	SMEARING_LIB
	"${HELPER_LIB}"
	"${IOTOOLS_LIB}"
	"${TOOLS_LIB}"
	OpenMP::OpenMP_CXX
)

add_executable(find_sizes
	src/find_sizes/main.cc
	src/bin/find_sizes.cc
)

target_link_libraries(find_sizes
	"${HELPER_LIB}"
)

install(TARGETS sublattice_algebra DESTINATION "${CMAKE_SOURCE_DIR}/lib")
install(TARGETS
	multilevel
	find_sizes
DESTINATION "${CMAKE_SOURCE_DIR}/bin")
