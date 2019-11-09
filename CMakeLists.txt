cmake_minimum_required(VERSION 3.10)
project(multilevel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

set(NC 2 CACHE STRING "number of colors for SU(NC)")
if(NOT NC MATCHES "^[23]$")
	message(FATAL_ERROR "NC must be 2 or 3")
endif()

add_library(SUN_O OBJECT IMPORTED)
if(NC EQUAL 2)
	set_property(TARGET SUN_O PROPERTY IMPORTED_OBJECTS
		"${SUN_PACKAGE}/io.o"
		"${SUN_PACKAGE}/ranlxs.o"
		"${SUN_PACKAGE}/ranlxd.o"
		"${SUN_PACKAGE}/ranlux.o"
		"${SUN_PACKAGE}/heatbath.o"
		"${SUN_PACKAGE}/smearing_techniques.o"
		"${LatticeTools_DIR}/lib/LinkPath.o"
		"${LatticeTools_DIR}/lib/MCSU2Inteferace.o"
	)
	set(SUN_PACKAGE "${SU2_PACKAGE}")
elseif(NC EQUAL 3)
	set_property(TARGET SUN_O PROPERTY IMPORTED_OBJECTS
#		... CL2QCD libs ...
		"${LatticeTools_DIR}/lib/LinkPath.o"
		"${LatticeTools_DIR}/lib/CL2QCDInterface.o"
	)
	set(SUN_PACKAGE "${SU3_PACKAGE}")
endif()
add_library(SUN_LIB $<TARGET_OBJECTS:SUN_O>)
set_target_properties(SUN_LIB PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	"${Tools_DIR}/include"
	"${LatticeTools_DIR}/include"
	"${SUN_PACKAGE}"
	include
)

find_library(HELPER_LIB helper_functions "${Tools_DIR}/lib")
find_library(IOTOOLS_LIB io_tools "${Tools_DIR}/lib")

add_library(mlvl
	src/sublattice_algebra.cc
	src/twolink_operator_functions.cc
	src/T_field.cc
	src/TwolinkComputer.cc
	src/TwolinkOperator.cc
	src/LevelDef.cc
	src/MultilevelConfig.cc
	src/MultilevelAnalyzer.cc
)
target_compile_definitions(mlvl
	PUBLIC __SUN_N__=${NC}
)

add_executable(multilevel
	src/parse_parameters.cc
	src/multilevel.cc
)
target_link_libraries(multilevel
	"${HELPER_LIB}"
	"${IOTOOLS_LIB}"
	SUN_LIB
	mlvl
)
