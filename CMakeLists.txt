cmake_minimum_required(VERSION 3.10)
project(multilevel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(NC 2 CACHE STRING "number of colors for SU(NC)")
if(NOT NC MATCHES "^[23]$")
	message(FATAL_ERROR "NC must be 2 or 3")
endif()

if(NC EQUAL 2)
	set(SUN_PACKAGE "${SU2_PACKAGE}")
	set(SUN_INTERFACE_OBJ "${SU2Interface_OBJ}")
elseif(NC EQUAL 3)
	set(SUN_PACKAGE "${SU3_PACKAGE}")
	set(SUN_INTERFACE_OBJ "${SU3Interface_OBJ}")
endif()

include_directories(
	"${Tools_DIR}/include"
	"${LatticeTools_DIR}/include"
	"${SUN_PACKAGE}"
	include
)

find_library(HELPER_LIB helper_functions "${Tools_DIR}/lib")
find_library(IOTOOS_LIB io_tools "${Tools_DIR}/lib")
add_library(SUN_O OBJECT IMPORTED)
set_property(TARGET SUN_O IMPORTED_OBJECTS
	"${SUN_PACKAGE}/smearing_techniques.o"
	"${LatticeTools_DIR}/lib/LinkPath.o"
	"${SUN_INTERFACE}"
)
add_library(SUN_LIB $<TARGET_OBJECTS:SUN_O>)
set_target_properties(SUN_O LINKER_LANGUAGE CXX)

add_library(mlvl
	"${HELPER_LIB}"
	SUN_LIB
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

add_library(multilevel
	"${IOTOOLS_LIB}"
	SUN_LIB
	parse_parameters.cc
	multilevel.cc
)
target_link_libraries(multilevel

