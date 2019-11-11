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
		"${LatticeTools_DIR}/lib/LinkPath.o"
		"${LatticeTools_DIR}/lib/MCSU2Interface.o"
		"${SU2_ROOT}/io.o"
		"${SU2_ROOT}/heatbath.o"
		"${SU2_ROOT}/smearing_techniques.o"
		"${SU2_ROOT}/ranlxs.o"
		"${SU2_ROOT}/ranlxd.o"
		"${SU2_ROOT}/ranlux.o"
	)
	set(SUN_ROOT "${SU2_ROOT}")
elseif(NC EQUAL 3)
	set_property(TARGET SUN_O PROPERTY IMPORTED_OBJECTS
		"${LatticeTools_DIR}/lib/LinkPath.o"
		"${LatticeTools_DIR}/lib/CL2QCDInterface.o"
		
	)
	set(SUN_ROOT "${SU3_ROOT}")
endif()
add_library(SUN_LIB $<TARGET_OBJECTS:SUN_O>)
set_target_properties(SUN_LIB PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	"${Tools_DIR}/include"
	"${LatticeTools_DIR}/include"
	"${SUN_ROOT}"
	include
)

find_library(HELPER_LIB helper_functions "${Tools_DIR}/lib")
find_library(IOTOOLS_LIB io_tools "${Tools_DIR}/lib")

add_library(mlvl
	src/MultilevelAnalyzer.cc
	src/MultilevelConfig.cc
	src/LevelDef.cc
	src/TwolinkOperator.cc
	src/TwolinkComputer.cc
	src/T_field.cc
	src/twolink_operator_functions.cc
	src/sublattice_algebra.cc
)
target_compile_definitions(mlvl
	PUBLIC __SUN_N__=${NC}
)

add_executable(multilevel
	src/multilevel.cc
	src/parse_parameters.cc
)
target_link_libraries(multilevel
	mlvl
	SUN_LIB
	"${HELPER_LIB}"
	"${IOTOOLS_LIB}"
)
