cmake_minimum_required(VERSION 3.10)
project(multilevel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(NC 2 CACHE STRING "number of colors for SU(NC)")

add_definitions(-D__SUN_N__=${NC})

if(NC EQUAL 2)
	find_library(SUN_INTERFACE_LIB MCSU2Gaugefield "${LatticeTools_ROOT}/lib")
	add_library(SU2_LIB 
		"${SU2_ROOT}/io.o"
		"${SU2_ROOT}/heatbath.o"
		"${SU2_ROOT}/ranlxs.o"
		"${SU2_ROOT}/ranlxd.o"
		"${SU2_ROOT}/ranlux.o"
	)
	set_target_properties(SU2_LIB PROPERTIES LINKER_LANGUAGE CXX)
	set(SUN_LIB SU2_LIB)
	set(SUN_ROOT "${SU2_ROOT}")
	
elseif(NC EQUAL 3)
	find_library(SUN_INTERFACE_LIB CL2QCDGaugefield "${LatticeTools_ROOT}/lib")
	
	find_library(PHYSICS_ALGORITHMS algorithms "${CL2QCD_ROOT}/build/physics/algorithms")
	find_library(PHYSICS_GAUGEFIELD gaugefield "${CL2QCD_ROOT}/build/physics/lattices")
	find_library(PHYSICS_PRNG prng "${CL2QCD_ROOT}/build/physics")
	find_library(META meta "${CL2QCD_ROOT}/build/meta")
	find_library(HARDWARE hardware "${CL2QCD_ROOT}/build/hardware")
	find_library(EXCEPTIONS exceptions "${CL2QCD_ROOT}/build/executables")
	find_library(HARDWARE_LATTICES hardware_lattices "${CL2QCD_ROOT}/build/hardware/lattices")
	find_library(ILDG_IO ildg_io "${CL2QCD_ROOT}/build/ildg_io")
	find_library(CONTRACTIONCODE_IO contractioncode_io "${CL2QCD_ROOT}/build/contractioncode_io")
	find_library(LOGGER logger "${CL2QCD_ROOT}/build/host_functionality")
	find_library(EINHARD einhard "${CL2QCD_ROOT}/build/einhard")
	find_library(HOST_FUNCTIONALITY host_functionality "${CL2QCD_ROOT}/build/host_functionality")
	find_library(GEOMETRY geometry "${CL2QCD_ROOT}/build/geometry")
	find_library(HARDWARE_CODE code "${CL2QCD_ROOT}/build/hardware/code")
	find_library(HARDWARE_TRANSFER transfer "${CL2QCD_ROOT}/build/hardware/transfer")
	find_library(HARDWARE_BUFFERS buffers "${CL2QCD_ROOT}/build/hardware/buffers")
	find_library(RANLUX ranlux "${CL2QCD_ROOT}/build/ranlux")
	find_library(CRYPTO crypto "${CL2QCD_ROOT}/build/crypto" NO_DEFAULT_PATH)
	find_library(ILDG_IO_LIME lime "${CL2QCD_ROOT}/build/ildg_io/lime" NO_DEFAULT_PATH)
	find_library(ILDG_IO_SOURCEFILEPARAMETERS sourcefileParameters "${CL2QCD_ROOT}/build/ildg_io/sourcefileParameters")

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
		"${PHYSICS_ALGORITHMS}"
		"${PHYSICS_GAUGEFIELD}"
		"${PHYSICS_PRNG}"
		"${META}"
		"${HARDWARE}"
		"${EXCEPTIONS}"
		"${HARDWARE_LATTICES}"
		"${ILDG_IO}"
		"${CONTRACTIONCODE_IO}"
		"${LOGGER}"
		"${EINHARD}"
		"${HOST_FUNCTIONALITY}"
		"${GEOMETRY}"
		"${HARDWARE_CODE}"
		"${HARDWARE_TRANSFER}"
		"${HARDWARE_BUFFERS}"
		"${RANLUX}"
		"${CRYPTO}"
		"${ILDG_IO_LIME}"
		"${ILDG_IO_SOURCEFILEPARAMETERS}"
#
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

add_executable(try_operator_on_timeslices
	src/bin/try_operator_on_timeslices.cc
)

target_link_libraries(try_operator_on_timeslices
	"${HELPER_LIB}"
)

install(TARGETS sublattice_algebra DESTINATION "${CMAKE_SOURCE_DIR}/lib")
install(TARGETS
	multilevel
	try_operator_on_timeslices
DESTINATION "${CMAKE_SOURCE_DIR}/bin")
