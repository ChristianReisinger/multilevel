cmake_minimum_required(VERSION 3.10)
project(multilevel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(NC 2 CACHE STRING "number of colors for SU(NC)")

if(NC EQUAL 2)
	find_library(SUN_INTERFACE_LIB MCSU2Interface "${LatticeTools_ROOT}/lib")
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
	find_library(SUN_INTERFACE_LIB CL2QCDInterface "${LatticeTools_ROOT}/lib")
	
	find_library(RANLUX ranlux "${CL2QCD_ROOT}/build/ranlux")
	find_library(META meta "${CL2QCD_ROOT}/build/meta")
	find_library(PHYSICS physics "${CL2QCD_ROOT}/build/physics")
	find_library(PHYSICS_ALGORITHMS algorithms "${CL2QCD_ROOT}/build/physics/algorithms")
	find_library(PHYSICS_ALGORITHMS_SOLVERS solvers "${CL2QCD_ROOT}/build/physics/algorithms/solvers")
	find_library(PHYSICS_OBSERVABLES_1 observables "${CL2QCD_ROOT}/build/physics/observables")
	find_library(PHYSICS_OBSERVABLES_2 gaugeObservables "${CL2QCD_ROOT}/build/physics/observables")
	find_library(PHYSICS_FERMIONMATRIX fermionmatrix "${CL2QCD_ROOT}/build/physics/fermionmatrix")
	find_library(PHYSICS_LATTICES_1 gaugefield "${CL2QCD_ROOT}/build/physics/lattices")
	find_library(PHYSICS_LATTICES_2 lattices "${CL2QCD_ROOT}/build/physics/lattices")
	find_library(PHYSICS prng "${CL2QCD_ROOT}/build/physics")
	find_library(KLEPSYDRA klepsydra "${CL2QCD_ROOT}/build")
	find_library(CRYPTO crypto "${CL2QCD_ROOT}/build/crypto")
	find_library(HOST_FUNCTIONALITY_1 logger "${CL2QCD_ROOT}/build/host_functionality")
	find_library(HOST_FUNCTIONALITY_2 host_functionality "${CL2QCD_ROOT}/build/host_functionality")
	find_library(EXECUTABLES_1 executables "${CL2QCD_ROOT}/build/executables")
	find_library(EXECUTABLES_2 exceptions "${CL2QCD_ROOT}/build/executables")
	find_library(ILDG_IO_LIME lime "${CL2QCD_ROOT}/build/ildg_io/lime")
	find_library(ILDG_IO_SOURCEFILEPARAMETERS sourcefileParameters "${CL2QCD_ROOT}/build/ildg_io/sourcefileParameters")
	find_library(ILDG_IO ildg_io "${CL2QCD_ROOT}/build/ildg_io")
	find_library(EINHARD einhard "${CL2QCD_ROOT}/build/einhard")
	find_library(OPTIMAL optimal "${CL2QCD_ROOT}/build")
	find_library(HARDWARE_BUFFERS buffers "${CL2QCD_ROOT}/build/hardware/buffers")
	find_library(HARDWARE_1 hardwareTestUtilities "${CL2QCD_ROOT}/build/hardware")
	find_library(HARDWARE_CODE_1 code "${CL2QCD_ROOT}/build/hardware/code")
	find_library(HARDWARE_CODE_2 kernelTestUtilities "${CL2QCD_ROOT}/build/hardware/code")
	find_library(HARDWARE_2 hardware "${CL2QCD_ROOT}/build/hardware")
	find_library(HARDWARE_LATTICES hardware_lattices "${CL2QCD_ROOT}/build/hardware/lattices")
	find_library(HARDWARE_TRANSFER transfer "${CL2QCD_ROOT}/build/hardware/transfer")
	find_library(GEOMETRY geometry "${CL2QCD_ROOT}/build/geometry")
	set(SUN_LIB
		"${RANLUX}"
		"${META}"
		"${PHYSICS}"
		"${PHYSICS_ALGORITHMS}"
		"${PHYSICS_ALGORITHMS_SOLVERS}"
		"${PHYSICS_OBSERVABLES_1}"
		"${PHYSICS_OBSERVABLES_2}"
		"${PHYSICS_FERMIONMATRIX}"
		"${PHYSICS_LATTICES_1}"
		"${PHYSICS_LATTICES_2}"
		"${PHYSICS}"
		"${KLEPSYDRA}"
		"${CRYPTO}"
		"${HOST_FUNCTIONALITY_1}"
		"${HOST_FUNCTIONALITY_2}"
		"${EXECUTABLES_1}"
		"${EXECUTABLES_2}"
		"${ILDG_IO_LIME}"
		"${ILDG_IO_SOURCEFILEPARAMETERS}"
		"${ILDG_IO}"
		"${EINHARD}"
		"${OPTIMAL}"
		"${HARDWARE_BUFFERS}"
		"${HARDWARE_1}"
		"${HARDWARE_CODE_1}"
		"${HARDWARE_CODE_2}"
		"${HARDWARE_2}"
		"${HARDWARE_LATTICES}"
		"${HARDWARE_TRANSFER}"
		"${GEOMETRY}"
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
find_library(LINKPATH_LIB LinkPath "${LatticeTools_ROOT}/lib")

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
	"${SUN_INTERFACE_LIB}"
	"${LINKPATH_LIB}"
	"${SUN_LIB}"
	SMEARING_LIB
	"${HELPER_LIB}"
	"${IOTOOLS_LIB}"
)

install(TARGETS multilevel DESTINATION "${CMAKE_SOURCE_DIR}/bin")
