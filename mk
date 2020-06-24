#!/bin/bash

depends_template="\
Tools_ROOT=
LatticeTools_ROOT=
SU2_ROOT=
SU3_ROOT=
CL2QCD_ROOT=
BOOST_ROOT:PATH=
CMAKE_MODULE_PATH="

if [ -s dependencies ]; then
	cmake_defines="$(cat dependencies | sed -r 's/^/-D/' | tr '\n' ' ')"
else
	echo -ne "$depends_template" > dependencies
	echo "dependencies not set .. creating template"
	exit
fi

mkdir -p build

(cd build

if [ "${1,,}" == su2 ]; then
	cmake $cmake_defines ..
elif [ "${1,,}" == su3 ]; then
	cmake -DNC=3 $cmake_defines ..
elif [ "$1" == clear ]; then
	rm -rf ../build/* ../bin/*
elif [ "$1" == clean ]; then
	make clean
	rm -rf ../bin/*
elif [ "$1" == "" ]; then
	make install
else
	echo "Invalid target '$1'"
fi

)
