#!/bin/bash

cmake_defines="$(cat dependencies | sed -r 's/^/-D/' | tr '\n' ' ')"
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