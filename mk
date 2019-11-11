#!/bin/bash

if [ "$1" == build ]; then
	(cd build
	cmake $(cat ../dependencies | sed -r 's/^/-D/' | tr '\n' ' ') ..)
elif [ "$1" == clear ]; then
	rm -rf build/* bin/*
elif [ "$1" == clean ]; then
	(cd build
	make clean)
elif [ "$1" == "" ]; then
	(cd build
	make)
else
	echo "Invalid target '$1'"
fi