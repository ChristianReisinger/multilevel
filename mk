#!/bin/bash

if [ "$1" == build ]; then
	(cd build
	cmake $(cat ../multilevel.ini | sed -r 's/^/-D/' | tr '\n' ' ') ..)
elif [ "$1" == clear ]; then
	rm -rf build/* bin/*
elif [ "$1" == clean ]; then
	(cd build
	make clean)
else
	(cd build
	make)
fi
