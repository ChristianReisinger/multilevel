#!/bin/bash

echo -e "\nWilson loop should be:"
"$ilp_program" su2_T20L20_b2.96.multilevel.1 20 20 0 <<< "$scr" 2> /dev/null

echo -e "\naverage:"
cat WL.test

rm su2_T20L20_b2.96.multilevel.1* WL.test
