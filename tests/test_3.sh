#!/bin/bash

#test correct averaging for operators with multiple definition on irregular timeslice pattern

R=4
Nc=2	#number of configs at level 1

../bin/multilevel -e test -b 2.96 -s 46421 -u 2,2 -w 20 20 $R 0 1,${Nc} comps_3.dat su2_T20L20_b2.96 1
rm su2_T20L20_b2.96.multilevel.1*

echo -e "\n__averages__"
for f in *.test; do
	echo "$f:"
	cat "$f"
	rm "$f"
	echo ""
done

echo "result should be .."
echo "at R = 3: WL = '(2*WL1 + WL2 + WL3) / 4'"
echo "at R = 4: WL = '(WL4 + WL5 + WL6) / 3'"
