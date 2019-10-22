#!/bin/bash

ilp_program="interactive_link_products"

echo -e "\nWilson loop should be:"
"$ilp_program" su2_T20L20_b2.96.multilevel.1 20 20 0 <<< "$scr" 2> /dev/null
rm su2_T20L20_b2.96.multilevel.1*

echo -e "\n__averages__"
for f in *.test; do
	echo "$f:"
	cat "$f"
	rm "$f"
	echo ""
done
