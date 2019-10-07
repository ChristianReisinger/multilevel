#!/bin/bash

#compare to multilevel Wilson loop output at <t>,<x>,<y>,<z> (in direction <dir>)
#--> in Wilson loop computing for ( ... ), do
#	if (vector<int>{t,x,y,z,i} == vector<int>{...}) { cout << ... curr_WL ... }

ilp_program="interactive_link_products"

R=4
Nc=2 #number of configs at level 1

if [ $# -ne 5 ]; then
	echo "Usage: $0 <t> <x> <y> <z> <dir>"
	exit
fi

../bin/multilevel -e test -b 2.96 -s 46421 -u 5,5 -w 20 20 $R 0 1,${Nc} comps_1.dat su2_T20L20_b2.96 1

t0=$1
x0=$2
y0=$3
z0=$4
dir=$5

tT=$(($t0+4))
t1=$(($t0+1))
xR=$x0
yR=$y0
zR=$z0

if [ $dir -eq 1 ]; then
	xR=$(($x0+$R))
	i=x
elif [ $dir -eq 2 ]; then
	yR=$(($y0+$R))
	i=y
elif [ $dir -eq 3 ]; then
	zR=$(($z0+$R))
	i=z
fi

rpath="$(for n in $(seq 1 $R); do echo "$i "; done; echo "end")"
n0_0="($t0,$x0,$y0,$z0)"
n1_0="($t1,$x0,$y0,$z0)"
nT_0="($tT,$x0,$y0,$z0)"
n0_R="($t0,$xR,$yR,$zR)"
n1_R="($t1,$xR,$yR,$zR)"

scr='
link S0
link ST

twolink T1_0
T1_0 zero
twolink T3_1
T3_1 zero

S0 path '$n0_0'
'"$rpath"'

ST path '$nT_0'
'$rpath'

'"$(for c in $(seq 1 $Nc); do echo '
load su2_T20L20_b2.96.multilevel.1.'$c' 0

link U1_0_0
U1_0_0 path '$n0_0'
t end

link U1_0_R
U1_0_R path '$n0_R'
t end

twolink cT1_0
cT1_0 =X U1_0_0 U1_0_R
del U1_0_0
del U1_0_R


link U3_1_0
U3_1_0 path '$n1_0'
t t t end

link U3_1_R
U3_1_R path '$n1_R'
t t t end

twolink cT3_1
cT3_1 =X U3_1_0 U3_1_R
del U3_1_0
del U3_1_R

T1_0 += cT1_0
del cT1_0

T3_1 += cT3_1
del cT3_1
'
done)"'

twolink T
T *= T1_0
T *= T3_1
T *re '"$(echo "1 / (${Nc}*${Nc})" | bc -l)"'

T WL S0 ST

exit
'

. ./test_fin.sh
