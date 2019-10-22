#!/bin/bash

#test correct computation on irregular timeslice pattern - 3 levels

#compare to multilevel Wilson loop output at <t>,<x>,<y>,<z> (in direction <dir>)
#--> in Wilson loop computing for ( ... ), do
#	if (vector<int>{t,x,y,z,i} == vector<int>{...}) { cout << ... curr_WL ... }

R=4
Nc=2	#number of configs at level 1
Ncc=2	#number of configs at level 2

if [ $# -ne 5 ]; then
	echo "Usage: $0 <t> <x> <y> <z> <dir>"
	exit
fi

../bin/multilevel -e test -b 2.96 -s 46421 -u 5,5,5 -w 20 20 $R 0 1,${Nc},${Ncc} comps_2.dat su2_T20L20_b2.96 1

t0=$1
x0=$2
y0=$3
z0=$4
dir=$5

t2=$(($t0+2))
t4=$(($t0+4))
t7=$(($t0+7))
tT=$(($t0+10))
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

rpath="$(for n in $(seq 1 $R); do echo -n "$i "; done; echo "end")"
n0_0="($t0,$x0,$y0,$z0)"
n0_R="($t0,$xR,$yR,$zR)"
n2_0="($t2,$x0,$y0,$z0)"
n2_R="($t2,$xR,$yR,$zR)"
n4_0="($t4,$x0,$y0,$z0)"
n4_R="($t4,$xR,$yR,$zR)"
n7_0="($t7,$x0,$y0,$z0)"
n7_R="($t7,$xR,$yR,$zR)"
nT_0="($tT,$x0,$y0,$z0)"

scr='
link S0
link ST

twolink T4_0
T4_0 zero
twolink T6_4
T6_4 zero

S0 path '$n0_0'
'"$rpath"'

ST path '$nT_0'
'$rpath'

'"$(for c in $(seq 1 $Nc); do echo '

twolink cT2_0
cT2_0 zero
twolink cT2_2
cT2_2 zero
twolink cT3_4
cT3_4 zero
twolink cT3_7
cT3_7 zero

'"$(for cc in $(seq 1 $Ncc); do echo '
load su2_T20L20_b2.96.multilevel.1.'$c'.'$cc' 0

twolink ccT2_0
twolink ccT2_2
twolink ccT3_4
twolink ccT3_7

link U2_0_0
link U2_0_R
link U2_2_0
link U2_2_R
link U3_4_0
link U3_4_R
link U3_7_0
link U3_7_R

U2_0_0 p '$n0_0'
t t end
U2_0_R p '$n0_R'
t t end
U2_2_0 p '$n2_0'
t t end
U2_2_R p '$n2_R'
t t end
U3_4_0 p '$n4_0'
t t t end
U3_4_R p '$n4_R'
t t t end
U3_7_0 p '$n7_0'
t t t end
U3_7_R p '$n7_R'
t t t end

ccT2_0 =X U2_0_0 U2_0_R
ccT2_2 =X U2_2_0 U2_2_R
ccT3_4 =X U3_4_0 U3_4_R
ccT3_7 =X U3_7_0 U3_7_R

del U2_0_0
del U2_0_R
del U2_2_0
del U2_2_R
del U3_4_0
del U3_4_R
del U3_7_0
del U3_7_R

cT2_0 += ccT2_0
cT2_2 += ccT2_2
cT3_4 += ccT3_4
cT3_7 += ccT3_7

del ccT2_0
del ccT2_2
del ccT3_4
del ccT3_7
'
done)"'

twolink cT4_0
cT4_0 *= cT2_0
cT4_0 *= cT2_2
del cT2_0
del cT2_2

twolink cT6_4
cT6_4 *= cT3_4
cT6_4 *= cT3_7
del cT3_4
del cT3_7


T4_0 += cT4_0
del cT4_0
T6_4 += cT6_4
del cT6_4
'
done)"'

twolink T
T *= T4_0
T *= T6_4
T *re '"$(echo "1 / (${Nc}*${Nc}*${Ncc}*${Ncc}*${Ncc}*${Ncc})" | bc -l)"'

T WL S0 ST

exit
'

. ./test_fin.sh
