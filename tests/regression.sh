#!/bin/bash

tolerance=0.0001

rm *.dat *.inp *.xyz *.rest

cp examples/${2}.* .

${1} < ${2}.inp
res=$?
if ((res != 0))
then
	exit 1
fi
head -10 ${2}.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > a.dat
head -10 reference/${2}.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > b.dat
check=$(paste a.dat b.dat | awk 'function abs(x) {return (x<0?-x:x);} {for (i=1; i<5; ++i) {j=i+4; if (abs($i-$j) > '${tolerance}') {printf("failed %.6f vs %.6f\n", $i, $j);} } }')

rm -f a.dat b.dat

if [ -n "${check}" ]
then
	echo "test on ${2} failed"
	exit 1
else
	echo "test on ${2} passed under tolerance ${tolerance}"
fi

exit 0
