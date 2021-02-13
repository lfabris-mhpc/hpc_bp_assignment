#!/bin/bash

function process() {
	rm *.dat *.inp *.xyz *.rest 
	
	cp examples/${ref}.* .
	
	../${exe} < ${ref}.inp
	head -10 ${ref}.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > a.dat
	head -10 reference/${ref}.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > b.dat
	check=$(paste a.dat b.dat | awk 'function abs(x) {return (x<0?-x:x);} {for (i=1; i<5; ++i) {j=i+4; if (abs($i-$j) > 0.01) {printf("failed %.6f vs %.6f\n", $i, $j);} } }')
	
	rm -f a.dat b.dat

	if [ -n "${check}" ]
	then
		echo "${check}"
		exit 1
	fi
}

export OMP_NUM_THREADS=$1
exe=ljmd.x

#compare with references
ref=argon_108
process

exit 0
