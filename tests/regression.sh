rm -f *.xyz *.dat

exe=ljmd.x

#dry run
../${exe}

#compare with references
cd ../../examples
../${exe} < argon_108.inp
head -10 argon_108.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > a.dat
head -10 ../reference/argon_108.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > b.dat
cmp a.dat b.dat || exit 1

../${exe} < argon_2916.inp
head -10 argon_2916.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }'> a.dat
head -10 ../reference/argon_2916.dat | awk '{ printf("%d %.6f %.6f %.6f\n", $1, $2, $3, $4); }' > b.dat
cmp a.dat b.dat || exit 1

rm -f a.dat b.dat