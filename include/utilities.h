#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <stdio.h>

int get_a_line(FILE *fp, char *buf);


double wallclock();


void azzero(double *d, const int n);


double pbc(double x, const double boxby2);


#endif
