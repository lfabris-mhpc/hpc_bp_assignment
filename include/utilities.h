#ifndef __UTILITIES_H__
#define __UTILITIES_H__

static int get_a_line(FILE *fp, char *buf);


static double wallclock();


static void azzero(double *d, const int n);


static double pbc(double x, const double boxby2);


#endif
