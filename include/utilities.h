#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "defs.h"
#include "input_output.h"


double wallclock();


void azzero(double *d, const int n);


double pbc(double x, const double boxby2);


#endif
