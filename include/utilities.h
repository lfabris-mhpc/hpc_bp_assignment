#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "defs.h"
#include "input_output.h"

double wallclock();

void azzero(double* d, const int n);

double pbc(double x, const double boxby2);

void init_displs_counts(const int natoms, const int nranks, int* displs, int* counts);

void init_displs_counts_even(const int natoms, const int nranks, int* displs, int* counts);

#endif
