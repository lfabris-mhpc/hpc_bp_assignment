#ifndef __ENGINE_H__
#define __ENGINE_H__

#include <math.h>
#include "defs.h"
#include "utilities.h"

void ekin(mdsys_t* sys);

void force(mdsys_t* sys);

void verlet_1(mdsys_t* sys);

void verlet_2(mdsys_t* sys);

void force_openmp(mdsys_t* sys); 

void force_openmp_nonew(mdsys_t* sys); 

#endif
