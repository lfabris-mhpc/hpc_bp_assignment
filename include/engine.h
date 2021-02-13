#ifndef __ENGINE_H__
#define __ENGINE_H__

#include <math.h>
#include "defs.h"
#include "utilities.h"

/* compute kinetic energy */
void ekin(mdsys_t* sys);

void force(mdsys_t* sys);

/* velocity verlet - first step before force calculation */
void verlet_1(mdsys_t* sys);

/* velocity verlet - second step must be called after force calculation */
void verlet_2(mdsys_t* sys);

#endif
