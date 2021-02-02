#ifndef __ENGINE_H__
#define __ENGINE_H__

#include <math.h>
#include "defs.h"
#include "utilities.h"

void ekin(mdsys_t *sys);

void force(mdsys_t *sys);

void velverlet(mdsys_t *sys);

#endif
