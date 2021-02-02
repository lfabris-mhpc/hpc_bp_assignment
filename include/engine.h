#ifndef __ENGINE_H__
#define __ENGINE_H__

#include "defs.h"

void ekin(mdsys_t* sys);

void force(mdsys_t* sys);

void velverlet(mdsys_t* sys);

#ifdef __MPI
void force_mpi(mdsys_t* sys, const int nranks, const int rank);

void velverlet_mpi(mdsys_t* sys, const int nranks, const int rank);
#endif

#endif
