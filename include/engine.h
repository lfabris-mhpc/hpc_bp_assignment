#ifndef __ENGINE_H__
#define __ENGINE_H__

#include <math.h>
#include "defs.h"
#include "utilities.h"

/* compute kinetic energy */
void ekin(mdsys_t* sys);

void force(mdsys_t* sys);

void force_pair(mdsys_t* sys, double* irx, double* iry, double* irz, double* jrx, double* jry,
                double* jrz, double* epot, double* ipfx, double* ipfy, double* ipfz);

void force_pair_symmetric(mdsys_t* sys, double* irx, double* iry, double* irz, double* jrx,
                          double* jry, double* jrz, double* epot, double* ipfx, double* ipfy,
                          double* ipfz, double* jpfx, double* jpfy, double* jpfz);

void force_mpi_basic(mdsys_t* sys);

void force_mpi_ibasic(mdsys_t* sys);

void force_mpi_ibasic_even(mdsys_t* sys);

void force_mpi_primitive(mdsys_t* sys);

void force_mpi_slice(mdsys_t* sys);

void force_mpi_ring(mdsys_t* sys);

void force_mpi_symmring(mdsys_t* sys);

void force_openmp_wredux(mdsys_t* sys);

/* velocity verlet - first step before force calculation */
void verlet_1(mdsys_t* sys);

/* velocity verlet - second step must be called after force calculation */
void verlet_2(mdsys_t* sys);

#endif
