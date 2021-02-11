#ifndef __ENGINE_H__
#define __ENGINE_H__

#include <math.h>
#include "defs.h"
#include "utilities.h"

/* compute kinetic energy */
void ekin(mdsys_t* sys);

void force(mdsys_t* sys);

__attribute__((always_inline)) inline void force_pair(mdsys_t* sys, double* irx, double* iry,
                                                      double* irz, double* jrx, double* jry,
                                                      double* jrz, double* epot, double* ipfx,
                                                      double* ipfy, double* ipfz) {
    /* get distance between particle i and j */
    const double rx = pbc(*irx - *jrx, 0.5 * sys->box);
    const double ry = pbc(*iry - *jry, 0.5 * sys->box);
    const double rz = pbc(*irz - *jrz, 0.5 * sys->box);
    const double r = sqrt(rx * rx + ry * ry + rz * rz);

    /* compute force and energy if within cutoff */
    if (r < sys->rcut) {
        const double ffac = -4.0 * sys->epsilon *
            (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

        *epot += 0.5 * 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

        *ipfx += rx / r * ffac;
        *ipfy += ry / r * ffac;
        *ipfz += rz / r * ffac;
    }
};

__attribute__((always_inline)) inline void
force_pair_symmetric(mdsys_t* sys, double* irx, double* iry, double* irz, double* jrx, double* jry,
                     double* jrz, double* epot, double* ipfx, double* ipfy, double* ipfz,
                     double* jpfx, double* jpfy, double* jpfz) {
    /* get distance between particle i and j */
    const double rx = pbc(*irx - *jrx, 0.5 * sys->box);
    const double ry = pbc(*iry - *jry, 0.5 * sys->box);
    const double rz = pbc(*irz - *jrz, 0.5 * sys->box);
    const double r = sqrt(rx * rx + ry * ry + rz * rz);

    /* compute force and energy if within cutoff */
    if (r < sys->rcut) {
        const double ffac = -4.0 * sys->epsilon *
            (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

        *epot += 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

        *ipfx += rx / r * ffac;
        *ipfy += ry / r * ffac;
        *ipfz += rz / r * ffac;

        *jpfx -= rx / r * ffac;
        *jpfy -= ry / r * ffac;
        *jpfz -= rz / r * ffac;
    }
};

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
