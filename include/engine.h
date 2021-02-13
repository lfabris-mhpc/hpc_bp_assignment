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
    double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;

    /* get distance between particle i and j */
    double rx = pbc(*irx - *jrx, 0.5 * sys->box);
    double ry = pbc(*iry - *jry, 0.5 * sys->box);
    double rz = pbc(*irz - *jrz, 0.5 * sys->box);
    double rsq = rx * rx + ry * ry + rz * rz;

    /* compute force and energy if within cutoff */
    if (rsq < rcsq) {
        double rinv = 1.0 / rsq;
        double r6 = rinv * rinv * rinv;

        double ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
        *epot += 0.5 * r6 * (c12 * r6 - c6);

        *ipfx += rx * ffac;
        *ipfy += ry * ffac;
        *ipfz += rz * ffac;
    }
};

__attribute__((always_inline)) inline void
force_pair_symmetric(mdsys_t* sys, double* irx, double* iry, double* irz, double* jrx, double* jry,
                     double* jrz, double* epot, double* ipfx, double* ipfy, double* ipfz,
                     double* jpfx, double* jpfy, double* jpfz) {
    double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;

    /* get distance between particle i and j */
    double rx = pbc(*irx - *jrx, 0.5 * sys->box);
    double ry = pbc(*iry - *jry, 0.5 * sys->box);
    double rz = pbc(*irz - *jrz, 0.5 * sys->box);
    double rsq = rx * rx + ry * ry + rz * rz;

    /* compute force and energy if within cutoff */
    if (rsq < rcsq) {
        double rinv = 1.0 / rsq;
        double r6 = rinv * rinv * rinv;

        double ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
        *epot += r6 * (c12 * r6 - c6);

        *ipfx += rx * ffac;
        *ipfy += ry * ffac;
        *ipfz += rz * ffac;

        *jpfx -= rx * ffac;
        *jpfy -= ry * ffac;
        *jpfz -= rz * ffac;
    }
};

void force_mpi_basic(mdsys_t* sys);

void force_mpi_ibasic(mdsys_t* sys);

void force_mpi_ibasic_even(mdsys_t* sys);

void force_mpi_primitive(mdsys_t* sys);

void force_mpi_slice(mdsys_t* sys);

void force_mpi_ring(mdsys_t* sys);

void force_mpi_symmring(mdsys_t* sys);

/* velocity verlet - first step before force calculation */
void verlet_1(mdsys_t* sys);

/* velocity verlet - second step must be called after force calculation */
void verlet_2(mdsys_t* sys);

#endif
