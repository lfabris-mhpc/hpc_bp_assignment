#ifndef __DEFS_H__
#define __DEFS_H__

#include <stdlib.h>

//#include <mpi.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms, nfi, nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    // MPI
    // MPI_Comm comm;
    void* comm;
    int nranks, rank;
    int *displs, *counts;
    double *pfx, *pfy, *pfz;
    double *srx, *sry, *srz;
    double *rrx, *rry, *rrz;
};
typedef struct _mdsys mdsys_t;

/* allocate memory */
int mdsys_init(mdsys_t* sys);

/* deallocate memory */
void mdsys_free(mdsys_t* sys);

#endif
