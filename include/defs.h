#ifndef __DEFS_H__
#define __DEFS_H__

#include <stdlib.h>

#if !defined(__cplusplus)
#include <mpi.h>
#else
struct ompi_communicator_t;
typedef struct ompi_communicator_t* MPI_Comm;
#endif

#define UNUSED(x) (void)(x)

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
    MPI_Comm comm;
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
