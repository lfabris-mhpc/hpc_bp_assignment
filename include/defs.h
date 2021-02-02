#ifndef __DEFS_H__
#define __DEFS_H__

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
#ifdef __MPI
	//partial forces
	double *px, *py, *pz;
#endif
};
typedef struct _mdsys mdsys_t;

#endif
