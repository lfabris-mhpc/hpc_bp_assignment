
#include <math.h>
#include <assert.h>
#include "engine.h"
#include "defs.h"
#include "utilities.h"

#ifdef __MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/* a few physical constants */
const double kboltz = 0.0019872067; /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* compute kinetic energy */
void ekin(mdsys_t* sys) {
    int i;

    sys->ekin = 0.0;
    for (i = 0; i < sys->natoms; ++i) {
        sys->ekin += 0.5 * mvsq2e * sys->mass *
            (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] + sys->vz[i] * sys->vz[i]);
    }
    sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}

/* compute forces */
void force(mdsys_t* sys) {
    double r, ffac;
    double rx, ry, rz;
    int i, j;

    /* zero energy and forces */
    sys->epot = 0.0;
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    for (i = 0; i < (sys->natoms); ++i) {
        for (j = 0; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
            if (i == j)
                continue;

            /* get distance between particle i and j */
            rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                sys->epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->fx[i] += rx / r * ffac;
                sys->fy[i] += ry / r * ffac;
                sys->fz[i] += rz / r * ffac;
            }
        }
    }
}

/* velocity verlet */
void velverlet(mdsys_t* sys) {
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt * sys->vx[i];
        sys->ry[i] += sys->dt * sys->vy[i];
        sys->rz[i] += sys->dt * sys->vz[i];
    }

    /* compute forces and potential energy */
    force(sys);

    /* second part: propagate velocities by another half step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}

#ifdef __MPI
void force_mpi(mdsys_t* sys, const int nranks, const int rank) {
    double r, ffac;
    double rx, ry, rz;
    int i, j;

    /* zero energy and forces */
    //sys->epot = 0.0;
	
	MPI_Request reqs[4];
	MPI_Status statuses[4];
	
	int ret = MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, reqs);
	assert(ret == MPI_SUCCESS);
	ret = MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, reqs + 1);
	assert(ret == MPI_SUCCESS);
	ret = MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, reqs + 2);
	assert(ret == MPI_SUCCESS);

	const int begin = rank * (sys->natoms / nranks);
	const int end = begin + sys->natoms /nranks + ((sys->natoms % nranks) != 0);
#ifndef NDEBUG
	//printf("rank %d will calc forces on [%d, %d)\n", rank, begin, end);
#endif
    /*
	azzero(sys->px + begin, end - begin);
    azzero(sys->py + begin, end - begin);
    azzero(sys->pz + begin, end - begin);
	*/
#ifdef _OPENMP
	double *px = sys->px;
	double *py = sys->py;
	double *pz = sys->pz;
#endif
	azzero(sys->px, sys->natoms);
    azzero(sys->py, sys->natoms);
    azzero(sys->pz, sys->natoms);
	azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);
	double epot = 0.0;

	ret = MPI_Waitall(3, reqs, statuses);
	assert(ret == MPI_SUCCESS);

#ifdef _OPENMP
	#pragma omp parallel for shared(sys, begin, end) private(i, j, r, ffac, rx, ry, rz) reduction(+: px[begin:end], py[begin:end], pz[begin:end], epot)
#endif
    for (i = begin; i < end; ++i) {
        for (j = 0; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
            if (i == j)
                continue;

            /* get distance between particle i and j */
            rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
    		rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->px[i] += rx / r * ffac;
                sys->py[i] += ry / r * ffac;
                sys->pz[i] += rz / r * ffac;
            }
        }
    }

	ret = MPI_Ireduce(sys->px, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, reqs);
	assert(ret == MPI_SUCCESS);
	ret = MPI_Ireduce(sys->py, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, reqs + 1);
	assert(ret == MPI_SUCCESS);
	ret = MPI_Ireduce(sys->pz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, reqs + 2);
	assert(ret == MPI_SUCCESS);
	ret = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, reqs + 3);
	assert(ret == MPI_SUCCESS);

	ret = MPI_Waitall(4, reqs, statuses);
	assert(ret == MPI_SUCCESS);
}

/* velocity verlet */
void velverlet_mpi(mdsys_t* sys, const int nranks, const int rank) {
    int i;

#ifdef __MPI
	if (!rank)
#endif
    /* first part: propagate velocities by half and positions by full step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt * sys->vx[i];
        sys->ry[i] += sys->dt * sys->vy[i];
        sys->rz[i] += sys->dt * sys->vz[i];
    }

    /* compute forces and potential energy */
    force_mpi(sys, nranks, rank);

#ifdef __MPI
	if (!rank)
#endif
    /* second part: propagate velocities by another half step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}
#endif
