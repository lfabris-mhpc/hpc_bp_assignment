#include <math.h>
#include <assert.h>

#include <omp.h>

#if !defined(__cplusplus)
#include <mpi.h>
#endif

#include "engine.h"
#include "defs.h"
#include "utilities.h"

/* a few physical constants */
const double kboltz = 0.0019872067; /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

void ekin(mdsys_t* sys) {
    int i;

    double mass_const = 0.5 * mvsq2e * sys->mass;

    sys->ekin = 0.0;
    for (i = 0; i < sys->natoms; ++i) {
        sys->ekin += mass_const *
            (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] + sys->vz[i] * sys->vz[i]);
    }
    sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}

void force(mdsys_t* sys) {
    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    int check = MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->comm);
    assert(check == MPI_SUCCESS);
    check = MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->comm);
    assert(check == MPI_SUCCESS);
    check = MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->comm);
    assert(check == MPI_SUCCESS);

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;
	
	double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
	double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
	double rcsq = sys->rcut * sys->rcut;

#pragma omp parallel for schedule(dynamic) shared(sys, c12, c6, rcsq, rx, ry, rz) reduction(+: epot, pfx[:sys->natoms], pfy[:sys->natoms], pfz[:sys->natoms])
    for (int i = 0; i < (sys->natoms - 1); i += sys->nranks) {
        int ii = i + sys->rank;
        if (ii < (sys->natoms - 1)) {
            for (int j = ii + 1; j < sys->natoms; ++j) {
				/* get distance between particle i and j */
				double prx = pbc(rx[ii] - rx[j], 0.5 * sys->box);
				double pry = pbc(ry[ii] - ry[j], 0.5 * sys->box);
				double prz = pbc(rz[ii] - rz[j], 0.5 * sys->box);
				double rsq = prx * prx + pry * pry + prz * prz;

				/* compute force and energy if within cutoff */
				if (rsq < rcsq) {
					double rinv = 1.0 / rsq;
					double r6 = rinv * rinv * rinv;

					double ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
					epot += r6 * (c12 * r6 - c6);

					pfx[ii] += prx * ffac;
					pfy[ii] += pry * ffac;
					pfz[ii] += prz * ffac;

					pfx[j] -= prx * ffac;
					pfy[j] -= pry * ffac;
					pfz[j] -= prz * ffac;
				}
            }
        }
    }

    check = MPI_Reduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm);
    assert(check == MPI_SUCCESS);
    check = MPI_Reduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm);
    assert(check == MPI_SUCCESS);
    check = MPI_Reduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm);
    assert(check == MPI_SUCCESS);
    check = MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm);
    assert(check == MPI_SUCCESS);

    UNUSED(check);
}

void verlet_1(mdsys_t* sys) {
    int i;

    double mass_time_const = 0.5 * sys->dt / mvsq2e / sys->mass;

    /* first part: propagate velocities by half and positions by full step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += mass_time_const * sys->fx[i];
        sys->vy[i] += mass_time_const * sys->fy[i];
        sys->vz[i] += mass_time_const * sys->fz[i];

        sys->rx[i] += sys->dt * sys->vx[i];
        sys->ry[i] += sys->dt * sys->vy[i];
        sys->rz[i] += sys->dt * sys->vz[i];
    }
}

void verlet_2(mdsys_t* sys) {
    int i;

    double mass_time_const = 0.5 * sys->dt / mvsq2e / sys->mass;

    /* second part: propagate velocities by another half step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += mass_time_const * sys->fx[i];
        sys->vy[i] += mass_time_const * sys->fy[i];
        sys->vz[i] += mass_time_const * sys->fz[i];
    }
}
