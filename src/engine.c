

#include "engine.h"

/* a few physical constants */
const double kboltz = 0.0019872067; /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* compute kinetic energy */
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

/* compute forces */
void force(mdsys_t* sys) {
    double rsq, rcsq, ffac;
    double rx, ry, rz;
    double c6, c12;
    int i, j;

    /* zero energy and forces */
    sys->epot = 0.0;
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    rcsq = sys->rcut * sys->rcut;

    for (i = 0; i < (sys->natoms) - 1; ++i) {
        for (j = i + 1; j < (sys->natoms); ++j) {
            /* get distance between particle i and j */
            rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            rsq = rx * rx + ry * ry + rz * rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double r6, rinv;
                rinv = 1.0 / rsq;
                r6 = rinv * rinv * rinv;

                ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;

                sys->epot += r6 * (c12 * r6 - c6);

                sys->fx[i] += rx * ffac;
                sys->fy[i] += ry * ffac;
                sys->fz[i] += rz * ffac;

                sys->fx[j] -= rx * ffac;
                sys->fy[j] -= ry * ffac;
                sys->fz[j] -= rz * ffac;
            }
        }
    }
}

/* velocity verlet */
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
