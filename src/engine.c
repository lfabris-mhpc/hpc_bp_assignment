#include <mpi.h>
#include <omp.h>

#include "engine.h"

/* a few physical constants */
const double kboltz = 0.0019872067; /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

void ekin(mdsys_t* sys) {
    int i;

    sys->ekin = 0.0;
    for (i = 0; i < sys->natoms; ++i) {
        sys->ekin += 0.5 * mvsq2e * sys->mass *
            (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] + sys->vz[i] * sys->vz[i]);
    }
    sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}

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

// must tweak
void force_from_slices(mdsys_t* sys, const int ibegin, const int iend, const int jbegin,
                       const int jend, double* irx, double* iry, double* irz, double* jrx,
                       double* jry, double* jrz, double* epot, double* pfx, double* pfy,
                       double* pfz) {
    for (int i = ibegin; i < iend; ++i) {
        for (int j = jbegin; j < iend; ++j) {
            /* get distance between particle i and j */
            const double rx = pbc(irx[i] - jrx[j], 0.5 * sys->box);
            const double ry = pbc(iry[i] - jry[j], 0.5 * sys->box);
            const double rz = pbc(irz[i] - jrz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                *epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                pfx[i - ibegin] += rx / r * ffac;
                pfy[i - ibegin] += ry / r * ffac;
                pfz[i - ibegin] += rz / r * ffac;
            }
        }
    }
}

void force_mpi_basic(mdsys_t* sys) {
    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    int check = MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);
    check = MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);
    check = MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);

    // this is unbalanced
    for (int i = 0; i < (sys->natoms - 1); i += sys->nranks) {
        int ii = i + sys->rank;
        if (ii >= (sys->natoms - 1))
            break;

        for (int j = ii + 1; j < sys->natoms; ++j) {
            /* get distance between particle i and j */
            const double rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->pfx[ii] += rx / r * ffac;
                sys->pfy[ii] += ry / r * ffac;
                sys->pfz[ii] += rz / r * ffac;

                sys->pfx[j] -= rx / r * ffac;
                sys->pfy[j] -= ry / r * ffac;
                sys->pfz[j] -= rz / r * ffac;
            }
        }
    }

    check =
        MPI_Reduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);
    check =
        MPI_Reduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);
    check =
        MPI_Reduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);
    check = MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm));
    assert(check == MPI_SUCCESS);
}

void force_mpi_ibasic(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check =
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    // this is unbalanced
    for (int i = 0; i < (sys->natoms - 1); i += sys->nranks) {
        int ii = i + sys->rank;
        if (ii >= (sys->natoms - 1))
            break;

        for (int j = ii + 1; j < sys->natoms; ++j) {
            /* get distance between particle i and j */
            const double rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->pfx[ii] += rx / r * ffac;
                sys->pfy[ii] += ry / r * ffac;
                sys->pfz[ii] += rz / r * ffac;

                sys->pfx[j] -= rx / r * ffac;
                sys->pfy[j] -= ry / r * ffac;
                sys->pfz[j] -= rz / r * ffac;
            }
        }
    }

    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm),
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);
}

void force_mpi_ibasic_even(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check =
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    const int ibegin = sys->displs[sys->rank];
    const int iend = ibegin + sys->counts[sys->rank];

    for (int i = ibegin; i < iend; ++i) {
        for (int j = i + 1; j < sys->natoms; ++j) {
            /* get distance between particle i and j */
            const double rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->pfx[i] += rx / r * ffac;
                sys->pfy[i] += ry / r * ffac;
                sys->pfz[i] += rz / r * ffac;

                sys->pfx[j] -= rx / r * ffac;
                sys->pfy[j] -= ry / r * ffac;
                sys->pfz[j] -= rz / r * ffac;
            }
        }
    }

    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm),
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);
}

void force_mpi_primitive(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check =
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    const int begin = sys->rank * (sys->natoms / sys->nranks);
    int end = begin + (sys->natoms / sys->nranks);
    if (end > sys->natoms) {
        end = sys->natoms;
    }
    // const int count = sys->counts[sys->rank];
    // const int begin = sys->displs[sys->rank];
    // const int end = begin + count;

    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    for (int i = begin; i < end; ++i) {
        for (int j = 0; j < (sys->natoms); ++j) {
            if (i == j)
                continue;

            const double rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->pfx[i] += rx / r * ffac;
                sys->pfy[i] += ry / r * ffac;
                sys->pfz[i] += rz / r * ffac;
            }
        }
    }

    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm),
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);
}

void force_mpi_slice(mdsys_t* sys) {
    const int count = sys->counts[sys->rank];
    const int begin = sys->displs[sys->rank];
    const int end = begin + count;

    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check =
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check =
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, count);
    azzero(sys->pfy, count);
    azzero(sys->pfz, count);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    for (int i = begin; i < end; ++i) {
        for (int j = 0; j < sys->natoms; ++j) {
            if (i == j)
                continue;

            /* get distance between particle i and j */
            const double rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->pfx[i - begin] += rx / r * ffac;
                sys->pfy[i - begin] += ry / r * ffac;
                sys->pfz[i - begin] += rz / r * ffac;
            }
        }
    }

    // reconstruct full arrays on master, reduce epot
    nreqs = 0;

    check = MPI_Igatherv(sys->pfx, count, MPI_DOUBLE, sys->fx, sys->counts, sys->displs, MPI_DOUBLE,
                         0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfy, count, MPI_DOUBLE, sys->fy, sys->counts, sys->displs, MPI_DOUBLE,
                         0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfz, count, MPI_DOUBLE, sys->fz, sys->counts, sys->displs, MPI_DOUBLE,
                         0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm),
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);
}

void force_mpi_ring(mdsys_t* sys) {
    const int count = sys->counts[sys->rank];

    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    // init the first send buffers with a scatterv
    int check = MPI_Iscatterv(sys->rx, sys->counts, sys->displs, MPI_DOUBLE, sys->srx, count,
                              MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->ry, sys->counts, sys->displs, MPI_DOUBLE, sys->sry, count,
                          MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->rz, sys->counts, sys->displs, MPI_DOUBLE, sys->srz, count,
                          MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    /*
    int check = MPI_Scatterv(sys->rx, sys->counts, sys->displs, MPI_DOUBLE, sys->srx, count,
    MPI_DOUBLE, 0, *((MPI_Comm*) sys->comm)); assert(check == MPI_SUCCESS); check =
    MPI_Scatterv(sys->ry, sys->counts, sys->displs, MPI_DOUBLE, sys->sry, count, MPI_DOUBLE, 0,
    *((MPI_Comm*) sys->comm)); assert(check
    == MPI_SUCCESS); check = MPI_Scatterv(sys->rz, sys->counts, sys->displs, MPI_DOUBLE, sys->srz,
    count, MPI_DOUBLE, 0, *((MPI_Comm*) sys->comm)); assert(check == MPI_SUCCESS); tcomm +=
    MPI_Wtime();
    */

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, count);
    azzero(sys->pfy, count);
    azzero(sys->pfz, count);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    if (sys->rank) {
        // copy the first sendbuf into own arrays
        for (int i = 0; i < count; ++i) {
            sys->rx[i] = sys->srx[i];
            sys->ry[i] = sys->sry[i];
            sys->rz[i] = sys->srz[i];
        }
    }
    const int nextrank = (sys->rank + 1) % sys->nranks;
    const int prevrank = (sys->rank + sys->nranks - 1) % sys->nranks;
    int recvcount = sys->counts[prevrank];
    int sendcount = count;
    // rank owner of the current slice being processed
    int other = sys->rank;

    for (int iter = 0; iter < sys->nranks - 1; ++iter) {
        int nreqs_iter = 0;
        MPI_Request reqs_iter[6];
        MPI_Status statuses_iter[6];

        check = MPI_Isend(sys->srx, sendcount, MPI_DOUBLE, nextrank, sys->rank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->sry, sendcount, MPI_DOUBLE, nextrank, sys->rank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->srz, sendcount, MPI_DOUBLE, nextrank, sys->rank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        check = MPI_Irecv(sys->rrx, recvcount, MPI_DOUBLE, prevrank, prevrank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rry, recvcount, MPI_DOUBLE, prevrank, prevrank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rrz, recvcount, MPI_DOUBLE, prevrank, prevrank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < sendcount; ++j) {
                if ((other == sys->rank) && (i == j))
                    continue;

                /* get distance between particle i and j */
                const double rx = pbc(sys->rx[i] - sys->srx[j], 0.5 * sys->box);
                const double ry = pbc(sys->ry[i] - sys->sry[j], 0.5 * sys->box);
                const double rz = pbc(sys->rz[i] - sys->srz[j], 0.5 * sys->box);
                const double r = sqrt(rx * rx + ry * ry + rz * rz);

                /* compute force and energy if within cutoff */
                if (r < sys->rcut) {
                    const double ffac = -4.0 * sys->epsilon *
                        (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                    epot += 0.5 * 4.0 * sys->epsilon *
                        (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                    sys->pfx[i] += rx / r * ffac;
                    sys->pfy[i] += ry / r * ffac;
                    sys->pfz[i] += rz / r * ffac;
                }
            }
        }

        check = MPI_Waitall(nreqs_iter, reqs_iter, statuses_iter);
        assert(check == MPI_SUCCESS);

        // swap sr* and rr*
        double* tmpp = sys->srx;
        sys->srx = sys->rrx;
        sys->rrx = tmpp;

        tmpp = sys->sry;
        sys->sry = sys->rry;
        sys->rry = tmpp;

        tmpp = sys->srz;
        sys->srz = sys->rrz;
        sys->rrz = tmpp;

        sendcount = recvcount;
        other = (other + sys->nranks - 1) % sys->nranks;
        recvcount = sys->counts[other];
    }

    // process last sendbuf
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < sendcount; ++j) {
            if ((other == sys->rank) && (i == j))
                continue;

            /* get distance between particle i and j */
            const double rx = pbc(sys->rx[i] - sys->srx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[i] - sys->sry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[i] - sys->srz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->pfx[i] += rx / r * ffac;
                sys->pfy[i] += ry / r * ffac;
                sys->pfz[i] += rz / r * ffac;
            }
        }
    }

    // reconstruct full arrays on master, reduce epot
    nreqs = 0;

    check = MPI_Igatherv(sys->pfx, count, MPI_DOUBLE, sys->fx, sys->counts, sys->displs, MPI_DOUBLE,
                         0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfy, count, MPI_DOUBLE, sys->fy, sys->counts, sys->displs, MPI_DOUBLE,
                         0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfz, count, MPI_DOUBLE, sys->fz, sys->counts, sys->displs, MPI_DOUBLE,
                         0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm),
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);
}

void force_mpi_symmring(mdsys_t* sys) {
    const int count = sys->counts[sys->rank];

    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    // init the first send buffers with a scatterv
    int check = MPI_Iscatterv(sys->rx, sys->counts, sys->displs, MPI_DOUBLE, sys->srx, count,
                              MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->ry, sys->counts, sys->displs, MPI_DOUBLE, sys->sry, count,
                          MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->rz, sys->counts, sys->displs, MPI_DOUBLE, sys->srz, count,
                          MPI_DOUBLE, 0, *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    if (sys->rank) {
        // copy the first sendbuf into own arrays
        for (int i = 0; i < count; ++i) {
            sys->rx[i] = sys->srx[i];
            sys->ry[i] = sys->sry[i];
            sys->rz[i] = sys->srz[i];
        }
    }
    const int nextrank = (sys->rank + 1) % sys->nranks;
    const int prevrank = (sys->rank + sys->nranks - 1) % sys->nranks;
    int recvcount = sys->counts[prevrank];
    int sendcount = count;
    // rank owner of the current slice being processed
    int other = sys->rank;

    const int ibegin = sys->displs[sys->rank];

    for (int iter = 0; iter < sys->nranks - 1; ++iter) {
        int nreqs_iter = 0;
        MPI_Request reqs_iter[6];
        MPI_Status statuses_iter[6];

        check = MPI_Isend(sys->srx, sendcount, MPI_DOUBLE, nextrank, sys->rank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->sry, sendcount, MPI_DOUBLE, nextrank, sys->rank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->srz, sendcount, MPI_DOUBLE, nextrank, sys->rank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        check = MPI_Irecv(sys->rrx, recvcount, MPI_DOUBLE, prevrank, prevrank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rry, recvcount, MPI_DOUBLE, prevrank, prevrank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rrz, recvcount, MPI_DOUBLE, prevrank, prevrank,
                          *((MPI_Comm*)sys->comm), reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        if (other >= sys->rank) {
            const int jbegin = sys->displs[other];

            for (int i = 0; i < count; ++i) {
                for (int j = 0; j < sendcount; ++j) {
                    if ((other == sys->rank) && (i >= j))
                        continue;

                    /* get distance between particle i and j */
                    const double rx = pbc(sys->rx[i] - sys->srx[j], 0.5 * sys->box);
                    const double ry = pbc(sys->ry[i] - sys->sry[j], 0.5 * sys->box);
                    const double rz = pbc(sys->rz[i] - sys->srz[j], 0.5 * sys->box);
                    const double r = sqrt(rx * rx + ry * ry + rz * rz);

                    /* compute force and energy if within cutoff */
                    if (r < sys->rcut) {
                        const double ffac = -4.0 * sys->epsilon *
                            (-12.0 * pow(sys->sigma / r, 12.0) / r +
                             6 * pow(sys->sigma / r, 6.0) / r);

                        epot += 4.0 * sys->epsilon *
                            (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                        sys->pfx[ibegin + i] += rx / r * ffac;
                        sys->pfy[ibegin + i] += ry / r * ffac;
                        sys->pfz[ibegin + i] += rz / r * ffac;

                        sys->pfx[jbegin + j] -= rx / r * ffac;
                        sys->pfy[jbegin + j] -= ry / r * ffac;
                        sys->pfz[jbegin + j] -= rz / r * ffac;
                    }
                }
            }
        }

        check = MPI_Waitall(nreqs_iter, reqs_iter, statuses_iter);
        assert(check == MPI_SUCCESS);

        // swap sr* and rr*
        double* tmpp = sys->srx;
        sys->srx = sys->rrx;
        sys->rrx = tmpp;

        tmpp = sys->sry;
        sys->sry = sys->rry;
        sys->rry = tmpp;

        tmpp = sys->srz;
        sys->srz = sys->rrz;
        sys->rrz = tmpp;

        sendcount = recvcount;
        other = (other + sys->nranks - 1) % sys->nranks;
        recvcount = sys->counts[other];
    }

    // process last sendbuf
    if (other >= sys->rank) {
        const int jbegin = sys->displs[other];

        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < sendcount; ++j) {
                if ((other == sys->rank) && (i >= j))
                    continue;

                /* get distance between particle i and j */
                const double rx = pbc(sys->rx[i] - sys->srx[j], 0.5 * sys->box);
                const double ry = pbc(sys->ry[i] - sys->sry[j], 0.5 * sys->box);
                const double rz = pbc(sys->rz[i] - sys->srz[j], 0.5 * sys->box);
                const double r = sqrt(rx * rx + ry * ry + rz * rz);

                /* compute force and energy if within cutoff */
                if (r < sys->rcut) {
                    const double ffac = -4.0 * sys->epsilon *
                        (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                    epot +=
                        4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                    sys->pfx[ibegin + i] += rx / r * ffac;
                    sys->pfy[ibegin + i] += ry / r * ffac;
                    sys->pfz[ibegin + i] += rz / r * ffac;

                    sys->pfx[jbegin + j] -= rx / r * ffac;
                    sys->pfy[jbegin + j] -= ry / r * ffac;
                    sys->pfz[jbegin + j] -= rz / r * ffac;
                }
            }
        }
    }

    // reconstruct full arrays on master, reduce epot
    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0,
                        *((MPI_Comm*)sys->comm), reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, *((MPI_Comm*)sys->comm),
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);
}

void force_openmp_wnewton(mdsys_t* sys) {
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    double epot = 0.0;
    double* fx = sys->fx;
    double* fy = sys->fy;
    double* fz = sys->fz;

#pragma omp parallel for schedule(dynamic) shared(sys) reduction(+: epot, fx[:sys->natoms], fy[:sys->natoms], fz[:sys->natoms])
    for (int i = 0; i < sys->natoms - 1; ++i) {
        for (int j = i + 1; j < sys->natoms; ++j) {
            // get distance between particle i and j
            const double rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            const double ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            const double rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            const double r = sqrt(rx * rx + ry * ry + rz * rz);

            // compute force and energy if within cutoff
            if (r < sys->rcut) {
                const double ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                epot += 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                fx[i] += rx / r * ffac;
                fy[i] += ry / r * ffac;
                fz[i] += rz / r * ffac;

                fx[j] -= rx / r * ffac;
                fy[j] -= ry / r * ffac;
                fz[j] -= rz / r * ffac;
            }
        }
    }

    sys->epot = epot;
}

void verlet_1(mdsys_t* sys) {
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
}

void verlet_2(mdsys_t* sys) {
    int i;

    /* second part: propagate velocities by another half step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}
