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
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;

    double epot = 0.0;
    double* fx = sys->fx;
    double* fy = sys->fy;
    double* fz = sys->fz;

#pragma omp parallel for schedule(dynamic) shared(sys, c12, c6, rcsq) reduction(+: epot, fx[:sys->natoms], fy[:sys->natoms], fz[:sys->natoms])
    for (int i = 0; i < sys->natoms - 1; ++i) {
        for (int j = i + 1; j < sys->natoms; ++j) {
            /* get distance between particle i and j */
            double rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            double ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            double rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            double rsq = rx * rx + ry * ry + rz * rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double rinv = 1.0 / rsq;
                double r6 = rinv * rinv * rinv;

                double ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;

                epot += r6 * (c12 * r6 - c6);

                fx[i] += rx * ffac;
                fy[i] += ry * ffac;
                fz[i] += rz * ffac;

                fx[j] -= rx * ffac;
                fy[j] -= ry * ffac;
                fz[j] -= rz * ffac;
            }
        }
    }

    sys->epot = epot;
}

void force_mpi_basic(mdsys_t* sys) {
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

#pragma omp parallel for schedule(dynamic) shared(sys, rx, ry, rz) reduction(+: epot, pfx[:sys->natoms], pfy[:sys->natoms], pfz[:sys->natoms])
    for (int i = 0; i < (sys->natoms - 1); i += sys->nranks) {
        int ii = i + sys->rank;
        if (ii < (sys->natoms - 1)) {
            for (int j = ii + 1; j < sys->natoms; ++j) {
                force_pair_symmetric(sys, rx + ii, ry + ii, rz + ii, rx + j, ry + j, rz + j, &epot,
                                     pfx + ii, pfy + ii, pfz + ii, pfx + j, pfy + j, pfz + j);
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

void force_mpi_ibasic(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check = MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;

#pragma omp parallel for schedule(dynamic) shared(sys, rx, ry, rz) reduction(+: epot, pfx[:sys->natoms], pfy[:sys->natoms], pfz[:sys->natoms])
    for (int i = 0; i < (sys->natoms - 1); i += sys->nranks) {
        int ii = i + sys->rank;
        if (ii >= (sys->natoms - 1)) {
            for (int j = ii + 1; j < sys->natoms; ++j) {
                force_pair_symmetric(sys, rx + ii, ry + ii, rz + ii, rx + j, ry + j, rz + j, &epot,
                                     pfx + ii, pfy + ii, pfz + ii, pfx + j, pfy + j, pfz + j);
            }
        }
    }

    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    UNUSED(check);
}

void force_mpi_ibasic_even(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check = MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    if (!sys->nfi) {
        init_displs_counts_even(sys->natoms, sys->nranks, sys->displs, sys->counts);
    }

    const int ibegin = sys->displs[sys->rank];
    const int iend = ibegin + sys->counts[sys->rank];

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;

#pragma omp parallel for schedule(dynamic) shared(sys, rx, ry, rz) reduction(+: epot, pfx[:sys->natoms], pfy[:sys->natoms], pfz[:sys->natoms])
    for (int i = ibegin; i < iend; ++i) {
        for (int j = i + 1; j < sys->natoms; ++j) {
            force_pair_symmetric(sys, rx + i, ry + i, rz + i, rx + j, ry + j, rz + j, &epot,
                                 pfx + i, pfy + i, pfz + i, pfx + j, pfy + j, pfz + j);
        }
    }

    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    UNUSED(check);
}

void force_mpi_primitive(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check = MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    const int ibegin = sys->displs[sys->rank];
    int iend = ibegin + sys->counts[sys->rank];

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, sys->natoms);
    azzero(sys->pfy, sys->natoms);
    azzero(sys->pfz, sys->natoms);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;

#pragma omp parallel for schedule(dynamic) shared(sys, rx, ry, rz) reduction(+: epot, pfx[:sys->natoms], pfy[:sys->natoms], pfz[:sys->natoms])
    for (int i = ibegin; i < iend; ++i) {
        for (int j = 0; j < (sys->natoms); ++j) {
            if (i == j)
                continue;
            force_pair(sys, rx + i, ry + i, rz + i, rx + j, ry + j, rz + j, &epot, pfx + i, pfy + i,
                       pfz + i);
        }
    }

    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    UNUSED(check);
}

void force_mpi_slice(mdsys_t* sys) {
    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    int check = MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    const int count = sys->counts[sys->rank];
    const int ibegin = sys->displs[sys->rank];
    const int iend = ibegin + count;

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->pfx, count);
    azzero(sys->pfy, count);
    azzero(sys->pfz, count);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;

#pragma omp parallel for schedule(dynamic) shared(sys, ibegin, iend, rx, ry, rz) reduction(+: epot, pfx[:count], pfy[:count], pfz[:count])
    for (int i = ibegin; i < iend; ++i) {
        for (int j = 0; j < sys->natoms; ++j) {
            if (i == j)
                continue;
            force_pair(sys, rx + i, ry + i, rz + i, rx + j, ry + j, rz + j, &epot, pfx + i - ibegin,
                       pfy + i - ibegin, pfz + i - ibegin);
        }
    }

    // reconstruct full arrays on master, reduce epot
    nreqs = 0;

    check = MPI_Igatherv(sys->pfx, count, MPI_DOUBLE, sys->fx, sys->counts, sys->displs, MPI_DOUBLE,
                         0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfy, count, MPI_DOUBLE, sys->fy, sys->counts, sys->displs, MPI_DOUBLE,
                         0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfz, count, MPI_DOUBLE, sys->fz, sys->counts, sys->displs, MPI_DOUBLE,
                         0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    UNUSED(check);
}

void force_mpi_ring(mdsys_t* sys) {
    const int count = sys->counts[sys->rank];

    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    // init the first send buffers with a scatterv
    int check = MPI_Iscatterv(sys->rx, sys->counts, sys->displs, MPI_DOUBLE, sys->srx, count,
                              MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->ry, sys->counts, sys->displs, MPI_DOUBLE, sys->sry, count,
                          MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->rz, sys->counts, sys->displs, MPI_DOUBLE, sys->srz, count,
                          MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    const int nextrank = (sys->rank + 1) % sys->nranks;
    const int prevrank = (sys->rank + sys->nranks - 1) % sys->nranks;
    int recvcount = sys->counts[prevrank];
    int sendcount = count;
    // rank owner of the current slice being processed
    int other = sys->rank;

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

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;

    double* srx = sys->srx;
    double* sry = sys->sry;
    double* srz = sys->srz;

    double* rrx = sys->rrx;
    double* rry = sys->rry;
    double* rrz = sys->rrz;

    for (int iter = 0; iter < sys->nranks - 1; ++iter) {
        int nreqs_iter = 0;
        MPI_Request reqs_iter[6];
        MPI_Status statuses_iter[6];

        check = MPI_Isend(sys->srx, sendcount, MPI_DOUBLE, nextrank, sys->rank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->sry, sendcount, MPI_DOUBLE, nextrank, sys->rank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->srz, sendcount, MPI_DOUBLE, nextrank, sys->rank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        check = MPI_Irecv(sys->rrx, recvcount, MPI_DOUBLE, prevrank, prevrank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rry, recvcount, MPI_DOUBLE, prevrank, prevrank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rrz, recvcount, MPI_DOUBLE, prevrank, prevrank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

#pragma omp parallel for schedule(dynamic) shared(sys, count, sendcount, rx, ry, rz, srx, sry, srz, rrx, rry, rrz) reduction(+: epot, pfx[:sendcount], pfy[:sendcount], pfz[:sendcount])
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < sendcount; ++j) {
                if ((other != sys->rank) || (i != j)) {
                    force_pair(sys, rx + i, ry + i, rz + i, srx + j, sry + j, srz + j, &epot,
                               pfx + i, pfy + i, pfz + i);
                }
            }
        }

        check = MPI_Waitall(nreqs_iter, reqs_iter, statuses_iter);
        assert(check == MPI_SUCCESS);

        // swap sr* and rr*
        double* tmpp = srx;
        srx = rrx;
        rrx = tmpp;

        tmpp = sry;
        sry = rry;
        rry = tmpp;

        tmpp = srz;
        srz = rrz;
        rrz = tmpp;

        sendcount = recvcount;
        other = (other + sys->nranks - 1) % sys->nranks;
        recvcount = sys->counts[other];
    }

// process last sendbuf
#pragma omp parallel for schedule(dynamic) shared(sys, count, sendcount, rx, ry, rz, srx, sry, srz, rrx, rry, rrz) reduction(+: epot, pfx[:sendcount], pfy[:sendcount], pfz[:sendcount])
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < sendcount; ++j) {
            if ((other != sys->rank) || (i != j)) {
                force_pair(sys, rx + i, ry + i, rz + i, srx + j, sry + j, srz + j, &epot, pfx + i,
                           pfy + i, pfz + i);
            }
        }
    }

    // reconstruct full arrays on master, reduce epot
    nreqs = 0;

    check = MPI_Igatherv(sys->pfx, count, MPI_DOUBLE, sys->fx, sys->counts, sys->displs, MPI_DOUBLE,
                         0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfy, count, MPI_DOUBLE, sys->fy, sys->counts, sys->displs, MPI_DOUBLE,
                         0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Igatherv(sys->pfz, count, MPI_DOUBLE, sys->fz, sys->counts, sys->displs, MPI_DOUBLE,
                         0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    UNUSED(check);
}

void force_mpi_symmring(mdsys_t* sys) {
    const int count = sys->counts[sys->rank];

    int nreqs = 0;
    MPI_Request reqs[4];
    MPI_Status statuses[4];

    // init the first send buffers with a scatterv
    int check = MPI_Iscatterv(sys->rx, sys->counts, sys->displs, MPI_DOUBLE, sys->srx, count,
                              MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->ry, sys->counts, sys->displs, MPI_DOUBLE, sys->sry, count,
                          MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Iscatterv(sys->rz, sys->counts, sys->displs, MPI_DOUBLE, sys->srz, count,
                          MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    const int nextrank = (sys->rank + 1) % sys->nranks;
    const int prevrank = (sys->rank + sys->nranks - 1) % sys->nranks;
    int recvcount = sys->counts[prevrank];
    int sendcount = count;
    // rank owner of the current slice being processed
    int other = sys->rank;

    const int ibegin = sys->displs[sys->rank];

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

    double* rx = sys->rx;
    double* ry = sys->ry;
    double* rz = sys->rz;

    double* pfx = sys->pfx;
    double* pfy = sys->pfy;
    double* pfz = sys->pfz;

    double* srx = sys->srx;
    double* sry = sys->sry;
    double* srz = sys->srz;

    double* rrx = sys->rrx;
    double* rry = sys->rry;
    double* rrz = sys->rrz;

    for (int iter = 0; iter < sys->nranks - 1; ++iter) {
        int nreqs_iter = 0;
        MPI_Request reqs_iter[6];
        MPI_Status statuses_iter[6];

        check = MPI_Isend(sys->srx, sendcount, MPI_DOUBLE, nextrank, sys->rank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->sry, sendcount, MPI_DOUBLE, nextrank, sys->rank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Isend(sys->srz, sendcount, MPI_DOUBLE, nextrank, sys->rank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        check = MPI_Irecv(sys->rrx, recvcount, MPI_DOUBLE, prevrank, prevrank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rry, recvcount, MPI_DOUBLE, prevrank, prevrank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);
        check = MPI_Irecv(sys->rrz, recvcount, MPI_DOUBLE, prevrank, prevrank, sys->comm,
                          reqs_iter + nreqs_iter++);
        assert(check == MPI_SUCCESS);

        if (other >= sys->rank) {
            const int jbegin = sys->displs[other];

#pragma omp parallel for schedule(dynamic) shared(sys, count, sendcount, rx, ry, rz, srx, sry, srz, rrx, rry, rrz) reduction(+: epot, pfx[:sendcount], pfy[:sendcount], pfz[:sendcount])
            for (int i = 0; i < count; ++i) {
                for (int j = 0; j < sendcount; ++j) {
                    if ((other == sys->rank) && (i >= j))
                        continue;
                    force_pair_symmetric(sys, rx + i, ry + i, rz + i, srx + j, sry + j, srz + j,
                                         &epot, pfx + ibegin + i, pfy + ibegin + i,
                                         pfz + ibegin + i, pfx + jbegin + j, pfy + jbegin + j,
                                         pfz + jbegin + j);
                }
            }
        }

        check = MPI_Waitall(nreqs_iter, reqs_iter, statuses_iter);
        assert(check == MPI_SUCCESS);

        // swap sr* and rr*
        double* tmpp = srx;
        srx = rrx;
        rrx = tmpp;

        tmpp = sry;
        sry = rry;
        rry = tmpp;

        tmpp = srz;
        srz = rrz;
        rrz = tmpp;

        sendcount = recvcount;
        other = (other + sys->nranks - 1) % sys->nranks;
        recvcount = sys->counts[other];
    }

    // process last sendbuf
    if (other >= sys->rank) {
        const int jbegin = sys->displs[other];

#pragma omp parallel for schedule(dynamic) shared(sys, count, sendcount, rx, ry, rz, srx, sry, srz, rrx, rry, rrz) reduction(+: epot, pfx[:sendcount], pfy[:sendcount], pfz[:sendcount])
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < sendcount; ++j) {
                if ((other == sys->rank) && (i >= j))
                    continue;
                force_pair_symmetric(sys, rx + i, ry + i, rz + i, srx + j, sry + j, srz + j, &epot,
                                     pfx + ibegin + i, pfy + ibegin + i, pfz + ibegin + i,
                                     pfx + jbegin + j, pfy + jbegin + j, pfz + jbegin + j);
            }
        }
    }

    // reconstruct full arrays on master, reduce epot
    nreqs = 0;

    check = MPI_Ireduce(sys->pfx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(sys->pfz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->comm,
                        reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ireduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
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
