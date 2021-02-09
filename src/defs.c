#include "defs.h"
#include "utilities.h"

int mdsys_init(mdsys_t* sys) {
    int ret = 0;

    sys->rx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->ry = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->ry == NULL;
    sys->rz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rz == NULL;

    sys->vx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->vx == NULL;
    sys->vy = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->vy == NULL;
    sys->vz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->vz == NULL;

    sys->fx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->fx == NULL;
    sys->fy = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->fy == NULL;
    sys->fz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->fz == NULL;

    // MPI
    sys->displs = (int*)malloc(sys->nranks * sizeof(int));
    ret += sys->displs == NULL;
    sys->counts = (int*)malloc(sys->nranks * sizeof(int));
    ret += sys->counts == NULL;

    // init_displs_counts(sys->natoms, sys->nranks, sys->displs, sys->counts);
    init_displs_counts_even(sys->natoms, sys->nranks, sys->displs, sys->counts);

    sys->pfx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->pfx == NULL;
    sys->pfy = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->pfy == NULL;
    sys->pfz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->pfz == NULL;

    sys->srx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->srx == NULL;
    sys->sry = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->sry == NULL;
    sys->srz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->srz == NULL;

    sys->rrx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rrx == NULL;
    sys->rry = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rry == NULL;
    sys->rrz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rrz == NULL;

    return ret;
}

void mdsys_free(mdsys_t* sys) {
    free(sys->rx);
    sys->rx = NULL;
    free(sys->ry);
    sys->ry = NULL;
    free(sys->rz);
    sys->rz = NULL;

    free(sys->vx);
    sys->vx = NULL;
    free(sys->vy);
    sys->vy = NULL;
    free(sys->vz);
    sys->vz = NULL;

    free(sys->fx);
    sys->fx = NULL;
    free(sys->fy);
    sys->fy = NULL;
    free(sys->fz);
    sys->fz = NULL;

    // MPI
    free(sys->displs);
    sys->displs = NULL;
    free(sys->counts);
    sys->counts = NULL;

    free(sys->pfx);
    sys->pfx = NULL;
    free(sys->pfy);
    sys->pfy = NULL;
    free(sys->pfz);
    sys->pfz = NULL;

    free(sys->srx);
    sys->srx = NULL;
    free(sys->sry);
    sys->sry = NULL;
    free(sys->srz);
    sys->srz = NULL;

    free(sys->rrx);
    sys->rrx = NULL;
    free(sys->rry);
    sys->ry = NULL;
    free(sys->rrz);
    sys->rrz = NULL;
}
