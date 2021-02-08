
#include "defs.h"

/* allocate memory */
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

    return ret;
}

int mdsys_free(mdsys_t* sys) {
    int ret = 0;

    free(sys->rx);
    ret += sys->rx != NULL;
    free(sys->ry);
    ret += sys->ry != NULL;
    free(sys->rz);
    ret += sys->rz != NULL;
    free(sys->vx);
    ret += sys->vx != NULL;
    free(sys->vy);
    ret += sys->vy != NULL;
    free(sys->vz);
    ret += sys->vz != NULL;
    free(sys->fx);
    ret += sys->fx != NULL;
    free(sys->fy);
    ret += sys->fy != NULL;
    free(sys->fz);
    ret += sys->fz != NULL;

    return ret;
}
