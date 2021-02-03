
#include "defs.h"

/* allocate memory */
int mdsys_init(mdsys_t* sys) {
    int ret = 0;

    sys->rx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->ry = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->rz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->vx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->vy = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->vz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->fx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->fy = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;
    sys->fz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->rx == NULL;

    return ret;
}
