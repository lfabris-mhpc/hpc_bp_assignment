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

    //init_displs_counts(sys->natoms, sys->nranks, sys->displs, sys->counts);

    sys->pfx = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->pfx == NULL;
    sys->pfy = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->pfy == NULL;
    sys->pfz = (double*)malloc(sys->natoms * sizeof(double));
    ret += sys->pfz == NULL;
	
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
}

int mdsys_synch(mdsys_t* sys) {
    int check;
    int nreqs = 0;
    MPI_Request reqs[9];
    MPI_Status statuses[9];

    // broadcasts of mdsys parameters
    // int natoms, nfi, nsteps;
    check = MPI_Ibcast(&(sys->natoms), 1, MPI_INT, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->nfi), 1, MPI_INT, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->nsteps), 1, MPI_INT, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    // double dt, mass, epsilon, sigma, box, rcut;
    check = MPI_Ibcast(&(sys->dt), 1, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->mass), 1, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->epsilon), 1, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->sigma), 1, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->box), 1, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys->rcut), 1, MPI_DOUBLE, 0, sys->comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    return check != MPI_SUCCESS;
}
