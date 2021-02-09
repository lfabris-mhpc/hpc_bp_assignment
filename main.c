/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include "engine.h"
#include "input_output.h"
#include "defs.h"
#include "utilities.h"

// force_mpi_basic
// force_mpi_ibasic
// force_mpi_ibasic_even
// force_mpi_primitive
// force_mpi_slice
// force_mpi_ring
// force_mpi_symmring
#define FORCE force_mpi_ibasic_even

int main(int argc, char** argv) {
    int check = MPI_Init(&argc, &argv);
    assert(check == MPI_SUCCESS);

    mdsys_t sys;
    sys.comm = MPI_COMM_WORLD;

    check = MPI_Comm_size(sys.comm, &sys.nranks);
    assert(check == MPI_SUCCESS);
    check = MPI_Comm_rank(sys.comm, &sys.rank);
    assert(check == MPI_SUCCESS);

    if (!sys.rank) {
        printf("LJMD version %3.1f\n", LJMD_VERSION);
    }

    double t_start = wallclock();

    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    int nprint;
    if (!sys.rank) {
        check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
        assert(check == 0);
    }

    // TODO: refactor in a mdsys_synch
    int nreqs = 0;
    MPI_Request reqs[9];
    MPI_Status statuses[9];
    // broadcasts of mdsys parameters
    // int natoms, nfi, nsteps;
    check = MPI_Ibcast(&(sys.natoms), 1, MPI_INT, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.nfi), 1, MPI_INT, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.nsteps), 1, MPI_INT, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    // double dt, mass, epsilon, sigma, box, rcut;
    check = MPI_Ibcast(&(sys.dt), 1, MPI_DOUBLE, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.mass), 1, MPI_DOUBLE, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.epsilon), 1, MPI_DOUBLE, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.sigma), 1, MPI_DOUBLE, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.box), 1, MPI_DOUBLE, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);
    check = MPI_Ibcast(&(sys.rcut), 1, MPI_DOUBLE, 0, sys.comm, reqs + nreqs++);
    assert(check == MPI_SUCCESS);

    check = MPI_Waitall(nreqs, reqs, statuses);
    assert(check == MPI_SUCCESS);

    check = mdsys_init(&sys);
    assert(check == 0);

    if (!sys.rank) {
        printf("displs: [");
        for (int i = 0; i < sys.nranks; ++i) {
            printf("%d, ", sys.displs[i]);
        }
        printf("]\n");
        printf("counts: [");
        for (int i = 0; i < sys.nranks; ++i) {
            printf("%d, ", sys.counts[i]);
        }
        printf("]\n");
    }

    if (!sys.rank) {
        check = read_restart(&sys, restfile);
        assert(check == 0);
    }

    /* initialize forces and energies.*/
    sys.nfi = 0;
    FORCE(&sys);
    ekin(&sys);

    FILE* erg = NULL;
    FILE* traj = NULL;
    if (!sys.rank) {
        erg = fopen(ergfile, "w");
        traj = fopen(trajfile, "w");

        printf("Startup time: %10.3fs\n", wallclock() - t_start);
        printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);

        /* reset timer */
        t_start = wallclock();
    }

    /**************************************************/
    /* main MD loop */
    for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        if (!sys.rank) {
            /* write output, if requested */
            if ((sys.nfi % nprint) == 0)
                output(&sys, erg, traj);

            /* propagate system and recompute energies */
            verlet_1(&sys);
        }

        FORCE(&sys);

        if (!sys.rank) {
            verlet_2(&sys);
            ekin(&sys);
        }
    }

    /**************************************************/
    if (!sys.rank) {
        /* clean up: close files, free memory */
        printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
        fclose(erg);
        fclose(traj);
    }

    MPI_Finalize();

    mdsys_free(&sys);

    return 0;
}
