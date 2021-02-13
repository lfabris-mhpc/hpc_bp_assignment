/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * MPI-enabled c version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <omp.h>
#include <mpi.h>

#include "engine.h"
#include "input_output.h"
#include "defs.h"
#include "utilities.h"

#if !defined(FORCE)
#define FORCE force_mpi_basic
#endif

int main(int argc, char** argv) {
    int mpi_thread_provided;
    // int check = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_provided);
    int check = MPI_Init(&argc, &argv);
    assert(check == MPI_SUCCESS);
    UNUSED(check);

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

    check = mdsys_synch(&sys);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    if (!sys.rank) {
        check = read_restart(&sys, restfile);
        assert(check == 0);
    }

    // initialize forces and energies.
    sys.nfi = 0;
    FORCE(&sys);
    ekin(&sys);

#ifdef VERBOSE
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
#endif

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
