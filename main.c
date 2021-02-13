/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "engine.h"
#include "input_output.h"
#include "defs.h"
#include "utilities.h"

#include <omp.h>

int main(int argc, char** argv) {
    int nprint, check;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *traj, *erg;
    mdsys_t sys;
    double t_start;

//    printf("LJMD version %3.1f\n", LJMD_VERSION);


    t_start = wallclock();

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    check = read_restart(&sys, restfile);
    assert(check == 0);

    // initialize forces and energies.
    sys.nfi = 0;

    force(&sys);

    ekin(&sys);

    erg = fopen(ergfile, "w");
    traj = fopen(trajfile, "w");

    printf("Startup time: %10.3fs\n", wallclock() - t_start);
    printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */
    for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        verlet_1(&sys);

        force(&sys);
       
        verlet_2(&sys);

        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
    fclose(erg);
    fclose(traj);

    check = mdsys_free(&sys);
    assert(check == 0);

    return 0;
}
