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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    int nprint, check;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *traj, *erg;
    mdsys_t sys;
    double t_start;
    int nthreads, tid;

#ifdef _OMP
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num();
    printf("Using nthreads=%d, tid=%d\n", nthreads, tid);
#endif

    printf("LJMD version %3.1f\n", LJMD_VERSION);


    t_start = wallclock();

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    check = read_restart(&sys, restfile);
    assert(check == 0);

#ifdef _OMP    
    //allocate memory
    #pragma omp parallel

    sys.rx=(double*)malloc(sys.natoms*sizeof(double));
    sys.ry=(double*)malloc(sys.natoms*sizeof(double));
    sys.rz=(double*)malloc(sys.natoms*sizeof(double));

    sys.vx=(double*)malloc(sys.natoms*sizeof(double));
    sys.vy=(double*)malloc(sys.natoms*sizeof(double));
    sys.vz=(double*)malloc(sys.natoms*sizeof(double));

    sys.fx=(double*)malloc(nthreads*sys.natoms*sizeof(double));
    sys.fy=(double*)malloc(nthreads*sys.natoms*sizeof(double));
    sys.fz=(double*)malloc(nthreads*sys.natoms*sizeof(double));
#endif

    // initialize forces and energies.
    sys.nfi = 0;
#ifdef _OPENMP
    force_openmp(&sys, nthreads, tid);
#else
    force(&sys);
#endif

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
#ifdef _OPENMP
        force_openmp(&sys, nthreads, tid);
#else
        force(&sys);
#endif
        verlet_2(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
    fclose(erg);
    fclose(traj);

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    return 0;
}
