/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
//#include <string.h>
//#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
//#include <math.h>
//#include <sys/time.h>

#include "engine.h"
#include "input_output.h"
#include "defs.h"
#include "utilities.h"
#include "engine.h"

#ifdef __MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
    int nprint = 1, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp = NULL, *traj = NULL, *erg = NULL;
    mdsys_t sys;
    double t_start;

#ifdef __MPI
	int ret = MPI_Init(&argc, &argv);
	assert(ret == MPI_SUCCESS);
	int rank, nranks;
	ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(ret == MPI_SUCCESS);
	ret = MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	assert(ret == MPI_SUCCESS);
	
	if (!rank)
#endif
    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

#ifdef __MPI
	if (!rank) {
#endif
    /* read input file */
    if (get_a_line(stdin, line))
        return 1;
    sys.natoms = atoi(line);
    if (get_a_line(stdin, line))
        return 1;
    sys.mass = atof(line);
    if (get_a_line(stdin, line))
        return 1;
    sys.epsilon = atof(line);
    if (get_a_line(stdin, line))
        return 1;
    sys.sigma = atof(line);
    if (get_a_line(stdin, line))
        return 1;
    sys.rcut = atof(line);
    if (get_a_line(stdin, line))
        return 1;
    sys.box = atof(line);
    if (get_a_line(stdin, restfile))
        return 1;
    if (get_a_line(stdin, trajfile))
        return 1;
    if (get_a_line(stdin, ergfile))
        return 1;
    if (get_a_line(stdin, line))
        return 1;
    sys.nsteps = atoi(line);
    if (get_a_line(stdin, line))
        return 1;
    sys.dt = atof(line);
    if (get_a_line(stdin, line))
        return 1;
    nprint = atoi(line);
#ifdef __MPI
	}
	
	//MPI datatype?
	//int natoms, nfi, nsteps;
    //double dt, mass, epsilon, sigma, box, rcut;
	ret = MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.nfi, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	ret = MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ret = MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    /* allocate memory */
    sys.rx = (double*)malloc(sys.natoms * sizeof(double));
    sys.ry = (double*)malloc(sys.natoms * sizeof(double));
    sys.rz = (double*)malloc(sys.natoms * sizeof(double));
    sys.vx = (double*)malloc(sys.natoms * sizeof(double));
    sys.vy = (double*)malloc(sys.natoms * sizeof(double));
    sys.vz = (double*)malloc(sys.natoms * sizeof(double));
    sys.fx = (double*)malloc(sys.natoms * sizeof(double));
    sys.fy = (double*)malloc(sys.natoms * sizeof(double));
    sys.fz = (double*)malloc(sys.natoms * sizeof(double));
#ifdef __MPI
    sys.px = (double*)malloc(sys.natoms * sizeof(double));
    sys.py = (double*)malloc(sys.natoms * sizeof(double));
    sys.pz = (double*)malloc(sys.natoms * sizeof(double));
#endif

#ifdef __MPI
	int restfile_ok = 0;
	if (!rank) {
#endif
    /* read restart */
    // REMEBER TO CHECK THIS VARIABLE, THAT IS USELESS
    int matched;
    fp = fopen(restfile, "r");
    if (fp) {
#ifdef __MPI
		restfile_ok = 1;
#endif
        for (i = 0; i < sys.natoms; ++i) {
            matched = fscanf(fp, "%lf%lf%lf", sys.rx + i, sys.ry + i, sys.rz + i);
            assert(matched == 3);
        }
        for (i = 0; i < sys.natoms; ++i) {
            matched = fscanf(fp, "%lf%lf%lf", sys.vx + i, sys.vy + i, sys.vz + i);
            assert(matched == 3);
        }
        fclose(fp);
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }
#ifdef __MPI
	}

	ret = MPI_Bcast(&restfile_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
	assert(ret == MPI_SUCCESS);

	MPI_Request reqs[3];
	MPI_Status statuses[3];
	if (restfile_ok) {
		ret = MPI_Ibcast(sys.rx, sys.natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, reqs);
		assert(ret == MPI_SUCCESS);
		ret = MPI_Ibcast(sys.ry, sys.natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, reqs + 1);
		assert(ret == MPI_SUCCESS);
		ret = MPI_Ibcast(sys.rz, sys.natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, reqs + 2);
		assert(ret == MPI_SUCCESS);
	}

	azzero(sys.px, sys.natoms);
	azzero(sys.py, sys.natoms);
	azzero(sys.pz, sys.natoms);
	
	if (restfile_ok) {
		ret = MPI_Waitall(3, reqs, statuses);
		assert(ret == MPI_SUCCESS);
	}
#endif

    /* initialize forces and energies.*/
    sys.nfi = 0;
#ifdef __MPI
	force_mpi(&sys, nranks, rank);
#else
    force(&sys);
#endif
	ekin(&sys);

#ifdef __MPI
	if (!rank) {
#endif
    erg = fopen(ergfile, "w");
    traj = fopen(trajfile, "w");
		
    printf("Startup time: %10.3fs\n", wallclock() - t_start);
    printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);
#ifdef __MPI
	}
#endif

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */
    for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {
#ifdef __MPI
	if (!rank) {
#endif
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);
#ifdef __MPI
	}
#endif

        /* propagate system and recompute energies */
#ifdef __MPI
		velverlet_mpi(&sys, nranks, rank);
#else
		velverlet(&sys);
#endif
		ekin(&sys);
    }
    /**************************************************/

#ifdef __MPI
	if (!rank) {
#endif
    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
#ifdef __MPI
	}
#endif
	
#ifdef __MPI
	MPI_Finalize();
#endif

	if (erg)
    	fclose(erg);
	if (traj)
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
#ifdef __MPI
    free(sys.px);
    free(sys.py);
    free(sys.pz);
#endif

    return 0;
}
