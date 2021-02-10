
#include <math.h>
#include <assert.h>
#include "engine.h"
#include "defs.h"
#include "utilities.h"

#include <omp.h>

/* a few physical constants */
const double kboltz = 0.0019872067; /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* compute kinetic energy */
void ekin(mdsys_t* sys) {
    int i;

    sys->ekin = 0.0;
    for (i = 0; i < sys->natoms; ++i) {
        sys->ekin += 0.5 * mvsq2e * sys->mass *
            (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] + sys->vz[i] * sys->vz[i]);
    }
    sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}

/* compute forces */
void force(mdsys_t* sys) {
    double r, ffac;
    double rx, ry, rz;
    int i, j;

    /* zero energy and forces */
    sys->epot = 0.0;
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);
    
    /*Third Newton law: eliminate if(i==j)*/
    for (i = 0; i < (sys->natoms)-1; ++i) {
        for (j = i+1; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
//            if (i == j)
//                continue;

            /* get distance between particle i and j */
            rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
            r = sqrt(rx * rx + ry * ry + rz * rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0 * sys->epsilon *
                    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

                sys->epot += 0.5 * 4.0 * sys->epsilon *
                    (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

                sys->fx[i] += rx / r * ffac;
                sys->fy[i] += ry / r * ffac;
                sys->fz[i] += rz / r * ffac;
	
                sys->fx[j] -= rx / r * ffac;
                sys->fy[j] -= ry / r * ffac;
                sys->fz[j] -= rz / r * ffac;
            }
        }
    }
}



/* velocity verlet */
void verlet_1(mdsys_t* sys) {
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt * sys->vx[i];
        sys->ry[i] += sys->dt * sys->vy[i];
        sys->rz[i] += sys->dt * sys->vz[i];
    }
}

void verlet_2(mdsys_t* sys) {
    int i;

    /* second part: propagate velocities by another half step */
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}


// compute forces using threads/
void force_openmp(mdsys_t* sys) {

    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    double epot = 0.0;
    double *fx = sys->fx;
    double *fy = sys->fy;
    double *fz = sys->fz;

	#pragma omp parallel shared(sys) reduction(+:epot, fx[:sys->natoms], fy[:sys->natoms], fz[:sys->natoms])
    	{
		const int nthreads = omp_get_num_threads();
		const int tid = omp_get_thread_num();	
    
	    //Third Newton law: eliminate if(i==j)/
	    for (int i = 0; i < (sys->natoms)-1; i += nthreads) {
		    int ii = i + tid;
		    if ( ii >= (sys->natoms) - 1)
			    break;

        	for (int j = ii + 1; j < (sys->natoms); ++j) {
            // particles have no interactions with themselves 
//            if (i == j)
//                continue;

	            // get distance between particle i and j 
        	    const double rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
	            const double ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
        	    const double rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
	            const double r = sqrt(rx * rx + ry * ry + rz * rz);

        	    // compute force and energy if within cutoff 
	            if (r < sys->rcut) {
        	        const double ffac = -4.0 * sys->epsilon *
                	    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

	                epot += 4.0 * sys->epsilon * (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

	                fx[ii] += rx / r * ffac;
        	        fy[ii] += ry / r * ffac;
	                fz[ii] += rz / r * ffac;

        	        fx[j] -= rx / r * ffac;
                	fy[j] -= ry / r * ffac;
	                fz[j] -= rz / r * ffac;

		    }
		}
	    }
	 
	}
	   
	sys->epot = epot;
}


