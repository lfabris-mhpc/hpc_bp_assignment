
#include <math.h>
#include <assert.h>
#include "engine.h"
#include "defs.h"
#include "utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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


#ifdef _OPENMP
// compute forces using threads/
void force_openmp(mdsys_t* sys, const int nthreads, const int tid) {
    double r, ffac;

    // zero energy and forces /
//    sys->epot = 0.0;

    double epot = 0.0;

	#pragma omp parallel reduction(+:epot)
    	{
		int i, j;
		int fromidx, toidx;

	    	double rx, ry, rz;
		double *fx, *fy, *fz;

//		int nthreads = omp_get_num_threads();

//		int tid = omp_get_thread_num();	

		fx = sys->fx + (tid * sys->natoms);
		fy = sys->fy + (tid * sys->natoms);
		fz = sys->fz + (tid * sys->natoms);

		azzero(fx, sys->natoms);
		azzero(fy, sys->natoms);
	 	azzero(fz, sys->natoms);

    
	    //Third Newton law: eliminate if(i==j)/
	    for (i = tid; i < (sys->natoms)-1; i += nthreads) {
        	for (j = i+1; j < (sys->natoms); ++j) {
            // particles have no interactions with themselves 
//            if (i == j)
//                continue;

	            // get distance between particle i and j 
        	    rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
	            ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
        	    rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
	            r = sqrt(rx * rx + ry * ry + rz * rz);

        	    // compute force and energy if within cutoff 
	            if (r < sys->rcut) {
        	        ffac = -4.0 * sys->epsilon *
                	    (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

	                sys->epot += 0.5 * 4.0 * sys->epsilon *
        	            (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

	                fx[i] += rx / r * ffac;
        	        fy[i] += ry / r * ffac;
	                fz[i] += rz / r * ffac;

        	        fx[j] -= rx / r * ffac;
                	fy[j] -= ry / r * ffac;
	                fz[j] -= rz / r * ffac;

		    }
		}
	    }

	#pragma omp barrier

	    i = 1 + (sys->natoms / nthreads);
	    fromidx = tid * i;
	    toidx = fromidx + i;

	    if (toidx > sys->natoms) toidx = sys->natoms;

	    for (i=1; i < nthreads; ++i ) {
		    int offs = i * sys->natoms;
		    for (j = fromidx; j < toidx; ++j) {
			    sys->fx[j] += sys->fx[offs+j];
			    sys->fy[j] += sys->fy[offs+j];
			    sys->fz[j] += sys->fz[offs+j];
		    }
	    }
	    
	}
	   
	sys->epot = epot;
}
#endif

