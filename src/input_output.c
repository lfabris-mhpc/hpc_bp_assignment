
#include "input_output.h"

/* append data to output. */
void output(mdsys_t* sys, FILE* erg, FILE* traj) {
    int i;

    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot,
           sys->ekin + sys->epot);
    fprintf(erg, "% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin,
            sys->epot, sys->ekin + sys->epot);
    fprintf(traj, "%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin + sys->epot);
    for (i = 0; i < sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
    }
}

/* read input file */
int read_input_file(mdsys_t* sys, FILE* fin, char* line, char* restfile, char* trajfile,
                    char* ergfile, int* nprint) {
    if (get_a_line(fin, line))
        return 1;
    sys->natoms = atoi(line);
    if (get_a_line(fin, line))
        return 1;
    sys->mass = atof(line);
    if (get_a_line(fin, line))
        return 1;
    sys->epsilon = atof(line);
    if (get_a_line(fin, line))
        return 1;
    sys->sigma = atof(line);
    if (get_a_line(fin, line))
        return 1;
    sys->rcut = atof(line);
    if (get_a_line(fin, line))
        return 1;
    sys->box = atof(line);
    if (get_a_line(fin, restfile))
        return 1;
    if (get_a_line(fin, trajfile))
        return 1;
    if (get_a_line(fin, ergfile))
        return 1;
    if (get_a_line(fin, line))
        return 1;
    sys->nsteps = atoi(line);
    if (get_a_line(fin, line))
        return 1;
    sys->dt = atof(line);
    if (get_a_line(fin, line))
        return 1;
    *nprint = atoi(line);
    return 0;
}

int read_restart(mdsys_t* sys, char* restfile) {
    /* read restart */
    int matched, i;
    FILE* fp;
    fp = fopen(restfile, "r");
    if (fp) {
        for (i = 0; i < sys->natoms; ++i) {
            matched = fscanf(fp, "%lf%lf%lf", sys->rx + i, sys->ry + i, sys->rz + i);
            assert(matched == 3);
        }
        for (i = 0; i < sys->natoms; ++i) {
            matched = fscanf(fp, "%lf%lf%lf", sys->vx + i, sys->vy + i, sys->vz + i);
            assert(matched == 3);
        }
        fclose(fp);
        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }
    return 0;
}

int write_restart(mdsys_t* sys, char* restfile) {
    /* write restart */
    FILE* fp;
    fp = fopen(restfile, "w");
    if (fp) {
        for (int i = 0; i < sys->natoms; ++i) {
            fprintf(fp, "%lf  %lf  %lf\n", sys->rx[i], sys->ry[i], sys->rz[i]);
        }
        for (int i = 0; i < sys->natoms; ++i) {
            fprintf(fp, "%lf  %lf  %lf\n", sys->vx[i], sys->vy[i], sys->vz[i]);
        }
        fclose(fp);
    } else {
        perror("cannot open restart file for writing");
        return 3;
    }
    return 0;
}

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE* fp, char* buf) {
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp, BLEN, fp)) {
        int i;

        ptr = strchr(tmp, '#');
        if (ptr)
            *ptr = '\0';
        i = strlen(tmp);
        --i;
        while (isspace(tmp[i])) {
            tmp[i] = '\0';
            --i;
        }
        ptr = tmp;
        while (isspace(*ptr)) {
            ++ptr;
        }
        i = strlen(ptr);
        strcpy(buf, tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}
