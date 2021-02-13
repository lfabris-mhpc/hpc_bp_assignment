
#include "utilities.h"

/* helper function: get current time in seconds since epoch */

double wallclock() {
    struct timeval t;
    gettimeofday(&t, 0);
    return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

/* helper function: zero out an array */
void azzero(double* d, const int n) {
    int i;
    for (i = 0; i < n; ++i) {
        d[i] = 0.0;
    }
}

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2) {
    while (x > boxby2)
        x -= 2.0 * boxby2;
    while (x < -boxby2)
        x += 2.0 * boxby2;
    return x;
}

/*
void init_displs_counts(const int natoms, const int nranks, int* displs, int* counts) {
    assert(displs);
    assert(counts);
    assert(natoms >= 0);
    assert(nranks > 0);

    const int block = natoms / nranks;
    const int rem = natoms % nranks;

    displs[0] = 0;
    counts[0] = block + (rem != 0);
    for (int i = 1; i < nranks; ++i) {
        displs[i] = displs[i - 1] + counts[i - 1];
        counts[i] = block + (i < rem);
    }
}

void init_displs_counts_even(const int natoms, const int nranks, int* displs, int* counts) {
    assert(displs);
    assert(counts);
    assert(natoms >= 0);
    assert(nranks > 0);

    const int target = natoms * (natoms - 1) / (2 * nranks) + 1;

    int iatom = 0;
    displs[0] = 0;
    for (int i = 1; i < nranks; ++i) {
        int work = 0;
        for (counts[i - 1] = 0; work < target && iatom < natoms; ++iatom) {
            work += natoms - iatom - 1;
            ++(counts[i - 1]);
        }
        displs[i] = displs[i - 1] + counts[i - 1];
    }
    counts[nranks - 1] = natoms - displs[nranks - 1];
}
*/