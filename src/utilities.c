
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

void init_displs_counts(const int natoms, const int nranks, int* displs, int* counts) {
    assert(displs);
    assert(counts);
    assert(natoms >= 0);
    assert(nranks > 0);
    assert(natoms >= nranks);

    const int block = natoms / nranks;
    const int rem = natoms % nranks;

    counts[nranks - 1] = block + (rem != 0);
    displs[nranks - 1] = natoms - counts[nranks - 1];

    for (int i = 1; i < nranks; ++i) {
        counts[nranks - 1 - i] = block + (i < rem);
        displs[nranks - 1 - i] = displs[nranks - i] - counts[nranks - 1 - i];
    }
    displs[0] = 0;
}

void init_displs_counts_even(const int natoms, const int nranks, int* displs, int* counts) {
    assert(displs);
    assert(counts);
    assert(natoms >= 0);
    assert(nranks > 0);
    assert(natoms >= nranks);

    const int target = natoms * (natoms - 1) / (2 * nranks) + 1;

    int atom = natoms;
    for (int i = nranks - 1; i > 0; --i) {
        counts[i] = 0;
        displs[i] = atom;

        int iworkload = 0;
        while (iworkload < target) {
            ++(counts[i]);
            --atom;
            iworkload += natoms - 1 - atom;
        }
        displs[i] = atom;
    }

    displs[0] = 0;
    counts[0] = displs[1];
}