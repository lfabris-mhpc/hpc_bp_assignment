#include <gtest/gtest.h>
#include <stdlib.h>

#include <mpi.h>

extern "C" {
#include "utilities.h"
#include "engine.h"
#include "input_output.h"
}

class MPITest : public ::testing::Test {
    protected:
    mdsys_t* sys;

    void SetUp() {
        constexpr int len = 2;

        sys = new mdsys_t;

        sys->comm = MPI_COMM_WORLD;
        MPI_Comm_size(sys->comm, &sys->nranks);
        MPI_Comm_rank(sys->comm, &sys->rank);

        std::cout << "nranks: " << sys->nranks << std::endl;
        std::cout << "rank: " << sys->rank << std::endl;

        sys->displs = new int[sys->nranks];
        sys->counts = new int[sys->nranks];

        // only 2 atoms
        /*
           sys->displs[0] = 0;
        sys->counts[0] = (sys->nranks > 1) ? 1 : len;

        if (sys->nranks > 1) {
            sys->displs[1] = 1;
            sys->counts[1] = 1;
        }
        for (int i = 2; i < sys->nranks; ++i) {
            sys->displs[i] = sys->displs[i - 1] + sys->counts[i - 1];
            sys->counts[i] = 0;
        }
        */

        sys->natoms = len;
        sys->mass = 1.0;
        // TODO: put real values
        sys->nfi = 1;
        sys->dt = 5.0;
        sys->epsilon = 0.2379;
        sys->sigma = 3.405;
        sys->box = 17.1580;
        sys->rcut = 8.5;
        sys->ekin = 0;
        sys->epot = 0;
        sys->temp = 0;

        init_displs_counts(sys->natoms, sys->nranks, sys->displs, sys->counts);

        sys->rx = new double[len];
        sys->ry = new double[len];
        sys->rz = new double[len];

        sys->vx = new double[len];
        sys->vy = new double[len];
        sys->vz = new double[len];

        sys->fx = new double[len];
        sys->fy = new double[len];
        sys->fz = new double[len];

        sys->pfx = new double[len];
        sys->pfy = new double[len];
        sys->pfz = new double[len];

        sys->srx = new double[len];
        sys->sry = new double[len];
        sys->srz = new double[len];

        sys->rrx = new double[len];
        sys->rry = new double[len];
        sys->rrz = new double[len];

        sys->rx[0] = 1.0;
        sys->ry[0] = 1.0;
        sys->rz[0] = 1.0;

        sys->vx[0] = -1.0;
        sys->vy[0] = 1.0;
        sys->vz[0] = 1.0;

        sys->fx[0] = 0.0;
        sys->fy[0] = 0.0;
        sys->fz[0] = 0.0;

        sys->rx[1] = -1.158;
        sys->ry[1] = 0.0;
        sys->rz[1] = 0.0;

        sys->vx[1] = -1.0;
        sys->vy[1] = 0.0;
        sys->vz[1] = 0.0;

        sys->fx[1] = 0.0;
        sys->fy[1] = 0.0;
        sys->fz[1] = 0.0;

        for (int i = 0; i < len; ++i) {
            sys->pfx[i] = 0.0;
            sys->pfy[i] = 0.0;
            sys->pfz[i] = 0.0;

            sys->srx[i] = 0.0;
            sys->sry[i] = 0.0;
            sys->srz[i] = 0.0;

            sys->rrx[i] = 0.0;
            sys->rry[i] = 0.0;
            sys->rrz[i] = 0.0;
        }
    }

    void TearDown() {
        delete[] sys->rx;
        delete[] sys->ry;
        delete[] sys->rz;

        delete[] sys->vx;
        delete[] sys->vy;
        delete[] sys->vz;

        delete[] sys->fx;
        delete[] sys->fy;
        delete[] sys->fz;

        delete[] sys->displs;
        delete[] sys->counts;

        delete[] sys->pfx;
        delete[] sys->pfy;
        delete[] sys->pfz;

        delete sys;
    }
};

TEST_F(MPITest, force_mpi_basic) {
    const double Fx_byhand{ 93.540791569 };
    const double Fy_byhand{ 43.346057261 };
    const double Fz_byhand{ 43.346057261 };

    MPI_Barrier(sys->comm);
    force_mpi_basic(sys);

    if (!sys->rank) {
        std::cout << "testing force_mpi_basic" << std::endl;

        ASSERT_NEAR(sys->fx[0], Fx_byhand, 0.005);
        ASSERT_NEAR(sys->fy[0], Fy_byhand, 0.005);
        ASSERT_NEAR(sys->fz[0], Fz_byhand, 0.005);

        ASSERT_NEAR(sys->fx[1], -Fx_byhand, 0.005);
        ASSERT_NEAR(sys->fy[1], -Fy_byhand, 0.005);
        ASSERT_NEAR(sys->fz[1], -Fz_byhand, 0.005);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    MPI_Init(&argc, &argv);

    int result = 0, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    result = RUN_ALL_TESTS();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return result;
}
