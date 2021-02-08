// -*- c++ -*-

#include <gtest/gtest.h>
#include <stdlib.h>

extern "C" {
#include "engine.h"
#include "input_output.h"
}

TEST(EngineEkin, ekin) {
    int nprint, i, check;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    mdsys_t sys;

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    check = read_restart(&sys, restfile);
    assert(check == 0);

    const double Ek_byhand{ 4780.114723067 };

    ekin(&sys);
    ASSERT_DOUBLE_EQ(sys.ekin, Ek_byhand);

    check = mdsys_free(&sys);
	assert(check == 0);
}

TEST(EngineVerlet_1, verlet_1) {
    int nprint, i, check;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    mdsys_t sys;

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    check = read_restart(&sys, restfile);
    assert(check == 0);

    // values computed by hand
    const double VERLET_1_byhand[6]{ 0.012246398,  0.005674874,   0.005674874,
                                     0.0024492796, 0.00113497486, 0.00113497486 };

    verlet_1(&sys);
    ASSERT_DOUBLE_EQ(sys.rx[0], VERLET_1_byhand[0]);
    ASSERT_DOUBLE_EQ(sys.ry[0], VERLET_1_byhand[1]);
    ASSERT_DOUBLE_EQ(sys.rz[0], VERLET_1_byhand[2]);

    ASSERT_DOUBLE_EQ(sys.vx[0], VERLET_1_byhand[3]);
    ASSERT_DOUBLE_EQ(sys.vy[0], VERLET_1_byhand[4]);
    ASSERT_DOUBLE_EQ(sys.vz[0], VERLET_1_byhand[5]);
}

TEST(EngineVerlet_2, verlet_2) {
    int nprint, i, check;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    mdsys_t sys;

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    check = read_restart(&sys, restfile);
    assert(check == 0);

    // values computed by hand
    const double VERLET_2_byhand[3]{ 0.0024492796, 0.00113497486, 0.00113497486 };

    verlet_2(&sys);

    ASSERT_DOUBLE_EQ(sys.vx[0], VERLET_2_byhand[0]);
    ASSERT_DOUBLE_EQ(sys.vy[0], VERLET_2_byhand[1]);
    ASSERT_DOUBLE_EQ(sys.vz[0], VERLET_2_byhand[2]);

    check = mdsys_free(&sys);
	assert(check == 0);
}

TEST(EngineForce, force) {
    int nprint, i, check;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    mdsys_t sys;

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(check == 0);

    check = mdsys_init(&sys);
    assert(check == 0);

    check = read_restart(&sys, restfile);
    assert(check == 0);

    const double Fx_byhand{ 93.540791569 };
    const double Fy_byhand{ 43.346057261 };
    const double Fz_byhand{ 43.346057261 };

    force(&sys);
    ASSERT_DOUBLE_EQ(sys.fx[0], Fx_byhand);
    ASSERT_DOUBLE_EQ(sys.fy[0], Fy_byhand);
    ASSERT_DOUBLE_EQ(sys.fz[0], Fz_byhand);
    ASSERT_DOUBLE_EQ(sys.fx[1], -Fx_byhand);
    ASSERT_DOUBLE_EQ(sys.fy[1], -Fy_byhand);
    ASSERT_DOUBLE_EQ(sys.fz[1], -Fz_byhand);

    check = mdsys_free(&sys);
	assert(check == 0);
}
