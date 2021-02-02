// -*- c++ -*-

#include <gtest/gtest.h>
#include <stdlib.h>

extern "C" {
#include <engine.h>
#include "input_output.h"
}


TEST(EngineVerlet_1, verlet_1)
{
	int nprint, i, check;
    	char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    	mdsys_t sys;

	check=read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
	assert(check==0);
	
	check = mdsys_init(&sys);
	assert(check==0);

	check = read_restart(&sys, restfile);
	assert(check==0);

	const double VERLET_1_byhand{rxo,ryo,rzo,vxo,vyo,vzo};

	verlet_1(&sys);
	ASSERT_DOUBLE_EQ(sys.velverlet_1->rx[0], VERLET_1_byhand[0]);	
	ASSERT_DOUBLE_EQ(sys.velverlet_1->ry[0], VERLET_1_byhand[1]);
	ASSERT_DOUBLE_EQ(sys.velverlet_1->rz[0], VERLET_1_byhand[2]);

	ASSERT_DOUBLE_EQ(sys.velverlet_1->vx[0], VERLET_1_byhand[3]);	
	ASSERT_DOUBLE_EQ(sys.velverlet_1->vy[0], VERLET_1_byhand[4]);	
	ASSERT_DOUBLE_EQ(sys.velverlet_1->vz[0], VERLET_1_byhand[5]);
	
}

TEST(EngineVerlet_2, verlet_2)
{
	int nprint, i, check;
    	char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    	mdsys_t sys;

	check=read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
	assert(check==0);
	
	check = mdsys_init(&sys);
	assert(check==0);

	check = read_restart(&sys, restfile);
	assert(check==0);

	const double VERLET_2_byhand{vxo,vyo,vzo};

	verlet_2(&sys);

	ASSERT_DOUBLE_EQ(sys.verlet_2->vx[0], VERLET_2_byhand[0]);	
	ASSERT_DOUBLE_EQ(sys.verlet_2->vy[0], VERLET_2_byhand[1]);	
	ASSERT_DOUBLE_EQ(sys.verlet_2->vz[0], VERLET_2_byhand[2]);
	
}
