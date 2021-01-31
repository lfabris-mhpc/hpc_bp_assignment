// -*- c++ -*-

#include <gtest/gtest.h>
#include <stdlib.h>

extern "C" {
#include <engine.h>
#include "input_output.h"
}


TEST(EngineEkin, ekin)
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
	
	const double Ek_byhand{4780.114723067};

	ekin(&sys);
	ASSERT_DOUBLE_EQ(sys.ekin, Ek_byhand);
	
}
