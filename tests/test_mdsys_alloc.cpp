// -*- c++ -*-

#include <gtest/gtest.h>
#include <stdlib.h>

extern "C" {
#include "defs.h"
#include "input_output.h"
}

TEST(DefsMdsysInit, mdsys_init) {
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    int check, nprint;
    int CHECK;
    mdsys_t sys;

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);

    CHECK = mdsys_init(&sys);
    assert(CHECK == 0);

    mdsys_free(&sys);
}

TEST(Input_outputRead, read_input_file) {
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    int check, nprint;
    int CHECK;
    mdsys_t sys;

    CHECK = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);
    assert(CHECK == 0);

    mdsys_free(&sys);
}

TEST(Input_outputRestart, read_restart) {
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    int check, nprint;
    int CHECK;
    mdsys_t sys;

    check = read_input_file(&sys, stdin, line, restfile, trajfile, ergfile, &nprint);

    CHECK = read_restart(&sys, restfile);
    assert(CHECK == 0);

    mdsys_free(&sys);
}
