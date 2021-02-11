#ifndef __INPUT_OUTPUT_H__
#define __INPUT_OUTPUT_H__

#include <assert.h>
#include <string.h>
#include "defs.h"
#include "utilities.h"

void output(mdsys_t* sys, FILE* erg, FILE* traj);

int read_input_file(mdsys_t* sys, FILE* fin, char* line, char* restfile, char* trajfile,
                    char* ergfile, int* nprint);

int read_restart(mdsys_t* sys, char* restfile);

int write_restart(mdsys_t* sys, char* restfile);

int get_a_line(FILE* fp, char* buf);

#endif
