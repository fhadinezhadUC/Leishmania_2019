/*
 * runtimeargs.h
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#ifndef RUNTIMEARGS_H_
#define RUNTIMEARGS_H_

#include "util.h"
#include <string>
#include <string.h>
#include <fstream>
#include <omp.h>
#include <algorithm>
#include <iostream>

using namespace std;

typedef struct s_runtime_args {
	string file,resFile,coreFile, outFile,format,metric,mBed,outMode,distMatOut, distMatIn, fileName;
	int threadNr, bootNr, word_size;
	double cutOff;
	char return_flag, help,verbose;
	bool aligned;
} runtime_args;

runtime_args get_runtime_args(int argc, char* argv[]);

#endif /* RUNTIMEARGS_H_ */
