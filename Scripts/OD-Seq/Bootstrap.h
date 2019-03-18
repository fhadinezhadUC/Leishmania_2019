/*
 * Bootstrap.h
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#ifndef BOOTSTRAP_H_
#define BOOTSTRAP_H_

using namespace std;

#include "Protein.h"
#include "runtimeargs.h"
#include <vector>
#include <math.h>

class Bootstrap {
private:
	vector<vector<unsigned int> > meanBoot;
	vector<vector<int> > meanBootInt;
	vector<vector<double> > meanBootDouble;
	vector<double> drawMean;
	vector<double> drawSd;
	double meanSd, mean;
public:
	Bootstrap(vector<Protein*>,vector<int>, runtime_args);
	Bootstrap(vector<Protein*>,vector<double>, runtime_args);
	~Bootstrap();
	double getMeanSd();
	double getMean();
};

#endif /* BOOTSTRAP_H_ */
