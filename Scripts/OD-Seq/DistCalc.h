/*
 * DistCalc.h
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#ifndef DISTCALC_H_
#define DISTCALC_H_

#include "Protein.h"
#include "runtimeargs.h"
#include "PairwiseAl.h"
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class DistCalc {
private:
	vector<vector<double> > distMat;
	vector<double> meanDistVec;
	vector<int> scoreVec;
	unsigned int seqNr, seedNr;
	vector<Protein*> seedSeq;
	vector<unsigned int> randSeq;
	int calcAffineDist(Protein*, Protein*);
	int calcComulativeDist(Protein*, Protein*);
	int calcLinearDist(Protein*, Protein*);
	int decodeBinDist(unsigned int long);
	int calcMeanDist();
public:
	DistCalc(vector<Protein*>, runtime_args);
	~DistCalc();
	vector<double> getMeanDistVec();
	vector<vector<double> > getDistMat();
	unsigned int getSeedNr();
};

#endif /* DISTCALC_H_ */
