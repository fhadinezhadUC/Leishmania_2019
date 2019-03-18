/*
 * DistMatReader.h
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#ifndef DISTMATREADER_H_
#define DISTMATREADER_H_

#include "Protein.h"
#include "util.h"
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

class DistMatReader {
private:
	vector<vector<double> >distMat;
	vector<double> meanDistVec;
	vector<Protein*> seqs;
public:
	DistMatReader(string);
	vector<double> getMeanDistVec();
	vector<Protein*> getSeqs();
};

#endif /* DISTMATREADER_H_ */
