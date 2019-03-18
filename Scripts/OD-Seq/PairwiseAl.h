/*
 * PairwiseAl.h
 *
 *  Created on: 9 Apr 2014
 *      Author: peter
 */

#ifndef PAIRWISEAL_H_
#define PAIRWISEAL_H_

#include "Protein.h"
#include "util.h"
#include <iostream>
#include <map>

using namespace std;

class PairwiseAl {
private:
	int score;
	vector<vector<int> > dynMat;
public:
	PairwiseAl(Protein* , Protein*);
	~PairwiseAl();
	unsigned int getScore();
};


#endif /* PAIRWISEAL_H_ */
