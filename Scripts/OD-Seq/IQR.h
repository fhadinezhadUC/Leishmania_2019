/*
 * IQR.h
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#ifndef IQR_H_
#define IQR_H_

#include <vector>
#include <algorithm>
#include "Protein.h"

using namespace std;


class IQR {
private:
	int median;
	int Q1;
	int Q3;
	int min;
	int max;
	int range;
public:
	IQR(vector<int> meanDistVec, vector<Protein*> seqs);
	IQR(vector<double> meanDistVec, vector<Protein*> seqs);
	~IQR();
	int getQ1();
	int getQ3();
	int getRange();
	int getMin();
	int getMax();
	int getMedian();
};


#endif /* IQR_H_ */
