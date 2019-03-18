/*
 * AliReader.h
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#ifndef ALIREADER_H_
#define ALIREADER_H_

#include "Protein.h"
#include "runtimeargs.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

class AliReader {
private:
	vector<Protein*> allSeq;
	bool checkFlag;
public:
	AliReader(runtime_args);
	~AliReader();
	vector<Protein*> getSeqs();
	bool getCheckFlag();
	bool checkSeqs();
};

#endif /* ALIREADER_H_ */
