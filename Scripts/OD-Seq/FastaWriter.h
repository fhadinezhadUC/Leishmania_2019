/*
 * FastaWriter.h
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#ifndef FASTAWRITER_H_
#define FASTAWRITER_H_

#include "AliReader.h"
#include "util.h"

using namespace std;

class FastaWriter{
private:
public:
	FastaWriter(string, AliReader*, double, bool);
};


#endif /* FASTAWRITER_H_ */
