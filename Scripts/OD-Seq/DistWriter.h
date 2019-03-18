/*
 * DistWriter.h
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#ifndef DISTWRITER_H_
#define DISTWRITER_H_

#include "AliReader.h"
#include "DistCalc.h"

using namespace std;

class DistWriter{
private:
public:
	DistWriter(AliReader*, DistCalc*, string, string);
};

#endif /* DISTWRITER_H_ */
