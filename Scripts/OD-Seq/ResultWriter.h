/*
 * ResultWriter.h
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#ifndef RESULTWRITER_H_
#define RESULTWRITER_H_

#include "DistMatReader.h"
#include "IQR.h"
#include "runtimeargs.h"
#include "AliReader.h"
#include "Bootstrap.h"

using namespace std;

class ResultWriter{
private:
public:
	ResultWriter();
	void write(runtime_args, AliReader*, Bootstrap*);
	void write(runtime_args, AliReader*, IQR*);
	void write(runtime_args, DistMatReader*, Bootstrap*);
	void write(runtime_args, DistMatReader*, IQR*);
	~ResultWriter();
	void closeFile();
};


#endif /* RESULTWRITER_H_ */
