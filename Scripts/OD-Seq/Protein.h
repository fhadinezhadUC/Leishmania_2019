/*
 * Protein.h
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include <vector>
#include <string>

using namespace std;

class Protein{
private:
	string name, sequence;
	vector<bool> binSeq;
	vector<unsigned long int> binRep;
	double outNr;
	vector<unsigned long int> calcBinRep(vector<bool>);
public:
	Protein(string);
	Protein(string, string);
	~Protein();
	string getName() , getSeq();
	vector<bool> getBinSeq();
	void setOutNr(double);
	double getOutNr();
	vector<unsigned long int> getBinRep();
};

#endif /* PROTEIN_H_ */
