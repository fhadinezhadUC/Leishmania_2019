#include "Protein.h"

using namespace std;


Protein::Protein(string nameIn){

	outNr = 0;
	name = nameIn;
	sequence = "";
}

Protein::Protein(string nameIn, string sequenceIn) {

	outNr = 0;
	name = nameIn;
	sequence = sequenceIn;
	binSeq.resize(sequence.length());
	for (unsigned i = 0; i < sequence.length(); ++i) {
		if (sequence.at(i) == '-' || sequence.at(i) == '.' || sequence.at(i) == '~') {
			binSeq.at(i) = 1;
		} else {
			binSeq.at(i) = 0;
		}
	}
	binRep = calcBinRep(binSeq);
}

vector<unsigned long int> Protein::calcBinRep(vector<bool> a){

	vector<unsigned long int> bin;
	unsigned short count = 0;
	unsigned long int tmp = 0;
	for(unsigned int i = 0 ; i < a.size(); i++){
		tmp = tmp << 1;
		if(a.at(i) == 1){
			tmp++;
		}
		count++;
		if(count > (sizeof(unsigned long int)*8)-2){
			count = 0;
			bin.push_back(tmp);
			tmp = 0;
		}
	}
	if(count>0){
		bin.push_back(tmp);
	}
	return bin;
}

Protein::~Protein() {

}

string Protein::getName() {
	return name;
}

string Protein::getSeq() {
	return sequence;
}

vector<bool> Protein::getBinSeq() {
	return binSeq;
}

void Protein::setOutNr(double out) {
	outNr = out;
}

double Protein::getOutNr() {
	return outNr;
}

vector<unsigned long int> Protein::getBinRep(){
	return binRep;
}
