/*
 * DistMatReader.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "DistMatReader.h"

using namespace std;

DistMatReader::DistMatReader(string path){
	ifstream fileIn(path.c_str());
	string line;
	string name;
	if(fileIn.is_open()){
		getline(fileIn,line);
		unsigned int seqNr = atoi(line.c_str());
		meanDistVec.resize(seqNr);
		while(getline(fileIn, line)){
			seqs.push_back(new Protein(line.substr(0,line.find_first_of(" "))));
			line = trim(line.substr(line.find_first_of(" "), line.size()));
			distMat.push_back(tokenizeToDouble(line, ' '));
		}
	} else {
		cout << "File bad" << endl;
	}
	for(unsigned int i = 0; i < distMat.size(); i++){
		meanDistVec.at(i) = 0;
	}
	for(unsigned int i = 0 ; i < distMat.size(); i++){
		for(unsigned int j = 0 ; j < distMat.at(i).size(); j++){
			meanDistVec.at(i) += distMat.at(i).at(j);
		}
		meanDistVec.at(i) = meanDistVec.at(i) / distMat.at(i).size();
	}
}

vector<double> DistMatReader::getMeanDistVec(){
	return meanDistVec;
}

vector<Protein*> DistMatReader::getSeqs(){
	return seqs;
}
