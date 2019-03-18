/*
 * FastaWriter.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "FastaWriter.h"

using namespace std;

FastaWriter::FastaWriter(string path, AliReader* ar, double cutOff, bool core){
	ofstream fastaFile(path.c_str());
	for(unsigned int i = 0; i < ar->getSeqs().size(); i++){
		if(core){
			if(cutOff > ar->getSeqs().at(i)->getOutNr()){
				fastaFile << ">" << ar->getSeqs().at(i)->getName() << endl;
				fastaFile << replaceAll(ar->getSeqs().at(i)->getSeq(), "-", "") << endl;
			}
		} else {
			if(cutOff < ar->getSeqs().at(i)->getOutNr()){
				fastaFile << ">" << ar->getSeqs().at(i)->getName() << endl;
				fastaFile << replaceAll(ar->getSeqs().at(i)->getSeq(), "-", "") << endl;
			}
		}
	}
	fastaFile.close();
}
