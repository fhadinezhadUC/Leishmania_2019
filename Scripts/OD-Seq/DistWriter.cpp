/*
 * DistWriter.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "DistWriter.h"

using namespace std;

DistWriter::DistWriter(AliReader* ar, DistCalc* dc, string path, string mBed){
	ofstream distFile(path.c_str());
	distFile << ar->getSeqs().size() << endl;
	for(unsigned int i = 0; i < ar->getSeqs().size(); i++){
		distFile << ar->getSeqs().at(i)->getName() << "\t";
		if(mBed.compare("mBed") == 0 ){
			for(unsigned int j = 0 ; j < dc->getDistMat().at(i).size(); j++){
				distFile << dc->getDistMat().at(i).at(j) << " ";
			}
		} else {
			for(unsigned int j = 0 ; j < ar->getSeqs().size(); j++){
				if(j < dc->getDistMat().at(i).size()){
					distFile << dc->getDistMat().at(i).at(j) << " ";
				} else {
					distFile << dc->getDistMat().at(j).at(i) << " ";
				}
			}
		}
		distFile << endl;
	}
	distFile.close();
}
