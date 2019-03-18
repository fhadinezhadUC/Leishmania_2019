/*
 * ResultWriter.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "ResultWriter.h"

using namespace std;

ResultWriter::ResultWriter(){

}

void ResultWriter::write(runtime_args args, AliReader* ar, Bootstrap* bs){

	ofstream resFile(args.resFile.c_str());
	resFile << "File name: " << args.fileName << endl << "Outlier score: " << args.cutOff << endl << "MBed: " << args.mBed << endl << "Metric mode: " << args.metric << endl;
	resFile << "Sequences: " << ar->getSeqs().size() << endl << "Length: " << ar->getSeqs().at(0)->getSeq().length() << endl;
	resFile << "Bootstrap mean: " << bs->getMean() << endl	<< "Bootstrap Sd: " << bs->getMeanSd() << endl;
	resFile << "SdNr:" << endl;
	for (unsigned int i = 0; i < ar->getSeqs().size(); i++) {
		resFile << ar->getSeqs().at(i)->getName() << ": " << ar->getSeqs().at(i)->getOutNr() << endl;
		if(args.verbose == 1){
			cout << ar->getSeqs().at(i)->getName() << ": " << ar->getSeqs().at(i)->getOutNr() << endl;
		}
	}
	resFile << "//" << endl;
	resFile.close();
	cout << "Results printed" << endl;
}

void ResultWriter::write(runtime_args args, AliReader* ar, IQR* iq){

	ofstream resFile(args.resFile.c_str());
	resFile << "File name: " << args.fileName << endl << "Outlier score: " << args.cutOff << endl << "MBed: " << args.mBed << endl << "Metric mode: " << args.metric << endl;
	resFile << "Sequences: " << ar->getSeqs().size() << endl << "Length: " << ar->getSeqs().at(0)->getSeq().length() << endl;
	resFile << "IQR median: " << iq->getMedian() << endl << "IQR range: " << iq->getRange() << endl;
	resFile << "SdNr:" << endl;
	for (unsigned int i = 0; i < ar->getSeqs().size(); i++) {
		resFile << ar->getSeqs().at(i)->getName() << ": " << ar->getSeqs().at(i)->getOutNr() << endl;

		if(args.verbose == 1){
			cout << ar->getSeqs().at(i)->getName() << ": " << ar->getSeqs().at(i)->getOutNr() << endl;
		}
	}
	resFile << "//" << endl;
	resFile.close();
	cout << "Results printed" << endl;
}

void ResultWriter::write(runtime_args args, DistMatReader* dc, Bootstrap* bs){

	ofstream resFile(args.resFile.c_str());
	resFile << "File name: " << args.fileName << endl << "Outlier score: " << args.cutOff << endl << "MBed: " << args.mBed << endl << "Metric mode: " << args.metric << endl << "Sequences: " << dc->getSeqs().size() << endl;
	resFile << "Bootstrap mean: " << bs->getMean() << endl	<< "Bootstrap Sd: " << bs->getMeanSd() << endl;
	resFile << "SdNr:" << endl;
	for (unsigned int i = 0; i < dc->getSeqs().size(); i++) {
		resFile << dc->getSeqs().at(i)->getName() << ": " << dc->getSeqs().at(i)->getOutNr() << endl;
		if(args.verbose == 1){
			cout << dc->getSeqs().at(i)->getName() << ": " << dc->getSeqs().at(i)->getOutNr() << endl;
		}
	}
	resFile << "//" << endl;
	resFile.close();
	cout << "Results printed" << endl;
}

void ResultWriter::write(runtime_args args, DistMatReader* dc, IQR* iq){

	ofstream resFile(args.resFile.c_str());
	resFile << "File name: " << args.fileName << endl << "Outlier score: " << args.cutOff << endl << "MBed: " << args.mBed << endl << "Metric mode: " << args.metric << endl << "Sequences: " << dc->getSeqs().size() << endl;
	resFile << "IQR median: " << iq->getMedian() << endl << "IQR range: " << iq->getRange() << endl;
	resFile << "SdNr:" << endl;
	for (unsigned int i = 0; i < dc->getSeqs().size(); i++) {
		resFile << dc->getSeqs().at(i)->getName() << ": " << dc->getSeqs().at(i)->getOutNr() << endl;
		if(args.verbose == 1){
			cout << dc->getSeqs().at(i)->getName() << ": " << dc->getSeqs().at(i)->getOutNr() << endl;
		}
	}
	resFile << "//" << endl;
	resFile.close();
	cout << "Results printed" << endl;
}

ResultWriter::~ResultWriter(){

}
