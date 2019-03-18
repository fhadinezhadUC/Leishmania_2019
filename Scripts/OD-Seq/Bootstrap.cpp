/*
 * Bootstrap.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "Bootstrap.h"

using namespace std;

Bootstrap::Bootstrap(vector<Protein*> seqs,vector<int> meanDistVec, runtime_args args) {

	meanBootInt.resize(args.bootNr);
	for(int i = 0 ; i < args.bootNr ; i++){
		meanBootInt.at(i).resize(meanDistVec.size());
	}
	drawMean.resize(args.bootNr);
	drawSd.resize(args.bootNr);
	meanSd = 0;
	mean = 0;
	//#pragma omp parallel for schedule(static)
	for (int i = 0; i < args.bootNr; i++) {
		if(args.verbose == 1){
			cout << "\rCalculating bootstrap replicate " << i+1 << " of " << args.bootNr << "              " << flush;
		}
		drawMean.at(i) = 0;
		drawSd.at(i) = 0;
		for (unsigned int j = 0; j < meanDistVec.size(); j++) {
			// Randomly pick the meanDist of a sequence
			meanBootInt.at(i).at(j) = meanDistVec.at(rand() % meanDistVec.size());
			drawMean.at(i) += meanBootInt.at(i).at(j);
		}
		drawMean.at(i) = (double) drawMean.at(i) / (double) meanDistVec.size();
		for( unsigned int j = 0 ; j < meanDistVec.size(); j++){
			drawSd.at(i) += pow((meanBootInt.at(i).at(j) - drawMean.at(i)), 2);
		}
		drawSd.at(i) = sqrt((double)drawSd.at(i) / (double) (meanDistVec.size() - 1));
		meanSd += drawSd.at(i);
		mean += drawMean.at(i);
	}
	meanSd = (meanSd / (double) args.bootNr);
	mean = (mean / (double) args.bootNr);
	if(args.verbose == 1){
		cout << endl << "Done" << endl;
	}
	for (unsigned int i = 0; i < seqs.size(); i++) {
		if( meanSd != 0){
			seqs.at(i)->setOutNr((meanDistVec.at(i) - mean) / meanSd);
		} else {
			seqs.at(i)->setOutNr(0);
		}
	}
}

Bootstrap::Bootstrap(vector<Protein*> seqs,vector<double> meanDistVec, runtime_args args) {

	meanBootDouble.resize(args.bootNr);
	for(int i = 0 ; i < args.bootNr ; i++){
		meanBootDouble.at(i).resize(meanDistVec.size());
	}
	drawMean.resize(args.bootNr);
	drawSd.resize(args.bootNr);
	meanSd = 0;
	mean = 0;
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < args.bootNr; i++) {
		if(args.verbose == 1){
			cout << "\rCalculating bootstrap replicate " << i+1 << " of " << args.bootNr << "              " << flush;
		}
		drawMean.at(i) = 0;
		drawSd.at(i) = 0;
		for (unsigned int j = 0; j < meanDistVec.size(); j++) {
			// Randomly pick the meanDist of a sequence
			meanBootDouble.at(i).at(j) = meanDistVec.at(rand() % meanDistVec.size());
			drawMean.at(i) += meanBootDouble.at(i).at(j);
		}
		drawMean.at(i) = (double) drawMean.at(i) / (double) meanDistVec.size();
		for( unsigned int j = 0 ; j < meanDistVec.size(); j++){
			drawSd.at(i) += pow((meanBootDouble.at(i).at(j) - drawMean.at(i)), 2);
		}
		drawSd.at(i) = sqrt((double)drawSd.at(i) / (double) (meanDistVec.size() - 1));
		meanSd += drawSd.at(i);
		mean += drawMean.at(i);
	}
	meanSd = (meanSd / (double) args.bootNr);
	mean = (mean / (double) args.bootNr);
	if(args.verbose == 1){
		cout << endl << "Done" << endl;
	}
	for (unsigned int i = 0; i < seqs.size(); i++) {
		if( meanSd != 0){
			seqs.at(i)->setOutNr((meanDistVec.at(i) - mean) / meanSd);
		} else {
			seqs.at(i)->setOutNr(0);
		}
	}
}

Bootstrap::~Bootstrap() {

}

double Bootstrap::getMeanSd() {
	return meanSd;
}

double Bootstrap::getMean() {
	return mean;
}


