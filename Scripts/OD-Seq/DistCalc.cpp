/*
 * DistCalc.cpp
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#include "DistCalc.h"

using namespace std;

DistCalc::DistCalc(vector<Protein*> seqs, runtime_args args) {



	seqNr = seqs.size();
	seedNr = seqNr;
	scoreVec.resize(seqs.at(0)->getSeq().size());
	vector<Protein*> tmpSeqs = seqs;
	if (args.mBed == "mBed") {
		seedNr = round(log(seqNr) / log(2));
		if (seedNr > seqNr) {
			seedNr = seqNr;
		} else if (seedNr < 2){
			seedNr = 2;
		}
		for (unsigned int i = 0; i < seedNr; i++) {
			unsigned int tmp = rand() % tmpSeqs.size();
			randSeq.push_back(tmp);
			seedSeq.push_back(tmpSeqs.at(tmp));
			tmpSeqs.erase(tmpSeqs.begin() + tmp);
		}
		distMat.resize(seqNr);
		for(unsigned int i = 0; i < seqNr; i++){
			distMat.at(i).resize(seedNr);
		}
		if(args.aligned){
			#pragma omp parallel for schedule(static)
			for (unsigned int i = 0; i < seqs.size(); i++) {
				if(args.verbose == 1){
					cout << "\rCalculating distance matrix for sequence " << i+1 << " of " << seqNr << "              " << flush;
				}
				for (unsigned int j = 0; j < seedNr; j++) {
					if(i == j){
						distMat.at(i).at(j) = 0;
					} else if(args.metric == "linear"){
						distMat.at(i).at(j) = calcLinearDist(seqs.at(i), seedSeq.at(j));
					} else if (args.metric == "affine"){
						distMat.at(i).at(j) = calcAffineDist(seqs.at(i), seedSeq.at(j));
					} else if (args.metric == "cumulative"){
						distMat.at(i).at(j) = calcComulativeDist(seqs.at(i), seedSeq.at(j));
					}
				}
			}
		} else {
			string tmpFilePath = random_string(32);
			ofstream tmpFile(tmpFilePath.c_str());
			for(unsigned int i = 0 ; i < seedSeq.size(); i++){
				tmpFile << ">" <<  seedSeq.at(i)->getName() << endl << seedSeq.at(i)->getSeq() << endl;
			}
			string cmd = "blastp -query " + args.file + " -subject " + tmpFilePath + " -outfmt 7 -evalue 100";
			string result = exec(cmd.c_str());
			remove(tmpFilePath.c_str());
			vector<string> hits = tokenizeToString(result, '\n');
			vector<double> scores;
			for(unsigned int i = 0 ; i < hits.size(); i++){
				if((hits.at(i).substr(0,3).compare("# B") != 0) && (hits.at(i).substr(0,3).compare("# F") != 0) && (hits.at(i).substr(0,3).compare("# Q") != 0) && (hits.at(i).substr(0,3).compare("# S") != 0)){
					if(hits.at(i).substr(0,3).compare("# 0") == 0){
						scores.push_back(0);
					} else if(hits.at(i).substr(0,1).compare("#") == 0) {
						double count = tokenizeToDouble(hits.at(i),' ').at(1);
						double tmpScore = 0;
						for(unsigned int j = i+1 ; j < i+count+1; j++){
							if(tokenizeToDouble(hits.at(j), '\t').at(11) > tmpScore){
								tmpScore = tokenizeToDouble(hits.at(j), '\t').at(11);
							}
						}
						scores.push_back(tmpScore);
						i=i+count;
					}
				}
			}
			for(unsigned int i = 0 ; i < distMat.size(); i++){
				for(unsigned int j = 0 ; j < distMat.at(i).size(); j++ ){
					distMat.at(i).at(j) = scores.at((i*seedNr)+j);
				}
			}
		}
	} else {
		if(args.aligned){
			distMat.resize(seqNr);
			for(unsigned int i = 0 ; i < seqNr; i++){
				distMat.at(i).resize(i+1);
			}
			#pragma omp parallel for schedule(dynamic)
			for (unsigned int i = 0; i < seqNr; i++) {
				distMat.at(i).at(i) = 0;
				if(args.verbose == 1){
					cout << "\rCalculating distance matrix for sequence " << i+1 << " of " << seqNr << "              " << flush;
				}
				for (unsigned int j = 0; j < i; j++) {
					if(args.metric == "linear"){
						distMat.at(i).at(j) = calcLinearDist(seqs.at(i), seqs.at(j));
					} else if (args.metric == "affine"){
						distMat.at(i).at(j) = calcAffineDist(seqs.at(i), seqs.at(j));
					} else if (args.metric == "cumulative"){
						distMat.at(i).at(j) = calcComulativeDist(seqs.at(i), seqs.at(j));
					}
				}
			}
		} else {
			distMat.resize(seqNr);
			for(unsigned int i = 0 ; i < seqNr; i++){
				distMat.at(i).resize(seqNr);
			}
			stringstream ss;
			ss << args.word_size;
			string word_size = ss.str();
			string cmd = "blastp -query " + args.file + " -subject " + args.file + " -outfmt 7 -word_size " + word_size;
			string result = exec(cmd.c_str());
			vector<string> hits = tokenizeToString(result, '\n');
			vector<double> scores;
			for(unsigned int i = 0 ; i < hits.size(); i++){
				if((hits.at(i).substr(0,3).compare("# B") != 0) && (hits.at(i).substr(0,3).compare("# F") != 0) && (hits.at(i).substr(0,3).compare("# Q") != 0) && (hits.at(i).substr(0,3).compare("# S") != 0)){
					if(hits.at(i).substr(0,3).compare("# 0") == 0){
						scores.push_back(0);
					} else if(hits.at(i).substr(0,1).compare("#") == 0) {
						double count = tokenizeToDouble(hits.at(i),' ').at(1);
						double tmpScore = 0;
						for(unsigned int j = i+1 ; j < i+count+1; j++){
							if(tokenizeToDouble(hits.at(j), '\t').at(11) > tmpScore){
								tmpScore = tokenizeToDouble(hits.at(j), '\t').at(11);
								}
							}
							scores.push_back(tmpScore);
							i=i+count;
						}
					}
				}
				for(unsigned int i = 0 ; i < distMat.size(); i++){
					for(unsigned int j = 0 ; j < distMat.at(i).size(); j++ ){
						distMat.at(i).at(j) = scores.at((i*distMat.size())+j);
					}
				}
			}
		}
	if(args.verbose == 1){
		cout << endl << "Done" << endl;
	}
	meanDistVec.resize(distMat.size());
	if(args.aligned){
		if(args.mBed == "mBed"){
			for (unsigned int i = 0; i < distMat.size(); i++) {
				double tmp = 0;
				for (unsigned int j = 0; j < distMat.at(i).size(); j++) {
					tmp += distMat.at(i).at(j);
				}
				meanDistVec.at(i) = tmp;
			}
		} else {
			for (unsigned int i = 0; i < distMat.size(); i++) {
				double tmp = 0;
				for (unsigned int j = 0; j < distMat.size(); j++) {
					if(i > j){
						tmp += distMat.at(i).at(j);
					} else {
						tmp += distMat.at(j).at(i);
					}
				}
				meanDistVec.at(i) = tmp;
			}
		}
	} else {
		for (unsigned int i = 0; i < distMat.size(); i++) {
			double tmp = 0;
			for (unsigned int j = 0; j < distMat.at(i).size(); j++) {
				tmp += distMat.at(i).at(j);
			}
			meanDistVec.at(i) = tmp;
		}
	}
}

DistCalc::~DistCalc() {

}

int DistCalc::calcAffineDist(Protein* seq1,Protein* seq2){
	unsigned int count = 0;
	if(seq1->getBinSeq().at(0) == seq2->getBinSeq().at(0)){
		scoreVec.at(0) = 0;
	} else {
		scoreVec.at(0) = 3;
	}
	for (unsigned int k = 1; k < seq1->getSeq().length(); k++) {
		if (seq1->getBinSeq().at(k) == seq2->getBinSeq().at(k)) {
			scoreVec.at(k) = 0;
		} else {
			if( scoreVec.at(k-1) == 0){
				scoreVec.at(k) = 3;
			} else {
				scoreVec.at(k) = 1;
			}
		}
	}
	for (unsigned int k = 0; k < seq1->getSeq().length(); k++) {
		count = count + scoreVec.at(k);
	}
	return count;
}

int DistCalc::calcLinearDist(Protein* seq1,Protein* seq2){

	unsigned int count = 0;
	for (unsigned int k = 0; k < seq1->getBinRep().size(); k++) {
		count = count + decodeBinDist((unsigned long int)seq1->getBinRep().at(k)^(unsigned long int)seq2->getBinRep().at(k));
	}
	return count;
}

int DistCalc::calcComulativeDist(Protein* seq1,Protein* seq2){
	unsigned int count = 0;
	if( seq1->getBinSeq().at(0) == seq2->getBinSeq().at(0)){
		scoreVec.at(0) = 0;
	} else {
		scoreVec.at(0) = 1;
	}
	for (unsigned int k = 1; k < seq1->getSeq().length(); k++) {
		if (seq1->getBinSeq().at(k) == seq2->getBinSeq().at(k)) {
			scoreVec.at(k) = 0;
		} else {
			if( scoreVec.at(k-1) == 0){
				scoreVec.at(k) = 1;
			} else {
				scoreVec.at(k) = scoreVec.at(k-1) + 1;
			}
		}
	}
	for (unsigned int k = 0; k < seq1->getSeq().length(); k++) {
		count = count + scoreVec.at(k);
	}
	return count;
}

int DistCalc::decodeBinDist(unsigned long int a){
	  unsigned long int c; // c accumulates the total bits set in v
	  for (c = 0; a; c++)
	    a &= a - 1; // clear the least significant bit set
	  return c;
}


vector<double> DistCalc::getMeanDistVec() {
	return meanDistVec;
}

vector<vector<double> > DistCalc::getDistMat(){
	return distMat;
}

unsigned int DistCalc::getSeedNr(){
	return seedNr;
}
