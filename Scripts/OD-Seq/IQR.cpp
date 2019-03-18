/*
 * IQR.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "IQR.h"

using namespace std;

IQR::IQR(vector<int> meanDistVec, vector<Protein*> seqs){

	unsigned int size = seqs.size();
	vector<int> tmpDistVec(meanDistVec);
	make_heap(tmpDistVec.begin(), tmpDistVec.end());
	sort_heap(tmpDistVec.begin(), tmpDistVec.end());
	min = tmpDistVec.at(0);
	max = tmpDistVec.at(tmpDistVec.size()-1);
	if(size%2 == 0){
		median = (tmpDistVec[size/2] + tmpDistVec[(size/2)-1]) / 2;
		if((size/2) % 2 == 0){
			Q1 = (tmpDistVec[size/4] + tmpDistVec[(size/4)-1]) / 2;
			Q3 = (tmpDistVec[(size/4)*3] + tmpDistVec[((size/4)*3)-1]) / 2;
		} else {
			Q1 = tmpDistVec[((size/2)-1)/2];
			Q3 = tmpDistVec[((size/2) + (size-1))/2];
		}
	} else {
		median = tmpDistVec[(size-1)/2];
		if(((size-1)/2) % 2 == 0){
			Q1 = (tmpDistVec[(size-1)/4] + tmpDistVec[((size-1)/4)-1]) /2;
			Q3 = (tmpDistVec[(size-1)/4*3] + tmpDistVec[((size-1)/4*3)+1]) /2;
		} else {
			Q1 = tmpDistVec[(((size-1)/2)-1)/2];
			Q3 = tmpDistVec[(size-1) - ((((size-1)/2)-1)/2)];
		}
	}
	range = Q3 - Q1;
	for(unsigned int i = 0 ; i < size; i++){
		if(meanDistVec.at(i) < Q1){
			if(range == 0){
				seqs.at(i)->setOutNr(99999);
			} else {
				seqs.at(i)->setOutNr((double)(Q1 - meanDistVec.at(i))/ (double)range);
			}
		} else if (meanDistVec.at(i) > Q3){
			if(range == 0){
				seqs.at(i)->setOutNr(99999);
			} else {
				seqs.at(i)->setOutNr((double)(meanDistVec.at(i) - Q3)/ (double)range);
			}
		} else {
			seqs.at(i)->setOutNr(0);
		}
	}
}

IQR::IQR(vector<double> meanDistVec, vector<Protein*> seqs){
	unsigned int size = seqs.size();
	vector<double> tmpDistVec(meanDistVec);
	make_heap(tmpDistVec.begin(), tmpDistVec.end());
	sort_heap(tmpDistVec.begin(), tmpDistVec.end());
	min = tmpDistVec.at(0);
	max = tmpDistVec.at(tmpDistVec.size()-1);
	if(size%2 == 0){
		median = (tmpDistVec[size/2] + tmpDistVec[(size/2)-1]) / 2;
		if((size/2) % 2 == 0){
			Q1 = (tmpDistVec[size/4] + tmpDistVec[(size/4)-1]) / 2;
			Q3 = (tmpDistVec[(size/4)*3] + tmpDistVec[((size/4)*3)-1]) / 2;
		} else {
			Q1 = tmpDistVec[((size/2)-1)/2];
			Q3 = tmpDistVec[((size/2) + (size-1))/2];
		}
	} else {
		median = tmpDistVec[(size-1)/2];
		if(((size-1)/2) % 2 == 0){
			Q1 = (tmpDistVec[(size-1)/4] + tmpDistVec[((size-1)/4)-1]) /2;
			Q3 = (tmpDistVec[(size-1)/4*3] + tmpDistVec[((size-1)/4*3)+1]) /2;
		} else {
			Q1 = tmpDistVec[(((size-1)/2)-1)/2];
			Q3 = tmpDistVec[(size-1) - ((((size-1)/2)-1)/2)];
		}
	}
	range = Q3 - Q1;
	for(unsigned int i = 0 ; i < size; i++){
		if(meanDistVec.at(i) < Q1){
			if(range == 0){
				seqs.at(i)->setOutNr(99999);
			} else {
				seqs.at(i)->setOutNr((Q1 - meanDistVec.at(i))/ range);
			}
		} else if (meanDistVec.at(i) > Q3){
			if(range == 0){
				seqs.at(i)->setOutNr(99999);
			} else {
				seqs.at(i)->setOutNr((meanDistVec.at(i) - Q3)/ range);
			}
		} else {
			seqs.at(i)->setOutNr(0);
		}
	}
}

IQR::~IQR(){

}

int IQR::getQ1(){
	return Q1;
}

int IQR::getQ3(){
	return Q3;
}

int IQR::getRange(){
	return range;
}

int IQR::getMin(){
	return min;
}

int IQR::getMax(){
	return max;
}

int IQR::getMedian(){
	return median;
}
