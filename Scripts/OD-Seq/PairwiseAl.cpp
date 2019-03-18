/*
 * PairwiseAl.cpp
 *
 *  Created on: 9 Apr 2014
 *      Author: peter
 */

#include "PairwiseAl.h"

using namespace std;

PairwiseAl::PairwiseAl(Protein* seq1, Protein* seq2){

	map<char, map<char,int> > blosum;
	vector<char> bloChar;
	bloChar.push_back('A');
	bloChar.push_back('R');
	bloChar.push_back('N');
	bloChar.push_back('D');
	bloChar.push_back('C');
	bloChar.push_back('Q');
	bloChar.push_back('E');
	bloChar.push_back('G');
	bloChar.push_back('H');
	bloChar.push_back('I');
	bloChar.push_back('L');
	bloChar.push_back('K');
	bloChar.push_back('M');
	bloChar.push_back('F');
	bloChar.push_back('P');
	bloChar.push_back('S');
	bloChar.push_back('T');
	bloChar.push_back('W');
	bloChar.push_back('Y');
	bloChar.push_back('V');
	bloChar.push_back('B');
	bloChar.push_back('Z');
	bloChar.push_back('X');
	bloChar.push_back('-');

	int bloNumAr[] = {4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,
	-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,
	-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,
	-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,
	 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,
	-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,
	-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
	 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,
	-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,
	-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,
	-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,
	-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,
	-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,
	-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,
	-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,
	 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,
	 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,
	-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,
	-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,
	 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,
	-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,
	-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
	 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,
	-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1};

	vector<int> bloNum (bloNumAr, bloNumAr + sizeof(bloNumAr) / sizeof(int));

	map<char,int> tmp;
	for(unsigned int i = 0 ; i < bloChar.size(); i++){
		tmp.clear();
		for(unsigned int j = 0 ; j < bloChar.size(); j++){
			tmp.insert(pair<char,int>(bloChar.at(j), bloNum.at(i * bloChar.size() + j)));
		}
		blosum.insert(pair<char,map<char,int> >(bloChar.at(i), tmp));
	}

	int GAP = -4;

	dynMat.resize(seq1->getSeq().size()+1);
	for(unsigned int i = 0 ; i < seq1->getSeq().size()+1; i++){
		dynMat.at(i).resize(seq2->getSeq().size()+1);
		for(unsigned int j = 0 ; j < seq2->getSeq().size()+1 ; j++){
			dynMat.at(i).at(j) = 0;
		}
	}

	for(unsigned int i = 0; i < seq1->getSeq().size()+1; i++){
		dynMat.at(i).at(0) = i * GAP;
	}
	for(unsigned int i = 1 ; i < seq2->getSeq().size()+1; i++){
		dynMat.at(0).at(i) = i * GAP;
	}
	for(unsigned int i = 1 ; i < seq1->getSeq().size()+1; i++){
		for(unsigned int j = 1 ; j < seq2->getSeq().size()+1; j++){
			dynMat.at(i).at(j) = max((dynMat.at(i-1).at(j-1) + blosum.find(seq1->getSeq().at(i-1))->second.find(seq2->getSeq().at(j-1))->second), (dynMat.at(i).at(j-1) + GAP), (dynMat.at(i-1).at(j) + GAP));
		}
	}

	score = dynMat.at(seq1->getSeq().size()).at(seq2->getSeq().size());
}

unsigned int PairwiseAl::getScore(){

	return score;
}

PairwiseAl::~PairwiseAl(){

}
