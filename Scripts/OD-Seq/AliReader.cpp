/*
 * AliReader.cpp
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#include "AliReader.h"

AliReader::AliReader(runtime_args args) {

	checkFlag = false;
	if (args.format.compare("auto") == 0){
		string line;
		ifstream file(args.file.c_str());
		if(file.is_open()){
			getline(file,line);
			if(line.substr(0,1).compare(">") == 0){
				args.format = "fasta";
				if( args.verbose == 1){
					cout << "Reading in " << args.format << endl;
				}
			} else if(line.substr(0,11).compare("# STOCKHOLM") == 0){
				args.format = "stockholm";
				if( args.verbose == 1){
					cout << "Reading in " << args.format << endl;
				}
			} else if(line.substr(0,7).compare("CLUSTAL") == 0) {
				args.format = "clustal";
				if( args.verbose == 1){
					cout << "Reading in " << args.format << endl;
				}
			} else if(line.substr(0,2).compare("!!") == 0){
				args.format = "msf";
				if( args.verbose == 1){
					cout << "Reading in " << args.format << endl;
				}
			} else if(line.substr(0,4).compare("#=SQ") == 0){
				args.format = "selex";
				if( args.verbose == 1){
					cout << "Reading in " << args.format << endl;
				}
			} else if(line.find_first_of(" ") != line.find_last_of(" ") && isInteger(line.substr(line.find_last_of(" ") +1,line.size())) && isInteger(line.substr(1,line.substr(1, line.size()).find_first_of(" "))) && line.substr(0,1).compare(" ") == 0){
				args.format = "phylip";
				if( args.verbose == 1){
					cout << "Reading in " << args.format << endl;
				}
			} else {
				cout << "Format not recognised: available formats, fasta, vienna, clustal, stockholm, selex, phylip, msf" << endl;
				checkFlag = true;
			}
		}
		file.close();
	}
	if(args.format.compare("fasta") == 0 || args.format.compare("vienna") == 0){
		if( args.verbose == 1 && args.format.compare("fasta") == 0){
			cout << "Reading alignment file in fasta format" << endl;
		} else if(args.verbose == 1 && args.format.compare("vienna") == 0){
			cout << "Reading alignment file in vienna format" << endl;
		}
		ifstream file(args.file.c_str());
		string line, name ,sequence = "";
		if (file.is_open()) {
			while (getline(file, line)) {
				if (line.at(0) == '>') {
					if(name != ""){
						allSeq.push_back(new Protein(name,sequence));
						sequence = "";
						name = "";
					}
					name = line.substr(1, line.length());
				} else {
					sequence = sequence + line;
				}
			}
			allSeq.push_back(new Protein(name,sequence));
		}
		file.close();
	} else if (args.format.compare("clustal") == 0){
		if( args.verbose == 1){
			cout << "Reading alignment file in clustal format" << endl;
		}
		ifstream file(args.file.c_str());
		string line;
		bool flag = false;
		vector<string> names, seqs;
		if(file.is_open()){
			unsigned int count = 0;
			getline(file,line);
			getline(file,line);
			getline(file,line);
			while(getline(file, line)){
				count++;
				if(line.at(0) == ' '){
					count = 0;
					flag = true;
					getline(file,line);
					getline(file,line);
				}
				if(!flag){
					names.push_back(line.substr(0,line.find(" ")));
					seqs.push_back(line.substr(line.find_last_of(" ")+1, line.size()));
				} else {
					seqs.at(count) = seqs.at(count) + line.substr(line.find_last_of(" ")+1, line.size());
				}
			}
		}
		file.close();
		for(unsigned int i = 0 ; i < names.size(); i++){
			allSeq.push_back(new Protein(names.at(i),seqs.at(i)));
		}
	} else if(args.format.compare("stockholm") == 0){
		if( args.verbose == 1){
			cout << "Reading alignment file in stockholm format" << endl;
		}
		ifstream file(args.file.c_str());
		string line;
		bool flag = false;
		vector<string> names, seqs;
		if(file.is_open()){
			unsigned int count = 0;
			getline(file,line);
			getline(file,line);
			while(getline(file, line)){
				if(line == "//"){
					break;
				}
				count++;
				if(line.empty()){
					count = 0;
					flag = true;
					getline(file,line);
				}
				if(!flag){
					names.push_back(line.substr(0,line.find(" ")));
					seqs.push_back(line.substr(line.find_last_of(" ")+1, line.size()));
				} else {
					seqs.at(count) = seqs.at(count) + line.substr(line.find_last_of(" ")+1, line.size());
				}
			}
		}
		file.close();
		for(unsigned int i = 0 ; i < names.size(); i++){
			allSeq.push_back(new Protein(names.at(i),seqs.at(i)));
		}
	} else if(args.format.compare("selex") == 0){
		if( args.verbose == 1){
			cout << "Reading alignment file in selex format" << endl;
		}
		ifstream file(args.file.c_str());
		string line;
		vector<string> names, seqs;
		bool flag = false;
		if(file.is_open()){
			unsigned int count = 0;
			unsigned int count2 = 0;
			while(getline(file, line)){
				if(!line.empty()){
					if(line.substr(0,4).compare("#=SQ") == 0){
						count++;
						names.push_back(line.substr(5,(line.substr(5, line.size())).find_first_of(" ")));
					} else {
						if(!flag){
							for(unsigned int i = 0 ; i < count; i++){
								seqs.push_back("");
							}
							flag = true;
						}
						seqs.at(count2) = seqs.at(count2) + line.substr(line.find_last_of(" ") +1, line.size());
						count2++;
						if(count2 == count){
							count2 = 0;
						}
					}
				}
			}
		}
		file.close();
		for(unsigned int i = 0 ; i < names.size(); i++){
			allSeq.push_back(new Protein(names.at(i),seqs.at(i)));
		}
	} else if(args.format.compare("phylip") == 0){
		if( args.verbose == 1){
			cout << "Reading alignment file in phylip format" << endl;
		}
		ifstream file(args.file.c_str());
		string line;
		vector<string> names, seqs;
		unsigned int count = 0;
		bool flag = false;
		if(file.is_open()){
			getline(file,line);
			unsigned int seqNr = atoi(line.substr(1,(line.substr(1,line.size())).find(" ")).c_str());
			while(getline(file, line)){
				if(line.empty()){
					flag = true;
					getline(file,line);
				}
				if(!flag){
					names.push_back(line.substr(0,10));
					seqs.push_back(line.substr(10,line.size()));
				} else {
					if(count == seqNr){
						count = 0;
					}
					seqs.at(count) = seqs.at(count) + line;
					count++;
				}
			}
		}
		file.close();
		for(unsigned int i = 0 ; i < names.size(); i++){
			allSeq.push_back(new Protein(names.at(i),seqs.at(i)));
		}
	} else if(args.format.compare("msf") == 0){
		if( args.verbose == 1){
			cout << "Reading alignment file in msf format" << endl;
		}
		ifstream file(args.file.c_str());
		string line;
		vector<string> names, seqs;
		unsigned int count = 0;
		if(file.is_open()){
			getline(file,line);
			getline(file,line);
			getline(file,line);
			while(getline(file, line)){
				if(line.substr(0,6).compare(" Name:") == 0){
					names.push_back(line.substr(7,line.substr(7, line.size()).find_first_of(" ")));
					seqs.push_back("");
				} else if(!line.empty() && !(line.compare("//") == 0) && !(line.at(0) == ' ')){
					if(count == names.size()){
						count = 0;
					}
					string tmp;
					line = line.substr(line.find_first_of(" "), line.size());
					remove_copy(line.begin(), line.end(), std::back_inserter(tmp), ' ');
					seqs.at(count) = seqs.at(count) + tmp;
					count++;
				}
			}
		}
		file.close();
		for(unsigned int i = 0 ; i < names.size(); i++){
			allSeq.push_back(new Protein(names.at(i),seqs.at(i)));
		}
	}
	if(!args.aligned){
		for(unsigned int i = 0; i < allSeq.size(); i++){
			replaceAll(allSeq.at(i)->getSeq(), "-", "");
		}
	} else {
		unsigned int alLength =	allSeq.at(0)->getSeq().size();
		for(unsigned int i = 0; i < allSeq.size(); i++){
			if(alLength != allSeq.at(i)->getSeq().size()){
				cout << "Sequences of different length in the input file" << endl;
				checkFlag = true;
				break;
			}
		}
	}
}

AliReader::~AliReader() {
	for (int i = (allSeq.size() - 1); i > -1; i--) {
		delete allSeq.at(i);
	}
}

vector<Protein*> AliReader::getSeqs() {
	return allSeq;
}

bool AliReader::getCheckFlag(){
	return checkFlag;
}
