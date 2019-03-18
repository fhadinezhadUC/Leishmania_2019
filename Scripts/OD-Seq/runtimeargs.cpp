/*
 * runtimeargs.cpp
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#include "runtimeargs.h"


runtime_args get_runtime_args(int argc, char* argv[]){

	runtime_args args;

	// default values
	args.file = "";
	if(!(args.file.compare("") == 0)){
		args.fileName = args.file.substr(args.file.find_last_of("/") + 1,args.file.length());
	}
	args.resFile = "";
	args.coreFile = "";
	args.outFile = "";
	args.distMatIn = "";
	args.format = "fasta";
	args.bootNr = 1000;
	args.metric = "linear";
	args.mBed = "mBed";
	args.outMode = "boot";
	args.distMatOut = "";
	args.threadNr = 1;
	args.cutOff = 2;
	args.aligned = true;
	args.word_size = 4;
	args.verbose = 0;
	args.help = 0;
	args.return_flag = 0;
	string score = "2";
	bool boot = false;
	bool metric = false;


	for(int i = 1; i < argc; i=i+2){
		if((0 == strcmp(argv[i], "-i")) || (0 == strcmp(argv[i], "--input"))){
			// input file
			args.file = argv[i+1];
			args.fileName = args.file.substr(args.file.find_last_of("/") + 1,args.file.length());
		} else if((0 == strcmp(argv[i], "-r")) || (0 == strcmp(argv[i], "--result"))){
			// result file
			args.resFile = argv[i+1];
		} else if((0 == strcmp(argv[i], "-o")) || (0 == strcmp(argv[i], "--outlier"))){
			// outlier file
			args.outFile = argv[i+1];
		} else if((0 == strcmp(argv[i], "-c")) || (0 == strcmp(argv[i], "--core"))){
			// core file
			args.coreFile = argv[i+1];
		} else if((0 == strcmp(argv[i], "-s")) || (0 == strcmp(argv[i], "--score"))){
			// outlier score
			score = argv[i+1];
		} else if((0 == strcmp(argv[i], "-f")) || (0 == strcmp(argv[i], "--format"))){
			// input format
			args.format = argv[i+1];
		} else if((0 == strcmp(argv[i], "--dist-in"))){
			// input distance matrix
			args.distMatIn = argv[i+1];
			args.fileName = args.distMatIn.substr(args.distMatIn.find_last_of("/") + 1,args.distMatIn.length());
		} else if((0 == strcmp(argv[i], "-a")) || (0 == strcmp(argv[i], "--analysis"))){
			// outlier analysis mode
			args.outMode = argv[i+1];
		} else if((0 == strcmp(argv[i], "--boot-rep"))){
			// number of bootstrap pseudoreplicates
			boot = true;
			args.bootNr = atoi(argv[i+1]);
		} else if((0 == strcmp(argv[i], "-m")) || (0 == strcmp(argv[i], "--metric"))){
			// distance metric
			metric = true;
			args.metric = argv[i+1];
		} else if((0 == strcmp(argv[i], "-t")) || (0 == strcmp(argv[i], "--threads"))){
			// number of threads
			args.threadNr = atoi(argv[i+1]);
		} else if((0 == strcmp(argv[i], "--full"))){
			// full distance matrix
			args.mBed = "full";
			i--;
		} else if((0 == strcmp(argv[i], "--dist-out"))){
			// matrix output
			args.distMatOut = argv[i+1];
		} else if((0 == strcmp(argv[i], "-ws"))){
			//blastp wordsize
			args.word_size = atoi(argv[i+1]);
		} else if((0 == strcmp(argv[i], "-h")) || (0 == strcmp(argv[i], "--help"))){
			// help page
			args.help = 1;
			i--;
		} else if((0 == strcmp(argv[i], "--unaligned"))){
			// aligned
			args.aligned = false;
			i--;
		}else if((0 == strcmp(argv[i], "-v")) || (0 == strcmp(argv[i], "--verbose"))){
			// verbose mode
			args.verbose = 1;
			i--;
		} else {
			cout << argv[i] << ": option not recognised" << endl;
			args.return_flag = -1;
		}
	}

	if(args.help == 0){
		// test input files
		if(!args.file.empty()){
			ifstream inFile(args.file.c_str());
			if (!inFile.is_open()) {
				cout << "Alignment input file not available" << endl;
				args.return_flag = -1;
			}
			inFile.close();
		}
		if(!args.distMatIn.empty()) {
			ifstream inFile(args.distMatIn.c_str());
			if (!inFile.is_open()) {
				cout << "Distance input file not available" << endl;
				args.return_flag = -1;
			}
			inFile.close();
		}
		// test output file
		if(args.resFile.empty() && args.coreFile.empty() && args.outFile.empty()){
			cout << "Result, outlier and core file not set. At least one of them is needed" << endl;
			args.return_flag = -1;
		} else {
			if(!args.outFile.empty()){
				ofstream outFile(args.outFile.c_str());
				if (!outFile.is_open()) {
					cout << "Outlier file not available" << endl;
					args.return_flag = -1;
				}
				outFile.close();
			}
			if(!args.coreFile.empty()){
				ofstream coreFile(args.coreFile.c_str());
				if (!coreFile.is_open()) {
					cout << "Core file not available" << endl;
					args.return_flag = -1;
				}
				coreFile.close();
			}
			if(!args.resFile.empty()){
				ofstream resFile(args.resFile.c_str());
				if (!resFile.is_open()) {
					cout << "Result file not available" << endl;
					args.return_flag = -1;
				}
				resFile.close();
			}
		}

		// test input format
		if(args.format.compare("fa") == 0){
			args.format = "fasta";
		} else if (args.format.compare("clu") == 0){
			args.format = " clustal";
		} else if (args.format.compare("st") == 0){
			args.format = "stockholm";
		} else if (args.format.compare("vie") == 0){
			args.format = "vienna";
		} else if (args.format.compare("phy") == 0){
			args.format = "phylip";
		}
		if(args.format.compare("auto") == 0 || args.format.compare("fasta") == 0 || args.format.compare("clustal") == 0 || args.format.compare("stockholm") == 0 || args.format.compare("vienna") == 0 || args.format.compare("msf") == 0 || args.format.compare("phylip") == 0 || args.format.compare("selex") == 0){
		} else {
			cout << "Input format invalid. Use fa[sta],clu[stal],st[ockholm],vie[nna],msf,phy[lip],selex (default: auto)" << endl;
			args.return_flag = -1;
		}

		// test threads
		if(args.threadNr < 1 || args.threadNr > omp_get_max_threads()){
			cout << "Using maximum number of threads: " << omp_get_max_threads() << endl;
			omp_set_num_threads(omp_get_max_threads());
		} else {
			omp_set_num_threads(args.threadNr);
		}

		// test bootstrap replicates
		if ((! isdigit(args.bootNr) == 0) || (args.bootNr < 1)) {
			cout << "Number of pseudoreplicates is not numeric or < 1" << endl;
			args.return_flag = -1;
		}

		//blastp wordsize
		if ((! isdigit(args.word_size) == 0) || (args.word_size < 2)) {
			cout << "Word size for blastp is either not an integer or < 2" << endl;
			args.return_flag = -1;
		}

		// test metric
		if(!((args.metric.compare("linear") == 0) || (args.metric.compare("affine") == 0) || (args.metric.compare("cumulative") == 0))){
			cout << "No valid metric mode: linear, affine, cumulative" << endl;
			args.return_flag = -1;
		}

		// test analysis
		if(args.outMode.compare("bootstrap") == 0){
			args.outMode = "boot";
		}
		if(!((args.outMode.compare("boot") == 0) || (args.outMode.compare("iqr") == 0))){
			cout << "No valid outlier analysis mode: boot or iqr" << endl;
			args.return_flag = -1;
		}

		// test distance matrix output
		if(!args.distMatOut.empty()){
			ofstream outDistFile(args.distMatOut.c_str());
			if (!outDistFile.is_open()) {
				cout << "Distance matrix output file not available" << endl;
				args.return_flag = -1;
			}
			outDistFile.close();
		}

		// test cutOff
		if(isInteger(score) || isFloat(score)){
			args.cutOff = atof(score.c_str());
			if(args.cutOff < 0){
				cout << "Score is negative" << endl;
				args.return_flag = -1;
			}
		} else {
			cout << "Score is not numeric" << endl;
			args.return_flag = -1;
		}

		if(!args.distMatIn.empty() && (!args.outFile.empty() || !args.coreFile.empty() || !args.file.empty() || (args.format.compare("auto") != 0) || (args.mBed.compare("normal") == 0) || !args.distMatOut.empty())){
			if(!args.outFile.empty()){
				cout << "Distance matrix input selected. Can't write outlier file." << endl;
			}
			if(!args.coreFile.empty()){
				cout << "Distance matrix input selected. Can't write core file." << endl;
			}
			if(!args.file.empty()){
				cout << "Distance matrix input and alignment input selected. Can't read both files." << endl;
			}
			if(args.format.compare("auto") != 0){
				cout << "Format and distance matrix input set. Can't apply format to distance matrix." << endl;
			}
			if(args.mBed.compare("normal") == 0){
				cout << "Full distance matrix mode and distance matrix input set. Reading distance matrix is always in full mode." << endl;
			}
			if(!args.distMatOut.empty()){
				cout << "Distance matrix input and output not possible together" << endl;
			}
			args.return_flag = -1;
		}

		if((args.outMode.compare("iqr") == 0) && boot){
			cout << "Number of bootstrap replicates not used in inter quantile range analysis." << endl;
			args.return_flag = -1;
		}

		if(metric && !args.distMatIn.empty()){
			cout << "Metric not used when distance matrix is read." << endl;
			args.return_flag = -1;
		}
	}
	return args;
}
