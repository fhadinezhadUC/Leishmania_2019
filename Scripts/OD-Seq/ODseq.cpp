/*
 * OutDet.cpp
 *
 *  Created on: 8 Apr 2014
 *      Author: peter
 */

#include "ODseq.h"

using namespace std;

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp (){
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main(int argc, char** argv) {

	timestamp_t t0 = get_timestamp();

	runtime_args args;
	args = get_runtime_args(argc, argv);
	if(args.help == 1){
		display_help();
		return 0;
	}
	if(args.return_flag != 0){
		cout << "Run cancelled due to incorrect input" << endl;
		return -1;
	}

	srand(time(NULL));

	if(args.distMatIn.compare("") == 0){
		// Multiple sequence alignment reader
		AliReader* ar = new AliReader(args);
		if(ar->getCheckFlag()){
			return 0;
		}
		if(ar->getSeqs().size() == 0){
			cout << "No sequences read. File is either empty or the format is not recognized" << endl;
			return 0;
		}
		cout << "Alignment read" << endl;
		// Object to calculate the distance matrix and mean distance vector
		DistCalc* dc = new DistCalc(ar->getSeqs(), args);
		cout << "Distance matrix calculated" << endl;
		if(args.distMatOut.compare("") != 0){
			DistWriter* dw = new DistWriter(ar, dc, args.distMatOut, args.mBed);
			cout << "Distance matrix written" << endl;
			delete dw;
		}

		// Bootstrapping the distribution of the distance means
		if(args.outMode.compare("boot") == 0){
			Bootstrap* bs = new Bootstrap(ar->getSeqs(),dc->getMeanDistVec(), args);
			cout << "Bootstrap calculated" << endl;
			if(!args.resFile.empty()){
				ResultWriter* ow = new ResultWriter();
				ow->write(args,ar, bs);
				delete ow;
			}
			if(!args.coreFile.empty()){
				FastaWriter* fw = new FastaWriter(args.coreFile, ar, args.cutOff, true);
				delete fw;
			}
			if(!args.outFile.empty()){
				FastaWriter* fw = new FastaWriter(args.outFile, ar, args.cutOff, false);
				delete fw;
			}
			// Destructor Bootstrap
			delete bs;
		} else {
			IQR* iq = new IQR(dc->getMeanDistVec(), ar->getSeqs());
			cout << "IQR calculated" << endl;
			if(!args.resFile.empty()){
				ResultWriter* ow = new ResultWriter();
				ow->write(args,ar, iq);
				delete ow;
			}
			if(!args.coreFile.empty()){
				FastaWriter* fw = new FastaWriter(args.coreFile, ar, args.cutOff, true);
				delete fw;
			}
			if(!args.outFile.empty()){
				FastaWriter* fw = new FastaWriter(args.outFile, ar, args.cutOff, false);
				delete fw;
			}
			delete iq;
		}
		// Destructor DistCalc
		delete dc;
		//Destructor AliReader
		delete ar;
	} else {
		DistMatReader* dr = new DistMatReader(args.distMatIn);
		cout << "Distance matrix read" << endl;
		if(args.outMode.compare("boot") == 0){
			Bootstrap* bs = new Bootstrap(dr->getSeqs(),dr->getMeanDistVec(), args);
			cout << "Bootstrap calculated" << endl;
			if(!args.resFile.empty()){
				ResultWriter* rw = new ResultWriter();
				rw->write(args, dr, bs);
				delete rw;
			}
			// Destructor Bootstrap
			delete bs;
		} else {
			IQR* iq = new IQR(dr->getMeanDistVec(), dr->getSeqs());
			cout << "IQR calculated" << endl;
			if(!args.resFile.empty()){
				ResultWriter* rw = new ResultWriter();
				rw->write(args, dr, iq);
				delete rw;
			}
			delete iq;
		}
		delete dr;
	}
	if(args.verbose == 1){
		cout << "Program terminated" << endl;
	}
	timestamp_t t1 = get_timestamp();
	double time = (t1 - t0) / 1000000.0L;
	cout << "File name: " << args.fileName << "\nMBed mode: " << args.mBed << "\nMetric: " << args.metric << "\nAligned: " << args.aligned << "\nTime: " << time << "\n";
	return 0;
}
