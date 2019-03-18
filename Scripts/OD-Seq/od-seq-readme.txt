To compile OD-seq on linux use the following command.

g++ -fopenmp -o OD-seq AliReader.cpp Bootstrap.cpp DistCalc.cpp DistMatReader.cpp DistWriter.cpp FastaWriter.cpp IQR.cpp ODseq.cpp PairwiseAl.cpp Protein.cpp ResultWriter.cpp runtimeargs.cpp util.cpp

./OD-seq --help will show a description of the available input parameters.

A first test run could be ./OD-sed -i <alignment file> -o <result file>

For further questions send me an email:

peter.jehl@ucdconnect.ie
