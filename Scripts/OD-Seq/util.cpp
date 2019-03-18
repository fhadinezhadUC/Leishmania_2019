/*
 * util.cpp
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#include "util.h"

bool isInteger(const std::string & s){
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;
   char * p ;
   strtol(s.c_str(), &p, 10) ;
   return (*p == 0) ;
}

bool isFloat(string myString) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f;
    return iss.eof() && !iss.fail();
}

void display_help(){
	cout << "OutDet 1.0" << endl << endl;
	cout << "Usage: OutDet [-hv] [-i <file>] [-r <file>] [-c <file>] [-o <file>] [-f {fa[sta],clu[stal],st[ockholm],vie[nna],msf,phy[lip],selex] [--dist-in <file>] [-a {boot,iqr}] [--boot-rep <n>] [-m {linear,affine,cumulative}] [-t <n>] [--full] [--matrix-out <file>] [-s <n>]" << endl << endl;
	cout << "Input options:" << endl;
	cout << "  -i, --input             Multiple sequence alignment file" << endl;
	cout << "  -f, --format            Multiple sequence alignment format (fa[sta],clu[stal],st[ockholm],vie[nna],msf,phy[lip],selex) (default auto)" << endl;
	cout << "  --dist-in               Distance matrix input (clustal omega format)" << endl << endl;
	cout << "Distance matrix options:" << endl;
	cout << "  -m, --metric            Distance metric (linear, affine, cumulative) (default linear" << endl;
	cout << "  -t, --threads           Number of threads used (0 for maximum available, default 1)" << endl;
	cout << "  --full                  Computes full distance matrix (default mBed)" << endl << endl;
	cout << "Analysis options:" << endl;
	cout << "  -a, --analysis          Analysis mode (boot[strap], iqr) (default bootstrap)" << endl;
	cout << "  --boot-rep              Number of bootstrap pseudoreplicates (default 1000)" << endl;
	cout << "  -s, --score             Sets the threshold for outliers in numbers of standard deviations (default 2)" << endl << endl;
	cout << "Output options:" << endl;
	cout << "  -r, --result            Result file containing the sequence names and score" << endl;
	cout << "  -o, --outlier           Detected outlier sequences in fasta format" << endl;
	cout << "  -c, --core              Detected core sequences in fasta format" << endl;
	cout << "  --dist-out              Writes calculated distance matrix to file" << endl;
	cout << "  -v, --verbose           Verbose output" << endl;
	cout << "  -h, --help              Prints this help and exit" << endl;
}

vector<double> tokenizeToDouble(const string str, const char& ch) {
    string next;
    vector<double> result;

    // For each character in the string
    for (string::const_iterator it = str.begin(); it != str.end(); it++) {
        // If we've hit the terminal character
        if (*it == ch) {
            // If we have some characters accumulated
            if (!next.empty()) {
                // Add them to the result vector
                result.push_back(atof(next.c_str()));
                next.clear();
            }
        } else {
            // Accumulate the next character into the sequence
            next += *it;
        }
    }
    if (!next.empty())
         result.push_back(atof(next.c_str()));
    return result;
}

vector<string> tokenizeToString(const string str, const char& ch) {
    string next;
    vector<string> result;

    // For each character in the string
    for (string::const_iterator it = str.begin(); it != str.end(); it++) {
        // If we've hit the terminal character
        if (*it == ch) {
            // If we have some characters accumulated
            if (!next.empty()) {
                // Add them to the result vector
                result.push_back(next.c_str());
                next.clear();
            }
        } else {
            // Accumulate the next character into the sequence
            next += *it;
        }
    }
    if (!next.empty())
         result.push_back(next.c_str());
    return result;
}

// trim from start
string ltrim(string s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
        return s;
}

// trim from end
string rtrim(string s) {
        s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
        return s;
}

// trim from both ends
string trim(string s) {
        return ltrim(rtrim(s));
}

string replaceAll(string str, string from, string to) {
	for(unsigned int i = 0; i < str.length(); i++){
		if(str.substr(i,1).compare(from) == 0){
			str.replace(i, from.length(), to);
			i += to.length()-1;
		}
	}
	return str;
}

int max( int a, int b, int c )
{
   int max = ( a < b ) ? b : a;
   return ( ( max < c ) ? c : max );
}

string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}


string random_string(int len){

	string s;
    static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < len; ++i) {
        s = s + alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s;
}
