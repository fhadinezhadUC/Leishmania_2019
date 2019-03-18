/*
 * util.h
 *
 *  Created on: 20 Mar 2014
 *      Author: peter
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

bool isInteger(const std::string& s);
bool isFloat(string myString);
void display_help();
vector<double> tokenizeToDouble(const string str, const char& ch);
vector<string> tokenizeToString(const string str, const char& ch);
string ltrim(string s);
string rtrim(string s);
string trim(string s);
string replaceAll(string str, string from, string to);
int max(int, int, int);
string exec(const char*);
string random_string(int);

#endif /* UTIL_H_ */
