#ifndef CARDIO_H
#define CARDIO_H 

#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

void diff(const urowvec&, irowvec&, int);

u32 cutST(const frowvec&, u32, urowvec&, urowvec&, bool&);

fmat merge_xstt(const fmat&, const urowvec&, const urowvec&, const u32&);

void _midf(fvec, fvec&, u32);

void midfilter1(fmat&,u32);

void ddencmp(string, string, vector<double>, double&, char&, bool&);

vector<double>  wdencmp(string, vector<double>&, vector<int>, string, int, vector<double>&, double&, char&, bool&);

vector<double> waveletfilter(vector<double>, int, string);

void pretreat(fmat&, u32&, u32&, float&, u32&);

fmat learn(int, fmat&);

#endif