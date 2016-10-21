#include <armadillo>
#include <iostream>
#include "cardio.h"
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{

    if( argc != 2)
    {
        cout << "usage:" << argv[0] << " <filename> "<<endl;
        return 0;
    }

	string filepath = argv[1];
/*	fstream ecgFile(argv[1]);
	stringstream ecgS;
	string ecgData;
	string line1, line2;
	getline(ecgFile, line1);

	getline(ecgFile, line2);*/
//	fmat testVal(ecgS.str());
//	cout << testVal.n_cols << endl;
//	testVal.print(cout, "");
/*	while (getline(ecgFile, line)) {
		ecgData += line + "\n";
	}*/

/*	cout << line1 << flush << ";" << line2 << endl;
	
	cout << val.n_cols << endl;*/
	fmat val;
	bool status = val.load(filepath, raw_ascii);
	if ( !status )
    {
        cout << "loading error!" << endl;
        exit(1);
    }
	u32 sfreq = 0, gain = 0,repeat = 0;
	float TS = 0.0;
	pretreat(val, sfreq, repeat, TS, gain);
	//val.save("result/preated.txt", raw_ascii);
	cout << "Pretreat finished..." << endl;
	fmat xyz = val.cols(12, 14);
	fstream f("result/test.txt", ios::out);
	stringstream s;
//	s << sfreq;
//	u32 sTmp;
//	s >> sTmp;
//	cout << sTmp << endl;
//	cout << s.str() << endl;
//	s.clear();
	xyz.print(f, "");
	xyz.print(s, "");
	string tmpXYZ = s.str();
	replace(tmpXYZ.begin(), tmpXYZ.end(), '\n', ';');
//	cout << tmpXYZ << endl;
	fmat tmp(tmpXYZ);
	cout << tmp.n_cols << endl;
	tmp.save("result/Pretreat.txt", raw_ascii);
	urowvec J, K;
	bool T_inv_flag = false;
	u32 STT_count = cutST(xyz.col(0).t(), sfreq, J, K, T_inv_flag);
	//J.save("result/J.txt", raw_ascii);
	//K.save("result/K.txt", raw_ascii);
	fmat x_stt = merge_xstt(xyz, J, K,STT_count);
	//cout << "STT_count" << STT_count << endl;
	//x_stt.save("result/x_stt.txt", raw_ascii);
	cout << "learning..." << endl;
	fmat WS = learn(sfreq, x_stt);
	cout << "learn finished,saving result..." << endl;
	WS.save("result/WS.txt", raw_ascii);
    cout << "WS.txt saved...!" << endl;
	return 0;

}