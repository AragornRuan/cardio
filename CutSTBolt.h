#ifndef CUTSTBOLT_H
#define CUTSTBOLT_H 

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <armadillo>

#include "cardio.h"
#include "Storm.h"
#include "json/json.h"

using namespace storm;
using namespace std;
using namespace arma;

class CutSTBolt : public Bolt {
public:
	void Initialize(Json::Value conf, Json::Value context) { }
		
	void Process(Tuple &tuple);
};

void CutSTBolt::Process(Tuple &tuple) {
	std::string filename = tuple.GetValues()[0].asString();
	if (filename.empty()) {
		return;
	}

	u32 sfreq = tuple.GetValues()[1].asUInt();
	std::string xyzStr = tuple.GetValues()[2].asString();

	Log("Converting xyz from string to fmat.");
	fmat xyz(xyzStr);
	Log("Converted xyz from string to fmat.");

	urowvec J, K;
	bool T_inv_flag = false;
	Log("Executing cutST function.");
	u32 STT_count = cutST(xyz.col(0).t(), sfreq, J, K, T_inv_flag);
	Log("Executed cutST function.");

	Log("Executing merge_xstt function.");
	fmat x_stt = merge_xstt(xyz, J, K, STT_count);
	Log("Executed merge_xstt function.");

	Log("Converting x_stt, the return value of merge_xstt, from fmat to string.");
	stringstream s;
	x_stt.print(s, "");
	std::string x_sttStr = s.str();
	replace(x_sttStr.begin(), x_sttStr.end(), '\n', ';');
	Log("Converted x_stt, the return value of merge_xstt, from fmat to string.");

	Log("Generating tuple.");
	Json::Value result;
	result.append(filename);
	result.append(sfreq);
	result.append(x_sttStr);
	Tuple t(result);
	Log("Generated tuple.");

	Log("Emitting tuple.");
	Emit(t);
	Log("Emitted tuple");

}

#endif