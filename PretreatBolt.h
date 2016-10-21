#ifndef PRETREATBOLT
#define PRETREATBOLT

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

class PretreatBolt : public Bolt
{
	public:
		void Initialize(Json::Value conf, Json::Value context) { }
		
		void Process(Tuple &tuple);
};

void PretreatBolt::Process(Tuple &tuple) {
	std::string filename = tuple.GetValues()[0].asString();
	if (filename.empty()) {
		return;
	}
	std::string data = tuple.GetValues()[1].asString();

	Log("Converting ECG string to fmat.");
	fmat val(data);
	Log("Converted ECG string to fmat.");

	u32 sfreq = 0, gain = 0,repeat = 0;
	float TS = 0.0;
	Log("Pretreating ECG fmat.");
	pretreat(val, sfreq, repeat, TS, gain);
	Log("Pretreated ECG fmat.");
	Log("Saving xyz from pretreated result.");
	fmat xyz = val.cols(12, 14);
	Log("Saved xyz from pretreated result.");

	Log("Converting xyz to string.");
	stringstream s;
	xyz.print(s, "");
	std::string xyzStr = s.str();
	replace(xyzStr.begin(), xyzStr.end(), '\n', ';');
	Log("Converted xyz to string.");

	Log("Generating tuple.");
	Json::Value result;
	result.append(filename);
	result.append(sfreq);
	result.append(xyzStr);

	Tuple t(result);
	Log("Generated tuple");

	Log("Emitting tuple");
	Emit(t);
	Log("Emitted tuple");
}

#endif
