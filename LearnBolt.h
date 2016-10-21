#ifndef LEARNBOLT_H
#define LEARNBOLT_H 

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

class LearnBolt : public Bolt
{
	public:
		void Initialize(Json::Value conf, Json::Value context) { }
		
		void Process(Tuple &tuple);
};

void LearnBolt::Process(Tuple &tuple) {
	std::string filename = tuple.GetValues()[0].asString();
	if (filename.empty()) {
		return;
	}

	u32 sfreq = tuple.GetValues()[1].asUInt();
	std::string x_sttStr = tuple.GetValues()[2].asString();

	Log("Converting x_stt from string to fmat.");
	fmat x_stt(x_sttStr);
	Log("Converted x_stt from string to fmat.");

	Log("Executing learn function.");
	fmat WS = learn(sfreq, x_stt);
	Log("Executed learn function.");

	Log("Converting WS, the return value of learn function, from fmat to string.");
	stringstream s;
	WS.print(s, "");
	std::string WSStr = s.str();
	Log("Converted WS, the return value of learn function, from fmat to string.");

	Log("Generating tuple.");
	Json::Value result;
	result.append(filename);
	result.append(WSStr);
	Tuple t(result);
	Log("Generated tuple.");

	Log("Emitting tuple.");
	Emit(t);
	Log("Emitted tuple.");
}

#endif