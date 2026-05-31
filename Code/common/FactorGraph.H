#ifndef _FACTOR_GRAPH
#define _FACTOR_GRAPH

#include <map>
#include <fstream>
using namespace std;

class SlimFactor;

class FactorGraph
{
	public:
		FactorGraph();
		~FactorGraph();

		int setFactor(SlimFactor*);
		int getFactorCnt();
		SlimFactor* getFactorAt(int);
		vector<SlimFactor*>& getAllFactors();

	private:
		vector<SlimFactor*> factorSet;
};
#endif
