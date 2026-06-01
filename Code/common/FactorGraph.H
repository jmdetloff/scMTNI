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

		int getFactorCnt();

		void addFactor(int varID);
		SlimFactor* getFactor(int varID);

	private:
		unordered_map<int, SlimFactor*> factors;
};
#endif
