#ifndef _SPECIES_DATA_MANAGER_
#define _SPECIES_DATA_MANAGER_

#include <map>
#include <string>
using namespace std;

class VariableManager;
class FactorGraph;
class PotentialManager;

class SpeciesDataManager
{
	public:
		int setVariableManager(VariableManager*);
		int createFactorGraph();
		int setPotentialManager(PotentialManager*);
		int setOutputLoc(const char*);
		int setMotifNetwork(const char*);
		VariableManager* getVariableManager();
		FactorGraph* getFactorGraph();
		PotentialManager* getPotentialManager();
		const char* getOutputLoc();
		unordered_map<int,unordered_map<int,double>*>& getMotifNetwork();
	private:
		VariableManager* varMgr;
		FactorGraph* fgraph;
		PotentialManager* potMgr;
		char outputLoc[1024];
		unordered_map<int,unordered_map<int,double>*> motifNetwork;
};
#endif
