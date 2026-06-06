#ifndef _SPECIES_DATA_MANAGER_
#define _SPECIES_DATA_MANAGER_

#include <map>
#include <string>
using namespace std;

class VariableManager;
class PotentialManager;

class SpeciesDataManager
{
	public:
		~SpeciesDataManager();

		int setVariableManager(VariableManager*);
		int setPotentialManager(PotentialManager*);
		int setOutputLoc(const char*);
		int setMotifNetwork(const char*);

		VariableManager* getVariableManager();
		PotentialManager* getPotentialManager();
		const char* getOutputLoc();
		unordered_map<int,unordered_map<int,double>*>& getMotifNetwork();

		void createFactors();
		SlimFactor* getFactor(int varID);

	private:
		VariableManager* varMgr;
		PotentialManager* potMgr;
		char outputLoc[1024];
		unordered_map<int,unordered_map<int,double>*> motifNetwork;
		unordered_map<int, SlimFactor*> factors;
};
#endif
