#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <string>
#include <vector>
#include "CommonTypes.H"
#include <utility>
using namespace std;

class VariableManager;
class PotentialManager;
class MetaMove;
class SpeciesDistance;
class MappedOrthogroupReader;

class MetaLearner
{
	public:
		MetaLearner();
		~MetaLearner();

		int init();
		int setInputFName(const char*);
		int setMaxFactorSize(int);
		int setMaxFactorSize_Approx(int);
		int setConvergenceThreshold(double);

		// This list limits the regulators that will be considered for all species.
		int setRegulatorList(const char*);

		// This list limits the targets that will be considered for all species.
		int setTargetList(const char*);

		int setSpeciesDistances(SpeciesDistance*);
		int setBeta1(double);
		int setBeta2(double);
		void start();
		int setINDEP();

	private:
		void initSpeciesData(string speciesName, string tableFileName, string outputLoc, string motifNetwork);
		int addVariable(string varName);
		void showModelParameters();
		void process_mem_usage(double& vm_usage, double& resident_set);
		double getScore();
		double getPriorDelta();
		void precomputeEmptyGraphPrior();
		double precomputePerSpeciesPrior(int specID, Variable* target, SpeciesDataManager* sdm, unordered_map<int, unordered_map<int, double>> &edgePresenceProb);
		vector<MetaMove> collectMoves(int currK);
		vector<MetaMove> collectMoves_INDEP(int currK);
		double makeMoves(vector<MetaMove>& moves);
		void getNewPLLScore(int, Variable* u, Variable* v, double&, double&);
		double getPLLScore(int, SlimFactor* sFactor, int&);
		double getEdgePrior(int tfID, int targetID, SpeciesDataManager*);
		int dumpAllGraphs(int);

		unordered_map<string,int> speciesNameIDMap;
		vector<SpeciesDataManager*> speciesDataSet;
		int maxFactorSize;
		int maxFactorSizeApprox;
		string inputFName;
		double convThreshold;
		unordered_map<string,vector<int>*> affectedVarPairs;
		unordered_map<int, unordered_map<int, double>> edgePriors;
		SpeciesDistance* speciesData;
		vector<int> targetList;
		vector<int> regulatorList;
		vector<string> variableList;
		vector<unordered_map<int, double>> varNeighborhoodPrior_PerSpecies;
		vector<unordered_map<int, unordered_map<int, double>>> edgePresenceProb_PerSpecies;
		double beta1;
		double beta2;
		bool INDEP;
		double initGlobalScore;
		double edgePresenceProb;
};
#endif
