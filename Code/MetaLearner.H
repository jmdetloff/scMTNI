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
		int collectMoves(int);
		int collectMoves_INDEP(int);
		double makeMoves(int& successMove);
		void attemptMove(MetaMove*);
		void getNewPLLScore(int, Variable* u, Variable* v, double&, double&);
		double getPLLScore(int, SlimFactor* sFactor, int&);
		double getEdgePrior(int tfID, int targetID, SpeciesDataManager*);
		int dumpAllGraphs(int);
		int initCondsetMap_Tree(SpeciesDistance::Species* node);
		vector<unordered_map<int,int> > initCondsetMap_Tree_backtrack(SpeciesDistance::Species* node, vector<unordered_map<int,int> > mycondition, int parentstatus, int ntransition);
		int getLeaves(SpeciesDistance::Species* node, map<string,int>& leaves);
		void readEvidenceTable(string fileName);

		vector<string> speciesIDNameMap;
		unordered_map<string,int> speciesNameIDMap;
		vector<SpeciesDataManager*> speciesDataSet;
		unordered_set<int> intermediatecells;
		int maxFactorSize;
		int maxFactorSizeApprox;
		string inputFName;
		double convThreshold;
		unordered_map<string,vector<int>*> affectedVarPairs;
		unordered_map<int, unordered_map<int, double>> edgePriors;
		vector<MetaMove*> moveSet;
		vector<vector<int>> condsetMap_Tree;
		SpeciesDistance* speciesData;
		vector<int> targetList;
		vector<int> regulatorList;
		vector<string> variableList;
		vector<unordered_map<int, double>> varNeighborhoodPrior_PerSpecies;
		vector<unordered_map<int, unordered_map<int, double>>> edgePresenceProb_PerSpecies;
		double beta1;
		double beta2;
		vector<string> inputTable;
		bool INDEP;
		double initGlobalScore;
		double edgePresenceProb;
};
#endif
