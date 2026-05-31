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
		int setRestrictedList(const char*);
		int setTrueGraph(const char* aGraphName);
		int setInputOGList(const char*);

		int setOrthogroupReader(MappedOrthogroupReader*);
		int setSpeciesDistances(SpeciesDistance*);
		int setBeta1(double);
		int setBeta2(double);
		int doOneFold();
		int start(int);
		int showModelParameters(int);
		int setINDEP();
		int setsplitGenes();

	private:		
		void process_mem_usage(double& vm_usage, double& resident_set);
		double getInitScore();
		double getScore();
		double getPriorDelta();
		double precomputeEmptyGraphPrior();
		double precomputePerSpeciesPrior(int specID, int targetID, Variable* target, SpeciesDataManager* sdm, vector<vector<double>> &edgePresenceProb,map<int,MappedOrthogroup*>& orthogroupSet);
		int initCondsetMap_Nopool();
		INTINTMAP* getConditionSet(int cind);
		int setinputVariableNames();
		int collectMoves_Orthogroups(int);
		int collectMoves_Orthogroups_INDEP(int);
		double makeMoves(int& successMove);

		int attemptMove(MetaMove*);
		int getNewPLLScore(int,Variable* u, Variable* v,double&,double&,int);
		double getPLLScore_Condition_Tracetrick(int,SlimFactor* sFactor,int&);
		double getEdgePrior_PerSpecies(int tfID, int targetID,SpeciesDataManager*);
		int genCondSetKey(INTINTMAP& condSet, string& aKey);
		int dumpAllGraphs(int,int);

		int generateData(int sampleCnt, int burnin,int foldid);
		int generateInitSample(INTDBLMAP& initialSample,gsl_rng* rndgen,map<int,Potential*>*);
		int writeEvidenceTo(ofstream* oFile,INTDBLMAP& sample);
		int initCondsetMap_Tree(SpeciesDistance::Species* node);
		vector<unordered_map<int,int> > initCondsetMap_Tree_backtrack(SpeciesDistance::Species* node, vector<unordered_map<int,int> > mycondition, int parentstatus, int ntransition);
		int getLeaves(SpeciesDistance::Species* node, map<string,int>& leaves);

		void readEvidenceTable(string fileName);

		vector<string> speciesIDNameMap;
		unordered_map<string,int> speciesNameIDMap;
		vector<SpeciesDataManager*> speciesDataSet;
		unordered_set<int> intermediatecells;
		map<string,PotentialManager*> pooledPotentials;
		int maxFactorSize;
		int maxFactorSizeApprox;
		string inputFName;
		double convThreshold;
		char restrictedFName[1024];

		unordered_map<string,vector<int>*> affectedOGPairs;
		vector<vector<double> > ogpairPrior;
		FactorGraph* fg;
		vector<MetaMove*> moveSet;
		map<int,INTINTMAP*> condsetMap;
		vector<vector<int>> condsetMap_Tree;
		map<int,INTDBLMAP*> validationPLLs;
		char trueGraphFName[1024];
		map<int,double> finalScores;
		map<string,int> subgraphVarSet;
		map<string,int> bottomClip;
		map<string,int> topClip;
		double globalVar;
		MappedOrthogroupReader* ogr;
		SpeciesDistance* speciesData;
		char orthoMapFName[1024];
		char speciesOrder[1024];
		vector<int> inputOGList;
		vector<int>inputRegulatorOGs;
		unordered_set<int> inputVariables;
		unordered_map<int,unordered_set<string>> inputVariableNames;
		vector<vector<double>> varNeighborhoodPrior_PerSpecies;
		vector<vector<vector<double>>> edgePresenceProb_PerSpecies;
		double beta1;
		double beta2;
		vector<string> inputTable;
		bool INDEP;
		bool splitGenes;
		double initGlobalScore;
		double edgePresenceProb;
};
#endif
