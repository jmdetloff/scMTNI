#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "Matrix.H"
#include "gsl/gsl_randist.h"
using namespace std;

class Error;
class Evidence;

//EMAP stores the map of variable ids and their evidence. This corresponds
//to one line in the datafile
typedef map<int,Evidence*> EMAP;
typedef map<int,Evidence*>::iterator EMAP_ITER;

typedef vector<EMAP*> EVIDENCE_SET;

class VariableManager;
class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();

		// This associates the variable information with the variable id
		int setVariableManager(VariableManager*);

		Error::ErrorCode loadEvidenceFromTable(vector<string>& inputLines);

		//Get the evidence set at the specified id
		EMAP* getEvidenceAt(int);

		int setFoldCnt(int);
		int setPreRandomizeSplit();
		int setPreRandomizeSplitSeed(int);
		int splitData(int);
		vector<int>& getTrainingSet();
        vector<int>& getTestSet();
		INTINTMAP& getValidationSet();
        Matrix* getData();
        VariableManager* getVariableManager();
    
	private:
		int populateRandIntegers(gsl_rng*, int*,int);

		EVIDENCE_SET evidenceSet;
		VariableManager* vMgr;
        vector<int> testIndex;
        vector<int> trainIndex;
		INTINTMAP validationIndex;
		int foldCnt;
		INTDBLMAP maxAbsValues;
		bool preRandomizeSplit;
		int randseed;
        Matrix* dataMat;
};
#endif
