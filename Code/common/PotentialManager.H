/* The goal of this class is to take all the MarkovBlankets and the factors and
 *  the data and create all potentials involving valid pairings MarkovBlankets 
 *  and factors. We create all putative potentials upfront because we want to
 *  iterate over the data once and update counts of configurations in one go*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include "Matrix.H"
#include "gsl/gsl_randist.h"

using namespace std;

class Potential;
class SlimFactor;
class Error;
class VariableManager;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();

        int deleteData();
		int setOutputDir(const char*);
		int init();
		int reset();
        int populatePotential(Potential* aPot);
        Error::ErrorCode loadEvidenceFromTable(vector<string>& inputLines);
        Error::ErrorCode loadEvidenceFromTable(vector<string>& inputTable,unordered_set<string>& inputVariableNames);
        VariableManager* getVariableManager();
        double computeMeanVarPseudoLikelihood_onefold(int);
        double computePotentialMBCovMean(SlimFactor* sFactor,int&);
        int computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondvar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect);
        int dumpVarMB_PairwiseFormat(SlimFactor* sFactor, ofstream& oFile, VSET& varSet);

	private:
		int estimateCovariance_Eff(int, int);
		double computeLL(int dim, double determinant);

		char outputDir[1024];
		gsl_matrix* ludecomp;
		gsl_permutation* perm;
        Matrix* data;
        Matrix* meanMat;
        Matrix* covMat;
        double testdataSize;
        VariableManager* vMgr;
};
#endif
