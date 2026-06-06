/* The goal of this class is to take all the MarkovBlankets and the factors and
 *  the data and create all potentials involving valid pairings MarkovBlankets 
 *  and factors. We create all putative potentials upfront because we want to
 *  iterate over the data once and update counts of configurations in one go*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include <vector>
#include <unordered_map>
#include "Matrix.H"

using namespace std;

class SlimFactor;
class VariableManager;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();

        void deleteData();
		int setOutputDir(const char*);
        void loadEvidenceFromTable(const vector<string>& inputTable, const vector<string>& variableList);
        VariableManager* getVariableManager();
        double computeMeanVarPseudoLikelihood(int);
        double computePotentialMBCovMean(SlimFactor* sFactor, int&);
        void computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondvar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect);
        void dumpVarMB(SlimFactor* sFactor, ofstream& oFile);

	private:
		int estimateCovariance(int, int);
		double computeLL(int dim, double determinant);

		char outputDir[1024];
        Matrix* data;
        Matrix* meanMat;
        Matrix* covMat;
        double testdataSize;
        VariableManager* vMgr;
};
#endif
