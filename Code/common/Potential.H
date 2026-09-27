#ifndef _POTENTIAL
#define _POTENTIAL
//A potential over continuous variables needs only store the mean and variance parameters.
//However to compute the pdf, we need to store the inverse of the covariance matrix and the
//determinant.
#include <fstream>
#include "cmath"
#define PI M_PI //3.1472 3.1415926
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#include "Matrix.H"
#include <unordered_map>

class Variable;

class Potential
{
	public:
		Potential();
		~Potential();

		// This enum is used to decide what role the Variable
		// plays in the potential. A Variable can be either from a Markov Blanket
		// or from the factor. This is required while converting the joint into a conditional
		// of the factor variables give Markov blanket variables
		
		typedef enum
		{
			FACTOR,
			MARKOV_BNKT
		} VariableRole;
		
		int setAssocVariable(Variable*,VariableRole);
		unordered_map<int,Variable*>& getAssocVariables();

		int potZeroInit();
		
		int copyMe(Potential**);

		int updateMean(int,double);
		int updateCovariance(int,int,double);
		int makeValidJPD(gsl_matrix*,gsl_permutation*);

		double generateSample(INTDBLMAP&, int vId,gsl_rng*);

		int initMBCovMean();
		double getCondVariance();
		double getCondBias();
		unordered_map<int,double>& getCondWeight();

	protected:
		unordered_map<int,Variable*> varSet;
		INTINTMAP factorVariables;
		INTINTMAP markovBlnktVariables;
		INTDBLMAP meanWrap;
		Matrix* covariance;
		Matrix* inverse;
		Matrix* mean;
		double determinant;
		double normFactor;

		//Each variable is mapped to an index in the covariance or mean matrix
		//This map stores that information
		int matrixInd;
		INTINTMAP vIDMatIndMap;
		INTINTMAP matIndvIDMap;

		double mbcondVar;
		double mbcondMean_Part;
		unordered_map<int,double> mbcondMean_Vect;
};
#endif
