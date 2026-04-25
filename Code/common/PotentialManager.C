/*
 *  scMTNI: single-cell Multi-Task learning Network Inference
 *   Copyright 2022 Shilu Zhang (szhang256@wisc.edu) and  Sushmita Roy (sroy@biostat.wisc.edu)
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, 
 *   including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 *   subject to the following conditions:
 *
 *   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 *   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *   */

#include <iostream>
#include <cstring>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"
#include <chrono>
#include "Utils.H"
#include "VariableManager.H"


using namespace std::chrono;

PotentialManager::PotentialManager()
{
	ludecomp=NULL;
	perm=NULL;
    data=NULL;
    meanMat=NULL;
    covMat=NULL;
    testdataSize=0;
}

PotentialManager::~PotentialManager()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
    if (data!=NULL)
    {
        delete data;
    }
		if (meanMat!=NULL)
    {
        delete meanMat;
    }
    if (covMat!=NULL)
    {
        delete covMat;
    }
}

int
PotentialManager::deleteData()
{
    if (data!=NULL)
    {
        delete data;
    }
    return 0;
}

int 
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
PotentialManager::setOutputDir(const char* aDirName)
{
	strcpy(outputDir,aDirName);
	return 0;
}

//added from EvidenceManager
Error::ErrorCode
PotentialManager::loadEvidenceFromTable(vector<string>& inputTable)
{
    /* Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
       to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.
       This method does what VariableManager::readVariables() and EvidenceManager::loadEvidenceFromFile_Continuous()
       do for .model and .data files.
    */

    /* For each sample, an EMAP (Evidence Map) is created that has an Evidence object per variable (gene) indexed by
       variable ID. These EMAPs are saved in evidenceSet, a vector of EMAP*.
    */
    
    // First, let's collect the values per sample. The table has them per gene...
    //vector<vector<double>> valuesPerSample;
    int nodeCount = inputTable.size();
    // How many samples do we have?
    vector<string> substrs = Utils::split(inputTable[0], '\t');
    //cout <<"Evidence " <<substrs[0] << endl;
    int sampleCount = substrs.size() - 1;  // not counting the gene name
    data=new Matrix(nodeCount,sampleCount);
    int nodeNum = 0;
    vMgr = new VariableManager;
    for (auto line : inputTable)
    {
        substrs = Utils::split(line, '\t');
        // set variable name:
        string nodeName = substrs[0];
        vMgr->readVariablesFromData(nodeName,nodeNum);
        for (int sample = 0; sample < sampleCount; sample++)
        {
            double varVal =stod(substrs[sample + 1]);
            //valuesPerSample[sample][nodeNum] =varVal;
            data->setValue(varVal,nodeNum,sample);
        }
        ++nodeNum;
    }
    testdataSize=sampleCount;
    cout <<"PotentialManager::loadEvidenceFromTable varCount=" << data->getRowCnt() << " sampleCount=" << data->getColCnt()  << endl;
    substrs.clear();
    return Error::SUCCESS;
}

Error::ErrorCode
PotentialManager::loadEvidenceFromTable(vector<string>& inputTable,unordered_set<string>& inputVariableNames)
{
    /* Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
       to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.
       This method does what VariableManager::readVariables() and EvidenceManager::loadEvidenceFromFile_Continuous()
       do for .model and .data files.
    */

    /* For each sample, an EMAP (Evidence Map) is created that has an Evidence object per variable (gene) indexed by
       variable ID. These EMAPs are saved in evidenceSet, a vector of EMAP*.
    */
    
    // First, let's collect the values per sample. The table has them per gene...
    //vector<vector<double>> valuesPerSample;
    int nodeCount = inputVariableNames.size();
    // How many samples do we have?
    vector<string> substrs = Utils::split(inputTable[0], '\t');
    //cout <<"Evidence " <<substrs[0] << endl;
    int sampleCount = substrs.size() - 1;  // not counting the gene name
    data=new Matrix(nodeCount,sampleCount);
    int nodeNum = 0;
    vMgr = new VariableManager;
    for (auto line : inputTable)
    {
        substrs = Utils::split(line, '\t');
        // set variable name:
        string nodeName = substrs[0];
        if(inputVariableNames.find(nodeName)!=inputVariableNames.end())
        {
            vMgr->readVariablesFromData(nodeName,nodeNum);
            for (int sample = 0; sample < sampleCount; sample++)
            {
                double varVal =stod(substrs[sample + 1]);
                //valuesPerSample[sample][nodeNum] =varVal;
                data->setValue(varVal,nodeNum,sample);
            }
            ++nodeNum;
        }/*else{
            cout << nodeName << " is not in GeneList" << endl;
        }*/
        
    }
    testdataSize=sampleCount;
    cout <<"PotentialManager::loadEvidenceFromTable varCount=" << data->getRowCnt() << " sampleCount=" << data->getColCnt()  << endl;
    substrs.clear();
    return Error::SUCCESS;
}

VariableManager*
PotentialManager::getVariableManager()
{
    return vMgr;
}

int
PotentialManager::init(int f)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_%d.txt",outputDir,f);
	sprintf(sdFName,"%s/gauss_std_%d.txt",outputDir,f);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		vector<int>& trainEvidSet=evMgr->getTestSet();
        estimateAllMeanCov(trainEvidSet);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::init()  // one fold only!
{
    int varCnt=data->getRowCnt();//number of variables
    meanMat=new Matrix(varCnt,1);
    meanMat->setAllValues(0);
    covMat=new Matrix(varCnt,varCnt);
    covMat->setAllValues(-1);
    return 0;
}

int
PotentialManager::reset()
{
    if (meanMat!=NULL)
    {
        delete meanMat;
    }
    if (covMat!=NULL)
    {
        delete covMat;
    }
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
	return 0;
}

int
PotentialManager::estimateAllMeanCov(vector<int>& trainEvidSet)
{
    int evidCnt=trainEvidSet.size();

	// First get the mean and then the variance
    // data is the data matrix which will have the variable by sample information
    EMAP* evidMap=evMgr->getEvidenceAt(trainEvidSet[0]);
    int varCnt=evidMap->size();
    if(data==NULL)
    {
        data=new Matrix(varCnt,trainEvidSet.size());
        meanMat=new Matrix(varCnt,1);
        meanMat->setAllValues(0);
        covMat=new Matrix(varCnt,varCnt);
        covMat->setAllValues(-1);
    }

    Matrix *data=evMgr->getData();
    cout << " evidCnt=" << evidCnt << " varCnt=" << varCnt ;
    for(int eIter=0;eIter<trainEvidSet.size();eIter++)
    {
        EMAP* evidMap = evMgr->getEvidenceAt(trainEvidSet[eIter]);
        for(int vId=0;vId<evidMap->size();vId++)
        {
            Evidence* evid=(*evidMap)[vId];
            double val=evid->getEvidVal();
            data->setValue(val,vId,trainEvidSet[eIter]);
        }
    }

    //Done copying. Now we can go over the rows of data and get the means
    for(int i=0;i<varCnt;i++)
    {
        double s=0;
        for(int j=0;j<data->getColCnt();j++)
        {
            double val=data->getValue(i,j);
            s=s+val;
        }
        double sampleSize=(double) data->getColCnt();
        meanMat->setValue(s/sampleSize,i,0);
        cout << "i=" << i << " mean=" <<s/sampleSize << endl;
    }

    return 0;
}

int
PotentialManager::estimateCovariance_Eff(int uId, int vId)
{
    double vmean=meanMat->getValue(vId,0);
    double umean=meanMat->getValue(uId,0);
    double ssd=data->vectorMultiply(vId,vmean,uId,umean);
    //Now estimate the variance
    double var=ssd/((double)(data->getColCnt()-1));
    covMat->setValue(var,uId,vId);
    covMat->setValue(var,vId,uId);
    return 0;
}

int
PotentialManager::readAllMeanCov(const char* mFName, const char* sdFName)
{

	ifstream mFile(mFName);
	ifstream sdFile(sdFName);
	char buffer[1024];
	while(mFile.good())
	{
		mFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId;
		double mean=0;
		while(tok!=NULL)
		{	
			if(tokCnt==0)
			{
				vId=atoi(tok);	
			}
			else if(tokCnt==1)
			{
				mean=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		//globalMean[vId]=mean;
        meanMat->setValue(mean,vId,0);
	}
	mFile.close();
	int lineNo=0;
	while(sdFile.good())
	{
		sdFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId=0;
		int uId=0;
		double covariance=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				uId=atoi(tok);
			}
			else if(tokCnt==1)
			{
				vId=atoi(tok);
			}
			else if(tokCnt==2)
			{
				covariance=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
        covMat->setValue(covariance,uId,vId);
        covMat->setValue(covariance,vId,uId);
		if(uId==0 && vId==0)
		{
			cout << "Found uId=0 vId=0 covar="<< covariance << " at lineno " << lineNo  << endl;
		}
		lineNo++;
	}
	sdFile.close();
	return 0;
}

//precompute covMat and meanMat:
// compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
double
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, int& status)//unordered_map<int,double>& mbcondMean_Vect)
{
    int mbsize=sFactor->mergedMB.size();
    if(mbsize==0)
    {
        return 0;
    }
    vector <int> ParChildID(sFactor->mergedMB.begin(),sFactor->mergedMB.end());
    ParChildID.push_back(sFactor->fId); //the last one is target
    //extract covariance
    Matrix* covariance=new Matrix(ParChildID.size(),ParChildID.size());
    for(int i=0;i<ParChildID.size();i++) {
        int vIter=ParChildID[i];  //index in covMat
        for(int j=i;j<ParChildID.size();j++)
        {
            int uIter=ParChildID[j];  //index in covMat
            double cval=covMat->getValue(vIter,uIter);
            covariance->setValue(cval,i,j);
            covariance->setValue(cval,j,i);
        }
    }
    int vId=sFactor->fId;
    double determinant=covariance->detMatrix();
    if(determinant<=0) //<1e-50)
    {
        //cout <<"Negative/Zero Determinant " << determinant << " for Target=" << sFactor->fId <<endl;
        status=-1;
        return 0;
    }
    double mbcondVar=covMat->getValue(vId,vId);
    
    Matrix* mbcov=covariance->getSubMatrix(0,0,mbsize,mbsize);
    Matrix* mbmargvar=covariance->getSubMatrix(ParChildID.size()-1,0,1,mbsize);
    double determinantP=mbcov->detMatrix();
    if(determinantP<=0) //<1e-50
    {
        //cout <<"Negative/Zero Determinant Parent " << determinantP << " for Target=" << sFactor->fId <<endl;
        status=-1;
        return 0;
    }

    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);

    for(int i=0;i<ParChildID.size()-1;i++)  //start from parent only, last one is target
    {
        int varID=ParChildID[i];
        double aVal=prod1->getValue(0,i);
        double bVal=mbmargvar->getValue(0,i);
        mbcondVar=mbcondVar-(aVal*bVal);
    }
    if(mbcondVar<0){
        status=-1;
        return 0;
    }
    double jointll1 = computeLL(ParChildID.size(), determinant);
    double jointll2 = computeLL(mbsize, determinantP);
    double pll = jointll1 - jointll2;
    if(isinf(pll)){
        status=-1;
        return 0;
    }
    
    ParChildID.clear();
    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    return pll;
}

int
PotentialManager::dumpVarMB_PairwiseFormat(SlimFactor* sFactor, ofstream& oFile, VSET& varSet)
{
    int mbsize=sFactor->mergedMB.size();
    if(mbsize==0)
    {
        return 0;
    }
    vector <int> ParChildID(sFactor->mergedMB.begin(),sFactor->mergedMB.end());
    ParChildID.push_back(sFactor->fId); //the last one is target
    //extract covariance
    Matrix* covariance=new Matrix(ParChildID.size(),ParChildID.size());
    for(int i=0;i<ParChildID.size();i++) {
        int vIter=ParChildID[i];  //index in covMat
        for(int j=i;j<ParChildID.size();j++)
        {
            int uIter=ParChildID[j];  //index in covMat
            double cval=covMat->getValue(vIter,uIter);
            covariance->setValue(cval,i,j);
            covariance->setValue(cval,j,i);
        }
    }

    Matrix* mbcov=covariance->getSubMatrix(0,0,mbsize,mbsize); //Matrix* mbcov=new Matrix(mbsize,mbsize);
    Matrix* mbmargvar=covariance->getSubMatrix(ParChildID.size()-1,0,1,mbsize); //Matrix* mbmargvar=new Matrix(1,mbsize); covariance(vId,j)

    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);

    for(int i=0;i<ParChildID.size()-1;i++)  //start from parent only, last one id=ParChildID.size()-1 is target
    {
        int varID=ParChildID[i];
        double aVal=prod1->getValue(0,i); // aVal is conditional weight
        oFile << varSet[varID]->getName()<< "\t"<< varSet[sFactor->fId]->getName() << "\t" << aVal << endl;
    }

    ParChildID.clear();
    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    return 0;
}

// compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
// for MetaLearner::showModelParameters() final parameters
int
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondVar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect)
{
    int mbsize=sFactor->mergedMB.size();
    if(mbsize==0)
    {
        return 0;
    }
    vector <int> ParChildID(sFactor->mergedMB.begin(),sFactor->mergedMB.end());
    ParChildID.push_back(sFactor->fId); //the last one is target
    //extract covariance
    Matrix* covariance=new Matrix(ParChildID.size(),ParChildID.size());
    for(int i=0;i<ParChildID.size();i++) {
        int vIter=ParChildID[i];  //index in covMat
        for(int j=i;j<ParChildID.size();j++)
        {
            int uIter=ParChildID[j];  //index in covMat
            double cval=covMat->getValue(vIter,uIter);
            covariance->setValue(cval,i,j);
            covariance->setValue(cval,j,i);
        }
    }
    int vId=sFactor->fId;
    mbcondVar=covMat->getValue(vId,vId);
    mbbias=meanMat->getValue(vId,0);
    
    Matrix* mbcov=covariance->getSubMatrix(0,0,mbsize,mbsize);
    Matrix* mbmargvar=covariance->getSubMatrix(ParChildID.size()-1,0,1,mbsize);
    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);

    for(int i=0;i<ParChildID.size()-1;i++)  //start from parent only, last one is target
    {
        int varID=ParChildID[i];
        double aVal=prod1->getValue(0,i);
        double bVal=mbmargvar->getValue(0,i);
        mbcondVar=mbcondVar-(aVal*bVal);
        mbcondMean_Vect[varID]=aVal;  // mbcondMean_Vect is conditional weight
        double cVal= meanMat->getValue(varID,0); //meanWrap[aIter->second];
        mbbias=mbbias-(cVal*aVal);
    }
    ParChildID.clear();
    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    return 0;
}

double
PotentialManager::computeLL(int dim, double determinant)
{
    double ll=testdataSize * (dim*log(2*PI)+log(determinant));
    double t=dim*(testdataSize-1);
    ll=(ll+t)*(-0.5);
    return ll;
}

int
PotentialManager::populatePotential(Potential* aPot)
{
    unordered_map<int,Variable*>& potVars=aPot->getAssocVariables();
    for(auto vIter=potVars.begin();vIter!=potVars.end(); vIter++)
    {
        double mean=meanMat->getValue(vIter->first,0);
        aPot->updateMean(vIter->first,mean);
        for(auto uIter=vIter;uIter!=potVars.end();uIter++)
        {
            double cval=covMat->getValue(vIter->first,uIter->first);
            if(cval==-1)
            {
                estimateCovariance_Eff(vIter->first,uIter->first);
                //cerr <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
                //exit(-1);
            }
            cval=covMat->getValue(vIter->first,uIter->first);
            aPot->updateCovariance(vIter->first,uIter->first,cval);
            aPot->updateCovariance(uIter->first,vIter->first,cval);
        }
    }
    aPot->makeValidJPD(ludecomp, perm);
    return 0;
}

double
PotentialManager::computeMeanVarPseudoLikelihood_onefold(int id) //SlimFactor* sFactor,VSET& varSet)
{
    // compute meanMat and Variance on covMat, PseudoLikelihood
    double vmean=data->RowMean(id);
    meanMat->setValue(vmean,id,0);
    double dataSetSize=data->getColCnt();
    double ssd=data->vectorMultiply(id,vmean,id,vmean);
    //double variance=covMat->getValue(id,id);
    //add 1e-10 to avoid singularity issues:
    double variance=ssd/(dataSetSize-1.0)+1e-10;  //(0.001+ssd)/((double)(data->getColCnt()-1))
    covMat->setValue(variance,id,id);
    double pll=-0.5*ssd/variance-0.5*log(2.0*PI*variance)*dataSetSize;
    //cout << "PotentialManager::computeMeanVarPseudoLikelihood_onefold id="<<id<<" variance=" <<variance << " mean=" << vmean << " ssd=" << ssd << " dataSetSize=" << dataSetSize <<" pll=" << pll << endl;
    for(int j=id+1;j<data->getRowCnt();j++)
    {
        estimateCovariance_Eff(id,j);
    }
    return pll;
}

double
PotentialManager::getPseudoLikelihood(SlimFactor* sFactor,vector<Variable*>& varSet, bool train)
{
	Potential* aPotFunc=new Potential;
	Variable* aVar=varSet[sFactor->fId];
	aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
	for(auto aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++)  //INTINTMAP_ITER
	{
		Variable* aVar=varSet[*aIter];
		aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	aPotFunc->potZeroInit();
    populatePotential(aPotFunc);
	//This function creates a submatrix of the covariance matrix and inverts it
	aPotFunc->initMBCovMean();
	vector<int>* dataSet=NULL;
	if(train)
	{
		dataSet=&(evMgr->getTrainingSet());
	}
	else
	{
		dataSet=&(evMgr->getTestSet());
	}
	INTDBLMAP subData;
	double pll=0;
	int thresholded=0;
    for(int dIter=0;dIter<dataSet->size();dIter++) //for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt((*dataSet)[dIter]); //dIter->first
		Evidence* evid=(*evidMap)[sFactor->fId];
		double val=evid->getEvidVal();
		subData[sFactor->fId]=val;
		for(auto vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end(); vIter++) //INTINTMAP_ITER
		{
			int vId=*vIter;
			Evidence* evid=(*evidMap)[vId];
			double val=evid->getEvidVal();
			subData[vId]=val;
		}
		double cll=aPotFunc->getCondPotValueFor(subData);

		if(cll<1e-50)
		{
			cll=1e-50;
			thresholded++;
		}
		pll=pll+log(cll);
	}
	subData.clear();
	delete aPotFunc;
	return pll;
}
