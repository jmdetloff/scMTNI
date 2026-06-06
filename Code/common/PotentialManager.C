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
#include "SlimFactor.H"
#include "PotentialManager.H"
#include <chrono>
#include "Utils.H"
#include "VariableManager.H"

using namespace std::chrono;

PotentialManager::PotentialManager()
{
    data=NULL;
    meanMat=NULL;
    covMat=NULL;
    testdataSize=0;
}

PotentialManager::~PotentialManager()
{
    if (data != NULL) {
        delete data;
    }
	if (meanMat != NULL) {
        delete meanMat;
    }
    if (covMat != NULL) {
        delete covMat;
    }
}

void
PotentialManager::deleteData()
{
    if (data != NULL) {
        delete data;
    }
}

int
PotentialManager::setOutputDir(const char* aDirName)
{
	strcpy(outputDir,aDirName);
	return 0;
}

void
PotentialManager::loadEvidenceFromTable(const vector<string>& inputTable, const vector<string>& variableList)
{
    // Reads gene names and expression levels from a tab-separated file. 
    // The first column is a gene name and remaining columns are expression levels.
    
    int varCount = variableList.size();

    // Count the columns and subtract 1 (for the gene name column) to get sample count.
    vector<string> substrs = Utils::split(inputTable[0], '\t');
    int sampleCount = substrs.size() - 1;

    data = new Matrix(varCount, sampleCount);
    vMgr = new VariableManager;

    for (auto line : inputTable) {
        substrs = Utils::split(line, '\t');
        string varName = substrs[0];

        vector<string>::const_iterator iter = find(variableList.begin(), variableList.end(), varName);

        // Variable list restricts which genes we include, so if it isn't present, skip it.
        if (iter == variableList.end()) {
            continue;
        }

        int varIndex = distance(variableList.begin(), iter);
        vMgr->addVariable(varName, varIndex);

        for (int sample = 0; sample < sampleCount; sample++) {
            double varVal = stod(substrs[sample + 1]);
            data->setValue(varVal, varIndex, sample);
        }
    }

    testdataSize = sampleCount;

    cout << "PotentialManager::loadEvidenceFromTable varCount=" << data->getRowCnt() << " sampleCount=" << data->getColCnt()  << endl;

    int varCnt = data->getRowCnt();
    meanMat = new Matrix(varCnt, 1);
    meanMat->setAllValues(0);
    covMat = new Matrix(varCnt, varCnt);
    covMat->setAllValues(-1);
}

VariableManager*
PotentialManager::getVariableManager()
{
    return vMgr;
}

int
PotentialManager::estimateCovariance(int uId, int vId)
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

double
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, int& status)
{
    int mbsize = sFactor->mergedMB.size();
    if (mbsize == 0) {
        return 0;
    }

    vector<int> parChildID(sFactor->mergedMB.begin(), sFactor->mergedMB.end());
    parChildID.push_back(sFactor->fId);

    //extract covariance
    Matrix* covariance = new Matrix(parChildID.size(), parChildID.size());
    for (int i = 0; i < parChildID.size(); i++) {
        int vID = parChildID[i];
        for (int j = i; j < parChildID.size(); j++) {
            int uID = parChildID[j];
            double cval = covMat->getValue(vID, uID);
            covariance->setValue(cval, i, j);
            covariance->setValue(cval, j, i);
        }
    }

    double determinant = covariance->detMatrix();
    if (determinant <= 0) {
        status = -1;
        return 0;
    }

    int vId = sFactor->fId;
    double mbcondVar = covMat->getValue(vId, vId);
    
    Matrix* mbcov = covariance->getSubMatrix(0, 0, mbsize, mbsize);
    Matrix* mbmargvar = covariance->getSubMatrix(parChildID.size() - 1, 0, 1, mbsize);

    double determinantP = mbcov->detMatrix();
    if (determinantP <= 0) {
        status = -1;
        return 0;
    }

    Matrix* covInv = mbcov->invMatrix();
    Matrix* prod1 = mbmargvar->multiplyMatrix(covInv);

    for (int i = 0; i < parChildID.size() - 1; i++) {
        int varID = parChildID[i];
        double aVal = prod1->getValue(0, i);
        double bVal = mbmargvar->getValue(0, i);
        mbcondVar = mbcondVar - (aVal * bVal);
    }

    if (mbcondVar < 0) {
        status=-1;
        return 0;
    }

    double jointll1 = computeLL(parChildID.size(), determinant);
    double jointll2 = computeLL(mbsize, determinantP);
    double pll = jointll1 - jointll2;

    if (isinf(pll)) {
        status=-1;
        return 0;
    }

    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;

    return pll;
}

void
PotentialManager::dumpVarMB(SlimFactor* sFactor, ofstream& oFile)
{
    int mbsize = sFactor->mergedMB.size();
    if (mbsize == 0) {
        return;
    }

    Variable* target = vMgr->getVariable(sFactor->fId);

    vector<int> parChildID(sFactor->mergedMB.begin(), sFactor->mergedMB.end());
    parChildID.push_back(sFactor->fId);

    //extract covariance
    Matrix* covariance = new Matrix(parChildID.size(), parChildID.size());
    for(int i = 0; i < parChildID.size(); i++) {
        int vID = parChildID[i];
        for(int j = i; j < parChildID.size(); j++) {
            int uID = parChildID[j];
            double cval = covMat->getValue(vID, uID);
            covariance->setValue(cval, i, j);
            covariance->setValue(cval, j, i);
        }
    }

    Matrix* mbcov = covariance->getSubMatrix(0, 0, mbsize, mbsize);
    Matrix* mbmargvar = covariance->getSubMatrix(parChildID.size() - 1, 0, 1, mbsize);

    Matrix* covInv = mbcov->invMatrix();
    Matrix* prod1 = mbmargvar->multiplyMatrix(covInv);

    // Start from parent only, last one id=ParChildID.size()-1 is target
    for(int i = 0; i < parChildID.size() - 1; i++) {
        Variable* parent = vMgr->getVariable(parChildID[i]);
        double conditionalWeight = prod1->getValue(0, i);
        oFile << parent->getName() << "\t"<< target->getName() << "\t" << conditionalWeight << endl;
    }

    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
}

void
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondVar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect)
{
    int mbsize = sFactor->mergedMB.size();
    if (mbsize == 0) {
        return;
    }

    vector <int> parChildID(sFactor->mergedMB.begin(), sFactor->mergedMB.end());
    parChildID.push_back(sFactor->fId);

    //extract covariance
    Matrix* covariance = new Matrix(parChildID.size(), parChildID.size());
    for(int i = 0; i < parChildID.size(); i++) {
        int vID = parChildID[i];
        for(int j = i; j < parChildID.size(); j++) {
            int uID = parChildID[j];
            double cval = covMat->getValue(vID, uID);
            covariance->setValue(cval, i, j);
            covariance->setValue(cval, j, i);
        }
    }

    int vId = sFactor->fId;
    mbcondVar = covMat->getValue(vId, vId);
    mbbias = meanMat->getValue(vId, 0);

    Matrix* mbcov = covariance->getSubMatrix(0, 0, mbsize, mbsize);
    Matrix* mbmargvar = covariance->getSubMatrix(parChildID.size() - 1, 0, 1, mbsize);
    Matrix* covInv = mbcov->invMatrix();
    Matrix* prod1 = mbmargvar->multiplyMatrix(covInv);

    // Start from parent only, last one is target
    for(int i = 0; i < parChildID.size() - 1; i++) {
        int varID = parChildID[i];
        double conditionalWeight = prod1->getValue(0, i);
        double bVal = mbmargvar->getValue(0, i);
        mbcondVar = mbcondVar - (conditionalWeight * bVal);
        mbcondMean_Vect[varID] = conditionalWeight;
        double cVal = meanMat->getValue(varID, 0);
        mbbias = mbbias - (cVal * conditionalWeight);
    }

    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
}

double
PotentialManager::computeLL(int dim, double determinant)
{
    double ll=testdataSize * (dim*log(2*PI)+log(determinant));
    double t=dim*(testdataSize-1);
    ll=(ll+t)*(-0.5);
    return ll;
}

double
PotentialManager::computeMeanVarPseudoLikelihood(int id)
{
    // compute meanMat and Variance on covMat, PseudoLikelihood
    double vmean=data->RowMean(id);
    meanMat->setValue(vmean,id,0);
    double dataSetSize=data->getColCnt();
    double ssd=data->vectorMultiply(id,vmean,id,vmean);
    //add 1e-10 to avoid singularity issues:
    double variance=ssd/(dataSetSize-1.0)+1e-10;  //(0.001+ssd)/((double)(data->getColCnt()-1))
    covMat->setValue(variance,id,id);
    double pll=-0.5*ssd/variance-0.5*log(2.0*PI*variance)*dataSetSize;
    for(int j=id+1;j<data->getRowCnt();j++)
    {
        estimateCovariance(id,j);
    }
    return pll;
}
