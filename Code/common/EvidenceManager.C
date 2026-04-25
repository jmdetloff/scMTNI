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

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "Utils.H"

EvidenceManager::EvidenceManager()
{
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=-1;
    dataMat = NULL;
}

EvidenceManager::~EvidenceManager()
{
    if (dataMat!=NULL)
    {
        delete dataMat;
    }
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

VariableManager*
EvidenceManager::getVariableManager()
{
    return vMgr;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromTable(vector<string>& inputTable)
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
    cout <<"EvidenceManager::loadEvidenceFromTable sampleCount=" << sampleCount << " nodeCount=" << nodeCount << endl;
	//valuesPerSample.reserve(sampleCount);
	/*for (int sample = 0; sample < sampleCount; sample++)
	{
		vector<double> sampleVec (nodeCount);
		valuesPerSample.push_back(sampleVec);
	}*/
    dataMat=new Matrix(nodeCount,sampleCount);
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
            dataMat->setValue(varVal,nodeNum,sample);
		}
		++nodeNum;
	}
    substrs.clear();
	return Error::SUCCESS;
}

EMAP* 
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

Matrix*
EvidenceManager::getData()
{
    return dataMat;
}

int 
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int 
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int 
EvidenceManager::setPreRandomizeSplitSeed(int seed)
{
	randseed=seed;
	return 0;
}

int 
EvidenceManager::splitData(int s)
{
     
	int testSetSize=(dataMat->getColCnt()-validationIndex.size())/foldCnt;  //(evidenceSet.size()-validationIndex.size())/foldCnt
    if(foldCnt==1)
    {
        // only testdata, no training or validation data
        return 0;
    }
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		if(randseed<0)
		{
			randseed = rand();
		}
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size());	
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if(validationIndex.find(eInd)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
            testIndex.push_back(eInd); //testIndex[eInd]=0;
		}
		else
		{
            trainIndex.push_back(eInd); //trainIndex[eInd]=0;
		}
		m++;
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}

vector<int>&
EvidenceManager::getTrainingSet()
{
    return trainIndex;
}

vector<int>&
EvidenceManager::getTestSet()
{
    return testIndex;
}
INTINTMAP&
EvidenceManager::getValidationSet()
{	
	return validationIndex;
}

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}
