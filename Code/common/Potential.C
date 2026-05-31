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
#include <math.h>

#include "Variable.H"
#include "Potential.H"

#include "gsl/gsl_randist.h"

Potential::Potential()
{
	mean=NULL;
	covariance=NULL;
	inverse=NULL;
	matrixInd=0;
	mbcondVar=0;
	mbcondMean_Part=0;
}

Potential::~Potential()
{
	varSet.clear();
	factorVariables.clear();
	markovBlnktVariables.clear();
	if(covariance!=NULL)
	{
		delete covariance;
	}
	if(inverse!=NULL)
	{
		delete inverse;
	}
	if(mean!=NULL)
	{
		delete mean;
	}
	meanWrap.clear();
	vIDMatIndMap.clear();
	matIndvIDMap.clear();
	mbcondMean_Vect.clear();
}

int 
Potential::setAssocVariable(Variable* var,Potential::VariableRole vRole)
{
    varSet[var->getID()]=var;
	vIDMatIndMap[var->getID()]=matrixInd;
	matIndvIDMap[matrixInd]=var->getID();
	matrixInd++;
	switch(vRole)
	{
		case Potential::FACTOR:
		{
			factorVariables[var->getID()]=0;
			break;
		}
		case Potential::MARKOV_BNKT:
		{
			markovBlnktVariables[var->getID()]=0;
			break;
		}
	}
	return 0;
}

unordered_map<int,Variable*>&
Potential::getAssocVariables()
{
	return varSet;
}

// initialize meanWrap[INTDBLMAP], covariance[matrix], mean[1d matrix]
int
Potential::potZeroInit()
{
	for(auto vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		meanWrap[vIter->first]=0;
	}
	int row=varSet.size();
	covariance=new Matrix(row,row);
	covariance->setAllValues(0);
	mean=new Matrix(row,1);
	mean->setAllValues(0);
	return 0;
}

int 
Potential::updateMean(int vID,double mVal)
{
	int mID=vIDMatIndMap[vID];
	meanWrap[vID]=mVal;
	mean->setValue(mVal,mID,0);
	return 0;
}

int 
Potential::updateCovariance(int vID,int uID,double sVal)
{
	int mID=vIDMatIndMap[vID];
	int nID=vIDMatIndMap[uID];
	covariance->setValue(sVal,mID,nID);
	return 0;
}

int
Potential::makeValidJPD(gsl_matrix* ludecomp, gsl_permutation* p)
{
	if(inverse!=NULL)
	{
		delete inverse;
	}

	inverse=covariance->invMatrix(ludecomp,p);
	determinant=covariance->detMatrix(ludecomp,p);
	if(determinant <0)
	{
		cout <<"Negative Determinant " << determinant << endl;
		covariance->showMatrix();
		for(INTINTMAP_ITER vIter=matIndvIDMap.begin();vIter!=matIndvIDMap.end();vIter++)
		{
			int vId=vIter->second;
			cout << vIter->first << " " << vId << " "<< varSet[vId]->getName() << endl;
		}
	}
	double n=((double)varSet.size())/2.0;
	normFactor=pow(2*PI,n)*determinant;
	normFactor=sqrt(normFactor);
	return 0;
}

double
Potential::generateSample(INTDBLMAP& jointConf, int vId,gsl_rng* r)
{
	if(jointConf.find(factorVariables.begin()->first)==jointConf.end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}
	double newmean=0;
	for(auto aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		if(jointConf.find(aIter->first)==jointConf.end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}
		double aval=jointConf[aIter->first];
		newmean=newmean+(aval*aIter->second);
	}
	newmean=newmean+mbcondMean_Part;
	double x=gsl_ran_gaussian(r,sqrt(mbcondVar));
	x=x+newmean;
	return x;
}


int 
Potential::copyMe(Potential** apot)
{
	*apot=new Potential;
	for(INTINTMAP_ITER vIter=matIndvIDMap.begin();vIter!=matIndvIDMap.end();vIter++)
	{
		Variable* v=varSet[vIter->second];
		if(factorVariables.find(vIter->second)!=factorVariables.end())
		{
			(*apot)->setAssocVariable(v,Potential::FACTOR);
		}
		else
		{
			(*apot)->setAssocVariable(v,Potential::MARKOV_BNKT);
		}
	}
	
	for(INTDBLMAP_ITER sdIter=meanWrap.begin();sdIter!=meanWrap.end();sdIter++)
	{
		(*apot)->updateMean(sdIter->first,sdIter->second);
	}
	for(INTDBLMAP_ITER uIter=meanWrap.begin();uIter!=meanWrap.end();uIter++)
	{
		int i=vIDMatIndMap[uIter->first];
		for(INTDBLMAP_ITER vIter=meanWrap.begin();vIter!=meanWrap.end();vIter++)
		{
			int j=vIDMatIndMap[vIter->first];
			double cv=covariance->getValue(i,j);
			(*apot)->updateCovariance(i,j,cv);
		}
	}
	return 0;
}

// compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
int 
Potential::initMBCovMean()
{
    int vId=factorVariables.begin()->first;
    int vIdmId=vIDMatIndMap[vId];
    mbcondVar=covariance->getValue(vIdmId,vIdmId);
    mbcondMean_Part=meanWrap[vId];
    //cout << "My target ID is vId=" << vId << " vIdmId=" << vIdmId << " mbcondVar=" << mbcondVar << " mbcondMean_Part="  << mbcondMean_Part <<" meanWrap.size=" << meanWrap.size() <<" covariance="  << endl;
    //covariance->showMatrix(cout);
    if(markovBlnktVariables.size()==0)
    {
        return 0;
    }
    Matrix* mbcov=new Matrix(markovBlnktVariables.size(),markovBlnktVariables.size());
    Matrix* mbmargvar=new Matrix(1,markovBlnktVariables.size());
    INTINTMAP localMatIDMap;
    for(INTDBLMAP_ITER uIter=meanWrap.begin();uIter!=meanWrap.end();uIter++)
    {
        int i=vIDMatIndMap[uIter->first];
        int inew=i;
        if(i>vIdmId)
        {
            inew--;
        }
        for(INTDBLMAP_ITER vIter=meanWrap.begin();vIter!=meanWrap.end();vIter++)
        {
            if(vIter->first==vId)
            {
                continue;
            }
            int j=vIDMatIndMap[vIter->first];
            double cv=covariance->getValue(i,j);
            if(j>vIdmId)
            {
                j--;
            }
            if(uIter->first==vId)
            {
                mbmargvar->setValue(cv,0,j);
            }
            else
            {
                mbcov->setValue(cv,inew,j);
            }
        }
        if(uIter->first!=vId)
        {
            localMatIDMap[inew]=uIter->first;
        }
    }
    //cout <<"mbcov:" <<endl;
    //mbcov->showMatrix(cout);
    //cout <<"mbmargvar:" <<endl;
    //mbmargvar->showMatrix(cout);
    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);
    //cout << "Potential::initMBCovMean get my conditional weight: " ;
    for(INTINTMAP_ITER aIter=localMatIDMap.begin();aIter!=localMatIDMap.end();aIter++)
    {
        double aVal=prod1->getValue(0,aIter->first);
        double bVal=mbmargvar->getValue(0,aIter->first);
        mbcondVar=mbcondVar-(aVal*bVal);
        mbcondMean_Vect[aIter->second]=aVal;
        double cVal=meanWrap[aIter->second];
        mbcondMean_Part=mbcondMean_Part-(cVal*aVal);
        //cout << " aIter->first=" << aIter->first<<" mbcondMean_Vect[" << aIter->second << "]=" << aVal << " ";
    }
    //cout <<endl;
    localMatIDMap.clear();
    if(mbcondVar<1e-10)
    {
        //mbcondVar=1e-10;
    }
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    
    return 0;
}

double
Potential::getCondVariance()
{
	return mbcondVar;
}

double
Potential::getCondBias()
{
	return mbcondMean_Part;
}

unordered_map<int,double>&
Potential::getCondWeight()
{
	return mbcondMean_Vect;
}
