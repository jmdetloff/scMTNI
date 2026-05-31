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
#include <string.h>
#include "Variable.H"
#include "Error.H"
#include "VariableManager.H"
#include "Potential.H"
#include "PotentialManager.H"
#include "SlimFactor.H"
#include "FactorGraph.H"
#include "SpeciesDataManager.H"

int
SpeciesDataManager::setVariableManager(VariableManager* aPtr)
{
	varMgr=aPtr;
	return 0;
}

int
SpeciesDataManager::createFactorGraph()
{
	vector<Variable*>& variableSet=varMgr->getVariableSet();
	
	fgraph=new FactorGraph;
	//Create the factors and the Markov blanket variables using the neighbours of each variable
    for(int globalFactorID=0;globalFactorID<variableSet.size();globalFactorID++)
	{
		SlimFactor* sFactor=new SlimFactor;
		sFactor->fId=globalFactorID;  //same as vIds
		fgraph->setFactor(sFactor);
	}
	return 0;
}

int 
SpeciesDataManager::setPotentialManager(PotentialManager* aPtr)
{
	potMgr=aPtr;
	return 0;
}

int
SpeciesDataManager::setOutputLoc(const char* aPtr)
{
	strcpy(outputLoc,aPtr);
	return 0;
}

int
SpeciesDataManager::setMotifNetwork(const char* aPtr)
{
	ifstream inFile(aPtr);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1024);
		if(strlen(buffer)<=0)
		{
			continue;
		}

		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string tfName;
		string tgtName;
		double score=0;

		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfName.append(tok);
			}
			else if(tokCnt==1)
			{
				tgtName.append(tok);
			}	
			else if(tokCnt==2)
			{
				score=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}

		int tfID=varMgr->getVarID(tfName.c_str());
		int tgtID=varMgr->getVarID(tgtName.c_str());
		if(tfID==-1 || tgtID==-1){
			continue;
		}
        unordered_map<int,double>* tgtSet=NULL;
		if(motifNetwork.find(tfID)==motifNetwork.end())
		{
			tgtSet=new unordered_map<int,double>;
			motifNetwork[tfID]=tgtSet;
		}
		else
		{
			tgtSet=motifNetwork[tfID];
		}
		(*tgtSet)[tgtID]=score;
		//cout << "TFvarID=" << tfID << " targetvarID=" << tgtID << " " << tfName << "->" << tgtName<<  "="<< score<<endl; 
	}
	cout << "motifNetwork.size() = " << motifNetwork.size() << endl;
	inFile.close();
	return 0;
}
	
VariableManager*
SpeciesDataManager::getVariableManager()
{
	return varMgr;
}

FactorGraph*
SpeciesDataManager::getFactorGraph()
{
	return fgraph;
}

PotentialManager*
SpeciesDataManager::getPotentialManager()
{
	return potMgr;
}

const char*
SpeciesDataManager::getOutputLoc()
{
	return outputLoc;
}

unordered_map<int,unordered_map<int,double>*>&
SpeciesDataManager::getMotifNetwork()
{
	return motifNetwork;
}
