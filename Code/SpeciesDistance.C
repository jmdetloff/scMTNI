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
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "SpeciesDistance.H"

SpeciesDistance::SpeciesDistance()
{
	root=NULL;
	proot=0.5;
}

SpeciesDistance::~SpeciesDistance()
{
    for(auto pIter=speciesSet.begin();pIter!=speciesSet.end();pIter++)
    {
        delete pIter->second;
    }
    speciesSet.clear();
    binEdgeKeyProbMap.clear();
    speciesNameIDMap.clear();
}

int
SpeciesDistance::setProot(double aVal)
{
    proot=aVal;
    return 0;
}

int
SpeciesDistance::setspeciesNameIDMap(unordered_map<string,int>& spNameIDMap){
    speciesNameIDMap=spNameIDMap;
}

int 
SpeciesDistance::readSpeciesTree(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string childSpeciesName;
		string parentSpeciesName;
		double p_gain=0;
		double p_maintain_edge=0;
        while(tok!=NULL)
        {
            if(tokCnt==0)
            {
                childSpeciesName.append(tok);
            }
            else if(tokCnt==1)
            {
                parentSpeciesName.append(tok);
            }
            else if(tokCnt==2)
            {
                p_gain=atof(tok);
            }
            else if(tokCnt==3)
            {
                p_maintain_edge=1-atof(tok);
            }
            tok=strtok(NULL,"\t");
            tokCnt++;
        }
		SpeciesDistance::Species* childSpecies=NULL;
		SpeciesDistance::Species* parentSpecies=NULL;
		if(speciesSet.find(childSpeciesName)==speciesSet.end())
		{
			childSpecies=new SpeciesDistance::Species;
			childSpecies->name.append(childSpeciesName.c_str());
			childSpecies->parent=NULL;
			speciesSet[childSpeciesName]=childSpecies;
		}
		else
		{
			childSpecies=speciesSet[childSpeciesName];
		}
		if(speciesSet.find(parentSpeciesName)==speciesSet.end())
		{
			parentSpecies=new SpeciesDistance::Species;
			parentSpecies->name.append(parentSpeciesName.c_str());
			parentSpecies->parent=NULL;
			speciesSet[parentSpeciesName]=parentSpecies;
		}
		else
		{
			parentSpecies=speciesSet[parentSpeciesName];
		}
		childSpecies->parent=parentSpecies;
		if(parentSpecies->parent==NULL)
		{
			root=parentSpecies;
		}
		//Probabilities are for the child node
		childSpecies->p_maintain_edge=p_maintain_edge;
		childSpecies->p_gain=p_gain;
        parentSpecies->children.push_back(childSpecies);
        cout << "("<<childSpeciesName << "|" << parentSpeciesName << ") p_maintain_edge=" << p_maintain_edge << " p_gain=" << p_gain << endl;
	}
    cout <<"speciesSet.size()=" <<speciesSet.size() << " Root is " << root->name << endl;
	inFile.close();
	return 0;
}

// Returns the probability of a child status given a parent status, based on the gain and maintain probabilities on the child species.
static double
getParentChildEdgeStatusProb(SpeciesDistance::Species *childSpecies, int parentStatus, int childStatus) {
    if (parentStatus == 0 && childStatus == 0)
    {
        return 1 - childSpecies->p_gain;
    }
    else if (parentStatus == 0 && childStatus == 1)
    {
        return childSpecies->p_gain;
    }
    else if (parentStatus == 1 && childStatus == 0)
    {
        return 1 - childSpecies->p_maintain_edge;
    }
    else
    {
        return childSpecies->p_maintain_edge;
    }
}

// Creates a unique integer key for a configuration of n edge statuses, by interpreting a present edge as a 1 and an absent
// edge as a 0 in an n digit binary number.
static int
createBinaryKeyWithEdgeStatus(vector<int>& edgeStatus) {
    int binaryKey = 0;
    int position = 0;
    for(int eIter=0; eIter < edgeStatus.size(); eIter++)
    {
        if(edgeStatus[eIter] != 0)
        {
            binaryKey += (int)pow(2, position);
        }
        position += 1;
    }
    return binaryKey;
}

// Returns the probability of an edge having a particular configuration of statuses across species.
double
SpeciesDistance::getEdgeStatusProb(vector<int>& edgeStatus)
{
    // If the probability is already cached, return it.
    int binary_key = createBinaryKeyWithEdgeStatus(edgeStatus);
    if(binEdgeKeyProbMap.find(binary_key) != binEdgeKeyProbMap.end())
    {
        return binEdgeKeyProbMap[binary_key];
    }

    // Start with the probability of the root status.
    double edgeStatusProb = (edgeStatus[0] == 0) ? 1 - proot : proot;

    // Multiply by the probabilities of the edge assignment in all children.
    for(int i = 0; i < root->children.size(); i++)
    {
        double childrenScore = getSubTreeProb(edgeStatus[0], root->children[i], edgeStatus);
        edgeStatusProb *= childrenScore;
    }

    binEdgeKeyProbMap[binary_key] = edgeStatusProb;

    return edgeStatusProb;
}

// Returns the probability of an edge having a particular configuration of statuses across a subtree of species.
double
SpeciesDistance::getSubTreeProb(bool parentStatus, Species* child, vector<int>& edgeStatus)
{
    int childStatus = edgeStatus[speciesNameIDMap[child->name]];

    double score = getParentChildEdgeStatusProb(child, parentStatus, childStatus);

    // If the child is a leaf, just return its own score
    if(child->children.empty())
    {
        return score;
    }

    // If it's not a leaf, multiply score by the scores of its children.
    for(int i = 0; i < child->children.size(); i++)
    {
        score *= getSubTreeProb(childStatus, child->children[i], edgeStatus);
    }

    return score;
}

SpeciesDistance::Species*
SpeciesDistance::getRoot()
{
	return root;
}
