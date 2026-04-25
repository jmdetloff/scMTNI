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

#include <cstring>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "Error.H"
#include "Utils.H"
#include "Variable.H"
#include "VariableManager.H"

Error::ErrorCode
VariableManager::readVariablesFromTable(vector<string>& inputTable)
{
	int nodeCount = inputTable.size();
	int nodeID = 0;
	for (auto line : inputTable)
	{
		vector<string> substrs = Utils::split(line, '\t');
		string nodeName = substrs[0];
		Variable* var = new Variable;
		var->setID(nodeID);
		var->setName(nodeName.c_str());
        variableSet.push_back(var);
		varNameIDMap[nodeName] = nodeID;
		++nodeID;
	}
	cout << "Read information about " << nodeCount << " variables" << endl;
	return Error::SUCCESS;
}

Error::ErrorCode
VariableManager::readVariablesFromData(string nodeName,int nodeID)
{
    Variable* var = new Variable;
    var->setID(nodeID);
    var->setName(nodeName.c_str());
    variableSet.push_back(var);
    varNameIDMap[nodeName] = nodeID;
    return Error::SUCCESS;
}

int
VariableManager::getVarID(const char* varName)
{
	string varKey(varName);
	if(varNameIDMap.find(varKey)==varNameIDMap.end())
	{
		return -1;
	}
	int vId=varNameIDMap[varKey];
	return vId;
}

vector<Variable*>&
VariableManager::getVariableSet()
{
	return variableSet;
}

Variable* 
VariableManager::getVariableAt(int vId)
{
	return variableSet[vId];
}
