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
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <utility>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

#include "FactorGraph.H"
#include "MetaMove.H"

#include "Utils.H"

#include "SpeciesDistance.H"
#include "SpeciesDataManager.H"
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "MetaLearner.H"
#include <chrono>
#include <unistd.h>
using namespace std::chrono;
using namespace std;

MetaLearner::MetaLearner()
{
    restrictedFName[0] = '\0';
    trueGraphFName[0] = '\0';
    convThreshold = 1e-4;
    beta1 = -0.9;
    beta2 = 4.0;
    INDEP = false;
    initGlobalScore = 0;
    splitGenes = false;
}

MetaLearner::~MetaLearner()
{
    inputVariables.clear();
    inputOGList.clear();
    inputRegulatorOGs.clear();
}

int MetaLearner::setInputFName(const char *aFName)
{
    // strcpy(inputFName,aFName);
    inputFName = aFName;
    return 0;
}

int MetaLearner::setMaxFactorSize(int aVal)
{
    maxFactorSize = aVal;
    return 0;
}

int MetaLearner::setMaxFactorSize_Approx(int aVal)
{
    maxFactorSizeApprox = aVal;
    return 0;
}

int MetaLearner::setINDEP()
{
    INDEP = true;
    return 0;
}

int MetaLearner::setsplitGenes()
{
    splitGenes = true;
    return 0;
}

int MetaLearner::setConvergenceThreshold(double aVal)
{
    convThreshold = aVal;
    return 0;
}

int MetaLearner::setRestrictedList(const char *aFName)
{
    strcpy(restrictedFName, aFName);
    ifstream inFile(restrictedFName);
    string buffer;
    while (inFile.good())
    {
        getline(inFile, buffer);
        if (buffer.length() <= 0)
        {
            continue;
        }
        int ogid = atoi(buffer.c_str());
        inputRegulatorOGs.push_back(ogid); // inputRegulatorOGs[ogid]=0;
        inputVariables.insert(ogid);
    }
    inFile.close();
    /*for(int i=0;i<inputRegulatorOGs.size();i++){
        cout <<"inputRegulatorOGs[" << i << "]=" << inputRegulatorOGs[i] << endl;
    }*/
    return 0;
}

int MetaLearner::setTrueGraph(const char *aGraphName)
{
    strcpy(trueGraphFName, aGraphName);
    return 0;
}

int MetaLearner::setInputOGList(const char *aFName)
{
    ifstream inFile(aFName);
    char buffer[1024];
    while (inFile.good())
    {
        inFile.getline(buffer, 1023);
        if (strlen(buffer) <= 0)
        {
            continue;
        }
        int og = atoi(buffer);
        inputOGList.push_back(og); // inputOGList[og]=0;
        inputVariables.insert(og);
    }
    inFile.close();
    /*for(int i=0;i<inputOGList.size();i++){
        cout <<"inputOGList[" << i << "]=" << inputOGList[i] << endl;
    }*/
    return 0;
}

int MetaLearner::setOrthogroupReader(MappedOrthogroupReader *aPtr)
{
    ogr = aPtr;
    return 0;
}

int MetaLearner::setSpeciesDistances(SpeciesDistance *aPtr)
{
    speciesData = aPtr;
    speciesData->setspeciesNameIDMap(speciesNameIDMap);
    return 0;
}

int MetaLearner::setBeta1(double b1)
{
    beta1 = b1;
    return 0;
}

int MetaLearner::setBeta2(double b2)
{
    beta2 = b2;
    return 0;
}

void MetaLearner::process_mem_usage(double &vm_usage, double &resident_set)
{
    vm_usage = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0 / 1024.0;                // MB
    resident_set = rss * page_size_kb;
}

int MetaLearner::init()
{
    ifstream inFile(inputFName);
    cout << "MetaLearner::init() read:" << inputFName << endl;
    string buffer;
    int datasetId = 0;
    if (splitGenes)
    {
        setinputVariableNames();
    }
    // set sparsity prior: edgePresenceProb
    map<int, MappedOrthogroup *> &orthogroupSet = ogr->getMappedOrthogroups();
    while (inFile.good())
    {
        getline(inFile, buffer);
        if (buffer.empty() || buffer.find("#") == 0)
        {
            continue;
        }
        vector<string> strs;

        // strip trailing newlines
        buffer.erase(buffer.find_last_not_of(" \n\r\t") + 1);

        strs = Utils::split(buffer, '\t');
        assert(strs.size() >= 3);
        string specName = strs[0];
        string datasetSuff = strs[1];
        string outputLoc = strs[2];
        // string regulatorFName = strs[3];
        // string targetFName = strs[4];
        string motifNetwork = strs[5];

        PotentialManager *potMgr = new PotentialManager;
        string tableFileName = datasetSuff;
        if (datasetSuff.find('.') == std::string::npos)
        {
            tableFileName = datasetSuff + ".table";
            std::cout << "no .table in input expression filename, add .table suffix:" << tableFileName << endl;
            
        }
        readEvidenceTable(tableFileName);
        Error::ErrorCode eCode;
        if (splitGenes)
        {
            eCode = potMgr->loadEvidenceFromTable(inputTable, inputVariableNames[datasetId]);
        }
        else
        {
            eCode = potMgr->loadEvidenceFromTable(inputTable);
        }
        if (eCode != Error::SUCCESS)
        {
            cout << Error::getErrorString(eCode) << endl;
            return -1;
        }
        VariableManager *varMgr = potMgr->getVariableManager();
        potMgr->setOutputDir(outputLoc.c_str());
        SpeciesDataManager *spMgr = new SpeciesDataManager;
        speciesDataSet.push_back(spMgr);
        spMgr->setVariableManager(varMgr);
        spMgr->setPotentialManager(potMgr);
        spMgr->createFactorGraph();
        spMgr->setOutputLoc(outputLoc.c_str());
        spMgr->setMotifNetwork(motifNetwork.c_str()); // add motif prior network in CVN for single cell
        speciesIDNameMap.push_back(specName);         // speciesIDNameMap[datasetId]=specName;
        speciesNameIDMap[specName] = datasetId;
        cout << datasetId << "=" << specName << " motifNetwork=" << motifNetwork << endl;

        // move one fold cross-validation here:
        potMgr->reset();
        potMgr->init();
        char foldOutputDirCmd[1024];
        sprintf(foldOutputDirCmd, "mkdir -p %s/fold0", outputLoc.c_str());
        system(foldOutputDirCmd);
        FactorGraph *condspecGraph = spMgr->getFactorGraph();
        // initialize factor graph with pll score:
        // sparsity prior for each target:
        vector<double> varNeighborhoodPrior;
        vector<Variable *> &varSet = varMgr->getVariableSet();
        vector<vector<double>> edgePresenceProb(condspecGraph->getFactorCnt(), vector<double>(condspecGraph->getFactorCnt(), 0));
        for (int f = 0; f < condspecGraph->getFactorCnt(); f++)
        {
            SlimFactor *sFactor = condspecGraph->getFactorAt(f);
            // double pll=getPLLScore_Condition_onefold((string&)eIter->first,sFactor);
            double pll = potMgr->computeMeanVarPseudoLikelihood_onefold(sFactor->fId);
            Variable *target = varSet[sFactor->fId];
            double priorScore = precomputePerSpeciesPrior(datasetId, sFactor->fId, target, spMgr, edgePresenceProb, orthogroupSet);
            varNeighborhoodPrior.push_back(priorScore);
            //cout << "cell=" << datasetId << " f=" << f << " sFactor->fId=" << sFactor->fId << " pll=" << pll << " priorScore=" << priorScore << endl;
            sFactor->mbScore = pll + priorScore;
            initGlobalScore = initGlobalScore + pll + priorScore;
            // sFactor->marginalLL=pll;
        }
        varNeighborhoodPrior_PerSpecies.push_back(varNeighborhoodPrior);
        edgePresenceProb_PerSpecies.push_back(edgePresenceProb);
        potMgr->deleteData();
        //cout << "cell=" << datasetId << " variable size=" << condspecGraph->getFactorCnt() << endl;
        datasetId += 1; //=datasetId*2;
        strs.clear();
        varNeighborhoodPrior.clear();
        for (int i = 0; i < edgePresenceProb.size(); i++)
        {
            edgePresenceProb[i].clear();
        }
        edgePresenceProb.clear();
    }
    inFile.close();
    inputTable.clear();
    cout << "-------------------------------------------------------------" << endl;
    return 0;
}

void MetaLearner::readEvidenceTable(string fileName)
{
    /*
     Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
     to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.
     This method just loads the meaningful lines into a vector that is then passed to VariableManager::readVariablesFromTable().
     It's done this way so as not to break encapsulation of private members of VariableManager.
     */
    cout << "readEvidenceTable:" << fileName << endl;
    ifstream inFile(fileName);
    assert(inFile);

    bool headerLine = true;
    string inputLine;
    inputTable.clear();

    // Read all the lines with meaningful data
    while (getline(inFile, inputLine))
    {
        if (inputLine.empty() || inputLine.find('#') == 0)
        {
            continue;
        }
        if (headerLine)
        {
            headerLine = false;
            continue;
        }
        inputTable.push_back(inputLine);
    }
}

void MetaLearner::start()
{
    if (strlen(trueGraphFName) > 0) {
        double globalScore = getScore();
        generateData(10000, 30000);
        showModelParameters();
        cout << "Final Score\t" << globalScore << endl;
        return;
    }

    auto start = high_resolution_clock::now();
    cout << "MetaLearner::start" << endl;
    cout << "inputRegulatorOGs.size() = " << inputRegulatorOGs.size() << " inputOGList.size() = " << inputOGList.size() << endl;

    double currGlobalScore = initGlobalScore;

    if (!INDEP) {
        precomputeEmptyGraphPrior();
        initCondsetMap_Tree(speciesData->getRoot());
    }

    int iter = 0;
    bool notConverged = true;

    // hardcode 100 iterations for now
    while (notConverged && iter < 100) {

        // collect the candidate edges
        if (INDEP) {
            collectMoves_Orthogroups_INDEP(maxFactorSizeApprox);
        } else {
            collectMoves_Orthogroups(maxFactorSizeApprox);
        }

        int successMove = 0;
        double diff = makeMoves(successMove);
        double priorChange = INDEP ? 0 : getPriorDelta();
        double newScore = currGlobalScore + diff;

        if (diff <= convThreshold) {
            notConverged = false;
        }

        currGlobalScore = newScore;

        double vm, rss;
        process_mem_usage(vm, rss);
        cout << "ITERATION " << iter << " newScore=" << newScore << " diffscore=" << diff << " priorChange=" << priorChange << " successMove=" << successMove << endl;

        iter++;
    }

    dumpAllGraphs(maxFactorSizeApprox);
    cout << "Final Score " << currGlobalScore << endl;
    showModelParameters();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "MetaLearner::start runtime: " << duration.count() << " ms" << endl;

    char scoreFName[1024];
    sprintf(scoreFName, "%s/scoreFile.txt", speciesDataSet[0]->getOutputLoc());

    ofstream sFile(scoreFName);
    sFile << currGlobalScore << endl;
    sFile.close();
}

double
MetaLearner::getScore()
{
    double gScore = 0;

    for (int i = 0; i < speciesDataSet.size(); i++)
    {
        FactorGraph *fg = speciesDataSet[i]->getFactorGraph();
        vector<SlimFactor *> &factorSet = fg->getAllFactors();
        for (int fIter = 0; fIter < factorSet.size(); fIter++)
        {
            SlimFactor *sFactor = factorSet[fIter];
            gScore = gScore + sFactor->mbScore;
        }
    }
    return gScore;
}

double
MetaLearner::getPriorDelta()
{
    double oldStructPrior = 0;
    double newStructPrior = 0;
    // cout << "MetaLearner::getPriorDelta() Update the OldpriorScore" << endl;
    // Need to consider the old contribution of the edges, delete that from the overall prior and add the new contribution
    for (auto edgeIter = affectedOGPairs.begin(); edgeIter != affectedOGPairs.end(); edgeIter++)
    {
        vector<string> keyid = Utils::split(edgeIter->first, '-');
        int regi = stoi(keyid[0]);    // edgeIter->first.first;
        int targeti = stoi(keyid[1]); // edgeIter->first.second;
        double aval = speciesData->getEdgeStatusProb(*(edgeIter->second));
        double edgePrior = log(aval);
        double oldEdgePrior = ogpairPrior[regi][targeti]; // oldEdgePrior=ogpairPrior[edgeIter->first];
        oldStructPrior += oldEdgePrior;
        newStructPrior += edgePrior;
        /*INTINTMAP* currEdgeStatus=edgeConditionMap[edgeIter->first];
        STRINTMAP* newEdgeStatus=edgeIter->second;
        for(STRINTMAP_ITER sIter=newEdgeStatus->begin();sIter!=newEdgeStatus->end();sIter++)
        {
            int specID=speciesNameIDMap[sIter->first];
            (*currEdgeStatus)[specID]=sIter->second;  //edgeConditionMap[edgeIter->first][specID]=sIter->second
        }*/
        ogpairPrior[regi][targeti] = edgePrior; // ogpairPrior[edgeIter->first]=edgePrior;
        // cout <<" ogpairPrior[" << regi<<"][" << targeti << "]=" <<ogpairPrior[regi][targeti] << endl;
        keyid.clear();
    }
    for (auto edgeIter = affectedOGPairs.begin(); edgeIter != affectedOGPairs.end(); edgeIter++)
    {
        edgeIter->second->clear();
        delete edgeIter->second;
    }
    affectedOGPairs.clear();
    double priorDelta = newStructPrior - oldStructPrior;
    return priorDelta;
}

double
MetaLearner::precomputeEmptyGraphPrior()
{
    vector<int> edgeStatus(speciesIDNameMap.size(), 0);
    double prior = speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior = log(prior);
    int edgeCnt = 0;
    for (int regi = 0; regi < inputRegulatorOGs.size(); regi++)
    {
        int ogno_tf = inputRegulatorOGs[regi];
        vector<double> targetPrior;
        for (int targeti = 0; targeti < inputOGList.size(); targeti++)
        {
            int ogno_tgt = inputOGList[targeti];
            if (ogno_tf == ogno_tgt)
            {
                targetPrior.push_back(0);
                continue;
            }
            targetPrior.push_back(logPrior);
            edgeCnt++;
        }
        ogpairPrior.push_back(targetPrior);
    }
    double emptyGraphPrior = edgeCnt * logPrior;
    return emptyGraphPrior;
}

// Shilu: compute varNeighborhoodPrior_PerSpecies[specID][targeti]: set to 0
// edgePresenceProb is the same for CVN for each edge per cell, so no need to store into a map
// targetID is variable ID
double
MetaLearner::precomputePerSpeciesPrior(int specID, int targetID, Variable *target, SpeciesDataManager *sdm, vector<vector<double>> &edgePresenceProb, map<int, MappedOrthogroup *> &orthogroupSet)
{
    // double initPrior=edgePresenceProb;
    string spec = speciesIDNameMap[specID];
    int ogno_tgt = ogr->getMappedOrthogroupID(target->getName().c_str(), spec.c_str());
    double NeighborhoodPrior = 0;
    for (int regi = 0; regi < inputRegulatorOGs.size(); regi++)
    {
        int ogno_tf = inputRegulatorOGs[regi];
        if (ogno_tf == ogno_tgt)
        {
            continue;
        }
        MappedOrthogroup *tfogrp = orthogroupSet[ogno_tf]; // oIter->second; four species
        vector<string> &tfgrpMembers = tfogrp->getOrthoMembers();
        VariableManager *vMgr = sdm->getVariableManager();
        int regID = vMgr->getVarID(tfgrpMembers[specID].c_str());
        if (regID == -1)
        {
            // cout << tfgrpMembers[specID]<< " not in " << spec << endl;
            continue;
        }
        double initPrior = getEdgePrior_PerSpecies(regID, targetID, sdm);
        edgePresenceProb[regID][targetID] = initPrior;
        NeighborhoodPrior += log(1 - initPrior);
        //cout << spec << " regi=" << regi << " TF_OGID=" << ogno_tf << " target_OGID=" << ogno_tgt << " regVarID=" << regID << " targetVarID=" << targetID << " " << tfgrpMembers[specID] << "->" << target->getName() << ": initPrior=" << initPrior << endl;
    }
    //cout << spec << "target=" << target->getName() << " NeighborhoodPrior=" << NeighborhoodPrior << endl;
    return NeighborhoodPrior;
}

int MetaLearner::initCondsetMap_Tree(SpeciesDistance::Species *node)
{
    // cout <<"rootnode=" << node->name << endl;
    //  this is root, it has 0/1 two options:
    int npos = pow(2, speciesIDNameMap.size());
    vector<unordered_map<int, int>> mycondition, mycondition1;
    int id = speciesNameIDMap[node->name];
    unordered_map<int, int> map0, map1;
    map0[id] = 0;
    mycondition.push_back(map0);
    map1[id] = 1;
    mycondition1.push_back(map1);
    map0.clear();
    map1.clear();
    if (!node->children.empty())
    {
        for (int i = 0; i < node->children.size(); i++)
        {
            mycondition = initCondsetMap_Tree_backtrack(node->children[i], mycondition, 0, 0);
            mycondition1 = initCondsetMap_Tree_backtrack(node->children[i], mycondition1, 1, 0);
        }
    }
    mycondition.insert(mycondition.end(), mycondition1.begin(), mycondition1.end());

    // copy mycondition to condsetMap_Tree
    // remove all 0s combination from mycondition:
    for (int i = 0; i < mycondition.size(); i++)
    {
        // cout<<"condition " << i << ": ";
        vector<int> cset(speciesIDNameMap.size(), 0);
        int zerocount = 0;
        for (auto sIter = mycondition[i].begin(); sIter != mycondition[i].end(); sIter++)
        {
            cset[sIter->first] = sIter->second; // cell id: sIter->first, edge status: sIter->second
            // cout << speciesIDNameMap[sIter->first] << "=" <<sIter->second <<" ";
            zerocount += 1 - sIter->second; // count number of 0s in the condition
        }
        if (zerocount < speciesIDNameMap.size())
        {
            condsetMap_Tree.push_back(cset);
        }
        mycondition[i].clear();
        cset.clear();
        // cout << endl;
    }
    // add cell-specific condition for intermediate cells (already included)
    for (auto it = intermediatecells.begin(); it != intermediatecells.end(); it++)
    {
        string species = speciesIDNameMap[*it];
        // cout << "add cell-specific condition for intermediate cell: " <<species << endl;
        //  only in this cell:
        vector<int> cset1; // cset0;
        for (int j = 0; j < speciesIDNameMap.size(); j++)
        {
            string specname = speciesIDNameMap[j];
            if (j == *it)
            {
                // cset0.push_back(0);
                cset1.push_back(1);
            }
            else
            {
                // cset0.push_back(1);
                cset1.push_back(0);
            }
        }
        // condsetMap_Tree.push_back(cset0);
        condsetMap_Tree.push_back(cset1);
        // cset0.clear();
        cset1.clear();
    }

    mycondition.clear();
    for (int i = 0; i < mycondition1.size(); i++)
    {
        mycondition1[i].clear();
    }
    mycondition1.clear();
    return 0;
}

// shilu: use backtracking to set up constrained conditions:
vector<unordered_map<int, int>>
MetaLearner::initCondsetMap_Tree_backtrack(SpeciesDistance::Species *node, vector<unordered_map<int, int>> mycondition, int parentstatus, int ntransition)
{

    int id = speciesNameIDMap[node->name];
    // cout <<"node=" << node->name << " id=" << id <<" parentstatus=" <<parentstatus<<" ntransition=" << ntransition << endl;
    //  intermediate cell:
    if ((!node->children.empty()) && node->parent != NULL)
    {
        // cout << "intermediate cell: " <<node->name << endl;
        intermediatecells.insert(id);
        if (ntransition == 0)
        {
            // not reach max transition yet in this branch, it has 0/1 two options:
            vector<unordered_map<int, int>> outcondition, cond1;
            for (int r = 0; r < mycondition.size(); r++)
            {
                unordered_map<int, int> temp1(mycondition[r]);
                temp1[id] = parentstatus;
                outcondition.push_back(temp1);
                unordered_map<int, int> temp2(mycondition[r]);
                temp2[id] = 1 - parentstatus;
                cond1.push_back(temp2);
                temp1.clear();
                temp2.clear();
            }
            if (!node->children.empty())
            {
                for (int i = 0; i < node->children.size(); i++)
                {
                    outcondition = initCondsetMap_Tree_backtrack(node->children[i], outcondition, parentstatus, ntransition);
                    cond1 = initCondsetMap_Tree_backtrack(node->children[i], cond1, 1 - parentstatus, ntransition + 1);
                }
            }
            outcondition.insert(outcondition.end(), cond1.begin(), cond1.end());

            for (int r = 0; r < mycondition.size(); r++)
            {
                mycondition[r].clear();
                cond1[r].clear();
            }
            cond1.clear();
            mycondition.clear();
            return outcondition;
        }
        else if (ntransition == 1)
        {
            // transition has occurred, it needs to be the same as parent:
            for (int r = 0; r < mycondition.size(); r++)
            {
                mycondition[r][id] = parentstatus;
            }
            if (!node->children.empty())
            {
                for (int i = 0; i < node->children.size(); i++)
                {
                    mycondition = initCondsetMap_Tree_backtrack(node->children[i], mycondition, parentstatus, ntransition);
                }
            }
            return mycondition;
        }
    }

    // reach leaves
    if (node->children.empty())
    {
        if (ntransition == 0)
        {
            // not reach max transition yet in this branch, it has 0/1 two options:
            vector<unordered_map<int, int>> outcondition;
            for (int r = 0; r < mycondition.size(); r++)
            {
                unordered_map<int, int> temp1(mycondition[r]);
                temp1[id] = parentstatus;
                outcondition.push_back(temp1);
                unordered_map<int, int> temp2(mycondition[r]);
                temp2[id] = 1 - parentstatus;
                outcondition.push_back(temp2);
                temp1.clear();
                temp2.clear();
                mycondition[r].clear();
            }
            mycondition.clear();
            return outcondition;
        }
        else if (ntransition == 1)
        {
            // transit twice, it needs to be the same as parent:
            for (int r = 0; r < mycondition.size(); r++)
            {
                mycondition[r][id] = parentstatus;
            }
            return mycondition;
        }
    }
}

// update by Shilu get all its children, not just leaves!
int MetaLearner::getLeaves(SpeciesDistance::Species *node, map<string, int> &leaves)
{
    if (node->children.empty())
    {
        leaves[node->name] = 0;
        return 0;
    }
    if (!node->children.empty())
    {
        for (int i = 0; i < node->children.size(); i++)
        {
            leaves[node->children[i]->name] = 0;
            getLeaves(node->children[i], leaves);
        }
    }
    return 0;
}

int MetaLearner::setinputVariableNames()
{
    //cout << "MetaLearner::setinputVariableNames()" << endl;
    map<int, MappedOrthogroup *> &orthogroupSet = ogr->getMappedOrthogroups();
    for (auto og = inputVariables.begin(); og != inputVariables.end(); og++)
    {
        MappedOrthogroup *targetogrp = orthogroupSet[*og];
        vector<string> &targetgrpMembers = targetogrp->getOrthoMembers();
        for (int specID = 0; specID < targetgrpMembers.size(); specID++)
        {
            if (strcmp(targetgrpMembers[specID].c_str(), "None") != 0)
            {
                inputVariableNames[specID].insert(targetgrpMembers[specID]);
                //cout << "og=" << *og << " specID=" << specID << " varName(inputVariableNames[specID])=" << targetgrpMembers[specID] << endl;
            }
        }
    }
    cout << "MetaLearner::setinputVariableNames number of variables: " << inputVariables.size() << endl;
    return 0;
}

int MetaLearner::collectMoves_Orthogroups(int currK)
{
    // auto start = high_resolution_clock::now();
    for (int i = 0; i < moveSet.size(); i++)
    {
        delete moveSet[i];
    }
    moveSet.clear();
    map<int, MappedOrthogroup *> &orthogroupSet = ogr->getMappedOrthogroups();

    // Now we will have a move for one orthogroup at a time
    for (int targeti = 0; targeti < inputOGList.size(); targeti++)
    {
        int targetOGIter = inputOGList[targeti];
        MappedOrthogroup *targetogrp = orthogroupSet[targetOGIter]; // oIter->second; four species
        // map<string,GeneMap*>& targetgrpMembers=targetogrp->getOrthoMembers();
        vector<string> &targetgrpMembers = targetogrp->getOrthoMembers();
        int n = speciesIDNameMap.size();

        // score the score of best regulator:
        vector<double> bestscore_PerSpecies(n, 0);
        vector<double> bestscoreImprovement_PerSpecies(n, 0);
        vector<int> besttarget_PerSpecies(n, -1);
        vector<int> besttf_PerSpecies(n, -1);
        // unordered_map<int,unordered_map<int,double>*> bestregWt_PerSpecies; only need for each species!
        int bestcsetid = -1;
        int bestregi = -1; // inputRegulatorOGs vector index
        // vector<double> bestScoreImprovementTF(n,0);
        double bestScoreImprovement_TF = 0;
        // The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        // term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.

        // for each regulator
        for (int regi = 0; regi < inputRegulatorOGs.size(); regi++)
        {
            int regOGIter = inputRegulatorOGs[regi];
            if (regOGIter == targetOGIter) // if(regOGIter->first==oIter->first)
            {
                continue;
            }
            MappedOrthogroup *tfogrp = orthogroupSet[regOGIter]; // orthogroupSet[regOGIter->first];
            // map<string,GeneMap*>& tfgrpMembers=tfogrp->getOrthoMembers();
            vector<string> &tfgrpMembers = tfogrp->getOrthoMembers();
            double oldpriorScore = ogpairPrior[regi][targeti]; // double oldpriorScore=ogpairPrior[ogPair];

            // only store current tf for each species!
            vector<double> score_PerSpecies(n, 0);
            vector<double> scoreImprovement_PerSpecies(n, 0);
            vector<int> target_PerSpecies(n, 0); // store variable index
            vector<int> tf_PerSpecies(n, 0);
            // map<int,INTDBLMAP*> regWt_PerSpecies;
            int nscoreImp = 0;

            // for each species:
            for (int specID = 0; specID < speciesIDNameMap.size(); specID++)
            {
                // unordered_map<int,double>* bestregWt_PerSpecies=new unordered_map<int,double>; no need to store weights
                // int specID=specIter->first;
                string spec = speciesIDNameMap[specID]; // specIter->second;
                SpeciesDataManager *sdm = speciesDataSet[specID];
                // map<string,int>& candidateregulators_Species=sdm->getRegulators();  //comment out
                FactorGraph *speciesGraph = sdm->getFactorGraph();
                VariableManager *vMgr = sdm->getVariableManager();
                vector<Variable *> &varSet = vMgr->getVariableSet();
                int targetID = vMgr->getVarID(targetgrpMembers[specID].c_str());
                int regID = vMgr->getVarID(tfgrpMembers[specID].c_str());
                //cout << "regi=" << regi << " regOGID=" << regOGIter << " targeti=" << targeti << " targetOGID=" << targetOGIter << " TFname=" << tfgrpMembers[specID] << " targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" << targetID << endl;
                if (targetID == -1 || regID == -1)
                {
                    //cout << "regulators or gene not present in " << spec << " TFname=" << tfgrpMembers[specID] << " targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" << targetID << endl;
                    continue;
                }
                /*
                //GeneMap* geneMap_Tgts=specIter->second;
                GeneMap* geneMap_Tgts=targetgrpMembers[spec]; //GeneMap* geneMap_Tgts=targetgrpMembers[specIter->first];
                map<string,map<string,STRINTMAP*>*>& speciesTargetSet=geneMap_Tgts->getGeneSet();
                GeneMap* geneMap_TFs=tfgrpMembers[spec]; //GeneMap* geneMap_TFs=tfgrpMembers[specIter->first];
                map<string,map<string,STRINTMAP*>*>& speciesTFSet=geneMap_TFs->getGeneSet();
                 */
                // If the species has multiple genes in it's list, we consider the member which will give the max improvement, that is
                // highest data likelihood. Similarly, we will extract the highest likelihood if there are multiple regulators
                double maxScore = -999999;
                double maxScoreImprovement = -999999;
                // Best here makes sense only if there are multiple TFs and targets
                int bestTarget = -1;
                int bestTF = -1;
                // speciesTargetSet.size()=1 and speciesTFSet.size()=1
                // for(map<string,map<string,STRINTMAP*>*>::iterator vIter=speciesTargetSet.begin();vIter!=speciesTargetSet.end();vIter++)
                //{
                // int targetID=vMgr->getVarID(vIter->first.c_str());
                Variable *target = varSet[targetID];
                SlimFactor *sFactor = speciesGraph->getFactorAt(targetID);
                // If the edge already exists in the MB of sFactor continue
                if (sFactor->mergedMB.find(regID) != sFactor->mergedMB.end())
                {
                    // cout <<"Edge exists in the MB of sFactor" << endl;
                    continue;
                }
                if (sFactor->mergedMB.size() >= currK)
                {
                    continue;
                }
                Variable *reg = varSet[regID];

                // Otherwise get the score delta of adding this reguluator in sFactor's MB.
                double scoreImprovement = 0;
                double newScore = 0;
                double newTargetScore = 0;
                getNewPLLScore(specID, reg, target, newScore, scoreImprovement, targetOGIter);
                //cout << "specID=" << specID << " regi=" << regi << " regOGID=" << regOGIter << " targeti=" << targeti << " targetOGID=" << targetOGIter << " TFname=" << tfgrpMembers[specID] << " targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" << targetID << " scoreImprovement=" << scoreImprovement << " newScore=" << newScore << endl;
                if (scoreImprovement > 0) //&& (scoreImprovement>maxScoreImprovement))
                {
                    bestTarget = targetID; // variable index
                    bestTF = regID;        // variable index
                    maxScoreImprovement = scoreImprovement;
                    maxScore = newScore;
                }
                else
                {
                    // regwt.clear();
                    continue;
                }
                // At this stage we are done with this species, and if maxDLL >0 we proceed with updating the information for this species
                scoreImprovement_PerSpecies[specID] = maxScoreImprovement;
                score_PerSpecies[specID] = maxScore;
                target_PerSpecies[specID] = bestTarget;
                tf_PerSpecies[specID] = bestTF;
                nscoreImp++;
                /*unordered_map<int,double>* regwtforspecies=new unordered_map<int,double>;  //copy of wts/regwt
                for(auto wIter=regwt.begin();wIter!=regwt.end();wIter++)
                {
                    (*regwtforspecies)[wIter->first]=wIter->second;
                }
                //regWt_PerSpecies[specID]=regwtforspecies;
                regwt.clear();*/
                // delete reg,sFactor,target,sdm,speciesGraph,vMgr;
                // varSet.clear();
            } // CVNvariant: species dataset end
            if (nscoreImp == 0)
            {
                continue;
                // cout <<"Skipping no hit regulator " <<  regOGIter << endl;
            }
            // Now we wish to see how good it would be add these edges in different conditions for all cell types:
            double bestImprovement = 0; // adding edge or not, which improvement is higher
            int csetid = -1;
            // cout <<"Compute for each condsetMap_Tree" << endl;
            for (int setIter = 0; setIter < condsetMap_Tree.size(); setIter++)
            {
                vector<int> &cset = condsetMap_Tree[setIter];
                vector<int> speciesEdgeStat(cset.size(), 0);
                int valid = 1;
                // compute score improvement+prior for each condition:
                double netImprovement = 0;
                for (int i = 0; i < cset.size(); i++)
                {
                    if (cset[i] == 1)
                    {
                        speciesEdgeStat[i] = 1;
                        if (scoreImprovement_PerSpecies[i] > 0)
                        {
                            netImprovement += scoreImprovement_PerSpecies[i];
                        }
                        else
                        {
                            valid = 0;
                            break;
                        }
                    }
                    else
                    {
                        speciesEdgeStat[i] = 0;
                    }
                }
                if (valid == 0)
                {
                    // This configuration is not valid
                    continue;
                }
                // compute the prior
                double ePrior = log(speciesData->getEdgeStatusProb(speciesEdgeStat));
                netImprovement = netImprovement + ePrior - oldpriorScore;

                if (netImprovement > bestImprovement)
                {
                    bestImprovement = netImprovement;
                    csetid = setIter;
                }
                // cout << "condition"<<setIter <<": ePrior=" << ePrior << " oldpriorScore=" << oldpriorScore << " netImprovementwPrior=" << netImprovement << endl;
                speciesEdgeStat.clear();
                // cset.clear();
            }
            // cout << "Best condition: "<<csetid << " bestImprovement=" << bestImprovement << endl;
            if (csetid == -1) // add nothing
            {
                continue;
            }
            if (bestImprovement > bestScoreImprovement_TF) // ePrior is added
            {
                // best condition:
                vector<int> &cset = condsetMap_Tree[csetid];
                bestscore_PerSpecies.clear();
                bestscoreImprovement_PerSpecies.clear();
                besttarget_PerSpecies.clear();
                besttf_PerSpecies.clear();
                for (int i = 0; i < cset.size(); i++)
                {
                    if (cset[i] == 0)
                    {
                        continue;
                    }
                    bestscore_PerSpecies[i] = score_PerSpecies[i];
                    bestscoreImprovement_PerSpecies[i] = scoreImprovement_PerSpecies[i];
                    besttarget_PerSpecies[i] = target_PerSpecies[i];
                    besttf_PerSpecies[i] = tf_PerSpecies[i];
                    // bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
                    // bestregWt_PerSpecies=regwtforspecies;
                }
                bestcsetid = csetid;
                bestScoreImprovement_TF = bestImprovement;
                bestregi = regi; // inputRegulatorOGs vector index
                //cout << "regi="<<regi << " regOGID=" << regOGIter << " targeti=" << targeti << " targetOGID=" << targetOGIter <<" best condition: "<<bestcsetid << " bestImprovement=" << bestImprovement << endl;
                
            }
            score_PerSpecies.clear();
            scoreImprovement_PerSpecies.clear();
            target_PerSpecies.clear();
            tf_PerSpecies.clear();
            // tfgrpMembers.clear();
            // delete tfogrp;
        } // all inputRegulatorOGs(regulators) end

        // add the best regulator for this target:
        if (bestcsetid == -1) // edge status cset[specID]
        {
            continue;
        }
        vector<int> &cset = condsetMap_Tree[bestcsetid];
        //cout << "best condition bestcsetid=" << bestcsetid << endl;
        for (int i = 0; i < cset.size(); i++)
        {
            if (cset[i] == 0)
            {
                continue;
            }
            MetaMove *move = new MetaMove;
            // int tfid=besttf_PerSpecies[i];
            // int tgtid=besttarget_PerSpecies[i];
            move->setSrcVertex(besttf_PerSpecies[i]); // TF variable id
            move->setTFID(bestregi);                  // inputRegulatorOGs vector index
            move->setConditionSetInd(i);
            move->setTargetVertex(besttarget_PerSpecies[i]); // target variable id
            move->setTargetID(targeti);                      // inputOGList vector index
            move->setTargetMBScore(bestscore_PerSpecies[i]);
            move->setScoreImprovement(bestscoreImprovement_PerSpecies[i]);
            moveSet.push_back(move);
            //cout << "Found edge for " << speciesIDNameMap[i] << " TFvarID=" << besttf_PerSpecies[i] << " TargetvarID=" << besttarget_PerSpecies[i] << " regOGidx=" << bestregi << " targetOGidx=" << targeti << " score improvement=" << bestscoreImprovement_PerSpecies[i] << endl;
            // bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
            // delete bestregWt_PerSpecies;
        } // species dataset end
        // double vm, rss;
        // process_mem_usage(vm, rss);
        // cout << "MetaLearner::collectMoves_Orthogroups VM: " << vm << " MB; RSS: " << rss << endl;
        /*for(auto wtIter=bestregWt_PerSpecies.begin();wtIter!=bestregWt_PerSpecies.end();wtIter++)
        {
            wtIter->second->clear();
            delete wtIter->second;
        }
        bestregWt_PerSpecies.clear();*/
        besttarget_PerSpecies.clear();
        bestscore_PerSpecies.clear();
        besttf_PerSpecies.clear();
        bestscoreImprovement_PerSpecies.clear();
        //cout <<"------------------------------------------------------------------------------------------------------" << endl;
        // targetgrpMembers.clear();
        // delete targetogrp;
        // bestScoreImprovementTF.clear();
    } // inputOGList ends
    // orthogroupSet.clear();
    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);
    // cout << "collectMoves_Orthogroups time: " << duration.count() << " microseconds" <<endl;
    return 0;
} // collectMoves_Orthogroups

int MetaLearner::collectMoves_Orthogroups_INDEP(int currK)
{
    cout << "MetaLearner::collectMoves_Orthogroups_INDEP" << endl;
    // auto start = high_resolution_clock::now();
    for (int i = 0; i < moveSet.size(); i++)
    {
        delete moveSet[i];
    }
    moveSet.clear();
    map<int, MappedOrthogroup *> &orthogroupSet = ogr->getMappedOrthogroups();

    // cout << "inputOGList.size() = " << inputOGList.size() << endl;
    // cout << "inputRegulatorOGs.size() = " << inputRegulatorOGs.size() << endl;
    // Now we will have a move for one orthogroup at a time
    for (int targeti = 0; targeti < inputOGList.size(); targeti++)
    {
        int targetOGIter = inputOGList[targeti];
        MappedOrthogroup *targetogrp = orthogroupSet[targetOGIter]; // oIter->second; four species
        // map<string,GeneMap*>& targetgrpMembers=targetogrp->getOrthoMembers();
        vector<string> &targetgrpMembers = targetogrp->getOrthoMembers();
        int n = speciesIDNameMap.size();

        // score the score of best regulator:
        vector<double> bestscore_PerSpecies(n, 0);
        vector<double> bestscoreImprovement_PerSpecies(n, 0);
        vector<int> besttarget_PerSpecies(n, -1);
        vector<int> besttf_PerSpecies(n, -1);
        // unordered_map<int,unordered_map<int,double>*> bestregWt_PerSpecies;
        int bestcsetid = -1;
        vector<double> bestScoreImprovementTF(n, 0);
        // The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        // term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.

        /*//only store current tf for each species!
        map<int,double> score_PerSpecies;
        map<int,double> scoreImprovement_PerSpecies;
        map<int,int> target_PerSpecies;
        map<int,int> tf_PerSpecies;
        map<int,INTDBLMAP*> regWt_PerSpecies;*/

        // we need to start from children to parents:loop from end of speciesIDNameMap to start of it!
        // for(map<int,string>::iterator specIter=speciesIDNameMap.begin();specIter!=speciesIDNameMap.end();specIter++)
        for (int specID = 0; specID < speciesIDNameMap.size(); specID++)
        // for(map<string,SpeciesDataManager*>::iterator specIter=speciesDataSet.begin();specIter!=speciesDataSet.end();specIter++)
        {
            // unordered_map<int,double>* bestregWt_PerSpecies=new unordered_map<int,double>; no need to store weights
            // int specID=specIter->first;
            string spec = speciesIDNameMap[specID]; // specIter->second;
            SpeciesDataManager *sdm = speciesDataSet[specID];
            // map<string,int>& candidateregulators_Species=sdm->getRegulators();  //comment out
            FactorGraph *speciesGraph = sdm->getFactorGraph();
            VariableManager *vMgr = sdm->getVariableManager();
            vector<Variable *> &varSet = vMgr->getVariableSet();
            int targetID = vMgr->getVarID(targetgrpMembers[specID].c_str());
            if (targetID == -1)
            {
                // cout << "gene not present in " << spec  << " targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" << targetID << endl;
                continue;
            }
            int bestregi = -1;
            // for each regulator
            // for(map<int,int>::iterator regOGIter=inputRegulatorOGs.begin();regOGIter!=inputRegulatorOGs.end();regOGIter++)
            for (int regi = 0; regi < inputRegulatorOGs.size(); regi++)
            {
                int regOGIter = inputRegulatorOGs[regi];
                if (regOGIter == targetOGIter) // if(regOGIter->first==oIter->first)
                {
                    continue;
                }
                MappedOrthogroup *tfogrp = orthogroupSet[regOGIter]; // orthogroupSet[regOGIter->first];
                // map<string,GeneMap*>& tfgrpMembers=tfogrp->getOrthoMembers();
                vector<string> &tfgrpMembers = tfogrp->getOrthoMembers();
                int regID = vMgr->getVarID(tfgrpMembers[specID].c_str());
                if (regID == -1)
                {
                    //cout << "regulators not present in " << spec << " TFname=" << tfgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" << targetID << endl;
                    continue;
                }
                // If the species has multiple genes in it's list, we consider the member which will give the max improvement, that is
                // highest data likelihood. Similarly, we will extract the highest likelihood if there are multiple regulators
                double maxScore = -999999;
                double maxScoreImprovement = -999999;
                // Best here makes sense only if there are multiple TFs and targets
                int bestTarget = -1;
                int bestTF = -1;
                // speciesTargetSet.size()=1 and speciesTFSet.size()=1
                // for(map<string,map<string,STRINTMAP*>*>::iterator vIter=speciesTargetSet.begin();vIter!=speciesTargetSet.end();vIter++)
                //{
                // int targetID=vMgr->getVarID(vIter->first.c_str());
                Variable *target = varSet[targetID];
                SlimFactor *sFactor = speciesGraph->getFactorAt(targetID);
                // If the edge already exists in the MB of sFactor continue
                if (sFactor->mergedMB.find(regID) != sFactor->mergedMB.end())
                {
                    // cout <<"Edge exists in the MB of sFactor" << endl;
                    continue;
                }
                if (sFactor->mergedMB.size() >= currK)
                {
                    continue;
                }
                Variable *reg = varSet[regID];

                // Otherwise get the score delta of adding this reguluator in sFactor's MB.
                double scoreImprovement = 0;
                double newScore = 0;
                double newTargetScore = 0;
                getNewPLLScore(specID, reg, target, newScore, scoreImprovement, targetOGIter);
                // cout << "specID=" << specID <<" regi=" << regi <<  " regOGID=" << regOGIter  <<" targeti=" << targeti << " targetOGID=" <<targetOGIter << " TFname=" << tfgrpMembers[specID] <<" targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" <<targetID << " scoreImprovement=" << scoreImprovement << " newScore=" << newScore << endl;
                // Now we wish to see how good it would be add these edges in different groups
                //  create my own edge status: add this one or not
                int regAdd = -1;
                if (scoreImprovement > 0) //&& (scoreImprovement>maxScoreImprovement))
                {
                    bestTarget = targetID; // variable index
                    bestTF = regID;        // variable index
                    maxScoreImprovement = scoreImprovement;
                    maxScore = newScore;
                    regAdd = 1;
                }
                else
                {
                    // regwt.clear();
                    continue;
                }
                /*unordered_map<int,double>* regwtforspecies=new unordered_map<int,double>;  //copy of wts/regwt
                for(auto wIter=regwt.begin();wIter!=regwt.end();wIter++)
                {
                    (*regwtforspecies)[wIter->first]=wIter->second;
                }
                //regWt_PerSpecies[specID]=regwtforspecies;
                regwt.clear();*/

                // add this one or not
                if (maxScoreImprovement > bestScoreImprovementTF[specID]) // ePrior is added
                {
                    bestscore_PerSpecies[specID] = maxScore;
                    bestscoreImprovement_PerSpecies[specID] = maxScoreImprovement;
                    besttarget_PerSpecies[specID] = bestTarget;
                    besttf_PerSpecies[specID] = bestTF;
                    // bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
                    // bestregWt_PerSpecies=regwtforspecies;
                    bestcsetid = regAdd;
                    bestScoreImprovementTF[specID] = maxScoreImprovement;
                    bestregi = regi;
                } /*else
                 {
                     regwtforspecies->clear();
                     delete regwtforspecies;
                 }*/
            }     // all inputRegulatorOGs(regulators) end

            // add the regulator: bestcsetid=1 add bestcsetid=0 continue
            if (bestcsetid != 1 || bestscoreImprovement_PerSpecies[specID] <= 0) // edge status cset[specID]
            {
                continue;
            }
            MetaMove *move = new MetaMove;
            // int tfid=besttf_PerSpecies[specID];
            // int tgtid=besttarget_PerSpecies[specID];
            move->setSrcVertex(besttf_PerSpecies[specID]);
            move->setTFID(bestregi);
            move->setConditionSetInd(specID);
            move->setTargetVertex(besttarget_PerSpecies[specID]);
            move->setTargetID(targeti);
            move->setTargetMBScore(bestscore_PerSpecies[specID]);
            move->setScoreImprovement(bestscoreImprovement_PerSpecies[specID]);
            moveSet.push_back(move);
            // bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
            // delete bestregWt_PerSpecies;

        } // species dataset end
        besttarget_PerSpecies.clear();
        bestscore_PerSpecies.clear();
        besttf_PerSpecies.clear();
        bestscoreImprovement_PerSpecies.clear();
        bestScoreImprovementTF.clear();
        /*for(auto wtIter=bestregWt_PerSpecies.begin();wtIter!=bestregWt_PerSpecies.end();wtIter++)
        {
            wtIter->second->clear();
            delete wtIter->second;
        }
        bestregWt_PerSpecies.clear();*/
    }
    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);
    // cout << "collectMoves_Orthogroups time: " << duration.count() << " microseconds" <<endl;
    return 0;
}

// currPrior is the species-specific prior, which we do not need
// priorScore is computed by not counted
// getNewPLLScore(specID,*cset,reg,target,newScore,scoreImprovement,oIter->first,regwt);
// shilu: removed INTINTMAP& conditionSet
// u is reg and v is target
// species-specific prior: sum_reg[log(p)]+sum_nonreg[log(1-p)]
// score=likelihood+species-specific prior:
int MetaLearner::getNewPLLScore(int cid, Variable *u, Variable *v, double &targetmbScore, double &scoreImprovement, int orthoGrpNo)
{
    SpeciesDataManager *sdm = speciesDataSet[cid];
    // PotentialManager* potMgr=sdm->getPotentialManager();
    // VSET& varSet=sdm->getVariableManager()->getVariableSet();
    FactorGraph *fg = sdm->getFactorGraph();
    // SlimFactor* sFactor=fg->getFactorAt(u->getID());  //regulator
    SlimFactor *dFactor = fg->getFactorAt(v->getID()); // target
    // map<int,double>* varNeighborhoodPrior=varNeighborhoodPrior_PerSpecies[speciesIDNameMap[cid]];
    vector<double> &varNeighborhoodPrior = varNeighborhoodPrior_PerSpecies[cid];
    vector<vector<double>> &edgePresenceProb = edgePresenceProb_PerSpecies[cid];
    // map<string,double>* edgePresenceProb=edgePresenceProb_PerSpecies[speciesIDNameMap[cid]];
    double currPrior = varNeighborhoodPrior[v->getID()]; // target
    bool toDel_d = true;
    
    /*if(dFactor->mergedMB.find(u->getID())!=dFactor->mergedMB.end())
    {
        toDel_d=false;
    }*/
    // already checked mergedMB before computing getNewPLLScore
    double plus = 0;
    double minus = 0;
    dFactor->mergedMB.insert(u->getID()); // Aug 23: dFactor->mergedMB[u->getID()]=0;
    int status = 0;
    for (auto mIter = dFactor->mergedMB.begin(); mIter != dFactor->mergedMB.end(); mIter++)
    {
        /*Variable* aVar=varSet[mIter->first];
        string regulatorKey(aVar->getName().c_str());
        regulatorKey.append("\t");
        regulatorKey.append(v->getName().c_str());*/
        double p = edgePresenceProb[*mIter][v->getID()]; //(*edgePresenceProb)[regulatorKey];
        // update by Shilu to avoid inf
        if (p == 0 or p == 1)
        {
            continue;
        }
        else
        {
            minus = minus + log(1 - p);
            plus = plus + log(p);
        }
        
    }
    double pll_d = getPLLScore_Condition_Tracetrick(cid, dFactor, status);
    //cout << "MetaLearner::getNewPLLScore cell=" << cid << " regID=" << u->getID() << " varID=" << v->getID()  << " currPrior="<< currPrior << " pll_d="<< pll_d << endl;
    if (status == -1)
    {
        scoreImprovement = -1;
        if (toDel_d)
        {
            auto dIter = dFactor->mergedMB.find(u->getID());
            dFactor->mergedMB.erase(dIter);
        }
        return 0;
    }
    currPrior = currPrior + plus - minus;
    //cout << "cell=" << cid << " regID=" << u->getID() << " varID=" << v->getID() << " oldsparsityPrior=" << varNeighborhoodPrior[v->getID()] << " sparsityPriorchanged=" << plus - minus << " currsparsityPrior=" << currPrior << " pll_d=" << pll_d << " newscore=" << pll_d + currPrior << " oldscore=" << dFactor->mbScore << endl;
    pll_d = pll_d + currPrior;
    /* comment by shilu
     double priorScore=0;
     double oldpriorScore=0;
     if(speciesData->getRoot()!=NULL)
     {
     priorScore=getEdgePrior(cid,u,v,orthoGrpNo);
     int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[cid].c_str());
     int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[cid].c_str());
     if(vogno!=-1)
     {
     char ogPairKey[256];
     sprintf(ogPairKey,"%d-%d",uogno,vogno);
     string ogPair(ogPairKey);
     oldpriorScore=ogpairPrior[ogPair];
     }
     }*/
    targetmbScore = pll_d;
    // scoreImprovement=0;
    // double dImpr=(targetmbScore+priorScore)-(dFactor->mbScore+oldpriorScore);
    // Don't include the prior. Just use the data likelihood improvement
    double dImpr = targetmbScore - dFactor->mbScore;
    //

    if (dImpr <= 0)
    {
        scoreImprovement = -1;
    }
    else
    {
        scoreImprovement = dImpr;
    }
    if (toDel_d)
    {
        auto dIter = dFactor->mergedMB.find(u->getID());
        dFactor->mergedMB.erase(dIter);
    }
    //cout << "getPLLScore_Condition_Tracetrick likelihood=" << pll_d << " dFactor->mbScore=" << dFactor->mbScore << " in " << speciesIDNameMap[cid] << " for regID=" << u->getID() << " varID=" << v->getID()  << " status=" << status << " toDel_d=" << toDel_d << endl;

    return 0;
}

// shilu: more efficient version
double
MetaLearner::getPLLScore_Condition_Tracetrick(int specID, SlimFactor *sFactor, int &status) // sFactor is target, unordered_map<int,double>& regWts
{
    PotentialManager *potMgr = speciesDataSet[specID]->getPotentialManager(); // grab managers
    double pll = potMgr->computePotentialMBCovMean(sFactor, status);
    return pll;
}

// prior probability: between regulator j and target k as a logistic function, remove motifweight
// tfID is variable ID for regulator, targetID is variable ID for target.
double
MetaLearner::getEdgePrior_PerSpecies(int tfID, int targetID, SpeciesDataManager *sdm)
{
    INTDBLMAP *regPriors = NULL;
    // double prior=1/(1+exp(-1*beta1));
    double motifweight = 0;
    // double chipweight=0;
    unordered_map<int, unordered_map<int, double> *> &motifNetwork = sdm->getMotifNetwork();
    if (motifNetwork.find(tfID) != motifNetwork.end())
    {
        unordered_map<int, double> *values = motifNetwork[tfID];
        if (values->find(targetID) != values->end())
        {
            motifweight = (*values)[targetID];
            //cout << "tfID=" << tfID << " targetID=" << targetID << " motifweight=" << motifweight << " prior=" << 1 / (1 + exp(-1 * (beta1 + motifweight * beta2))) << endl;
        }
    }
    double fwt = motifweight * beta2;
    double prior = 1 / (1 + exp(-1 * (beta1 + fwt)));
    if (prior < 1e-6)
    {
        prior = 1e-6;
    }
    if (prior == 1)
    {
        prior = 1 - 1e-6;
    }
    //cout << "tfID=" << tfID << " targetID=" << targetID << " motifweight=" << motifweight << " beta1=" << beta1 << " beta2=" << beta2<< " prior=" << prior<< endl;
    return prior;
}

double
MetaLearner::makeMoves(int &successMove)
{
    /*map<int,INTINTMAP*> affectedVariables;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        int cind=speciesNameIDMap[gIter->first];
        INTINTMAP* csVars=new INTINTMAP;
        affectedVariables[cind]=csVars;
    }*/
    // int successMove=0;
    double netScoreDelta = 0;
    // int net1Move=0;
    // int net2Move=0;
    for (int m = 0; m < moveSet.size(); m++)
    {
        MetaMove *move = moveSet[m];
        if (attemptMove(move) == 0) // if(attemptMove(move,affectedVariables)==0)
        {
            successMove++;
            netScoreDelta = netScoreDelta + move->getScoreImprovement();
        }
        else
        {
            cout << "Move not valid " << endl;
        }
    }
    // double vm, rss;
    // process_mem_usage(vm, rss);
    // cout << "MetaLearner::makeMoves VM: " << vm << " MB; RSS: " << rss << endl;
    /*for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
    {
        cIter->second->clear();
        delete cIter->second;
    }
    affectedVariables.clear();*/
    // cout <<"Total successful moves " << successMove << " out of total " << moveSet.size() << " with net score improvement " << netScoreDelta << endl;
    return netScoreDelta;
}

// affectedVars not needed, each target is added only once
int MetaLearner::attemptMove(MetaMove *move)
{
    int specID = move->getConditionSetInd();
    SpeciesDataManager *sdm = speciesDataSet[specID];
    // vector<Variable*>& varSet=sdm->getVariableManager()->getVariableSet();
    // Variable* u=varSet[move->getSrcVertex()];
    // Variable* v=varSet[move->getTargetVertex()];
    int regi = move->getTFID();        // inputRegulatorOGs vector index
    int targeti = move->getTargetID(); // inputOGList vector index
    FactorGraph *csGraph = sdm->getFactorGraph();
    SlimFactor *dFactor = csGraph->getFactorAt(move->getTargetVertex());
    dFactor->mergedMB.insert(move->getSrcVertex()); // Aug 23: dFactor->mergedMB[move->getSrcVertex()]=0;
    dFactor->mbScore = move->getTargetMBScore();

    vector<int> *newEdgeStatus; // STRINTMAP* newEdgeStatus=NULL;
    char ogpair[20];
    sprintf(ogpair, "%d-%d", regi, targeti);
    string ogpairKey(ogpair);
    // cout <<"Attempting move " << ogpairKey << " in species " << speciesIDNameMap[specID] << endl;
    if (affectedOGPairs.find(ogpairKey) == affectedOGPairs.end())
    {
        newEdgeStatus = new vector<int>(speciesIDNameMap.size(), 0); // new STRINTMAP;
        affectedOGPairs[ogpairKey] = newEdgeStatus;
    }
    else
    {
        newEdgeStatus = affectedOGPairs[ogpairKey];
    }
    // string& species=speciesIDNameMap[specID];
    (*newEdgeStatus)[specID] = 1;
    return 0;
}

int MetaLearner::dumpAllGraphs(int currK)
{
    cout << "MetaLearner::dumpAllGraphs" << endl;
    char aFName[1024];
    for (int i = 0; i < speciesDataSet.size(); i++)
    {
        SpeciesDataManager *sdm = speciesDataSet[i];
        const char *dirname = sdm->getOutputLoc();
        sprintf(aFName, "%s/var_mb_pw_k%d.txt", dirname, currK);
        ofstream oFile(aFName);
        FactorGraph *fg = sdm->getFactorGraph();
        vector<Variable *> &varSet = sdm->getVariableManager()->getVariableSet();
        vector<SlimFactor *> &factorSet = fg->getAllFactors();
        PotentialManager *potMgr = sdm->getPotentialManager();
        for (int aIter = 0; aIter < factorSet.size(); aIter++)
        {
            SlimFactor *sFactor = factorSet[aIter];
            potMgr->dumpVarMB_PairwiseFormat(sFactor, oFile, varSet);
        }
        oFile.close();
    }
    return 0;
}

int MetaLearner::showModelParameters()
{
    char aFName[1024];
    for (int i = 0; i < speciesDataSet.size(); i++)
    {
        SpeciesDataManager *sdm = speciesDataSet[i];
        const char *dirname = sdm->getOutputLoc();
        sprintf(aFName, "%s/modelparams.txt", dirname);
        ofstream oFile(aFName);
        FactorGraph *fg = sdm->getFactorGraph();
        PotentialManager *potMgr = sdm->getPotentialManager();
        vector<SlimFactor *> &slimFactorSet = fg->getAllFactors();
        vector<Variable *> &varSet = sdm->getVariableManager()->getVariableSet();
        for (auto sIter = slimFactorSet.begin(); sIter != slimFactorSet.end(); sIter++)
        {
            SlimFactor *sFactor = *sIter;
            double mbcondvar = 0;
            double mbbias = 0;
            unordered_map<int, double> mbwt;
            potMgr->computePotentialMBCovMean(sFactor, mbcondvar, mbbias, mbwt);
            Variable *var = varSet[sFactor->fId];
            oFile << "Var=" << var->getName() << "\tWt=-1" << "\tCondVar=" << mbcondvar << "\tCondBias=" << mbbias << "\tCondWt=";
            for (auto dIter = mbwt.begin(); dIter != mbwt.end(); dIter++)
            {
                if (dIter != mbwt.begin())
                {
                    oFile << ",";
                }
                Variable *mbVar = varSet[dIter->first];
                oFile << mbVar->getName() << "=" << dIter->second;
            }
            oFile << endl;
            mbwt.clear();
        }
        oFile.close();
    }
    return 0;
}

int MetaLearner::generateData(int sampleCnt, int burnin)
{
    map<int, ofstream *> newdataFiles;
    gsl_rng *rndgen = gsl_rng_alloc(gsl_rng_default);
    map<int, map<int, Potential *> *> potSet;

    for (int i = 0; i < speciesDataSet.size(); i++)
    {
        // int gid=speciesNameIDMap[gIter->first];
        FactorGraph *fg = speciesDataSet[i]->getFactorGraph();
        VariableManager *varMgr = speciesDataSet[i]->getVariableManager();
        vector<Variable *> &varSet = varMgr->getVariableSet();
        map<int, Potential *> *pSet = new map<int, Potential *>;
        potSet[i] = pSet;
        PotentialManager *potMgr = speciesDataSet[i]->getPotentialManager();
        for (map<string, int>::iterator vIter = subgraphVarSet.begin(); vIter != subgraphVarSet.end(); vIter++)
        {
            int vId = vIter->second;
            SlimFactor *sFactor = fg->getFactorAt(vId);
            Potential *sPot = new Potential;
            (*pSet)[vId] = sPot;
            sPot->setAssocVariable(varSet[sFactor->fId], Potential::FACTOR);
            for (auto mIter = sFactor->mergedMB.begin(); mIter != sFactor->mergedMB.end(); mIter++)
            {
                Variable *aVar = varSet[*mIter];
                sPot->setAssocVariable(aVar, Potential::MARKOV_BNKT);
            }
            sPot->potZeroInit();
            potMgr->populatePotential(sPot);
            sPot->initMBCovMean();
        }
    }

    for (int i = 0; i < speciesDataSet.size(); i++)
    {
        SpeciesDataManager *sdm = speciesDataSet[i];
        char fileName[1024];
        sprintf(fileName, "%s/newsamples.txt", sdm->getOutputLoc());
        ofstream *oFile = new ofstream(fileName);
        newdataFiles[i] = oFile;
        // Write the model file
        sprintf(fileName, "%s/newmodel.txt", sdm->getOutputLoc());
        ofstream mFile(fileName);
        mFile << "NodeCnt\t" << subgraphVarSet.size() << endl;
        mFile << "ContinuousNodes";
        int nodeId = 0;
        for (map<string, int>::iterator vIter = subgraphVarSet.begin(); vIter != subgraphVarSet.end(); vIter++)
        {
            mFile << "\t" << nodeId;
            nodeId++;
        }
        mFile << endl;
        mFile << "DiscreteNodes" << endl;
        nodeId = 0;
        for (map<string, int>::iterator vIter = subgraphVarSet.begin(); vIter != subgraphVarSet.end(); vIter++)
        {
            mFile << "NodeName=" << vIter->first << "\tNodeID=" << nodeId << "\tParents="
                  << "\tChildren="
                  << "\tValues=0,1,2" << endl;
            nodeId++;
        }
        mFile.close();
    }

    for (map<int, map<int, Potential *> *>::iterator psIter = potSet.begin(); psIter != potSet.end(); psIter++)
    {
        int did = psIter->first;
        INTDBLMAP initialSample;
        INTDBLMAP saveSample;
        int iter = 0;
        map<int, Potential *> *pSet = psIter->second;
        generateInitSample(initialSample, rndgen, pSet);

        while (iter < (sampleCnt + burnin))
        {
            int selectFirst = 0;
            // We just need to align localmodelid with the outputfile to which the data is generated
            if (iter >= burnin)
            {
                if (iter == burnin)
                {
                    cout << "Starting to write" << endl;
                }
                writeEvidenceTo(newdataFiles[did], initialSample);
            }
            for (map<string, int>::iterator vIter = subgraphVarSet.begin(); vIter != subgraphVarSet.end(); vIter++)
            {
                saveSample[vIter->second] = initialSample[vIter->second];
                FactorGraph *fg = speciesDataSet[psIter->first]->getFactorGraph();
                SlimFactor *sFactor = fg->getFactorAt(vIter->second);
                Potential *sPot = (*pSet)[vIter->second];
                double sample_s = sPot->generateSample(initialSample, sFactor->fId, rndgen);
                int siter = 0;
                while ((sample_s < -160 || sample_s > 160) && (siter < 50))
                {
                    sample_s = sPot->generateSample(initialSample, sFactor->fId, rndgen);
                    siter++;
                }
                if (sample_s < -160)
                {
                    sample_s = -160;
                }
                else if (sample_s > 160)
                {
                    sample_s = 160;
                }
                if ((isnan(sample_s)) || (isinf(sample_s)))
                {
                    cout << "Found Nan/Inf for " << sFactor->fId << "=" << sample_s << " Neighbors: ";
                    for (auto mbIter = sFactor->mergedMB.begin(); mbIter != sFactor->mergedMB.end(); mbIter++)
                    {
                        cout << " " << *mbIter << "=" << initialSample[*mbIter];
                    }
                    cout << endl;
                }
                initialSample[sFactor->fId] = sample_s;
            }
            iter++;
        }
    }
    gsl_rng_free(rndgen);
    cout << "Bottom clips " << endl;
    for (map<string, int>::iterator sIter = bottomClip.begin(); sIter != bottomClip.end(); sIter++)
    {
        cout << sIter->first.c_str() << " " << sIter->second << endl;
    }
    cout << "Top clips " << endl;
    for (map<string, int>::iterator sIter = topClip.begin(); sIter != topClip.end(); sIter++)
    {
        cout << sIter->first.c_str() << " " << sIter->second << endl;
    }
    return 0;
}

int MetaLearner::generateInitSample(INTDBLMAP &initialSample, gsl_rng *rndgen, map<int, Potential *> *pSet)
{
    double globalMean = 0;
    globalVar = 0;
    for (map<string, int>::iterator vIter = subgraphVarSet.begin(); vIter != subgraphVarSet.end(); vIter++)
    {
        int vId = vIter->second;
        Potential *sPot = (*pSet)[vId];
        cout << vId << "\t" << sPot->getCondBias() << "\t" << sPot->getCondVariance() << endl;
        globalMean = globalMean + sPot->getCondBias();
        globalVar = globalVar + sPot->getCondVariance();
    }
    double norm = (double)(subgraphVarSet.size() * speciesDataSet.size());
    globalMean = globalMean / norm;
    globalVar = globalVar / norm;
    for (map<string, int>::iterator vIter = subgraphVarSet.begin(); vIter != subgraphVarSet.end(); vIter++)
    {
        double sampleval = gsl_ran_gaussian(rndgen, sqrt(globalVar));
        sampleval = sampleval + globalMean;
        initialSample[vIter->second] = sampleval;
    }
    return 0;
}

int MetaLearner::writeEvidenceTo(ofstream *oFile, INTDBLMAP &sample)
{
    int nodeId = 0;
    for (map<string, int>::iterator dIter = subgraphVarSet.begin(); dIter != subgraphVarSet.end(); dIter++)
    {
        if (dIter != subgraphVarSet.begin())
        {
            (*oFile) << "\t";
        }
        double val = sample[dIter->second];
        if (val < -160)
        {
            val = -160;
            if (bottomClip.find(dIter->first) == bottomClip.end())
            {
                bottomClip[dIter->first] = 1;
            }
            else
            {
                bottomClip[dIter->first] = bottomClip[dIter->first] + 1;
            }
        }
        else if (val > 160)
        {
            val = 160;
            if (topClip.find(dIter->first) == topClip.end())
            {
                topClip[dIter->first] = 1;
            }
            else
            {
                topClip[dIter->first] = topClip[dIter->first] + 1;
            }
        }
        (*oFile) << nodeId << "=[" << exp(val) << "]";
        nodeId++;
    }
    (*oFile) << endl;
    return 0;
}
