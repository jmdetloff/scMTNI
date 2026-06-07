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
#include <math.h>
#include <string>
#include <vector>
#include <utility>
#include <chrono>
#include <unistd.h>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"
#include "MetaMove.H"
#include "Utils.H"
#include "SpeciesDistance.H"
#include "SpeciesDataManager.H"
#include "MetaLearner.H"

using namespace std::chrono;
using namespace std;

MetaLearner::MetaLearner()
{
    convThreshold = 1e-4;
    beta1 = -0.9;
    beta2 = 4.0;
    INDEP = false;
    initGlobalScore = 0;
}

MetaLearner::~MetaLearner()
{
    variableList.clear();
    targetList.clear();
    regulatorList.clear();
}

int MetaLearner::setInputFName(const char *aFName)
{
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

int MetaLearner::setConvergenceThreshold(double aVal)
{
    convThreshold = aVal;
    return 0;
}

int MetaLearner::setRegulatorList(const char *aFName)
{
    ifstream inFile(aFName);
    string buffer;
    while (inFile.good()) {
        getline(inFile, buffer);
        if (buffer.length() <= 0) {
            continue;
        }
        int varIndex = addVariable(buffer);
        regulatorList.push_back(varIndex);
    }
    inFile.close();
    return 0;
}

int MetaLearner::setTargetList(const char *aFName)
{
    ifstream inFile(aFName);
    string buffer;
    while (inFile.good()) {
        getline(inFile, buffer);
        if (buffer.length() <= 0) {
            continue;
        }
        int varIndex = addVariable(buffer);
        targetList.push_back(varIndex);
    }
    inFile.close();
    return 0;
}

int MetaLearner::addVariable(string varName)
{
    vector<string>::iterator iter = find(variableList.begin(), variableList.end(), varName);

    // If the variable already exists, just return its index.
    if (iter != variableList.end()) {
        return distance(variableList.begin(), iter);
    }

    // Otherwise append it.
    variableList.push_back(varName);
    return variableList.size() - 1;
}

int MetaLearner::setSpeciesDistances(SpeciesDistance *aPtr)
{
    speciesData = aPtr;
    speciesData->setSpeciesNameIDMap(speciesNameIDMap);
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
    cout << "MetaLearner::init() read:" << inputFName << endl;

    ifstream inFile(inputFName);
    string buffer;

    while (inFile.good()) {
        getline(inFile, buffer);

        if (buffer.empty() || buffer.find("#") == 0) {
            continue;
        }

        // strip trailing newlines
        buffer.erase(buffer.find_last_not_of(" \n\r\t") + 1);

        vector<string> strs = Utils::split(buffer, '\t');
        initSpeciesData(strs[0], strs[1], strs[2], strs[3]);
    }
    inFile.close();
    cout << "-------------------------------------------------------------" << endl;
    return 0;
}

void MetaLearner::initSpeciesData(string speciesName, string tableFileName, string outputLoc, string motifNetwork)
{
    PotentialManager *potMgr = new PotentialManager;
    potMgr->loadEvidenceFromTable(tableFileName, variableList);

    SpeciesDataManager *spMgr = new SpeciesDataManager;
    spMgr->setPotentialManager(potMgr);
    spMgr->setOutputLoc(outputLoc.c_str());
    spMgr->setMotifNetwork(motifNetwork.c_str());
    speciesDataSet.push_back(spMgr);

    int datasetID = speciesDataSet.size() - 1;
    speciesNameIDMap[speciesName] = datasetID;

    cout << datasetID << "=" << speciesName << " motifNetwork=" << motifNetwork << endl;

    unordered_map<int, double> varNeighborhoodPrior;
    unordered_map<int, unordered_map<int, double>> edgePresenceProb;

    VariableManager *varMgr = potMgr->getVariableManager();
    vector<Variable*>& variableSet = varMgr->getVariableSet();
    for (int i = 0; i < variableSet.size(); i++) {
        Variable *target = variableSet[i];
        SlimFactor *sFactor = spMgr->getFactor(target->getID());
        double pll = potMgr->computeUnivariateLL(sFactor->fId);
        double priorScore = precomputePerSpeciesPrior(datasetID, target, spMgr, edgePresenceProb);
        varNeighborhoodPrior[sFactor->fId] = priorScore;
        sFactor->mbScore = pll + priorScore;
        initGlobalScore = initGlobalScore + pll + priorScore;
    }

    varNeighborhoodPrior_PerSpecies.push_back(varNeighborhoodPrior);
    edgePresenceProb_PerSpecies.push_back(edgePresenceProb);
}

void MetaLearner::start()
{
    auto start = high_resolution_clock::now();
    cout << "MetaLearner::start" << endl;
    cout << "Regulator Count: " << regulatorList.size() << " Target Count: " << targetList.size() << endl;

    double currGlobalScore = initGlobalScore;

    if (!INDEP) {
        precomputeEmptyGraphPrior();
        speciesData->createConditionSets();
    }

    int iter = 0;
    bool notConverged = true;

    // hardcode 100 iterations for now
    while (notConverged && iter < 100) {

        // collect the candidate edges
        vector<MetaMove> moves = INDEP ? collectMoves_INDEP(maxFactorSizeApprox) : collectMoves(maxFactorSizeApprox);

        double diff = makeMoves(moves);
        double priorChange = INDEP ? 0 : getPriorDelta();
        double newScore = currGlobalScore + diff;

        if (diff <= convThreshold) {
            notConverged = false;
        }

        currGlobalScore = newScore;

        double vm, rss;
        process_mem_usage(vm, rss);
        cout << "ITERATION " << iter << " newScore=" << newScore << " diffscore=" << diff << " priorChange=" << priorChange << " successfulMoves=" << moves.size() << endl;

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

double MetaLearner::getScore()
{
    double gScore = 0;

    for (int i = 0; i < speciesDataSet.size(); i++) {
        SpeciesDataManager *speciesDataManager = speciesDataSet[i];
        VariableManager *varMgr = speciesDataManager->getVariableManager();
        vector<Variable*>& variableSet = varMgr->getVariableSet();
        for (int j = 0; j < variableSet.size(); j++) {
            Variable *var = variableSet[j];
            SlimFactor *sFactor = speciesDataManager->getFactor(var->getID());
            gScore = gScore + sFactor->mbScore;
        }
    }
    return gScore;
}

double MetaLearner::getPriorDelta()
{
    double oldStructPrior = 0;
    double newStructPrior = 0;
    // Need to consider the old contribution of the edges, delete that from the overall prior and add the new contribution
    for (auto edgeIter = affectedVarPairs.begin(); edgeIter != affectedVarPairs.end(); edgeIter++)
    {
        vector<string> keyid = Utils::split(edgeIter->first, '-');
        int regulatorID = stoi(keyid[0]);
        int targetID = stoi(keyid[1]);
        double aval = speciesData->getEdgeStatusProb(*(edgeIter->second));
        double edgePrior = log(aval);
        double oldEdgePrior = edgePriors[targetID][regulatorID];
        oldStructPrior += oldEdgePrior;
        newStructPrior += edgePrior;
        edgePriors[targetID][regulatorID] = edgePrior;
    }
    for (auto edgeIter = affectedVarPairs.begin(); edgeIter != affectedVarPairs.end(); edgeIter++)
    {
        delete edgeIter->second;
    }
    affectedVarPairs.clear();
    double priorDelta = newStructPrior - oldStructPrior;
    return priorDelta;
}

void MetaLearner::precomputeEmptyGraphPrior()
{
    vector<int> edgeStatus(speciesNameIDMap.size(), 0);
    double prior = speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior = log(prior);
    for (int i = 0; i < targetList.size(); i++) {
        int targetID = targetList[i];
        unordered_map<int, double> targetPriors;
        for (int j = 0; j < regulatorList.size(); j++) {
            int regulatorID = regulatorList[j];
            if (regulatorID == targetID) {
                targetPriors[regulatorID] = 0;
            } else {
                targetPriors[regulatorID] = logPrior;
            }    
        }
        edgePriors[targetID] = targetPriors;
    }
}

double MetaLearner::precomputePerSpeciesPrior(int specID, Variable *target, SpeciesDataManager *sdm, unordered_map<int, unordered_map<int, double>> &edgePresenceProb)
{
    VariableManager *vMgr = sdm->getVariableManager();
    int targetID = target->getID();
    double neighborhoodPrior = 0;
    for (int i = 0; i < regulatorList.size(); i++) {
        int regID = regulatorList[i];
        Variable *regulatorVar = vMgr->getVariable(regID);

        // Confirm that the regulator exists in this species' dataset, and that it isnt the target.
        if (regulatorVar == nullptr || regID == targetID) {
            continue;
        }

        double initPrior = getEdgePrior(regID, targetID, sdm);
        edgePresenceProb[regID][targetID] = initPrior;
        neighborhoodPrior += log(1 - initPrior);
    }
    return neighborhoodPrior;
}

vector<MetaMove> MetaLearner::collectMoves(int currK)
{
    vector<MetaMove> moveSet;

    int numSpecies = speciesNameIDMap.size();
    vector<vector<int>>& conditionSets = speciesData->getConditionSets();

    for (int i = 0; i < targetList.size(); i++) {
        int targetID = targetList[i];

        unordered_map<int, double> oldPriors = edgePriors[targetID];

        // score the score of best regulator:
        vector<double> bestscore_PerSpecies(numSpecies, 0);
        vector<double> bestscoreImprovement_PerSpecies(numSpecies, 0);
        vector<int> besttarget_PerSpecies(numSpecies, -1);
        vector<int> besttf_PerSpecies(numSpecies, -1);
        int bestcsetid = -1;
        int bestRegulatorID = -1;
        double bestScoreImprovement_TF = 0;

        // The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        // term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.

        for (int j = 0; j < regulatorList.size(); j++) {
            int regulatorID = regulatorList[j];

            if (targetID == regulatorID) {
                continue;
            }

            double oldPrior = oldPriors[regulatorID];

            vector<double> score_PerSpecies(numSpecies, 0);
            vector<double> scoreImprovement_PerSpecies(numSpecies, 0);
            int nscoreImp = 0;

            // for each species:
            for (int specID = 0; specID < speciesNameIDMap.size(); specID++)
            {
                SpeciesDataManager *sdm = speciesDataSet[specID];
                VariableManager *vMgr = sdm->getVariableManager();
                Variable *target = vMgr->getVariable(targetID);
                Variable *regulator = vMgr->getVariable(regulatorID);

                // Ensure that the target and regulator are both present in this species' dataset.
                if (target == nullptr || regulator == nullptr) {
                    continue;
                }

                SlimFactor *sFactor = sdm->getFactor(targetID);

                // If the edge already exists in the MB of sFactor continue
                if (sFactor->mergedMB.find(regulatorID) != sFactor->mergedMB.end()) {
                    continue;
                }

                // If the target already has the max num edges, continue.
                if (sFactor->mergedMB.size() >= currK) {
                    continue;
                }

                // Otherwise get the score delta of adding this regulator in sFactor's MB.
                double scoreImprovement = 0;
                double newScore = 0;
                getNewPLLScore(specID, regulator, target, newScore, scoreImprovement);

                // If adding the edge wouldn't improve score, we dont need to consider it.
                if (scoreImprovement <= 0) {
                    continue;
                }

                scoreImprovement_PerSpecies[specID] = scoreImprovement;
                score_PerSpecies[specID] = newScore;
                nscoreImp++;
            }

            if (nscoreImp == 0) {
                continue;
            }

            // Now we wish to see how good it would be add these edges in different conditions for all cell types:
            double bestImprovement = 0;
            int csetid = -1;

            for (int setIter = 0; setIter < conditionSets.size(); setIter++)
            {
                vector<int> &cset = conditionSets[setIter];
                vector<int> speciesEdgeStat(cset.size(), 0);
                int valid = 1;

                // compute score improvement+prior for each condition:
                double netImprovement = 0;
                for (int i = 0; i < cset.size(); i++) {
                    if (cset[i] == 1) {
                        speciesEdgeStat[i] = 1;
                        if (scoreImprovement_PerSpecies[i] <= 0) {
                            valid = 0;
                            break;
                        }
                        netImprovement += scoreImprovement_PerSpecies[i];
                    } else {
                        speciesEdgeStat[i] = 0;
                    }
                }

                // Confirm configuration is valid
                if (valid == 0) {
                    continue;
                }

                // compute the prior
                double ePrior = log(speciesData->getEdgeStatusProb(speciesEdgeStat));
                netImprovement = netImprovement + ePrior - oldPrior;

                if (netImprovement > bestImprovement) {
                    bestImprovement = netImprovement;
                    csetid = setIter;
                }
            }

            // If we couldn't find a score improving configuration of edges, continue.
            if (csetid == -1) {
                continue;
            }

            if (bestImprovement > bestScoreImprovement_TF) {
                vector<int> &cset = conditionSets[csetid];
                bestscore_PerSpecies.clear();
                bestscoreImprovement_PerSpecies.clear();
                besttarget_PerSpecies.clear();
                besttf_PerSpecies.clear();
                for (int i = 0; i < cset.size(); i++) {
                    if (cset[i] == 0) {
                        continue;
                    }
                    bestscore_PerSpecies[i] = score_PerSpecies[i];
                    bestscoreImprovement_PerSpecies[i] = scoreImprovement_PerSpecies[i];
                    besttarget_PerSpecies[i] = targetID;
                    besttf_PerSpecies[i] = regulatorID;
                }
                bestcsetid = csetid;
                bestScoreImprovement_TF = bestImprovement;
                bestRegulatorID = regulatorID;
            }
        }

        // If we didnt find a score improving regulator for this target, continue.
        if (bestcsetid == -1) {
            continue;
        }

        vector<int> &cset = conditionSets[bestcsetid];

        for (int i = 0; i < cset.size(); i++) {
            if (cset[i] == 0) {
                continue;
            }
            MetaMove move;
            move.setSrcVertex(besttf_PerSpecies[i]);
            move.setTFID(bestRegulatorID);
            move.setConditionSetInd(i);
            move.setTargetVertex(besttarget_PerSpecies[i]);
            move.setTargetID(targetID);
            move.setTargetMBScore(bestscore_PerSpecies[i]);
            move.setScoreImprovement(bestscoreImprovement_PerSpecies[i]);
            moveSet.push_back(move);
        }
    }

    return moveSet;
}

vector<MetaMove> MetaLearner::collectMoves_INDEP(int currK)
{
    vector<MetaMove> moveSet;

    int numSpecies = speciesNameIDMap.size();

    // Now we will have a move for one orthogroup at a time
    for (int i = 0; i < targetList.size(); i++) {
        int targetID = targetList[i];

        // score the score of best regulator:
        vector<double> bestscore_PerSpecies(numSpecies, 0);
        vector<double> bestscoreImprovement_PerSpecies(numSpecies, 0);
        vector<int> besttarget_PerSpecies(numSpecies, -1);
        vector<int> besttf_PerSpecies(numSpecies, -1);
        int bestcsetid = -1;
        vector<double> bestScoreImprovementTF(numSpecies, 0);

        // The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        // term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.

        for (int specID = 0; specID < speciesNameIDMap.size(); specID++) {
            SpeciesDataManager *sdm = speciesDataSet[specID];
            VariableManager *vMgr = sdm->getVariableManager();

            // Ensure that the target is present in this species' dataset.
            Variable *target = vMgr->getVariable(targetID);
            if (target == nullptr) {
                continue;
            }

            int bestRegulatorID = -1;

            for (int j = 0; j < regulatorList.size(); j++) {
                int regulatorID = regulatorList[j];

                if (targetID == regulatorID) {
                    continue;
                }

                // Ensure regulator exists in this species' dataset.
                Variable *regulator = vMgr->getVariable(regulatorID);
                if (regulator == nullptr) {
                    continue;
                }

                SlimFactor *sFactor = sdm->getFactor(targetID);

                // If the edge already exists in the MB of sFactor continue
                if (sFactor->mergedMB.find(regulatorID) != sFactor->mergedMB.end()) {
                    continue;
                }

                // If the target already has the max num edges, continue.
                if (sFactor->mergedMB.size() >= currK) {
                    continue;
                }

                double scoreImprovement = 0;
                double newScore = 0;
                getNewPLLScore(specID, regulator, target, newScore, scoreImprovement);

                if (scoreImprovement <= 0 || scoreImprovement <= bestScoreImprovementTF[specID]) {
                    continue;
                }

                bestscore_PerSpecies[specID] = newScore;
                bestscoreImprovement_PerSpecies[specID] = scoreImprovement;
                besttarget_PerSpecies[specID] = targetID;
                besttf_PerSpecies[specID] = regulatorID;
                bestcsetid = 1;
                bestScoreImprovementTF[specID] = scoreImprovement;
                bestRegulatorID = regulatorID;
            }

            if (bestcsetid != 1 || bestscoreImprovement_PerSpecies[specID] <= 0) {
                continue;
            }

            MetaMove move;
            move.setSrcVertex(besttf_PerSpecies[specID]);
            move.setTFID(bestRegulatorID);
            move.setConditionSetInd(specID);
            move.setTargetVertex(besttarget_PerSpecies[specID]);
            move.setTargetID(targetID);
            move.setTargetMBScore(bestscore_PerSpecies[specID]);
            move.setScoreImprovement(bestscoreImprovement_PerSpecies[specID]);
            moveSet.push_back(move);
        }
    }

    return moveSet;
}

// u is reg and v is target
// species-specific prior: sum_reg[log(p)]+sum_nonreg[log(1-p)]
// score=likelihood+species-specific prior:
void MetaLearner::getNewPLLScore(int cid, Variable *u, Variable *v, double &targetmbScore, double &scoreImprovement)
{
    SpeciesDataManager *sdm = speciesDataSet[cid];
    SlimFactor *dFactor = sdm->getFactor(v->getID()); // target
    unordered_map<int, double> &varNeighborhoodPrior = varNeighborhoodPrior_PerSpecies[cid];
    unordered_map<int, unordered_map<int, double>> &edgePresenceProb = edgePresenceProb_PerSpecies[cid];
    double currPrior = varNeighborhoodPrior[v->getID()]; // target
    bool toDel_d = true;
    // already checked mergedMB before computing getNewPLLScore
    double plus = 0;
    double minus = 0;
    dFactor->mergedMB.insert(u->getID());
    int status = 0;
    for (auto mIter = dFactor->mergedMB.begin(); mIter != dFactor->mergedMB.end(); mIter++)
    {
        double p = edgePresenceProb[*mIter][v->getID()];
        if (p == 0 || p == 1) {
            continue;
        }
        minus = minus + log(1 - p);
        plus = plus + log(p);
    }
    double pll_d = getPLLScore(cid, dFactor, status);
    if (status == -1)
    {
        scoreImprovement = -1;
        if (toDel_d)
        {
            auto dIter = dFactor->mergedMB.find(u->getID());
            dFactor->mergedMB.erase(dIter);
        }
        return;
    }
    currPrior = currPrior + plus - minus;
    pll_d = pll_d + currPrior;
    targetmbScore = pll_d;
    // Don't include the prior. Just use the data likelihood improvement
    double dImpr = targetmbScore - dFactor->mbScore;
    scoreImprovement = (dImpr <= 0) ? -1 : dImpr;
    if (toDel_d) {
        auto dIter = dFactor->mergedMB.find(u->getID());
        dFactor->mergedMB.erase(dIter);
    }
}

double MetaLearner::getPLLScore(int specID, SlimFactor *sFactor, int &status)
{
    PotentialManager *potMgr = speciesDataSet[specID]->getPotentialManager();
    double pll = potMgr->computeConditionalLL(sFactor, status);
    return pll;
}


double MetaLearner::getEdgePrior(int tfID, int targetID, SpeciesDataManager *sdm)
{
    double motifweight = 0;

    // Pull weight from the motif network, if it exists.
    unordered_map<int, unordered_map<int, double> *> &motifNetwork = sdm->getMotifNetwork();
    if (motifNetwork.find(tfID) != motifNetwork.end()) {
        unordered_map<int, double> *values = motifNetwork[tfID];
        if (values->find(targetID) != values->end()) {
            motifweight = (*values)[targetID];
        }
    }

    // prior probability between regulator j and target k as a logistic function
    double fwt = motifweight * beta2;
    double prior = 1 / (1 + exp(-1 * (beta1 + fwt)));

    // Don't allow prior to be 1 or 0.
    prior = (prior < 1e-6) ? 1e-6 : prior;
    prior = (prior == 1) ? 1 - 1e-6 : prior;

    return prior;
}

double MetaLearner::makeMoves(vector<MetaMove>& moveSet)
{
    double netScoreDelta = 0;
    for (int m = 0; m < moveSet.size(); m++) {
        MetaMove& move = moveSet[m];

        int specID = move.getConditionSetInd();
        int regulatorID = move.getTFID();
        int targetID = move.getTargetID();
        SpeciesDataManager *sdm = speciesDataSet[specID];
        SlimFactor *dFactor = sdm->getFactor(move.getTargetVertex());
        dFactor->mergedMB.insert(move.getSrcVertex());
        dFactor->mbScore = move.getTargetMBScore();

        char varPair[20];
        sprintf(varPair, "%d-%d", regulatorID, targetID);
        string varPairKey(varPair);

        vector<int> *newEdgeStatus;
        if (affectedVarPairs.find(varPairKey) == affectedVarPairs.end()) {
            newEdgeStatus = new vector<int>(speciesNameIDMap.size(), 0);
            affectedVarPairs[varPairKey] = newEdgeStatus;
        } else {
            newEdgeStatus = affectedVarPairs[varPairKey];
        }
        (*newEdgeStatus)[specID] = 1;

        netScoreDelta += move.getScoreImprovement();
    }
    return netScoreDelta;
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
        PotentialManager *potMgr = sdm->getPotentialManager();
        VariableManager *vMgr = sdm->getVariableManager();
        vector<Variable*>& variableSet = vMgr->getVariableSet();
        for (int j = 0; j < variableSet.size(); j++) {
            Variable *target = variableSet[j];
            SlimFactor *sFactor = sdm->getFactor(target->getID());
            potMgr->dumpVarMB(sFactor, oFile);
        }
        oFile.close();
    }
    return 0;
}

void MetaLearner::showModelParameters()
{
    char aFName[1024];
    for (int i = 0; i < speciesDataSet.size(); i++)
    {
        SpeciesDataManager *sdm = speciesDataSet[i];
        VariableManager *vMgr = sdm->getVariableManager();
        PotentialManager *potMgr = sdm->getPotentialManager();

        const char *dirname = sdm->getOutputLoc();
        sprintf(aFName, "%s/modelparams.txt", dirname);
        ofstream oFile(aFName);

        vector<Variable*>& variableSet = vMgr->getVariableSet();
        for (int j = 0; j < variableSet.size(); j++) {
            Variable *var = variableSet[j];
            SlimFactor *sFactor = sdm->getFactor(var->getID());

            double mbcondvar = 0;
            double mbbias = 0;
            unordered_map<int, double> mbwt;
            potMgr->computePotentialMBCovMean(sFactor, mbcondvar, mbbias, mbwt);

            oFile << "Var=" << var->getName() << "\tWt=-1" << "\tCondVar=" << mbcondvar << "\tCondBias=" << mbbias << "\tCondWt=";

            for (auto dIter = mbwt.begin(); dIter != mbwt.end(); dIter++) {
                if (dIter != mbwt.begin()) {
                    oFile << ",";
                }
                Variable *mbVar = vMgr->getVariable(dIter->first);
                oFile << mbVar->getName() << "=" << dIter->second;
            }
            oFile << endl;
        }

        oFile.close();
    }
}
