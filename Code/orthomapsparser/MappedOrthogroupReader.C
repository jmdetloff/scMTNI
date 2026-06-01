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
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"

//This function reads the orthogroup mapping produced in any one of the PROPER, FILLED and UNFILLED formats.
//We are essentialy interested in the first two columns
int 
MappedOrthogroupReader::readFile(const char* aFName)
{
    cout << "MappedOrthogroupReader::readFile - " << aFName << endl;
    ifstream inFile(aFName);
    char* buffer=NULL;
    int bufflen=0;
    string buffstr;
    int lineCnt=0;  // FIX ME - should require header to have initial #
    int oldogid=-2;
    MappedOrthogroup* ogrp=NULL;
    while(inFile.good())
    {
        getline(inFile,buffstr);
        if(lineCnt==0)
        {
            lineCnt++;
            continue;
        }
        if(bufflen<=buffstr.length())
        {
            if(buffer!=NULL)
            {
                delete[] buffer;
            }
            bufflen=buffstr.length()+1;
            buffer=new char[bufflen];
        }
        strcpy(buffer,buffstr.c_str());
        //Get the orthogroup and the genes
        int ogid=-1;
        string geneSet;
        char* tok=strtok(buffer,"\t");
        //cout << "  " << tok << endl;
        int tokCnt=0;
        while(tok!=NULL)
        {
            if(tokCnt==0)
            {
                char* pos=strchr(tok,'_');
                if(pos==NULL)
                {
                    cout <<"Bad OG name " << tok<<  endl;
                    exit(0);
                }
                *pos='\0';
                ogid=atoi(tok+2);
                //cout << lineCnt<< " oldogid=" << oldogid <<"ogid=" << ogid << endl;
                if(oldogid!=ogid)
                {
                    ogrp=new MappedOrthogroup;
                    ogrp->setID(ogid);
                    orthogroupSet[ogrp->getID()]=ogrp;
                    oldogid=ogid;
                }
                //ogrp->incrCnt();
            }
            else if(tokCnt==1)
            {
                addMembers(tok,ogrp);
            }
            tokCnt++;
            tok=strtok(NULL,"\t");
        }
        lineCnt++;
    }
    inFile.close();
    generateGeneOrthoMap();
    return 0;
}

int
MappedOrthogroupReader::readSpeciesMapping(const char* aFName)
{
    cout << "MappedOrthogroupReader::readSpeciesMapping() --- " <<aFName <<endl;
    ifstream inFile(aFName);
    char buffer[1024];
    int lineCnt=0;
    while(inFile.good())
    {
        inFile.getline(buffer,1023);
        if(strlen(buffer)<=0)
        {
            continue;
        }
        string speciesName(buffer);
        speciesIDNameMap[lineCnt]=speciesName;
        //cout << "speciesIDNameMap["<<lineCnt << "]=" << speciesName << endl;
        lineCnt++;
    }
    inFile.close();
    return 0;
}

int 
MappedOrthogroupReader::getMappedOrthogroupID(const char* geneName, const char* species)
{
    string speciesKey(species);
    string geneKey(geneName);
    int orthoID=(*geneOrthoMap[speciesKey])[geneKey];
    return orthoID;
}

map<int,MappedOrthogroup*>& 
MappedOrthogroupReader::getMappedOrthogroups()
{
    return orthogroupSet;
}

MappedOrthogroup*
MappedOrthogroupReader::getMappedOrthogroup(const char* geneName, const char* species)
{
    int orthoID=getMappedOrthogroupID(geneName,species);
    MappedOrthogroup* og=orthogroupSet[orthoID];
    return og;
}

int
MappedOrthogroupReader::addMembers(char* abuffer, MappedOrthogroup* ogrp)
{
    bool debug = false;  // can set this to true in the debugger to display useful stuff
    char* tok=abuffer;
    int tokCnt=0;
    while(tok!=NULL)
    {
        char* end=strchr(tok,',');
        if(end!=NULL)
        {
            *end='\0';
        }
        if(speciesIDNameMap.find(tokCnt)==speciesIDNameMap.end())
        {
            cout <<"No species with id " << tokCnt << endl;
            exit(0);
        }
        string geneName(tok);
        if(strcmp(geneName.c_str(),"NONE")!=0)
        {
            ogrp->addMembers(geneName);
        }
        tokCnt++;
        if(end!=NULL)
        {
            tok=end+1;
        }
        else
        {
            tok=end;
        }
    }
    return 0;
}

// generate map<string,STRINTMAP*> geneOrthoMap;
// (*geneOrthoMap[cell])[genename]=OGID
int 
MappedOrthogroupReader::generateGeneOrthoMap()
{
    for(map<int,MappedOrthogroup*>::iterator oIter=orthogroupSet.begin();oIter!=orthogroupSet.end();oIter++)
    {
        MappedOrthogroup* og=oIter->second;
        vector<string>& members=og->getOrthoMembers();
        for(int sIter=0;sIter<members.size();sIter++) //map<string,string>::iterator
        {
            string speciesName=speciesIDNameMap[sIter];
            STRINTMAP* allgenesForSpecies=NULL;
            if(geneOrthoMap.find(speciesName)==geneOrthoMap.end())
            {
                allgenesForSpecies=new STRINTMAP;
                geneOrthoMap[speciesName]=allgenesForSpecies;
            }
            else
            {
                allgenesForSpecies=geneOrthoMap[speciesName];
            }
            (*allgenesForSpecies)[members[sIter]]=oIter->first;
        }
    }
    return 0;
}
