#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <unordered_map>
#include <string>
#include <vector>

using namespace std;

class Error;
class Variable;

class VariableManager
{
	public:
		int getVarID(const char*);
		vector<Variable*>& getVariableSet();
		Variable* getVariableAt(int);
		Error::ErrorCode readVariablesFromTable(vector<string>&);
        Error::ErrorCode readVariablesFromData(string,int);

	private:
		unordered_map<string,int> varNameIDMap;
		vector<Variable*> variableSet;
};
#endif
