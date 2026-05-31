#ifndef _ORTHOGROUP_
#define _ORTHOGROUP_
#include <string>
#include <map>
#include <vector>
using namespace std;

class MappedOrthogroup
{
	public:
		int addMembers(string&);
		int setID(int);
		int getID();
		vector<string>& getOrthoMembers();

	private:
		vector<string> orthoMembers;
		int oid;
};
#endif
