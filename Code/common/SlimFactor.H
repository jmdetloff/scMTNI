#ifndef _SLIM_FACTOR
#define _SLIM_FACTOR
#include <iostream>
#include<unordered_set>
#include<unordered_map>
using namespace std;

#include "CommonTypes.H"

class SlimFactor
{
	public:
		SlimFactor();
		~SlimFactor();

		int fId;
		double mbScore;
        unordered_set<int> mergedMB;
};
#endif
