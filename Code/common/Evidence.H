#ifndef _EVIDENCE
#define _EVIDENCE
//This is the evidence for a particular variable
//This can take the place of the likelihood

#include <fstream>
using namespace std;
#include "CommonTypes.H"

class Evidence
{
	public:
		double getEvidVal();

	private:
		double contValue;
};

#endif
