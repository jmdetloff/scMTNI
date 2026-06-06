#ifndef _SPECIES_DATA_
#define _SPECIES_DATA_

#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

/*
SpeciesDistance encapsulates a lineage tree of cell species and the relationships
between edge probabilities across the networks inferred for the nodes of this tree.
*/
class SpeciesDistance
{
	public:

		struct Species
		{
			string name;
			Species* parent;
			vector<Species*> children;
			double p_maintain_edge;
			double p_gain;
		};

		SpeciesDistance();
		~SpeciesDistance();

		Species* getRoot();

		// Sets the probability of an edge existing in the root species.
		void setProot(double aVal);

		// Expects a path to a file describing the cell lineage, in the format
		// parent	child	gain_rate	loss_rate
		void readSpeciesTree(const char*);

		// Associates each species in the tree with an id.
		void setspeciesNameIDMap(unordered_map<string,int>& spNameIDMap);

		/*
		Returns the probability of an edge having a particular configuration of
		presence/absence in the networks of the species of this tree. edgeStatus should
		be ordered such that the status at each index maps to the species whose id, as
		set by setspeciesNameIDMap, is that index.
		*/
		double getEdgeStatusProb(vector<int>& edgeStatus);

	private:
		unordered_map<string,Species*> speciesSet;
		unordered_map<string,int> speciesNameIDMap;
		unordered_map<int,double> binEdgeKeyProbMap;
		Species* root;
		double proot;
		double getSubTreeProb(bool, Species*, vector<int>&);
};
#endif
