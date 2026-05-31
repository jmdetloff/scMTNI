/*The goal of this class is to describe a graph edit move we make to generate a new graph.
 * We will make a move on a set of variables, typically a pair of variables, but may be more
 * if we are doing exact search. We need to store the type of move and the score
 * because of this new move. The move could be an edge addition, edge deletion, edge reversal,
 * Markov blanket addition, which is equal to addition of a set of edges */
#ifndef _MOVE_
#define _MOVE_

class MetaMove
{
	public:
		int setScoreImprovement(double);
		int setTargetMBScore(double);
		int setSrcVertex(int);
		int setTargetVertex(int);
		int setConditionSetInd(int);
		int setTFID(int);
		int setTargetID(int);
		int getTFID();
		int getTargetID();
		int getSrcVertex();
		int getTargetVertex();
		int getConditionSetInd();
		double getTargetMBScore();
		double getScoreImprovement();

	private:
		double scoreDelta;
		double targetMBScore;
		int src;
		int target;
		int regi;
		int targeti;
		int conditionSetInd;
};

#endif
