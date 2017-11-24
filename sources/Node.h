

const int Nstate = 4;

#include <vector>
#include "Random.h"

class Node	{

	public:

	// returns the name of the node (the name of the taxon if leaf node)
	string GetNodeName()	{
		return nodename;
	}

	// allows you to change the age of the node
	void SetTime(double intime)	{
		time = intime;
	}

	// returns the age of the node
	double GetTime()	{
		return time;
	}

	// returns the length of the branch (the difference between the times of the two nodes delimiting the branch)
	double GetBranchLength()	{
		if (isRoot())	{
			return 0;
		}
		return up->time - time;
	}

	bool isRoot()	{
		return (up == 0);
	}

	bool isLeaf()	{
		return (left == 0);
	}

	Node(Node* inleft, Node* inright, string innodename, double intime)	{

		left = inleft;
		right = inright;
		up = 0;
		if (left)	{
			left->up = this;
		}
		if (right)	{
			right->up = this;
		}
		nodename = innodename;
		time = intime;
		branchname = "";

		condl = new double[Nstate];
		for (int i=0; i<Nstate; i++)	{
			condl[i] = 0;
		}
	}

	Node(Node* inleft, Node* inright, string innodename, string inbranchname)	{
		left = inleft;
		right = inright;
		up = 0;
		if (left)	{
			left->up = this;
		}
		if (right)	{
			right->up = this;
		}
		nodename = innodename;
		branchname = inbranchname;
		time = atof(nodename.c_str());

		condl = new double[Nstate];
		for (int i=0; i<Nstate; i++)	{
			condl[i] = 0;
		}

	}

	void Backup(Node* from)	{
		up = from->up;
		left = from->left;
		right = from->right;
		time = from->time;
	}

	Node* up;
	Node* left;
	Node* right;

	double* condl;
	string	nodename;
	string branchname;
	double time;

};
