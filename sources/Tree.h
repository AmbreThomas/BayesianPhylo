
#include "Node.h"

class Tree	{

	public:

	// constructor based on a file
	Tree(string filename);

	// construct a random tree, with n leaves, with names given by an array of strings <namelist>
	// draw from the coalescent
	Tree(int n, string* namelist);

	// read the tree from a stream
	void ReadFromStream(istream& is);

	// write the tree into a stream
	void WriteToStream(ostream& os);

	// return a pointer to the root node
	Node* GetRoot()	{
		return root;
	}

	// makes a backup of the current tree topology and branch lengths
	void Backup();

	// restore tree such as was saved upon the last backup
	void Restore();

	// -----------------
	// TO BE IMPLEMENTED
	// -----------------
	
	// return the total length of the tree (sum of lengths of all branches)
	// this should be calculated recursively
	// you can look at the GetSize() function for an example of a recursive function
	double GetLength()	{
		return 0;
	}

	// -----------------
	// ALREADY IMPLEMENTED
	// -----------------

	// returns the size of the tree (the number of leaves)
	// the size is calculated recursively (using RecursiveGetSize, declared below)
	int GetSize()	{
		return RecursiveGetSize(root);
	}

	int RecursiveGetSize(Node* node)	{
		if (! node->left)	{
			return 1;
		}
		return RecursiveGetSize(node->left) + RecursiveGetSize(node->right);
	}

	// METROPOLIS HASTINGS METHODS

	// -----------------
	// TO BE IMPLEMENTED
	// -----------------
	
	// choose a node at random (other than leaf or root)
	// move the age of that node (within the allowed time interval)
	void ProposeTimeMove(double tuning);

	// -----------------
	// ALREADY IMPLEMENTED
	// -----------------

	// choose a node at random (should not be root nor a direct child of root)
	// perform a SPR (Subtree Pruning and Regrafting) move
	// i.e. regrafting the subtree pending from the node
	// on any branch among those that are available at the same branching time
	void ProposeSPRMove();

	// choose an internal node (other than root or leaf) uniformly at random
	// return a pointer to that node
	Node* ChooseInternalNodeAtRandom();

	private:

	void RecursiveWriteToStream(Node* from, ostream& os);

	Node* ParseGroup(string s);

	void SetRoot(Node* inroot)	{
		root = inroot;
	}

	// utilities for the SPR move
	Node* ChooseSPRNodeAtRandom();
	void RecursiveGetTargets(Node* node, vector<Node*>& list, double time, Node* cousin);


	// data structures
	// a vector contaning the list of all nodes of the tree
	vector<Node*> nodelist;
	// a backup
	vector<Node*> bknodelist;
	// the root node
	Node* root;
};

