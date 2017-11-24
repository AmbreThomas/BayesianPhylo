
#include "Tree.h"

// -----------------
// TO BE IMPLEMENTED
// -----------------

void Tree::ProposeTimeMove(double tuning)	{

	// choose an internal node at random
	// move its age within allowed interval
}

// -----------------
// ALREADY IMPLEMENTED
// -----------------

Tree::Tree(string filename)	{

	ifstream is(filename.c_str());
	if (!is)	{
		cerr << "cannot find file : " << filename << '\n';
		exit(1);
	}
	ReadFromStream(is);
}

Tree::Tree(int n, string* namelist)	{

	cerr << "start with : " << n << '\n';
	// start with n lineages
	Node** node = new Node*[n];
	for (int i=0; i<n; i++)	{
		node[i] = new Node(0,0,namelist[i],0);
		nodelist.push_back(node[i]);
	}

	double time = 0;
	for (int k=n; k>1; k--)	{

		// draw exponential waiting time
		double t = 2 * Random::sExpo() / k / (k-1);
		time += t;

		int i1 = (int) (k*Random::Uniform());
		int i2 = (int) ((k-1)*Random::Uniform());
		if (i2 >= i1)	{
			i2++;
		}

		// lineage maxi will be merged with mini
		int mini = (i1 < i2) ? i1 : i2;
		int maxi = (i1 > i2) ? i1 : i2;

		// but first, swap lineage maxi and last lineage of the table
		// so that, after the merge, the first k-1 entries of the table correspond to the k-1 remaining lineages
		if (maxi < k-1)	{

			Node* tmpnode = node[maxi];
			node[maxi] = node[k-1];
			node[k-1] = tmpnode;
		}
			
		Node* newnode = new Node(node[mini],node[k-1],"",time);
		node[mini] = newnode;
		nodelist.push_back(newnode);
	}

	// renormalize times
	for (int i=0; i<2*n-1; i++)	{
		nodelist[i]->SetTime(nodelist[i]->GetTime() / time);
	}
	SetRoot(node[0]);

	delete[] node;

	cerr << "check size:\n";
	cerr << GetSize() << '\n';
}

void Tree::WriteToStream(ostream& os)	{

	RecursiveWriteToStream(GetRoot(),os);
	os << ";\n";
}

void Tree::RecursiveWriteToStream(Node* from, ostream& os)	{

	if (from->isLeaf())	{
		os << from->GetNodeName();
	}
	else	{
		os << '(';
		RecursiveWriteToStream(from->left,os);
		os << ',';
		RecursiveWriteToStream(from->right,os);
		os << ')';
	}
	if (! from->isRoot())	{
		os << ':';
		os << from->up->GetTime() - from->GetTime();
	} 
}

void Tree::ReadFromStream(istream& is)	{

	string expr = "";
	int cont = 1;
	while (cont)	{
		string s;
		is >> s;
		unsigned int k = 0;
		while ((k<s.length()) && (s[k] != ';'))	k++;
		expr += s.substr(0,k);
		cont = (! is.eof()) && (k==s.length());
	}
	SetRoot(ParseGroup(expr));
}

Node* Tree::ParseGroup(string input)	{

	try	{

		// parse input as (body)nodeval:branchval

		string body = "";
		unsigned int k = 0;
		if (input[0] == '(')	{
			int brack = 1;
			k = 1;
			while ((k<input.length()) && brack)	{
				char c = input[k];
				if (c == '(') brack++;
				if (c == ')') brack--;
				k++;
			}
			if (brack)	{
				cout << "in parse group: too many (\n";
				cout << input << '\n';
				throw;
			}
			body = input.substr(1,k-2);
		}

		int b = k;
		while ((k<input.length()) && (input[k]!=':')) k++;
		string nodeval = input.substr(b,k-b);

		string branchval = "";
		if (k<input.length())	{
			branchval = input.substr(k+1, input.length() - k);
		}

		// parse input as a list of strings separated by ','
		string ns1,ns2;
		int n = input.size();
		k = 0;
		int brack = 0;
		b = 0;
		while (k<n)	{
			char c = input[k];
			if (c == '(') brack++;
			if (c == ')') brack--;
			if ((!brack) && (c == ','))	{
				ns1 = (string) (input.substr(b,k-b));
				b = k+1;
				ns2 = (string) (input.substr(b,n-b));
			}

			if (brack<0)	{
				cout << "in parse list : too many )\n";
				cout << input.substr(0,k) << '\n';
				cout << input << '\n';
				throw;
			}
			k++;
		}
		Node* node1 = ParseGroup(ns1);
		Node* node2 = ParseGroup(ns2);
		// make a new node and a new branch
		Node* node = new Node(node1,node2,nodeval,branchval);
		nodelist.push_back(node);

		return node;
	}
	catch(...)	{
		cout << "exit in parse group\n";
		exit(1);
	}
}

void Tree::Backup()	{

	if (!bknodelist.size())	{
		for (size_t i=0; i<nodelist.size(); i++)	{
			Node* node = new Node(0,0,"",0);
			bknodelist.push_back(node);
		}
	}
	for (size_t i=0; i<nodelist.size(); i++)	{
		bknodelist[i]->Backup(nodelist[i]);
	}
}

void Tree::Restore()	{

	for (size_t i=0; i<nodelist.size(); i++)	{
		nodelist[i]->Backup(bknodelist[i]);
	}

}

Node* Tree::ChooseInternalNodeAtRandom()	{

	int n = GetSize() - 2;
	int choose = 1 + ((int) (n * Random::Uniform()));
	int k = 0;
	while (choose)	{

		while ((k < nodelist.size()) && ((nodelist[k]->isRoot()) || (nodelist[k]->isLeaf())))	{
			k++;
		}
		if (k == nodelist.size())	{
			cerr << "error in choose node: overflow\n";
			exit(1);
		}

		choose--;

		if (choose)	{
			k++;
		}
	}

	return nodelist[k];
}

Node* Tree::ChooseSPRNodeAtRandom()	{

	int n = 2 * GetSize() - 4;
	int choose = 1 + (int) (n * Random::Uniform());
	int k = 0;
	while (choose)	{

		while ((k < nodelist.size()) && ((nodelist[k]->isRoot()) || (nodelist[k]->up->isRoot())))	{
			k++;
		}
		if (k == nodelist.size())	{
			cerr << "error in choose node: overflow\n";
			exit(1);
		}

		choose--;

		if (choose)	{
			k++;
		}
	}

	return nodelist[k];
}

void Tree::RecursiveGetTargets(Node* node, vector<Node*>& list, double time, Node* cousin)	{

	if (! node->isRoot())	{
		if ((node->GetTime() < time) && (node->up->GetTime() > time) && (node != cousin))	{
			list.push_back(node);
		}
	}
	if (! node->isLeaf())	{
		RecursiveGetTargets(node->left,list,time,cousin);
		RecursiveGetTargets(node->right,list,time,cousin);
	}
}

void Tree::ProposeSPRMove()	{

	int sizebefore = GetSize();

	Node* node = ChooseSPRNodeAtRandom();
	Node* father = node->up;
	double time = father->GetTime();
	Node* grandfather = node->up->up;
	Node* cousin = (node == father->left) ? father->right : father->left;

	if (father == grandfather->left)	{
		grandfather->left = cousin;
	}
	else	{
		grandfather->right = cousin;
	}
	cousin->up = grandfather;

	vector<Node*> targets;
	RecursiveGetTargets(GetRoot(),targets,time,cousin);
	if (! targets.size())	{
		cerr << "error: could not find target\n";
		exit(1);
	}


	int n = (int) (targets.size() * Random::Uniform());
	Node* newcousin = targets[n]; 
	Node* newgrandfather = newcousin->up;

	if (newgrandfather->left == newcousin)	{
		newgrandfather->left = father;
	}
	else	{
		newgrandfather->right = father;
	}
	father->up = newgrandfather;

	if (father->left == node)	{
		father->right = newcousin;
	}
	else	{
		father->left = newcousin;
	}
	newcousin->up = father;

	int sizeafter = GetSize();
	if (sizeafter != sizebefore)	{
		cerr << "error: spr move changed tree size\n";
		cerr << "before : " << sizebefore << '\n';
		cerr << "after : " << sizeafter << '\n';
		cerr << targets.size() << '\n';
		cerr << newcousin->GetTime() << '\t' << father->GetTime() << '\t' << newgrandfather->GetTime() << '\n';
		exit(1);
	}
}

