#ifndef __DATA_H__ 
 #define __DATA_H__
// this class implements a container for a genetic dataset
// which is an alignment of <Ntaxa> DNA sequences, each of length <Nsite>

#include "Random.h"
#include <map>

class Data	{

	public:

	// constructor: create a dataset, reading the sequences from the file <infile>
	Data(string infile);

	// returns the number of sites (number of nucleotides) of the aligned sequences
	int GetNsite();

	// returns the number of species (number of taxa)
	int GetNtaxa();

	// returns the nucleotide state for taxon <tax> and site <site>
	// 0 : A
	// 1 : C
	// 2 : G
	// 3 : T
	// -1: unknown state
	int GetState(string tax, int site);

	// return a pointer to the list of names of the taxa
	// useful for creating a tree object whose leaves correspond to the taxa of the dataset
	string* GetNameList(); 

	private:

	map<string,int*> data;
	int Nsite;
	int Ntaxa;
	string*namelist;

};
#endif
