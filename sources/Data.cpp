
#include "Data.h"

int Data::GetNsite()	{
	return Nsite;
}

int Data::GetNtaxa()	{
	return Ntaxa;
}

int Data::GetState(string tax, int site)	{
	return data[tax][site];
}

string* Data::GetNameList()	{
	return namelist;
}

Data::Data(string infile)	{

	ifstream is(infile.c_str());
	is >> Ntaxa >> Nsite;
	namelist = new string[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		string tax;
		string seq;
		is >> tax >> seq;
		namelist[i] = tax;
		int* tmp = new int[Nsite];
		if (seq.size() != Nsite)	{
			cerr << "error: non matching sequence length\n";
			exit(1);
		}
		for (int i=0; i<Nsite; i++)	{
			char c = seq[i];
			if ((c == '-') || (c == '?'))	{
				tmp[i] = -1;
			}
			else if (c == 'A')	{
				tmp[i] = 0;
			}
			else if (c == 'C')	{
				tmp[i] = 1;
			}
			else if (c == 'G')	{
				tmp[i] = 2;
			}
			else if (c == 'T')	{
				tmp[i] = 3;
			}
			else	{
				cerr << "error: does not recognize character in nucleotide sequence\n";
				exit(1);
			}
		}
		data[tax] = tmp;
	}
}
