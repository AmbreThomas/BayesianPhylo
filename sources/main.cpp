
#include "Sampler.h"

int main(int argc, char* argv[])	{

	// name of the datafile
	string datafile = argv[1];

	// name of the MCMC run
	string name = argv[2];

	// create the MCMC sampler
	Sampler* sampler = new Sampler(datafile);

	// run the sampler forever...
	sampler->Run(name);
}
