
#include "Sampler.h"

// -----------------
// TO BE IMPLEMENTED
// -----------------

double Sampler::GetLogPrior()	{
	return 0;
}

double Sampler::GetLogLikelihood()	{
	return 0;
}

double Sampler::SiteLogLikelihood(int site)	{
	return 0;
}

void Sampler::Pruning(int site, Node* node)	{
}

int Sampler::RateMove(double tuning)	{
	return 0;
}

int Sampler::TimeMove(double tuning)	{
	return 0;
}

int Sampler::TopoMove()	{
	return 0;
}

void Sampler::Cycle()	{
}


