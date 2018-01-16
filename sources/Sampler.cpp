#include "Random.h"
#include "Sampler.h"
#include <cmath>
#include "Tree.h"

// -----------------
// TO BE IMPLEMENTED
// -----------------

double Sampler::GetLogPrior()	{
	return -rate;
}

double Sampler::GetLogLikelihood()	{
	double logLikelihood = 0;
	for(int i = 0 ; i < data->GetNsite() ; i++){
		logLikelihood += SiteLogLikelihood(i);
	}
	return logLikelihood;
}

double Sampler::SiteLogLikelihood(int site)	{
	RecursiveConditionnalSiteLikelihood(site, tree->GetRoot());
	double* likelihood_root = tree->GetRoot()->condl;
	double likelihood = 0;
	for(int i = 0 ; i < 4 ; i++){
		likelihood += likelihood_root[i];
	}
	return log(likelihood/4);
}

void Sampler::RecursiveConditionnalSiteLikelihood(int site, Node* node){
	double* likelihood = node->condl;
	if(! node->left){
		int state = data->GetState(node->GetNodeName(),site);
		if(state==-1){
			for(int i = 0 ; i < 4 ; i++){
				likelihood[i] = 1;
			}
		}
		else{
			for(int i = 0 ; i < 4 ; i++){
				if(i!=state){
				likelihood[i] = 0;
			}else{
				likelihood[i] = 1;
				}
			}
		}
	}else{
		RecursiveConditionnalSiteLikelihood(site, node->left);
		RecursiveConditionnalSiteLikelihood(site, node->right);
		double* likelihood_left = node->left->condl;
		double* likelihood_right = node->right->condl;
		double length_left = node->left->GetBranchLength();
		double length_right = node->right->GetBranchLength();
		double sum_left = 0;
		double sum_right = 0;
		double exp_left = exp(-rate*length_left);
		double exp_right = exp(-rate*length_right);
		for(int i = 0 ; i < 4 ; i++){
			sum_left += likelihood_left[i];
		}
		for(int i = 0 ; i < 4 ; i++){
			sum_right += likelihood_right[i];
		}
		
		for(int i = 0 ; i < 4 ; i++){
			likelihood[i] = ((1-exp_left)*sum_left/4+exp_left*likelihood_left[i]) * ((1-exp_right)*sum_right/4+exp_right*likelihood_right[i]);
		}
	}
}

int Sampler::RateMove(double tuning)	{
	double backupRate = this->rate;
	this->rate = abs(Random::sNormal() * tuning + this->rate);
	double newLogLikelihood = GetLogLikelihood();
	if(newLogLikelihood-logLikelihoodBackup < 0 && newLogLikelihood-logLikelihoodBackup < log(Random::Uniform()) && logLikelihoodBackup != 0){
		this->rate = backupRate;
		return 0;
	}else{
		logLikelihoodBackup = newLogLikelihood;
		return 1;
	}
}

int Sampler::TimeMove(double tuning)	{
	tree->Backup();
	tree->ProposeTimeMove(tuning);
	double newLogLikelihood = GetLogLikelihood();
	if(newLogLikelihood-logLikelihoodBackup < 0 && newLogLikelihood-logLikelihoodBackup < log(Random::Uniform()) && logLikelihoodBackup != 0){
		tree->Restore();
		return 0;
	}else{
		logLikelihoodBackup = newLogLikelihood;
		return 1;
	}
}

int Sampler::TopoMove()	{
	tree->Backup();
	tree->ProposeSPRMove();
	double newLogLikelihood = GetLogLikelihood();
	if(newLogLikelihood-logLikelihoodBackup < 0 && newLogLikelihood-logLikelihoodBackup < log(Random::Uniform()) && logLikelihoodBackup != 0){
		tree->Restore();
		return 0;
	}else{
		logLikelihoodBackup = newLogLikelihood;
		return 1;
	}
}

void Sampler::Cycle()	{
	acceptedRateMove += RateMove(0.02);
	acceptedTimeMove += TimeMove(0.04);
	acceptedTopoMove += TopoMove();
}


