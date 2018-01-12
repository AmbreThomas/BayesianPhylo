
#include "Sampler.h"

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
	RecursiveConditionnalSiteLikelihood(site, tree->root);
	double* likelihood_root = root->condl;
	double likelihood = 0;
	for(int i = 0 ; i < 4 ; i++){
		likelihood += likelihood_root[i];
	}
	return log(likelihood/4);
}

void Sampler::RecursiveConditionnalSiteLikelihood(int site, Node* node){
	double likelihood [4] = {0,0,0,0};
	if(! node->left){
		int state = data->GetState(node->GetGetNodeName,site);
		if(state==-1){
			for(int i = 0 ; i < 4 ; i++){
				likelihood[i] = 1;
			}
		}
		else{
			likelihood[state] = 1;
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
	node->condl = likelihood;
}

int Sampler::RateMove(double tuning)	{
	return 0;
}

int Sampler::TimeMove(double tuning)	{
	tree.ProposeTimeMove(tuning);
	return 0;
}

int Sampler::TopoMove()	{
	return 0;
}

void Sampler::Cycle()	{
}


