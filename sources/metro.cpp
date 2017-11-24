#include "Random.h"

#include <stdio.h>

// Donnees
int N = 10;
double t = 5.0;

double logposterior(double u){
	return -(t+1)*u + N*log(t*u);
}

int main()	{  
	int Ncycle = 10000;
	double delta = 5.0;
	
	double Nacc = 0;
  
	// Valeur initiale
	double u = 1.0;
  
	// Metropolis
	for(int i = 0 ; i < Ncycle ; i++){
		double bku = u;
		
		double logpost1 = logposterior(u);
		
		u = u + delta * (Random::Uniform() - 0.5);
		if(u < 0)
			u = -u;
			
		double logpost2 = logposterior(u);
		
		double logalpha = logpost2 - logpost1;
	
		if(log(Random::Uniform()) < logalpha){
			Nacc++;
			printf("%f\n",u);
		}else{
			u = bku;
		}
	}
	printf("Taux d'acceptation : %f\n",Nacc/Ncycle);
}
