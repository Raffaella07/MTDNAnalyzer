#include "TRandom.h"
#include "TMath.h"
#include <iostream>
#include <utility>
#include <math.h>
#include <time.h>


int collision(int , float );

int main(){
	int i;
	float sigma=0.872;
	TRandom* rdm = new TRandom();
	int PU=200, tries=10000;
	float mean;
	int vertices[tries];
	for(i=0;i<tries; i++){
	if(rdm->Uniform(1) < 1) vertices[i]=collision(200,sigma);
	else vertices[i]=200;
	}
	std::cout <<" the mean of number of vertices for " << tries << " LHC collisions at " << PU << "pileup is " << TMath::Mean(tries,vertices) << "#pm" << TMath::RMS(tries,vertices)/sqrt(tries) << std::endl;
	

}




int collision(int PU, float sigma){
	int i;
	int counter=0;
	TRandom* rdm = new TRandom();
	float value[PU];
	Double_t vertex;
	float min, max;
	float seed;
	//time(&seed);
	rdm->SetSeed();
	for(i=0;i<PU;i++){
		value[i] = rdm->Gaus(0.,5.);
		
	}
	
	vertex=value[rdm->Integer(200)];
	std::cout <<"higgs vertex " << counter << " probable vertices " << std::endl; 

	
	for(i=0;i<PU;i++){
		if(abs(value[i]-vertex)<sigma) counter++;
		
	}

	std::cout << "With a " << sigma << "resolution we have " << counter << " probable vertices " << std::endl; 
	return counter;
}
