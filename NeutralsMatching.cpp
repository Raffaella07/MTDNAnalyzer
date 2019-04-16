
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <utility>  
#include <iostream>
#include <vector>
#include "TH2.h"
#include <string>
#include "TString.h"
#include "TGraph.h"
#include "TLegend.h"
#include "NeutralsClass.h"
#include "NeutralsClass.cpp"




int main(int argc, char **argv){


	TFile* file = TFile::Open(argv[1]);
	std::string outname = argv[2];
	std::cout << outname << std::endl;
	TTree* gentree = (TTree*)file->Get("gen_tree");
	TTree* tracktree= (TTree*)file->Get("tracks_tree");
	TTree* mcttree = (TTree*)file->Get("mct_tree");
	TTree* candtree= (TTree*)file->Get("cand_tree");
	NeutralsClass track,gen,mct,cand;
	track.Init(tracktree);
	gen.Init(gentree);
	mct.Init(mcttree);
	cand.Init(candtree);

	int i, j,k;

	for(i=0;i<gen.fChain->GetEntries();i++){
		gen.fChain->GetEntry(i);
		mct.fChain->GetEntry(i);
		cand.fChain->GetEntry(i);
		std::cout << "_________EVENT" << gen.event << std::endl;
		for(int g =0; g<gen.pt->size();g++){
		//	for(j=0;j<mct.fChain->GetEntries();j++){
		//		mct.fChain->GetEntry(j);}
				for(int p=0; p<mct.pt->size();p++){
				//	for(k=0;k<cand.fChain->GetEntries();k++){
				//		cand.fChain->GetEntry(k);}
						for(int c =0; c<cand.pt->size();c++){
							if(gen.event == mct.event && mct.event == cand.event){
								if(gen.label->at(g) == mct.label->at(p)){
									if(mct.label->at(p) == cand.label->at(c) && mct.plabel->at(p) == cand.plabel->at(c)){
										std:: cout <<"cand" << c <<  "____      pt            phi             eta  " << std::endl;
										std::cout << "______________" << cand.pt->at(c) << "      " << cand.phi->at(c) << "      " << cand.eta->at(c) << std::endl;
										if(mct.PId->at(p)== 4)	std::cout << "sim photon____" << mct.pt->at(p) << "      " << mct.phi->at(p) << "     " << mct.eta->at(p) << std::endl;
										if(mct.PId->at(p)== 2)	std::cout << "sim electron__"  << mct.pt->at(p) << "     " << mct.phi->at(p) << "      " << mct.eta->at(p) << std::endl;
										std::cout << "gen photon____" << gen.pt->at(g) << "      " << gen.phi->at(g) << "      " << gen.eta->at(g) << std::endl;



									}


								}	
							}	
						}}
		}}
}
