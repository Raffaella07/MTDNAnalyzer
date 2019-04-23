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
#include "TEllipse.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLegend.h"
#include "NeutralsClass.h"
#include "NeutralsClass.cpp"
struct tr{
	std::vector <float> *pt=0;
	std::vector <float> *In_phi=0;
	std::vector <float> *In_eta=0;
	std::vector <float> *In_x=0;
	std::vector <float> *In_y=0;
	std::vector <float> *In_z=0;
	std::vector <float> *ECAL_phi=0;
	std::vector <float> *ECAL_eta=0;
	std::vector <float> *ECAL_x=0;
	std::vector <float> *ECAL_y=0;
	std::vector <float> *ECAL_z=0;
	std::vector <unsigned short int > *nhits=0;
	std::vector <float> *pt_in_phi=0;
	std::vector <float> *pt_in_eta=0;
	std::vector <float> *pt_out_phi=0;
	std::vector <float> *pt_out_eta=0;
	 


};

struct cnd{
		
	std::vector <float> *phi=0;
	std::vector <float> *eta=0;


}


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
	int i,j;
	struct tr tracks,matched_tracks;
	struct cnd matched_cand;


	int counter=0;
	TH1D * tracks_to_cand = new TH1D("","",r_bin,min_r,max_r);

	for(i=0;i<cand.fChain()->GetEntries();i++){

	
	track.fChain->GetEntry(i);
	cand.fChain->GetEntry(i);

	 tracks.pt=track.Track_pt;
         tracks.In_phi=track.Track_In_phi;
         tracks.In_eta=track.Track_In_eta;
         tracks.In_x=track.Track_In_x;
         tracks.In_y=track.Track_In_y;
         tracks.In_z=track.Track_In_z;
         tracks.ECAL_phi=track.Track_ECAL_phi;
         tracks.ECAL_eta=track.Track_ECAL_eta;
         tracks.ECAL_x=track.Track_ECAL_x;
         tracks.ECAL_y=track.Track_ECAL_y;
	 tracks.ECAL_z=track.Track_ECAL_z;
         tracks.nhits=track.Track_nhits;
         tracks.pt_in_phi=track.Track_pt_Iphi;
         tracks.pt_in_eta=track.Track_pt_Ieta;
         tracks.pt_out_phi=track.Track_pt_Ephi;
         tracks.pt_out_eta=track.Track_pt_Eeta;

	for(j=0; j<cand.pt->size();j++){
		float dr_min = 1e6;
		for(k=0;k<track.pt->size();k++){
		float cand_eta=-1;
		float cand_phi=-1;
		float track_eta=-1;
		float track_phi=-1;
		
		float dr = sqrt(pow((cand.ecaleta->at(j)-tracks.ECAL_eta->at(k)),2)+pow((cand.ecalphi->at(j)-tracks.ECAL_phi->at(k)),2))
			if(dr< dr_min){
			dr_min = dr;
			plotter();
			counter++;
			
			float cand_eta=cand.ecaleta->at(j);
			float cand_phi=cand.ecalphi->at(j);
			float track_eta=tracks.ECAL_eta->at(k);
			float track_phi=tracks.ECAL_phi->at(k);
				}

			}
		matched_tracks.ECAL_eta->push_back(track_eta);
		matched_tracks.ECAL_phi->push_back(track_phi);
		matched_cand.eta->push_back(cand_eta);
		matched_cand.phi->push_back(cand_phi);



	}


}



}


void plotter(){

	int i;
	TEllipse * ECAL= new TEllipse(0,0,129,129);	
	TEllipse * MTD= new TEllipse(0,0,120,120);	
	TBox* ECAL[cand.eta->size()];
	TArrow* ECALtrack[cand.eta->size];
	TArrow* Innertrack[cand.eta->size];
	
	for(i=0; i<track.eta->size();i++){

	ECAL[i]=new TBox(129*cos(cand.phi->at(i)),129*sin(cand.phi->at(i)),(129+15)*cos(cand.phi->at(i))+10*sin(cand.phi->at(i),(129+15)*sin(cand.phi->at(i)));


	}





}
