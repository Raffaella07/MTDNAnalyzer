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


float matcher(float,float,std::vector<float>*,std::vector <float>* ,std::pair <float,float>*, std::pair <float,float>* );

void SavePlot (char * , TH1D *, const char*, bool );

void Slicer(std::string ,int ,float, float ,std::string ,TH2D *,std::string );

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
	int i,j,k;
	float dr,label=-1;
	int bin_dr,bin_pt;
	float dr_min, dr_max,pt_min,pt_max;
	std::string PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events/VertexPointer");
	system(("mkdir "+PLOTPATH).c_str());

	bin_dr=15;
	bin_pt=15;
	dr_min=0;
	dr_max=1;
	pt_min=0;
	pt_max=10;
	
	TH1D * proj_full[4];
	TH2D* dr_photon_clus = new TH2D("","dr among clusters and matched mct photons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* dr_electron_clus = new TH2D("","dr among clusters and matched mct electrons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* photon_dr_clusMTD = new TH2D("","dr among clusters and MTD hit for mct photons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* ele_dr_clusMTD = new TH2D("","dr among clusters and Mtd hit for mct electrons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);

	for(i=0;i<cand.fChain->GetEntries();i++){
		cand.fChain->GetEntry(i);
		mct.fChain->GetEntry(i);
		for(j=0;j<mct.pt->size();j++){
			std::pair <float,float> matched_mct;
			std::pair <float,float> matched_clus;
			std::pair <float,float> matched_MTDhit;
			dr = matcher(mct.eta->at(j),mct.phi->at(j),cand.clus_eta, cand.clus_phi,&matched_mct,&matched_clus);
			if(dr != -1){
				label = mct.PId->at(j);
				if (label == 4 ) dr_photon_clus->Fill(dr,mct.pt->at(j));
				else if(label == 2) dr_electron_clus->Fill(dr,mct.pt->at(j));
			}
			dr = matcher(matched_clus.first,matched_clus.second,cand.MTDeta, cand.MTDphi,&matched_mct,&matched_MTDhit);
			if(dr !=-1){
				if (label == 4 ) photon_dr_clusMTD->Fill(dr,mct.pt->at(j));
				else if(label == 2) ele_dr_clusMTD->Fill(dr,mct.pt->at(j));

			}






		}

	}
	system(("mkdir"+PLOTPATH+"/controlplots").c_str());
	Slicer(PLOTPATH,bin_pt,pt_min,pt_max,"dr",photon_dr_clusMTD,"Sphoton_dr_clusMTD");
	Slicer(PLOTPATH,bin_pt,pt_min,pt_max,"dr",ele_dr_clusMTD,"Sele_dr_clusMTD");
	proj_full[0]=dr_photon_clus->ProjectionX("1",0,bin_pt-1);
	proj_full[1]=dr_electron_clus->ProjectionX("2",0,bin_pt-1);
        proj_full[2]=photon_dr_clusMTD->ProjectionX("3",0,bin_pt-1);
	proj_full[3]=ele_dr_clusMTD->ProjectionX("4",0,bin_pt-1);
	dr_photon_clus->GetXaxis()->SetTitle("dr");
	dr_electron_clus->GetXaxis()->SetTitle("dr");
	photon_dr_clusMTD->GetXaxis()->SetTitle("dr");
	ele_dr_clusMTD->GetXaxis()->SetTitle("dr");
	dr_photon_clus->GetYaxis()->SetTitle("entries");
	dr_electron_clus->GetYaxis()->SetTitle("entries");
	photon_dr_clusMTD->GetYaxis()->SetTitle("entries");
	ele_dr_clusMTD->GetYaxis()->SetTitle("entries");
	SavePlot("dr among matched mct photons and clusters",proj_full[0],(PLOTPATH+"/dr_photon_clus.pdf").c_str(),false);
	SavePlot("dr among matched mct electrons and clusters",proj_full[1],(PLOTPATH+"/dr_electron_clus.pdf").c_str(),false);
	SavePlot("dr among MTD hit and cluster for mct photons",proj_full[2],(PLOTPATH+"/photon_dr_clusMTD.pdf").c_str(),true);
	SavePlot("dr among MTD hit and cluster for mct electrons",proj_full[3],(PLOTPATH+"/ele_dr_clusMTD.pdf").c_str(),true);





}


float  matcher(float eta1,float phi1,std::vector<float>* eta2,std::vector <float>* phi2,std::pair <float,float>* pair1, std::pair <float,float>* pair2 ){

	int k;
	float dr_min = 1e6;
	for(k=0;k<phi2->size();k++){
		float dr = sqrt(pow((eta1-eta2->at(k)),2)+pow((phi1-phi2->at(k)),2));
		if(dr_min>dr){

			dr_min=dr;
			*pair1= std::make_pair(eta1,phi1);			
			*pair2= std::make_pair(eta2->at(k),phi2->at(k));
		}

	}


	std::cout << "dr_min" << dr_min << std::endl;
	if(dr_min <100 )return dr_min;		
	else return -1;



}

void SavePlot (char * titlestring, TH1D * histo, const char * filename, bool log=false){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	if (log) canvas->SetLogy();
	histo->Draw("hist");
	canvas->SaveAs(filename);

	canvas->Clear();
}

void Slicer(std::string PLOTPATH,int bin,float min, float max,std::string xaxis,TH2D *hist2D,std::string filename){

	int i;
	float x[bin],mean[bin],RMS[bin], *v=0;
	TH1D * proj[bin];
	for(i=0; i<bin;i++){
		x[i]=min+(max-min)/bin*i;
		char range[15];
		int n;
		n=sprintf(range,"_%.1f;%.1f_",x[i],x[i]+(max-min)/bin);
		proj[i]=hist2D->ProjectionX(("nhits in"+std::string(range)+xaxis+"range").c_str(),i,i+1);
		proj[i]->GetXaxis()->SetTitle(xaxis.c_str());
		proj[i]->GetYaxis()->SetTitle("entries");
		mean[i]=proj[i]->GetMean();
		RMS[i]=proj[i]->GetMeanError();

		SavePlot ("", proj[i],(PLOTPATH+"/controlplots/"+filename+std::string(range)+".pdf").c_str() ,false);


	}
	TGraphErrors* graph = new TGraphErrors(bin,x,mean,v,RMS);
	TCanvas* canvas = new TCanvas("","",600,550);
	canvas->Clear();
	graph->SetMarkerStyle(8);
	graph->GetYaxis()->SetRangeUser(0,1);
	graph->Draw("AP");
	canvas->SaveAs((PLOTPATH+"/"+filename+".pdf").c_str());
	delete canvas;
	delete graph;
	for(i=0;i<bin;i++){
		delete proj[i];
	}
}
