
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


void SavePlot (char *, TH1D *, const char * , bool );


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


	std::string  PLOTPATH;

		
	TH1D * mct_unconv_pt = new TH1D("","sim unconverted photons pt",10,0,10);
	TH1D * mct_electrons_pt = new TH1D("","sim electrons from conv pt",10,0,10);
	TH1D * mct_unconv_eta = new TH1D("","sim unconverted photons pt",10,-1.5,1-5);
	TH1D * mct_electrons_eta = new TH1D("","sim electrons from conv pt",10,-1.5,1.5);
	TH1D * cand_unconv_pt = new TH1D("","candidate unconverted photons pt",10,0,10);
	TH1D * cand_electrons_pt = new TH1D("","cand electrons from conv pt",10,0,10);
	TH1D * cand_unconv_eta = new TH1D("","cand unconverted photons pt",10,-1.5,1-5);
	TH1D * cand_electrons_eta = new TH1D("","cand electrons from conv pt",10,-1.5,1.5);
	TH2D * mct_eta_convRadius = new TH2D("","2d hist",15,-3,3,15,0,120);
	int i, j,k,p,c,cn;
	



	PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events");
	system (("mkdir "+ PLOTPATH).c_str());
	for(i=0;i<mct.fChain->GetEntries();i++){
		mct.fChain->GetEntry(i);
		cand.fChain->GetEntry(i);
		gen.fChain->GetEntry(i);
		
		std::cout << "event__" <<mct.event << std::endl;
		if(mct.pt->size()>0){
		for(p =0; p<mct.pt->size();p++){
			

			if(mct.PId->at(p)==4){ 
			mct_unconv_pt->Fill(mct.pt->at(p));
			mct_unconv_eta->Fill(mct.eta->at(p));
				
		//		for(j=0;j<cand.fChain->GetEntries();j++){
		//			cand.fChain->GetEntry(j);
					if(cand.event == mct.event && cand.pt->size()>0){
					for(c =0; c<cand.pt->size();c++){
						if(cand.plabel->at(c)==mct.plabel->at(p) && cand.label->at(c)==mct.label->at(p)){
						
						cand_unconv_pt->Fill(cand.pt->at(c));
						cand_unconv_eta->Fill(cand.eta->at(c));

						}

}

}
//}

}
			else if(mct.PId->at(p)==2){
			mct_eta_convRadius->Fill(gen.phi->at(mct.label->at(p)),mct.convRadius->at(p));
			mct_electrons_pt->Fill(mct.pt->at(p));
			mct_electrons_eta->Fill(mct.eta->at(p));
	//			for(j=0;j<cand.fChain->GetEntries();j++){
	//				cand.fChain->GetEntry(j);
					if(cand.event == mct.event && cand.pt->size()>0){
					for(cn =0; cn<cand.pt->size();cn++){
//					std::cout << "here index" << cn <<  " <" << cand.pt->size()<<  std::endl;
						if(cand.plabel->at(cn) == mct.plabel->at(p) && cand.label->at(cn)==mct.label->at(p)){
						
						cand_electrons_pt->Fill(cand.pt->at(cn));
						cand_electrons_eta->Fill(cand.eta->at(cn));

						}
				
}

}
//}


	}
}
}


}
	system (("mkdir "+ PLOTPATH+"/RecoCheck").c_str());
	SavePlot("sim unconv photons  pt", mct_unconv_pt,(PLOTPATH+"/RecoCheck/sim_unconv_pt.pdf").c_str(),false);
	SavePlot("sim unconv photons  eta", mct_unconv_eta,(PLOTPATH+"/RecoCheck/sim_unconv_eta.pdf").c_str(),false);
	SavePlot("sim electrons  pt", mct_electrons_pt,(PLOTPATH+"/RecoCheck/sim_electrons_pt.pdf").c_str(),false);
	SavePlot("sim electrons eta", mct_electrons_eta,(PLOTPATH+"/RecoCheck/sim_electrons_eta.pdf").c_str(),false);
	SavePlot("cand unconv photons  pt", cand_unconv_pt,(PLOTPATH+"/RecoCheck/cand_unconv_pt.pdf").c_str(),false);
	SavePlot("cand unconv photons  eta", cand_unconv_eta,(PLOTPATH+"/RecoCheck/cand_unconv_eta.pdf").c_str(),false);
	SavePlot("cand electrons  pt", cand_electrons_pt,(PLOTPATH+"/RecoCheck/cand_electrons_pt.pdf").c_str(),false);
	SavePlot("cand electrons eta", cand_electrons_eta,(PLOTPATH+"/RecoCheck/cand_electrons_eta.pdf").c_str(),false);
	TCanvas * super_pt = new  TCanvas ("","superposition of pt distributions",600, 550);
	TLegend* leg = new TLegend();
	mct_electrons_pt->GetXaxis()->SetTitle("pt(GeV/c)");
	mct_electrons_pt->GetYaxis()->SetTitle("N normalized");
	mct_electrons_pt->Divide(mct_electrons_pt);
	mct_electrons_pt->SetLineColor(kRed);
	mct_electrons_pt->Draw();
	//cand_electrons_pt->Draw("same");
	//leg->AddEntry(mct_electrons_pt,"simulated electrons");
	leg->AddEntry(cand_electrons_pt,"candidate electrons");
	leg->Draw();
	super_pt->SaveAs((PLOTPATH+"/RecoCheck/super_pt.pdf").c_str());
	TCanvas * his2D = new TCanvas("","",600,550);
	mct_eta_convRadius->Draw("COLZ");
	his2D->SaveAs((PLOTPATH+"/RecoCheck/convRadius_eta.pdf").c_str());
	
}





void SavePlot (char * titlestring, TH1D * histo, const char * filename, bool log=false){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	histo->Draw("hist");
	canvas->SaveAs(filename);

	if (log) canvas->SetLogy();
	canvas->Clear();
}
