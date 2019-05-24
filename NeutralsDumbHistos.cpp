#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <utility>  
#include <iostream>
#include <vector>
#include <math.h>
#include "TH2.h"
#include <string>
#include "TString.h"
#include "TEllipse.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TText.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "NeutralsClass.h"
#include "NeutralsClass.cpp"
struct track {
	std::vector <float> *pt=0;
	std::vector <int> *q=0;
	std::vector <float> *phi=0;
	std::vector <float> *eta=0;
};

struct ele {
	float pt=-99;
	int q=-99;
	float phi=-99;
	float eta=-99;
	bool IsFromConversion=false;
	float vertex=0;
	float mct_truth_pt;	
};
struct pho {
	float pt=-99;
	float eta=-99;
	float phi=-99;
	int nlegs=-99;
	float convRadius=-99;
	float particleId=-99;
	struct ele ele1;
	struct ele ele2;

};

void SavePlot (char *, TH1F *, const char *, bool);
void MaterialBudget(TH1F*,TH1F*,std::string );
float TracksMatching(struct track,std::vector<struct ele> ,float, std::string,TH1F*,TH1F*);
void setStyle();

TGraph* rad_lenght(TH1F* hist);


int main (int argc, char **argv){
	TFile* file = TFile::Open(argv[1]);
	std::string outname = argv[2];
	std::cout << outname << std::endl;
	TTree* gentree = (TTree*)file->Get("gen_tree");
	TTree* tracktree= (TTree*)file->Get("tracks_tree");
	TTree* mcttree = (TTree*)file->Get("mct_tree");
	TTree* candtree= (TTree*)file->Get("cand_tree");
	TTree* dumbtree= (TTree*)file->Get("dumb_tree");
	NeutralsClass track,gen,mct,cand;
	track.Init(tracktree);
	gen.Init(gentree);
	mct.Init(mcttree);
	cand.Init(candtree);
	/*std::string outname = argv[2];
	  std::cout << outname << std::endl;
	  TTree* neutree = (TTree*)file->Get("neu_tree");
	  TTree* tracktree= (TTree*)file->Get("tracks_tree");
	  */

	setStyle();
	int i,j,k,neu_event,track_event, counter,temp,TotMatching=0;
	std::vector <float>* particleId=0;
	std::vector <float>*  pt=0, *eta=0,*phi=0;
	std::vector <float>* ele1_pt=0, *ele2_pt=0,*ele1_phi=0, *ele2_phi=0, *match_ele1=0,*match_ele2=0, *convRadius=0,*gen_DR=0;
	std::vector<int>* nlegs=0;
	struct track Eventetracks;
	struct pho photon;
	struct ele electron;
	std::string  PLOTPATH;
	std::vector<struct pho> Eventphotons;
	std::vector<struct ele> Eventelectrons;

	tracktree->SetBranchAddress("event",&track_event);
	tracktree->GetEntry(23);
	std::cout << "________________trackevent________________--" << track_event << std::endl; 

	candtree->SetBranchAddress("pt",&pt);
	candtree->SetBranchAddress("eta",&eta);
	/*	neutree->SetBranchAddress("mct_ele1_pt",&ele1_pt);
		neutree->SetBranchAddress("mct_ele2_pt",&ele2_pt);
		neutree->SetBranchAddress("mct_ele1_phi",&ele1_phi);
		neutree->SetBranchAddress("mct_ele2_phi",&ele2_phi);*/
	//mcttree->SetBranchAddress("nlegs",&nlegs);
	//	mcttree->SetBranchAddress("convRadius",&convRadius);
	//	mcttree->SetBranchAddress("event",&neu_event);
	candtree->SetBranchAddress("phi",&phi);
	tracktree->SetBranchAddress("Track_pt",&Eventetracks.pt);
	tracktree->SetBranchAddress("Track_In_phi",&Eventetracks.phi);
	//	mcttree->SetBranchAddress("PId",&particleId);
	tracktree->SetBranchAddress("Track_charge",&Eventetracks.q);
	//	neutree->SetBranchAddress("mct_eles_dr",&gen_DR);
	//	neutree->SetBranchAddress("ele1_match",&match_ele1);
	//	neutree->SetBranchAddress("ele2_match",&match_ele2);
	tracktree->SetBranchAddress("Track_In_eta",&Eventetracks.eta);


	PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events");
	system (("mkdir "+ PLOTPATH).c_str());

	TH1F * histo_pt =new TH1F("photons pt","photons pt",10,0,10);
	TH1F * histo_eta =new TH1F("photons eta","photons eta",30,-2,2);
	TH1F * histo_rad =new TH1F("photons conversion radius tot","photons conversion radius(for unconv r=0)tot ",100,0,150);
	TH1F * histoconv_pt =new TH1F(" converted photons pt ","barrel (|#eta| < 1.4) conv photons pt (R < 1.20 m)",10,0,10);
	TH1F * histoconv_eta =new TH1F("converted photons eta ","barrel (|#eta| < 1.4) convertes photons eta (R < 1.10 m)",30,-2,2);
	TH1F * histoconv_rad =new TH1F("photons conversion radius","photons conversion radius(for unconv r=0)",100,0,150);
	TH1F * histocut_eta =new TH1F("mtd converted photons eta ","barrel (|#eta| < 1.4) mtd converted photons eta (1.10 m <R < 1.22 m)",30,-2,2);
	TH1F * pt_ele =new TH1F("pt distris","comparison among generated electrons and candidates pt distribu ",10,0,10);
	//	TH1F * pt_ele2 =new TH1F("pt generated electron2 from conv ","pt generated electron2 from conv ",10,0,10);
	TH1F * cand_ele =new TH1F("pt distris","comparison among generated electrons and candidates pt distributions ",10,0,10);
	//	TH1F * cand_ele2 =new TH1F("pt  electron2 from conv ","pt generated electron2 from conv ",10,0,10);

	TH1F * matched_t =new TH1F("pr of matched tracks to  electron1 from conv ","pt matched tracsk to electron1 from conv",10,0,10);
	TH1F * matched_e =new TH1F("pt distributions "," comparison among  matched tracks and candidates pt distributions ",10,0,10);
	TH2F * tracker =new TH2F("tracker radiograpy","tracker radiography(|#eta|< 1.4)",150,-150,150,150,-150,150);// radius in cm?
	//	TH2F * radius_DR =new TH2F("convRadius gen_DR correlation","convRadius gen_DR correlation",50,0,150,10,0,0.1);// radius in cm?

	for(i=0;i<cand.fChain->GetEntries();i++){
		cand.fChain->GetEntry(i);
		mct.fChain->GetEntry(i);
		gen.fChain->GetEntry(i);
		track.fChain->GetEntry(i);
		std::cout << "HEREEEEEE" << std::endl;	
		//		if (i==0) temp = event-1;

		std::cout << "neu EVENT____" << neu_event << std::endl;
		for (j =0; j < pt->size(); j++){	
			//	if(particleId->at(j) ==4 ){}
			std::cout << "label______" << cand.label->at(j) << std::endl;
			if(cand.label->at(j)<10 && cand.label->at(j)>0){

				std::cout << "inside if___labels__________________________________"<< mct.PId->at(cand.label->at(j)) <<  std::endl;
				if (mct.PId->at(cand.label->at(j)) == 4) {
					photon.pt = pt->at(j);
					photon.eta = eta->at(j);
					photon.phi = phi->at(j);
					photon.nlegs = mct.nlegs->at(cand.label->at(j));
					photon.convRadius = mct.convRadius->at(cand.label->at(j));
					/*	photon.ele1.pt = ele1_pt->at(j);
						photon.ele1.phi = ele1_phi->at(j);
						photon.ele2.pt = ele2_pt->at(j);
						photon.ele2.phi = ele2_phi->at(j);
						*/
					Eventphotons.push_back(photon);
				}
				else if (mct.PId->at(cand.label->at(j)) ==2 ){

					electron.pt = pt->at(j);
					electron.eta = eta->at(j);
					electron.phi = phi->at(j);
					electron.q = 1;
					std::cout << "inside if______"<< mct.PId->at(cand.label->at(j)) <<  std::endl;
					/*if (match_ele1->at(j) ==1 || match_ele2->at(j) == 1){}*/ electron.IsFromConversion = true;
					electron.vertex = mct.convRadius->at(cand.label->at(j));
					electron.mct_truth_pt=mct.pt->at(cand.label->at(j)) ;
					Eventelectrons.push_back(electron);
				}
			}
			//	if(DeltaElPhoton(etrack,mct_phi,particleId,10)>0){
			//	mismatch++;
			//	std::cout <<"Photon-Electron track deltaR @ ECAL" << DeltaElPhoton(etrack,mct_phi,particleId,10) << std::endl;

			//	}
			//
			//		std::cout << "	" <<  mct_convRadius*cos(mct_phi) << "   #######" << std::endl;
			//			
			//	radius_DR->Fill(photon.convRadius, gen_DR);
		}

		std::cout << "inside for k ______" <<  mct.convRadius->size() << std::endl;
		for(k=0;k<gen.convRadius->size();k++){
				
				histo_pt->Fill(gen.pt->at(k));
				histo_eta->Fill(gen.eta->at(k));
				histo_rad->Fill(gen.convRadius->at(k));
					if (gen.eta->at(k) < 1.4  &&  gen.eta->at(k) > -1.4){ 
						if (gen.convRadius->at(k) < 110 )
						{
							histoconv_pt->Fill(gen.pt->at(k));
							histoconv_eta->Fill(gen.eta->at(k));
							histoconv_rad->Fill(gen.convRadius->at(k));
						}

						else if( (gen.convRadius->at(k) > 114 && gen.convRadius->at(k) < 118) )
						{
							histocut_eta->Fill(gen.eta->at(k));
							////		histocut_pt>Fill(pt->at(k));
						}

						tracker->Fill(gen.convRadius->at(k)*cos(gen.phi->at(k)),gen.convRadius->at(k)*sin(gen.phi->at(k)));
					}	
				
			}
		
		for (j=0; j < Eventelectrons.size(); j++){
			if (Eventelectrons[j].vertex < 100 && Eventelectrons[j].IsFromConversion == true ){

				cand_ele->Fill(Eventelectrons[j].pt);

				pt_ele->Fill(Eventelectrons[j].mct_truth_pt);
			}
		}



		/*			for(k=0; k <tracktree->GetEntries(); k++ ){
					tracktree->GetEntry(k);		
					if (track_event == neu_event){ 	Eventetracks.push_back(etrack);		
					}*/

		/*		else if (particleId = -3){

				Eventetracks.push_back(etrack);				

				}*/

		//	}
		if (pt->size() < 4 ) TotMatching += TracksMatching(Eventetracks,Eventelectrons,0.3,PLOTPATH,matched_t,matched_e);

		Eventphotons.clear();
		Eventelectrons.clear();
		Eventetracks.pt->clear();
		Eventetracks.phi->clear();
		Eventetracks.eta->clear();
		Eventetracks.q->clear();
	}
	std::cout << "Number of converted photons matched to tracks" << TotMatching << std::endl;
	histoconv_pt->GetXaxis()->SetTitle("pt(Gev/c)");
	histoconv_pt->GetYaxis()->SetTitle("N_{conv}/N_{tot}");
	histoconv_pt->SetLineColor(kRed);

	histoconv_eta->GetXaxis()->SetTitle("#eta");
	histoconv_eta->GetYaxis()->SetTitle("N_{conv}/N_{tot}");
	histoconv_eta->SetLineColor(kRed);
	histoconv_rad->GetXaxis()->SetTitle("conversion radius(cm)");
	histoconv_rad->GetYaxis()->SetTitle("N_{conv}(normalized)");
	histoconv_rad->SetLineColor(kRed);
	system(("mkdir "+PLOTPATH+"/DumbHistoCheck").c_str());
	SavePlot("Photon pt", histo_pt,(PLOTPATH+"/DumbHistoCheck/_pt.pdf").c_str(),false);
	SavePlot("Photon eta", histo_eta,(PLOTPATH+"/DumbHistoCheck/_eta.pdf").c_str(),false);
	SavePlot("Converted ph per pt", histoconv_pt,(PLOTPATH+"/DumbHistoCheck/c_pt.pdf").c_str(),false);
	SavePlot("Converted ph per eta", histoconv_eta,(PLOTPATH+"/DumbHistoCheck/c_eta.pdf").c_str(),false);
	SavePlot("photons  per convradius", histo_rad,(PLOTPATH+"/DumbHistoCheck/_rad.pdf").c_str(),false);
	SavePlot("Converted ph per converted radius", histoconv_rad,(PLOTPATH+"/DumbHistoCheck/c_rad.pdf").c_str(),false);

	histoconv_pt->Divide(histo_pt);
	histoconv_eta->Divide(histo_eta);
	histocut_eta->Divide(histo_eta);
	histoconv_rad->Scale(1/histoconv_rad->Integral());

	histoconv_eta->SetAxisRange(-1.4,1.4,"X");	
	histocut_eta->SetAxisRange(-1.4,1.4,"X");	
	MaterialBudget(histoconv_eta,histocut_eta,PLOTPATH);
	SavePlot("Fraction of converted ph per pt", histoconv_pt,(PLOTPATH+"/ratio_pt.pdf").c_str(),true);
	SavePlot("Fraction of converted ph per eta", histoconv_eta,(PLOTPATH+"/ratio_eta.pdf").c_str(),true);
	SavePlot("Fraction of mtd converted ph per eta", histocut_eta,(PLOTPATH+"/ratio_mtd_eta.pdf").c_str(),true);
	SavePlot("Fraction of converted ph per convRadius", histoconv_rad,(PLOTPATH+"/ratio_rad.pdf").c_str(),true);
	SavePlot("pt distribution for generated from conversion electron1", pt_ele,(PLOTPATH+"/gen_ele_pt.pdf").c_str(),true);
	//SavePlot("pt distribution for generated from conversion electron2", pt_ele2,(PLOTPATH+"/gen_ele2_pt.pdf").c_str(),true);
	SavePlot("pt distribution for candidate  conversion electron1", cand_ele,(PLOTPATH+"/cand_ele_pt.pdf").c_str(),true);
	//SavePlot("pt distribution for candidate from conversion electron2", cand_ele2,(PLOTPATH+"/cand_ele2_pt.pdf").c_str(),true);
	SavePlot("matched electrons  pt  distribution ", matched_e,(PLOTPATH+"/matched_ele_pt.pdf").c_str(),true);
	SavePlot("matched  tracks pt distribution ", matched_t,(PLOTPATH+"/matched_tracks_pt.pdf").c_str(),true);
	TCanvas * cand_conversions_pt  =new TCanvas("pt distribution ","comparison among electron candidates and matching generated particles", 650,550); 
	pt_ele->GetXaxis()->SetTitle("pt (Gev/c)");
	pt_ele->GetYaxis()->SetTitle("N");
	pt_ele->SetLineColor(kBlue);
	cand_ele->SetLineColor(kRed);
	pt_ele->Draw("");
	cand_ele->Draw("same");
	TLegend * legend1 = new TLegend();
	legend1->AddEntry(pt_ele,"generated electrons","l");
	legend1->AddEntry(cand_ele,"candidate electrons","l");
	legend1->Draw();
	TCanvas * tracks_conversions_pt  =new TCanvas("pt distributions ","comparison among electron candidates and matching tracks pt", 650,550); 
	matched_e->GetXaxis()->SetTitle("pt (Gev/c)");
	matched_e->GetYaxis()->SetTitle("N");
	matched_e->SetLineColor(kBlue);
	matched_e->Draw();
	matched_t->SetLineColor(kOrange);
	matched_t->Draw("same");
	TLegend * legend = new TLegend();
	legend->AddEntry(matched_e,"candidate electrons","l");
	legend->AddEntry(matched_t,"electron tracks","l");
	legend->Draw();
	cand_conversions_pt->SaveAs((PLOTPATH+"/cand_gen_conversions_pt.pdf").c_str());
	tracks_conversions_pt->SaveAs((PLOTPATH+"/cand_tracks_conversions_pt.pdf").c_str());
	TCanvas * phototrack  =new TCanvas("tracker radiography","tracker radioography", 650,550); 
	tracker->SetMarkerStyle(8);
	tracker->SetMarkerSize(.2);
	tracker->GetXaxis()->SetTitle("X (cm)");
	tracker->GetYaxis()->SetTitle("Y (cm)");
	tracker->Draw();
	phototrack->SaveAs((PLOTPATH+"/phototracker.pdf").c_str());


	delete 	histo_pt;
	delete 	histo_eta;
	delete 	histo_rad;
	delete 	histoconv_pt;
	delete 	histoconv_eta;
	delete 	histoconv_rad;
	delete 	histocut_eta;
	delete 	pt_ele;
	delete matched_t;
	delete matched_e;	
	//	pt_ele2 =new TH1F("pt generated electron2 from conv ","pt generated electron2 from conv ",10,0,10);
	delete cand_ele;
	//	cand_ele2 =new TH1F("pt  electron2 from conv ","pt generated electron2 from conv ",10,0,10);

	delete tracker; // radius in cm?
	//	TH2F * radius_DR =new TH2F("convRadius gen_DR correlation","convRadius gen_DR correlation",50,0,150,10,0,0.1);// radius in cm?
	/*	TCanvas * rad_DR  =new TCanvas("radius dr corr"," radius dr corr", 650,550); 
		radius_DR->GetXaxis()->SetTitle("convRadius (cm)");
		radius_DR->GetYaxis()->SetTitle("gen_DR ");
		radius_DR->Draw("COLZ");
		rad_DR->SaveAs((PLOTPATH+"/RAD_DR.pdf").c_str());*/
}

void SavePlot (char * titlestring, TH1F * histo, const char * filename, bool log=false){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	histo->Draw("hist");
	canvas->SaveAs(filename);

	if (log) canvas->SetLogy();
	canvas->Clear();
}


void MaterialBudget(TH1F* hist1,TH1F* hist2, std::string PATH ){

	TGraph* tracker_budget = rad_lenght(hist1);
	TGraph* mtd_budget = rad_lenght(hist2);

	//mtd_budget->Print();
	TCanvas* budget = new TCanvas("material budget","material budget", 650,550); 
	tracker_budget->SetLineColor(kRed);
	tracker_budget->SetFillColor(kRed);
	mtd_budget->SetLineColor(kBlue);
	mtd_budget->GetYaxis()->SetRangeUser(0,0.75);
	mtd_budget->GetXaxis()->SetTitle("#eta");
	mtd_budget->GetYaxis()->SetTitle("X/X_{0}");
	mtd_budget->Draw("ALP");
	tracker_budget->Draw("LP");

	TLegend * budget_legend = new TLegend(0.35,0.8,0.65,0.9);
	budget_legend->AddEntry(tracker_budget, "tracker (R < 1.10 m) ", "l");
	budget_legend->AddEntry(mtd_budget, "mtd (1.10 m < R < 1.20 m) ", "l");
	//	std::cout << "_______________HERE____________" << std::endl;
	budget_legend->Draw();
	budget->SaveAs((PATH+"/materialbudget.pdf").c_str());

}
TGraph* rad_lenght(TH1F* hist){

	int i,counter=0;

	for (i=0;i<hist->GetNbinsX();i++){
		if(hist->GetBinCenter(i)< 1.41 &&  hist->GetBinCenter(i) > -1.41){
			counter ++;

		}
	}
	float radlenght[counter],eta[counter];
	counter =0;
	for (i=0;i<hist->GetNbinsX(); i++){ 
		if(hist->GetBinCenter(i)< 1.41 &&  hist->GetBinCenter(i) > -1.41){
			counter++;
			radlenght[counter-1]=-log(1-hist->GetBinContent(i));
			eta[counter-1]=hist->GetBinCenter(i);
		}	
	}
	std::cout << "##### NBIN " << hist->GetNbinsX() << std::endl;
	TGraph * MBudget = new TGraph(counter,eta,radlenght);

	return MBudget;
}

float TracksMatching(struct track tracks, std::vector<struct ele> electrons ,  float deltaRmax,std::string PLOTPATH,TH1F* matched_t,TH1F* matched_e){
	//	float radius =  129;
	int match=0;
	//	float cand_X,cand_Y, e_X, e_Y;
	int i,j;



	std::cout << "_____nphotons______ntracks_______" << electrons.size() << "________________" <<tracks.pt->size() << std::endl;
	for(i=0; i<electrons.size();i++){
		for(j=0; j < tracks.pt->size();j++){
			std:: cout << "track pt____________" << tracks.pt->at(j) << std::endl;
			if(electrons[i].vertex < 100 && electrons[j].IsFromConversion== true){
				//	cand_X = ElectronProp(photons[i].convRadius,photons[i].phi,photons[i].ele1.phi).first;
				//	cand_Y = ElectronProp(photons[i].convRadius,photons[i].phi,photons[i].ele1.phi).second;
				//	e_X = ElectronProp(0,0,tracks[j].phi).first;
				//	e_Y = ElectronProp(0,0,tracks[j].phi).second;
				float deltaR = sqrt(pow((electrons[i].phi-tracks.phi->at(j)),2)+pow((electrons[i].eta-tracks.eta->at(j)),2));
				std::cout << "DeltaR = " << deltaR << std::endl;
				if (deltaR < deltaRmax) {
					match ++;
					std::cout << "converted electron-track matching succeded"<<std::endl;
					std::cout << "DeltaR = " << deltaR << std::endl;
					matched_t->Fill(tracks.pt->at(j));
					matched_e->Fill(electrons[i].pt);
				}


			}	
		}
	}


	if (match > 0) return match;
	else return 0;
}
void setStyle() {


	// set the TStyle
	TStyle* style = new TStyle("DrawBaseStyle", "");
	style->SetCanvasColor(0);
	style->SetPadColor(0);
	style->SetFrameFillColor(0);
	style->SetStatColor(0);
	style->SetOptStat(0);
	style->SetOptFit(0);
	style->SetTitleFillColor(0);
	style->SetCanvasBorderMode(0);
	style->SetPadBorderMode(0);
	style->SetFrameBorderMode(0);
	style->SetPadBottomMargin(0.12);
	style->SetPadLeftMargin(0.12);
	style->cd();
	// For the canvas:
	style->SetCanvasBorderMode(0);
	style->SetCanvasColor(kWhite);
	style->SetCanvasDefH(600); //Height of canvas
	style->SetCanvasDefW(600); //Width of canvas
	style->SetCanvasDefX(0); //POsition on screen
	style->SetCanvasDefY(0);
	// For the Pad:
	style->SetPadBorderMode(0);
	style->SetPadColor(kWhite);
	style->SetPadGridX(false);
	style->SetPadGridY(false);
	style->SetGridColor(0);
	style->SetGridStyle(3);
	style->SetGridWidth(1);
	// For the frame:
	style->SetFrameBorderMode(0);
	style->SetFrameBorderSize(1);
	style->SetFrameFillColor(0);
	style->SetFrameFillStyle(0);
	style->SetFrameLineColor(1);
	style->SetFrameLineStyle(1);
	style->SetFrameLineWidth(1);
	// Margins:
	style->SetPadTopMargin(0.10);
	style->SetPadBottomMargin(0.14);//0.13);
	style->SetPadLeftMargin(0.16);//0.16);
	style->SetPadRightMargin(0.1);//0.02);
	// For the Global title:
	style->SetOptTitle(0);
	style->SetTitleFont(42);
	style->SetTitleColor(1);
	style->SetTitleTextColor(1);
	style->SetTitleFillColor(10);
	style->SetTitleFontSize(0.05);
	// For the axis titles:
	style->SetTitleColor(1, "XYZ");
	style->SetTitleFont(42, "XYZ");
	style->SetTitleSize(0.05, "XYZ");
	style->SetTitleXOffset(1.15);//0.9);
	style->SetTitleYOffset(1.5); // => 1.15 if exponents
	// For the axis labels:
	style->SetLabelColor(1, "XYZ");
	style->SetLabelFont(42, "XYZ");
	style->SetLabelOffset(0.007, "XYZ");
	style->SetLabelSize(0.045, "XYZ");
	// For the axis:
	style->SetAxisColor(1, "XYZ");
	style->SetStripDecimals(kTRUE);
	style->SetTickLength(0.03, "XYZ");
	style->SetNdivisions(510, "XYZ");
	style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
	style->SetPadTickY(1);
	// for histograms:
	style->SetHistLineColor(1);
	// for the pallete
	Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(5, stops, red, green, blue, 100);
	style->SetNumberContours(100);

	style->cd();

}

/*float TracksMatching(std::vector<struct track> tracks, std::vector<struct pho> photons ,  float deltaRmax,std::string PLOTPATH){
//	float radius =  129;
int match=0;
//	float cand_X,cand_Y, e_X, e_Y;
int i,j;

TH1F * matched_e1 =new TH1F("pr of matched tracks to  electron1 from conv ","pt matched tracsk to electron1 from conv",10,0,10);
TH1F * matched_e2 =new TH1F("pt of matched tracke to electron2 from conv ","pt of matched tracks to  generated electron2 from conv ",10,0,10);

std::cout << "_____nphotons______ntracks_______" << photons.size() << "________________" <<tracks.size() << std::endl;
for(i=0; i<photons.size();i++){
for(j=0; j < tracks.size();j++){
std:: cout << "track pt____________" << tracks[j].pt << std::endl;
//	cand_X = ElectronProp(photons[i].convRadius,photons[i].phi,photons[i].ele1.phi).first;
//	cand_Y = ElectronProp(photons[i].convRadius,photons[i].phi,photons[i].ele1.phi).second;
//	e_X = ElectronProp(0,0,tracks[j].phi).first;
//	e_Y = ElectronProp(0,0,tracks[j].phi).second;
float deltaR = sqrt(pow((photons[i].ele1.phi-tracks[j].phi),2)+pow((photons[i].ele1.eta-tracks[j].eta),2));
std::cout << "DeltaR = " << deltaR << std::endl;
if (deltaR < deltaRmax) {
match ++;
std::cout << "converted electron-track matching succeded"<<std::endl;
std::cout << "DeltaR = " << deltaR << std::endl;
matched_e1->Fill(tracks[j].pt);
}
else {

//	cand_X = ElectronProp(photons[i].convRadius,photons[i].phi,photons[i].ele2.phi).first;
//	cand_Y = ElectronProp(photons[i].convRadius,photons[i].phi,photons[i].ele2.phi).second;

float deltaR = sqrt(pow((photons[i].ele1.phi-tracks[j].phi),2)+pow((photons[i].ele1.eta-tracks[j].eta),2));
if (deltaR < deltaRmax) {
match ++;
std::cout << "converted electron-track matching succeded"<<std::endl;
std::cout << "DeltaR = " << deltaR << std::endl;
matched_e2->Fill(tracks[j].pt);
}					
}



}
}

SavePlot("matched 1  tracks pt  distribution ", matched_e1,(PLOTPATH+"/tracks_ele1_pt.pdf").c_str(),true);
SavePlot("matched  2 tracks pt distribution ", matched_e2,(PLOTPATH+"/tracks_ele2_pt.pdf").c_str(),true);
if (match > 0) return match;
else return 0;
}*/
