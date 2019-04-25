
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
#include "TStyle.h"
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

struct photon{

	float convRadius;
	float eta;
	float phi;


};

void SavePlot (char *, TH1D *, const char * , bool );

void Slicer(std::string, int bin,float min, float max,std::string xaxis,TH2D *hist2D,std::string filename);
void Trackplotter(std::string PLOTPATH,struct tr track,struct photon p1,struct photon p2,int event);
void setStyle();

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
	struct tr tracks;
	struct photon p1,p2;
	float min_r, max_r;
	float min_pt, max_pt;
	int r_bin,pt_bin,bin_hits;
	float hits_min,hits_max;
	std::string  PLOTPATH;
	r_bin=20;
	pt_bin=10;
	bin_hits=20;
	min_r=0;
	min_pt=0;
	hits_min=0;
	max_r=15;
	max_pt=10;
	hits_max=20;

	setStyle();

	TH1D * convRad = new TH1D("","converted photon convradius",r_bin,min_r,max_r);
	TH1D * tracks_convRad = new TH1D("","photon matched tracks  convradius",r_bin,min_r,max_r);
	TH1D * tracks_x_event = new TH1D("","tracks x event ",10,0,10);

	TH2D * nhits_pt = new TH2D("","nhits vs track  pt",pt_bin,min_pt,max_pt,bin_hits,hits_min,hits_max);
	TH2D * nhits_convRad = new TH2D("","nhits vs convRadius",r_bin,min_r,max_r,bin_hits,hits_min,hits_max);
	TH2D * eta_convRad = new TH2D("","eta vs convRadius",20,0,130,10,-1.4,1.4);


	PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events");
	system(("mkdir"+PLOTPATH+"/tracks_analysis").c_str());

	PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events/tracks_analysis");
	system(("mkdir "+PLOTPATH+"/tracks").c_str());

	int counter =0;
	for(i=0;i<track.fChain->GetEntries();i++){

		track.fChain->GetEntry(i);
		gen.fChain->GetEntry(i);

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

		p1.convRadius=gen.convRadius->at(0);
		p1.phi=gen.phi->at(0);
		p1.eta=gen.eta->at(0);
		p2.convRadius=gen.convRadius->at(1);
		p2.phi=gen.phi->at(1);
		p2.eta=gen.eta->at(1);



		std::cout << "event____" << track.event <<  std::endl;
		if (track.Track_pt->size()>0 && track.event%20==0)	Trackplotter(PLOTPATH+"/tracks",tracks,p1,p2,track.event);
		tracks_x_event->Fill(tracks.pt->size());
		//	if(gen.pt->size()>0){
		//	for(int j; j<gen.pt->size();j++){
		if(gen.convRadius->at(0) >0. && (p1.eta<(1+0.4) && p1.eta>-(1+0.4)) && p1.convRadius<100){

			convRad->Fill(p1.convRadius);
			convRad->Fill(p1.convRadius);
			eta_convRad->Fill(p1.convRadius,p1.eta);
		}		

		if( gen.convRadius->at(1) >0. &&  (p2.eta<(1+.4) && p2.eta>-(1+.4)) &&  p2.convRadius<100){

			convRad->Fill(p2.convRadius);
			convRad->Fill(p2.convRadius);
			eta_convRad->Fill(p2.convRadius,p2.eta);
		}		
		//	}
		//	}
		if(tracks.pt->size()>0){
			for(int j=0; j<tracks.pt->size();j++){
				std::cout << "angles    " << tracks.In_phi->at(j) << "__   _"<<p1.phi <<"   difference   " <<  abs(tracks.In_phi->at(j)-p1.phi)<< " < " << M_PI/2 << std::endl; 
				if(abs(tracks.In_phi->at(j)-p1.phi)< M_PI/2 /*abs(sqrt(pow(tracks.In_x->at(it),2)+pow(tracks.In_y->at(it),2))-gen.convRadius->at(0))<20*/){
					if(gen.convRadius->at(0) >0. && abs(tracks.In_eta->at(j))< 1.4 && gen.convRadius->at(0)<100){
						tracks_convRad->Fill(gen.convRadius->at(0));
						counter++;	
						if (gen.convRadius->at(0)>50)std::cout << "check_0  " << gen.convRadius->at(0) << "__" << tracks.nhits->at(j) << std::endl;
						nhits_convRad->Fill(gen.convRadius->at(0),tracks.nhits->at(j));

					}

				}else {

					if(gen.convRadius->at(1) >0. && abs(tracks.In_eta->at(j))< 1.4 && gen.convRadius->at(1)<100){

						tracks_convRad->Fill(gen.convRadius->at(1));
						counter++;
						if (gen.convRadius->at(1)>50)std::cout << "check_1  " << gen.convRadius->at(1) << "__" << tracks.nhits->at(j) << std::endl;
						nhits_convRad->Fill(gen.convRadius->at(1),tracks.nhits->at(j));

					}


				}		



				nhits_pt->Fill(tracks.pt->at(j),tracks.nhits->at(j));

			}


		}



	}
	std::cout << "entries radius " << counter << std::endl;
	system(("mkdir "+ PLOTPATH+"/controlplots").c_str());
	Slicer(PLOTPATH,pt_bin,min_pt,max_pt,"pt (Gev/c)",nhits_pt,"nhits_pt");
	Slicer(PLOTPATH,r_bin,min_r,max_r," nhits",nhits_convRad,"nhits_radius");
	Slicer(PLOTPATH,r_bin,min_r,max_r,"nhits",eta_convRad,"nhits_eta");
	tracks_convRad->GetXaxis()->SetTitle("convRadius(cm)");
	tracks_convRad->GetYaxis()->SetTitle("N_{tracks}");
	convRad->GetXaxis()->SetTitle("convRadius(cm)");
	convRad->GetYaxis()->SetTitle("N_{electrons}");
	SavePlot ("tracks in bin of convradius", tracks_convRad,(PLOTPATH+"/tracks_radius.pdf").c_str() ,false);
	tracks_convRad->Divide(convRad);
	tracks_convRad->GetXaxis()->SetTitle("convRadius(cm)");
	tracks_convRad->GetYaxis()->SetTitle("N_{tracks}/N_{electrons}");
	tracks_convRad->SetTitle("Ratio of tracks over electrons in bins of convRadius");
	tracks_x_event->GetXaxis()->SetTitle("N_{tracks}");
	tracks_x_event->GetYaxis()->SetTitle("entries");
	tracks_x_event->SetTitle("Number of tracks per event");
	for(i=0;i<r_bin;i++){
	if(tracks_convRad->GetBinContent(i)!=0)
	tracks_convRad->SetBinError(i,sqrt((tracks_convRad->GetBinContent(i))/pow(convRad->GetBinContent(i),2)+pow(tracks_convRad->GetBinContent(i),2)/pow(convRad->GetBinContent(i),3)));

	}
	SavePlot ("tracks/electrons in bin of convradius", tracks_convRad,(PLOTPATH+"/ratiotracks_radius.pdf").c_str() ,false);
	SavePlot ("tracks x event", tracks_x_event,(PLOTPATH+"/tracks_event.pdf").c_str() ,false);
	convRad->Scale(1/convRad->Integral());
	SavePlot ("electrons in bin of convradius",convRad,(PLOTPATH+"/electrons_radius.pdf").c_str() ,false);


}

void SavePlot (char * titlestring, TH1D * histo, const char * filename, bool log=false){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	if (log) canvas->SetLogy();
	histo->SetLineWidth(3);
	histo->SetMarkerStyle(8);
	histo->Draw("PE1");
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
		proj[i]=hist2D->ProjectionY(("nhits in"+std::string(range)+xaxis+"range").c_str(),i,i+1);
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
	graph->GetYaxis()->SetRangeUser(0,20);
	graph->Draw("AP");
	canvas->SaveAs((PLOTPATH+"/"+filename+".pdf").c_str());
	delete canvas;
	delete graph;
	for(i=0;i<bin;i++){
		delete proj[i];
	}
}

void Trackplotter(std::string PLOTPATH,struct tr track,struct photon p1,struct photon p2,int event){
	TEllipse * ECAL= new TEllipse(0,0,129,129);	
	TEllipse * MTD= new TEllipse(0,0,120,120);	
	TMarker* inpoint[track.pt->size()];
	TMarker* calopoint[track.pt->size()];
	TMarker* photon1, photon2;
	TArrow* indir[track.pt->size()];
	TLegend* legend= new TLegend();
	Color_t colors[4]={kBlue-6,kRed-6,kOrange+3,kGreen+3};
	TText * pt[track.pt->size()];
	int i;
	for(i=0;i<track.pt->size();i++){
		if(sqrt(pow(track.ECAL_x->at(i),2)+pow(track.ECAL_y->at(i),2))>127){
			indir[i] = new TArrow();
			inpoint[i]= new TMarker(track.In_x->at(i),track.In_y->at(i),5);
			inpoint[i]->SetMarkerColor(colors[i]);
			indir[i]->SetLineColor(colors[i]);
			indir[i]->SetFillColor(colors[i]);
			legend->AddEntry(inpoint[i],("track n"+std::to_string(i)).c_str(),"p");
			pt[i]=new TText();
			pt[i]->SetTextColor(colors[i]);
			pt[i]->SetTextSize(0.04);
			calopoint[i]= new TMarker(track.ECAL_x->at(i),track.ECAL_y->at(i),4);
			calopoint[i]->SetMarkerColor(colors[i]);
		}
	}
	TH2D * plotter = new TH2D("","",10,-140,140,10,-140,140);
	TCanvas * track_plot = new TCanvas(".","trackplotter",600,550);
	ECAL->SetLineColor(kGray+2);
	ECAL->SetLineWidth(4);
	MTD->SetLineColor(kGray+2);
	MTD->SetLineWidth(4);
	plotter->GetXaxis()->SetTitle("X(cm)");
	plotter->GetYaxis()->SetTitle("Y(cm)");
	plotter->Draw();
	ECAL->Draw("same");
	MTD->Draw("same");
	float Y=60;
	for(i=0;i<track.pt->size();i++){
		if(sqrt(pow(track.ECAL_x->at(i),2)+pow(track.ECAL_y->at(i),2))>127){
			std::cout << "Hereeee" << std::endl;
			inpoint[i]->Draw("same");
			char range[10];
			int n;
			n=sprintf(range,"%.2f",track.pt->at(i));
			pt[i]->DrawText(20,Y,("pt = "+std::string(range)+"Gev/c").c_str());
			indir[i]->DrawArrow(track.In_x->at(i),track.In_y->at(i),track.In_x->at(i)+15*cos(track.pt_in_phi->at(i)),track.In_y->at(i)+15* sin(track.pt_in_phi->at(i)),0.01,"|>");
			calopoint[i]->Draw("same");
			Y-=15;
		}
	}
	legend->Draw();
	track_plot->SaveAs((PLOTPATH+"/event"+std::to_string(event)+".pdf").c_str());
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
	style->SetPadRightMargin(0.04);//0.02);
	// For the Global title:
	style->SetOptTitle(1);
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





