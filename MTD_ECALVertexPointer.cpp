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
#include "TStyle.h"
#include "TMarker.h"
#include "TText.h"
#include "TF1.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLegend.h"
#include "NeutralsClass.h"
#include "NeutralsClass.cpp"

struct pos{

    float x;
    float y;
    float z;
    float phi;
    float eta;



};

float matcher(float,float,std::vector<float>*,std::vector <float>* ,std::vector <float>* ,std::vector <float>* ,std::vector <float>* ,struct pos*, std::pair <float,float>* );

void SavePlot (char * , TH1D *, const char*, bool );

void Slicer(std::string ,int ,float, float ,std::string ,TH2D *,std::string );

void setStyle();

float VertexDistance(TVector3* ,float , float , float , float , bool );

void Plotter(float ,float , float, float , TVector3 * ,float , float );

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
	setStyle();

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
	TH1D* photon_vert_dr = new TH1D("","distance of the real vertex wrt to the ECAL_MTD line",25,0,10);

	for(i=0;i<cand.fChain->GetEntries();i++){
		cand.fChain->GetEntry(i);
		mct.fChain->GetEntry(i);
		
		for(j=0;j<mct.pt->size();j++){
			std::pair <float,float>  matched_mct;
			std::pair <float,float>  pair_clus;
			struct pos matched_clus;
			struct pos matched_MTDhit;
			
			TVector3 *vertex=new TVector3();
			vertex->SetX(mct.primary_x->at(j));
			vertex->SetY(mct.primary_y->at(j));
			vertex->SetZ(mct.primary_z->at(j));
			dr = matcher(mct.eta->at(j),mct.phi->at(j),cand.ecaleta, cand.ecalphi,cand.ecalx,cand.ecaly,cand.ecalz,&matched_clus,&matched_mct);
				
			if(dr != -1){
				label = mct.PId->at(j);
				if (label == 4 ) dr_photon_clus->Fill(dr,mct.pt->at(j));
				else if(label == 2) dr_electron_clus->Fill(dr,mct.pt->at(j));
			}
			dr = matcher(matched_clus.eta,matched_clus.phi,cand.MTDeta, cand.MTDphi, cand.MTDx, cand.MTDy, cand.MTDz,&matched_MTDhit,&pair_clus);
			std::cout  << "_____" << sqrt(matched_clus.x*matched_clus.x+matched_clus.y*matched_clus.y) << std::endl;
			std::cout  << "_____" << sqrt(matched_MTDhit.x*matched_MTDhit.x+matched_MTDhit.y*matched_MTDhit.y) << std::endl;
			if(dr !=-1){
				if (label == 4 ){
				 photon_dr_clusMTD->Fill(dr,mct.pt->at(j));
				std::cout << "hereee" << std::endl;
					if(dr < 0.01){
						
						float dxy = VertexDistance(vertex,matched_MTDhit.x,matched_clus.x,matched_MTDhit.y,matched_clus.y,false);	
						float dzy = VertexDistance(vertex,matched_MTDhit.z,matched_clus.z,matched_MTDhit.y,matched_clus.y,true);	
						if(dxy != -1 && dzy != -1){	
						std::cout << "displacements" << dxy << "____" << dzy << std::endl;
						photon_vert_dr->Fill(sqrt(dxy*dxy + dzy*dzy));									
						}
					}
				}
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
	proj_full[0]->GetXaxis()->SetTitle("dr");
	proj_full[1]->GetXaxis()->SetTitle("dr");
	proj_full[2]->GetXaxis()->SetTitle("dr");
	proj_full[3]->GetXaxis()->SetTitle("dr");
	proj_full[0]->GetYaxis()->SetTitle("entries");
	proj_full[1]->GetYaxis()->SetTitle("entries");
	proj_full[2]->GetYaxis()->SetTitle("entries");
	proj_full[3]->GetYaxis()->SetTitle("entries");
	SavePlot("dr among matched mct photons and clusters",proj_full[0],(PLOTPATH+"/dr_photon_clus.pdf").c_str(),false);
	SavePlot("dr among matched mct electrons and clusters",proj_full[1],(PLOTPATH+"/dr_electron_clus.pdf").c_str(),false);
	SavePlot("dr among MTD hit and cluster for mct photons",proj_full[2],(PLOTPATH+"/photon_dr_clusMTD.pdf").c_str(),false);
	SavePlot("dr among MTD hit and cluster for mct electrons",proj_full[3],(PLOTPATH+"/ele_dr_clusMTD.pdf").c_str(),false);
	SavePlot("distance between real vertex and MTD_ECAL line",photon_vert_dr,(PLOTPATH+"/photon_vert_dr.pdf").c_str(),false);





}


float  matcher(float eta1,float phi1,std::vector<float>* eta2,std::vector <float>* phi2,std::vector <float>* x2,std::vector <float>* y2,std::vector <float>* z2,struct pos* pos, std::pair <float,float>* pair ){

	int k;
	float dr_min = 1e6;
	std::cout << "in matcherrr" << std::endl;
	for(k=0;k<phi2->size();k++){
		float dr = sqrt(pow((eta1-eta2->at(k)),2)+pow((phi1-phi2->at(k)),2));
		if(dr_min>dr){

			dr_min=dr;
			*pair= std::make_pair(eta1,phi1);			
			(*pos).x= x2->at(k);
			(*pos).y= y2->at(k);
			(*pos).z= z2->at(k);
			(*pos).phi= phi2->at(k);
			(*pos).eta= eta2->at(k);
	
	std::cout << "dr_min" << dr_min << std::endl;
		}

	}


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
	graph->GetXaxis()->SetTitle("pt(Gev/c)");
	graph->GetYaxis()->SetTitle("mean dr");
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

float VertexDistance(TVector3 *vert,float x1, float x2, float y1, float y2, bool longitudinal = false){
	float m1;
	float c1,c2;
	float x_int, y_int;
	float dr=-1;
	float x_v;
	
	m1 = (y2-y1)/(x2-x1);
	c1 = y1-m1*x1;
	if(longitudinal) x_v = vert->z();
	else x_v = vert->x();
	c2= vert->y()+x_v/m1;
	x_int = m1*(c2-c1)/(m1*m1+1);
	y_int = c2-(c2-c1)/(m1*m1+1);
	
	std::cout << "m1" << m1 <<"__" <<  c1 << std::endl;
	std::cout << "coord_x" << x1 <<"__" <<  x2 << std::endl;
	std::cout << "coord_y" << y1 << "__"<< y2 << std::endl;
	std::cout << "coord" << x_int << vert->x() << std::endl;
	if((y1 >0 && m1*x1 >0) || (y1< 0 && m1*x1<0)){
	if(longitudinal==true ){

	dr = sqrt((y_int-vert->y())*(y_int-vert->y())+(x_int-vert->z())*(x_int-vert->z()));
	std::cout << "dr " << dr << std::endl;
	return dr ;

	}else{

	dr = sqrt((y_int-vert->y())*(y_int-vert->y())+(x_int-vert->x())*(x_int-vert->x()));
	if(m1>10 && m1<30 )Plotter(x2,y2,x1,y1,vert,m1,c1);
	std::cout << "dr " << dr << std::endl;
	return dr ;

	}		

	}else return dr;


}

void Plotter(float ecal_x,float ecal_y, float mtd_x, float mtd_y, TVector3 * vert,float m1, float c1){
	TEllipse * ECAL= new TEllipse(0,0,129,129);	
	TEllipse * MTD= new TEllipse(0,0,117,117);	
	TMarker * ECALpoint;
	TMarker* MTDpoint;
	TMarker* vertex;
	TF1* MTD_ECAL;
	ECALpoint= new TMarker(ecal_x,ecal_y,kCircle);
	MTDpoint= new TMarker(mtd_x,mtd_y,kCircle);
	vertex= new TMarker(vert->x(),vert->y(),kMultiply);
	MTD_ECAL = new TF1("pointer","x*[0]+[1]",-120,120);
	MTD_ECAL->SetParameter(0,m1);
	MTD_ECAL->SetParameter(1,c1);	
	TH2D * plotter = new TH2D("","",10,-140,140,10,-140,140);
	TCanvas * vertex_plot = new TCanvas(".","trackplotter",600,550);
	ECAL->SetLineColor(kGray+2);
	ECAL->SetLineWidth(4);
	MTD->SetLineColor(kGray+2);
	MTD->SetLineWidth(4);
	plotter->GetXaxis()->SetTitle("X(cm)");
	plotter->GetYaxis()->SetTitle("Y(cm)");
	plotter->Draw();
	ECAL->Draw("same");
	MTD->Draw("same");
	vertex->Draw("same");
	ECALpoint->Draw("same");
	MTDpoint->Draw("same");
	MTD_ECAL->Draw("same");
	vertex_plot->SaveAs((" ~/CMSSW_10_4_0_mtd5_prova/plots/testAll_clusevents/VertexPointer/vertex_"+std::to_string(mtd_x)+".pdf").c_str());




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

