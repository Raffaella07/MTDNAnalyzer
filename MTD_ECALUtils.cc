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




struct pos{

    float x;
    float y;
    float z;
    float phi;
    float eta;
    float pt;



};

float  matcher(float eta1,float phi1,std::vector<float>* eta2,std::vector <float>* phi2,std::vector <float>* x2,std::vector <float>* y2,std::vector <float>* z2,std::vector <float>* pt2,struct pos* pos, std::pair <float,float>* pair ){

	int k;
	float dr_min = 1e6;
	//std::cout << "in matcherrr" << std::endl;
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
		if(k<pt2->size())	(*pos).pt= pt2->at(k);
	
//	std::cout << "dr_min" << dr_min  << std::endl;
		}

	}


	if(dr_min <100 )return dr_min;		
	else return -1;



}

void SavePlot (char * titlestring, TH1D * histo, const char * filename, bool log=false, TF1* fit= NULL){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	if (log) canvas->SetLogy();
	histo->GetXaxis()->SetMaxDigits(2);
	histo->Draw("hist");
	if(fit !=NULL){
	TPaveText* fitpar = new TPaveText(.6,.65,.95,.85,"NDC");
	fitpar->SetTextFont(62);
	fitpar->AddText("Gaussian fit parameters");
	fitpar->AddLine(0,.8,1,.8);
	char strng[15];
	int n = sprintf(strng, "%.2E #pm %.2E",fit->GetParameter(0),fit->GetParError(0));
	
	fitpar->AddText(("Constant  "+std::string(strng)).c_str());
	
	n = sprintf(strng, "%.2E #pm %.2E",fit->GetParameter(1),fit->GetParError(1));
	fitpar->AddText(("Mean  "+std::string(strng)).c_str());
	n = sprintf(strng, "%.2E #pm %.2E",fit->GetParameter(2),fit->GetParError(2));
	fitpar->AddText(("#sigma  "+std::string(strng)).c_str());
	fit->SetLineColor(kRed);
	fit->Draw("SAME");
	fitpar->Draw("SAME");

	}
	canvas->SaveAs(filename);
	canvas->Clear();
}
void SavePlot2D (char * titlestring, TH2D * histo, const char * filename, bool log=false){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	TH2D * plotter =new TH2D("-","",20,-20,20,20,0.0,0.4);
	if (log) canvas->SetLogy();
	plotter->Draw();
	histo->Draw("sameCOLZ");
	canvas->SaveAs(filename);

	canvas->Clear();
}

void Slicer(std::string PLOTPATH,int bin,float min, float max,std::string xaxis,TH2D *hist2D,std::string filename){

	int i;
	float x[bin],mean[bin],RMS[bin], *v=0;
	TH1D * proj[bin];
	for(i=0; i<bin;i++){
		x[i]=min +(max-min)/bin*i;
		char range[15]="";
		int n;
		n=sprintf(range,"%.3f;%.3f",x[i],x[i]+(max-min)/bin);
		proj[i]=hist2D->ProjectionX(("nhits in"+std::to_string(i)+"range").c_str(),i,i+1);
		proj[i]->GetXaxis()->SetTitle(xaxis.c_str());
		proj[i]->GetYaxis()->SetTitle("entries");
		mean[i]=proj[i]->GetMean();
		RMS[i]=proj[i]->GetMeanError();

		SavePlot ("", proj[i],(PLOTPATH+"/controlplots/"+filename+std::string(range)+"_.pdf").c_str() ,false);


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

void Plotter(float ecal_x,float ecal_y, float mtd_x, float mtd_y, TVector3 * vert,float m1, float c1,std::string PLOTPATH=""){
	TEllipse * ECAL= new TEllipse(0,0,129,129);	
	TEllipse * MTD= new TEllipse(0,0,117,117);	
	TMarker * ECALpoint;
	TMarker* MTDpoint;
	TMarker* vertex;
	TF1* MTD_ECAL;
	ECALpoint= new TMarker(ecal_x,ecal_y,kFullCircle);
	MTDpoint= new TMarker(mtd_x,mtd_y,kCircle);
	vertex= new TMarker(vert->x(),vert->y(),kFullDoubleDiamond);
	MTD_ECAL = new TF1("pointer","x*[0]+[1]",-120,120);
	MTD_ECAL->SetParameter(0,m1);
	MTD_ECAL->SetParameter(1,c1);	
	MTD_ECAL->SetLineColor(kAzure+1);	
	MTD_ECAL->SetLineWidth(5);	
	MTD_ECAL->SetLineStyle(9);	
	TH2D * plotter = new TH2D("","",10,-160,160,10,-160,160);
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
	MTD_ECAL->DrawF1(-30,30,"same");
	vertex_plot->SaveAs((PLOTPATH+"/vertex_"+std::to_string(mtd_x)+".pdf").c_str());




}
float VertexDistance(TVector3 *vert,float x1, float x2, float y1, float y2, bool longitudinal = false,std::string PLOTPATH=""){
	float m1;
	float c1,c2;
	float dr=-1;
	
	m1 = (y2-y1)/(x2-x1);
	c1 = y1-m1*x1;
	
	
	std::cout << "m1" << m1 <<"__" <<  c1 << std::endl;
	std::cout << "coord_MTDx" << x1 <<"ECALX__" <<  x2 << std::endl;
	std::cout << "coord_MTDy" << y1 << "ECALY__"<< y2 << std::endl;
	if((y1 >0 && m1*x1 >0) || (y1< 0 && m1*x1<0)){
	if(longitudinal==true ){

	dr = c1-vert->y();
	std::cout << "dr " << dr << std::endl;
	return dr ;

	}else{

	dr = -c1/m1-vert->x();
	if(m1>10 && m1<30 &&  sqrt(x2*x2+y2*y2)>116)Plotter(x2,y2,x1,y1,vert,m1,c1,PLOTPATH);
	std::cout << "dr " << dr << std::endl;
	return dr ;

	}		

	}else return dr;


}

void AngularCoordinates(TH1D **hist, std::pair <float,float> matched_mct, struct pos centroid, TVector3* vertex){



	hist[0]->Fill(centroid.eta-matched_mct.first);
	hist[1]->Fill(centroid.phi-matched_mct.second);
	TVector3 temp;
	TVector3 vert;
	vert.SetXYZ(vertex->x(),vertex->y(),vertex->z());
	temp.SetMagThetaPhi(137-vert.Mag(),2*atan(exp(-matched_mct.first)),matched_mct.second);
	temp+=vert;
	std::cout << "________________temp mag" << temp.Mag() << std::endl;
	hist[2]-> Fill(centroid.eta-temp.Eta());
	hist[3]-> Fill(centroid.phi-temp.Phi());



	

}


TF1* fitgaus(TH1D* hist){

	TF1* fit= new TF1("gaus","gaus",hist->GetMean()-hist->GetNbinsX()/4*hist->GetXaxis()->GetBinWidth(1),hist->GetMean()+hist->GetNbinsX()/4*hist->GetXaxis()->GetBinWidth(1));
	fit->SetParLimits(0,hist->GetMaximum()-50,hist->GetMaximum()+50);
	fit->SetParameter(1,hist->GetMean());
	fit->SetParameter(2,hist->GetRMS());
	hist->Fit("gaus","R0EM");
	return fit;



}


TF1* N_gausFit(TH1D* hist, int n_gaus){
	int i;
	double par[n_gaus*3];
	float max_range=-9999, min_range=1e6;
	TF1* fits[n_gaus];
	TF1* global_fit;
	std::string global_formula= "";
	std::pair <float,float> range[n_gaus*2];
	for(i=0;i<n_gaus;i++){
		if(i==0) range[i]=std::make_pair(hist->GetMean()-hist->GetNbinsX()/4*hist->GetXaxis()->GetBinWidth(1),hist->GetMean()+hist->GetNbinsX()/4*hist->GetXaxis()->GetBinWidth(1));
		else range[i]=std::make_pair(hist->GetMean()-hist->GetNbinsX()/2*hist->GetXaxis()->GetBinWidth(1),hist->GetMean()+hist->GetNbinsX()/2*hist->GetXaxis()->GetBinWidth(1));

		fits[i]=new TF1(("gaus"+std::to_string(i)).c_str(),"gaus",range[i].first,range[i].second);
		if(range[i].first<min_range) min_range=range[i].first;
		if(range[i].second>max_range) max_range=range[i].second;
		fits[i]->SetParameter(0,hist->GetMaximum());
		fits[i]->SetParameter(1,hist->GetMean());
		fits[i]->SetParameter(2,hist->GetRMS());
		if(i==0)fits[i]->SetParLimits(0,hist->GetMaximum()-3,hist->GetMaximum()+(int)hist->GetMaximum()/7);
		else fits[i]->SetParLimits(0,0,hist->GetMaximum()/2);
		std::cout << "___" << i << "___fit name___" << fits[i]->GetName() << std::endl; 
		hist->Fit(("gaus"+std::to_string(i)).c_str(),"0REM");
		fits[i]->GetParameters(&par[i*3]);
		global_formula+="gaus("+std::to_string(i*3)+")";
		if(i != n_gaus-1) global_formula+="+";
	}	
	global_fit= new TF1("global_fit",(global_formula).c_str(),min_range+5*hist->GetBinWidth(1),max_range-3*hist->GetBinWidth(1));
	global_fit->SetParameters(par);
	for(i=0;i<n_gaus*3;i++){
	if(i%3 != 1)	global_fit->SetParLimits(i,0,1000);

	}
	global_fit->SetParLimits(0,hist->GetMaximum()-(int)hist->GetMaximum()/10,hist->GetMaximum()+(int)hist->GetMaximum()/10);
	hist->Fit("global_fit","0REM");
	return global_fit;


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

