#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <utility>  
#include <iostream>
#include <vector>
#include "TH2.h"
#include "THStack.h"
#include <string>
#include "TString.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TAttLine.h"
#include "NeutralsClass.h"
#include "NeutralsClass.cpp"


void SavePlot (char * , TH1D *, const char*, bool );

TH1D*  projections(int, TH2* , std::string,std::string, std::string ,Color_t, std::string);

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


	const int Nsig=5; 
	int i,j,k,cn,cm,NbinX,NbinY,idx,ind;
	float Xmin,Xmax,Ymin,Ymax;
	std::vector<float> *particleId=0,*mct_pt=0,*convRadius=0;
	const char * signatures[Nsig];
	Color_t color[Nsig]={kRed,kBlue,kCyan,kOrange+7,kGreen+3};
	std::string PLOTPATH;

	std::cout << "labels   " << std::endl;

	NbinX = 15;
	NbinY= 10;
	Xmin = 0;
	Xmax = 100;
	Ymin = 0;
	Ymax = 10;
	signatures[0]= "e+e-";
	signatures[1]= "e_same";
	signatures[2]= "e_single";
	signatures[3]= "nomatch";
//	signatures[4]= "singlepm";
//	signatures[5]= "UnID";
	signatures[4]= "tot";


	PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events");
	std::string   pt_min[NbinY];
	std::string  pt_max[NbinY];
	system(("mkdir"+ PLOTPATH+"/RecoEfficiency").c_str());
	THStack * stack_signatures[NbinY];
	THStack * stack_sig_all= new THStack("","");
	//	TH2F * signature_pt =new TH2F("TH2","Ratio of recognized reco signatures in bins of photon pt for convRdius < 50 cm",8,0,7,10,0,11);
	TH2D * histo_sig[Nsig];
	histo_sig[0]=new TH2D("histo complete conversion","histo complete conversion",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
	histo_sig[1]=new TH2D("histo single electron","histo single electron",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
	histo_sig[2]=new TH2D("histo single photon","histo single photon",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
	histo_sig[3]=new TH2D("histo double photon","histo double photon",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
	histo_sig[4]=new TH2D("histo ele gamma","histo ele gamma",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
//	histo_sig[5]=new TH2D("histo mix","histo mix",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
//	histo_sig[6]=new TH2D("histo tot","histo tot",NbinX, Xmin,Xmax, NbinY,Ymin,Ymax);
	TH1D * proj[Nsig][NbinY];
	TH1D * h_all[Nsig];
	for(i=0;i<mct.fChain->GetEntries();i++){
		mct.fChain->GetEntry(i);
		cand.fChain->GetEntry(i);
		gen.fChain->GetEntry(i);
		float temp=0;
		float electrons[2];
		std::cout <<"____"<< cand.event << std::endl;
		for(j=0;j<mct.pt->size();j++){
			if(temp != mct.label->at(j) )temp = mct.label->at(j);
			for(k=0;k<mct.pt->size();k++){
				if(mct.label->at(k)==temp && mct.convRadius->at(k)<Xmax && mct.plabel->at(j)!= mct.plabel->at(k)){
	
				histo_sig[4]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
				electrons[0]=mct.plabel->at(k);
				electrons[1]=mct.plabel->at(j);
					std::cout << "mctlabel_mctlabel "<<electrons[0] <<"____"<< electrons[1] <<  std::endl;
				
					std::cout << "candsize "<< cand.plabel->size() <<  std::endl;
				for(idx=0; idx< cand.plabel->size();idx++){
					if(cand.plabel->at(idx)== electrons[0]){
					for(ind=0; ind< cand.pt->size();ind++){
					std::cout << "candlabel_candlabel "<<cand.plabel->at(ind) <<"____"<< cand.plabel->at(idx) <<  std::endl;
					if(cand.plabel->at(ind)== electrons[1] && idx != ind){
					 histo_sig[0]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					
					std::cout << "plabels " <<  std::endl;
					}
					else if(cand.plabel->at(ind) == electrons[0]) histo_sig[1]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					else if(cand.plabel->at(ind)==-1) histo_sig[2]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					break;
					}
				}
				
				else if(cand.plabel->at(idx)== electrons[1]){
					for(ind=0; ind< cand.pt->size();ind++){
					if(cand.plabel->at(ind)== electrons[0] && idx != ind){
					 histo_sig[0]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					std::cout << "plabels " <<  std::endl;
					}
					else if(cand.plabel->at(ind) == electrons[1]) histo_sig[1]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					else if (cand.plabel->at(ind==-1)) histo_sig[2]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					break;
					}
				}else if (cand.plabel->at(ind)==-1 && cand.plabel->at(idx)==-1){
					 histo_sig[3]->Fill(mct.convRadius->at(k),gen.pt->at(mct.label->at(k)));
					break;
						}

						}
			break;		}
				}
			}
		}
	
				
		/*{{	for(j=0;j<particleId->size();j++){
				if (particleId->at(j)==4){
					particleId->erase(particleId->begin()+j); //aggiungi tutti i vettori che usi a questa operazione
					convRadius->erase(convRadius->begin()+j); //aggiungi tutti i vettori che usi a questa operazione
					mct_pt->erase(mct_pt->begin()+j); //aggiungi tutti i vettori che usi a questa operazionei
					break;
				}
			}

			if(convRadius->at(1) < Xmax && mct_pt->at(1)== mct_pt->at(0) ){
				std::cout << "labl " << particleId->at(0) << "    " << particleId->at(1) << std::endl;
				histo_sig[6]->Fill(convRadius->at(1),mct_pt->at(1));

				if (particleId->at(0) == 2 )
				{
					if ( particleId->at(1)==2)histo_sig[0]->Fill(convRadius->at(1),mct_pt->at(1));
					else if ( particleId->at(1)==4) histo_sig[1]->Fill(convRadius->at(1),mct_pt->at(1));
					else if (particleId->at(1)==-1) histo_sig[4]->Fill(convRadius->at(1),mct_pt->at(1));
				}
				else if(particleId->at(0) == 4){
					if (particleId->at(1) == 4)  histo_sig[2]->Fill(convRadius->at(1),mct_pt->at(1));
					else if (particleId->at(1) == -1)  histo_sig[3]->Fill(convRadius->at(1),mct_pt->at(1));
				}		
				else if (particleId->at(1) == 2 )
				{
					if ( particleId->at(0)==2)histo_sig[0]->Fill(convRadius->at(1),mct_pt->at(1));
					else if ( particleId->at(0)==4) histo_sig[1]->Fill(convRadius->at(1),mct_pt->at(1));
					else if (particleId->at(0)==-1) histo_sig[4]->Fill(convRadius->at(1),mct_pt->at(1));
				}
				else if(particleId->at(1) == 4){
					if (particleId->at(0) == 4)  histo_sig[2]->Fill(convRadius->at(1),mct_pt->at(1));
					else if (particleId->at(0) == -1)  histo_sig[3]->Fill(convRadius->at(1),mct_pt->at(1));
				}		
				else  histo_sig[5]->Fill(convRadius->at(1),mct_pt->at(1));

			}
		}		
	}*/

	for(i=0;i<Nsig; i++){
		system(("mkdir "+ PLOTPATH+"/RecoEfficiency/"+signatures[i]).c_str());
		h_all[i]=histo_sig[i]->ProjectionX((std::string("Ratio of identified reco signature ")+signatures[i]+std::string(", convRadius < 80cm")).c_str());
		histo_sig[i]->GetXaxis()->SetTitle("convRadius(cm)");
		histo_sig[i]->GetYaxis()->SetTitle("pt (Gev/c)");

		for(j=0;j<NbinY;j++){
			float min_pt = (Ymax-Ymin)/NbinY *j;
			std::string  min = std::to_string((int)min_pt);
			std::string max = std::to_string( (int)(min_pt + (Ymax-Ymin)/NbinY));
			pt_min[j]=min;
			pt_max[j]=max;
			proj[i][j]= projections(j,histo_sig[i],min,max,signatures[i],color[i],PLOTPATH+"/RecoEfficiency");


		}
	}


	TCanvas * stack_canvas[NbinY];
	TLegend * legend = new TLegend(0.55,0.7,0.9,0.9);
	legend->AddEntry(proj[0][0], "e^{+} e^{-}","l");
	legend->AddEntry(proj[1][0], "same e^{#pm}","l");
	legend->AddEntry(proj[2][0],  "single e^{#pm}","l");
	legend->AddEntry(proj[3][0], "no matched","l");
//	legend->AddEntry(proj[4][0], "e^{+}","l");
//	legend->AddEntry(proj[5][0], "unrecognized","l");
	for(i=0;i<NbinY;i++){
		stack_signatures[i]= new THStack("","");
		for(j=0;j<Nsig-1;j++){
			proj[j][i]->Divide(proj[4][i]);

		proj[j][i]->GetYaxis()->SetRangeUser(0,0.8);
			stack_signatures[i]->Add(proj[j][i]);
		}

		stack_canvas[i] = new TCanvas("",("stacked signaturesVsconv radius, pt ["+pt_min[i]+","+pt_max[i]+"]").c_str(),600,550);
		stack_signatures[i]->Draw("nostack");
		stack_signatures[i]->SetTitle(("stacked signatures Vs convRadius, pt["+pt_min[i]+","+pt_max[i]+"](GeV/c)").c_str());
		stack_signatures[i]->GetXaxis()->SetTitle("convRadius(cm)");
		stack_signatures[i]->GetYaxis()->SetTitle("ratio of converted photons");
		legend->Draw();
		stack_canvas[i]->SaveAs((PLOTPATH+"/RecoEfficiency/stacked_signatures_ pt_"+pt_min[i]+"_"+pt_max[i]+".pdf").c_str());


	}





	/*	for (j=0;j< Nsig; j++){
		int counter =0;
		for(i=0;i<8;i++){
		if (signature_pt->GetBinContent(j,i)==0)
		{
		counter++;

		signature_pt->SetBinContent(signatures[j],0.00000001);
		std::cout << "empty bin " << signatures[j] << std::endl;


		}
		signature_pt->GetXaxis()->SetTitle("single photon reco signatures in bins of convRadius");
		signature_pt->GetYaxis()->SetTitle("pt (Gev/c)");
		TH1D * projection[NbinY];
		TCanvas * projs = new TCanvas ("projs","projs",600,550);
		for(i=0;i<NbinY;i++){
		float min_pt = (Ymax-Ymin)/NbinY *i;
		std::string  min = std::to_string(min_pt);
		std::string max = std::to_string( min_pt + (Ymax-Ymin)/NbinY);
		projection[i] = signature_pt->ProjectionX(("Ratio of recognized reco signatures for pt in ["+min +";"+max +" ], in convRadius bin").c_str(),i,i+1); 
		projection[i]->GetYaxis()->SetTitle("ratio of reco candidates");
		projs ->Clear();

		projs->SaveAs((PLOTPATH+"/RecoEfficiency/signs_pt_"+min+"to_"+max+"convR_50.pdf").c_str());
		projection->Reset();


		}




		}


		} */ 


	h_all[0]->Divide(h_all[4]);
	h_all[0]->SetLineWidth(6);
	h_all[0]->SetLineColor(color[0]);
	h_all[1]->Divide(h_all[4]);
	h_all[1]->SetLineWidth(6);
	h_all[1]->SetLineColor(color[1]);
	h_all[2]->Divide(h_all[4]);
	h_all[2]->SetLineWidth(6);
	h_all[2]->SetLineColor(color[2]);
	h_all[3]->Divide(h_all[4]);
	h_all[3]->SetLineWidth(6);
	h_all[3]->SetLineColor(color[3]);
	h_all[3]->GetYaxis()->SetRangeUser(0,0.8);

/*	h_all[4]->Divide(h_all[6]);
	h_all[5]->Divide(h_all[6]);
	h_all[4]->SetLineWidth(6);
	h_all[4]->SetLineColor(color[4]);
	h_all[5]->SetLineWidth(6);
	h_all[5]->SetLineColor(color[5]);*/

	stack_sig_all->Add(h_all[0]);
	stack_sig_all->Add(h_all[1]);
	stack_sig_all->Add(h_all[2]);
	stack_sig_all->Add(h_all[3]);
	//stack_sig_all->Add(h_all[4]);
	//stack_sig_all->Add(h_all[5]);
	SavePlot("histo complete conversion", h_all[0],(PLOTPATH+"/RecoEfficiency/histo_ee.pdf").c_str(),false);
	SavePlot("histo same electron", h_all[1],(PLOTPATH+"/RecoEfficiency/histo_singlee.pdf").c_str(),false);
	SavePlot("histo single electron", h_all[2],(PLOTPATH+"/RecoEfficiency/histo_singlegamma.pdf").c_str(),false);
	SavePlot("histo no match", h_all[3],(PLOTPATH+"/RecoEfficiency/histo_gammagamma.pdf").c_str(),false);
	SavePlot("histo tot conversion", h_all[4],(PLOTPATH+"/RecoEfficiency/histo_egamma.pdf").c_str(),false);
//	SavePlot("histo UnID", h_all[5],(PLOTPATH+"/RecoEfficiency/histo_mix.pdf").c_str(),false);
	TCanvas * stackAll_canvas= new TCanvas("stacked signatures","stacked signatures", 600, 550);
	stack_sig_all->Draw("nostack");
	stack_sig_all->SetTitle("stacked signatures");
	stack_sig_all->GetXaxis()->SetTitle("convRadius(cm)");
	stack_sig_all->GetYaxis()->SetTitle("ratio of converted photons");
	legend->Draw();
	stackAll_canvas->SaveAs((PLOTPATH+"/RecoEfficiency/stacked_signatures.pdf").c_str());	

}
void SavePlot (char * titlestring, TH1D * histo, const char * filename, bool log=false){

	TCanvas* canvas = new TCanvas(titlestring,titlestring,600,550);
	histo->Draw("hist");
	canvas->SaveAs(filename);

	if (log) canvas->SetLogy();
	canvas->Clear();
}


TH1D* projections(int i,TH2* histo, std::string min, std::string max, std::string dir,Color_t color,std::string PLOTPATH){


	TH1D* projection;
	TCanvas * projs = new TCanvas ("projs","projs",600,550);
	projection = histo->ProjectionX(("Ratio of recognized reco signature"+dir+"  for pt in ["+min +";"+max +"], in convRadius bin").c_str(),i,i+1); 
	std::cout << "limits      " << min << "  " << max << std::endl;
	projection->GetYaxis()->SetTitle("ratio of reco candidates");
	projection->SetLineColor(color);
	projection->SetLineWidth(6);
	projs ->Clear();
	projs->SetTitle(("Number of recognized reco signature"+dir+"  for pt in ["+min +";"+max +"], in convRadius bin").c_str());
	projection->Draw();
	projs->SaveAs((PLOTPATH+"/"+dir+"/signs_pt_"+min+"to_"+max+"convR<80.pdf").c_str());
	delete projs;



	return projection;


}
