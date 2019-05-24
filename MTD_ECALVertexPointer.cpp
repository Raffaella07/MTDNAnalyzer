//#include "TTree.h"
//#include "TMath.h"
//#include "TFile.h"
//#include "TH1.h"
//#include "TCanvas.h"
//#include <utility>  
//#include <iostream>
//#include <vector>
//#include "TH2.h"
//#include <string>
//#include "TString.h"
//#include "TEllipse.h"
//#include "TArrow.h"
//#include "TStyle.h"
//#include "TMarker.h"
//#include "TText.h"
//#include "TF1.h"
//#include "TVector3.h"
//#include "TGraphErrors.h"
//#include "TGraph.h"
//#include "TPaveText.h"
//#include "TLegend.h"
//#include "NeutralsClass.h"
//#include "NeutralsClass.cpp"
#include "MTD_ECALUtils.cc"




int main(int argc, char **argv){

	//--grab and initialize trees
	TFile* file = TFile::Open(argv[1]);
	std::string outname = argv[2];
	int nonoise = atoi(argv[3]);
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

	//--indices definition
	int i,j,k;
	float dr,label=-1;
	//--bin variables
	int bin_dr,bin_pt,bin_eta;
	float dr_min, dr_max,pt_min,pt_max,eta_min,eta_max;
	//--filenames
	std::string ang_filename[4]={"gen_clus_pho_deta.pdf","gen_clus_pho_dphi.pdf","gen_clus_gpho_eta.pdf", "gen_clus_gpho_phi.pdf"};



	Double_t ybin[2];
	//--Path to plots
	std::string PLOTPATH=std::string(" ~/CMSSW_10_4_0_mtd5_prova/plots/"+outname+"events/VertexPointer");
	system(("mkdir "+PLOTPATH).c_str());
	setStyle();
	//--set bin numbers and ranges
	bin_dr=50;
	bin_pt=10;
	dr_min=0;
	dr_max=0.4;
	pt_min=50;
	pt_max=60;
	bin_eta=10;
	eta_min=0;
	eta_max=3;
	ybin[0]=0.0;
	ybin[1]=0.02;
	
	//--Histos and fit definitions
	TH1D * proj_full[7];
//	TF1* narrow_gaus = new TF1("narrowgaus","gaus",-25,25);
//	TF1* wide_gaus = new TF1("widegaus","gaus",-15,15);


/*	Double_t par[6];*/
	TF1* doublegaus_fit;// =new TF1 ("doublegaus", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)",-25,25);
//	doublegaus_fit->SetParLimits(2,0,1000000);
//	doublegaus_fit->SetParLimits(5,0,1000000);
	


	
	TH1D* angular_dis[4];
	TF1* angular_fit[4];	
	angular_dis[0]= new TH1D("#eta_{clus}-#eta_{gen}","#eta_{clus}-#eta_{gen}",40,-0.1,0.1);
	angular_dis[1]= new TH1D("#phi_{clus}-#phi_{gen}","#phi_{clus}-#phi_{gen}",40,-0.05,0.05);
	angular_dis[2]= new TH1D("#eta_{clus}-#eta_{geom}","#eta_{clus}-#eta_{geom}",40,-0.03,0.03);
	angular_dis[3]= new TH1D("#phi_{clus}-#phi_{geom}","#phi_{clus}-#phi_{geom}",40,-0.03,0.03);
	TH1D* slices[2];
	TH2D* dr_photon_clus = new TH2D("","dr among clusters and matched mct photons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* dr_electron_clus = new TH2D("","dr among clusters and matched mct electrons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* photon_dr_clusMTD = new TH2D("","dr among clusters and MTD hit for mct photons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* ele_dr_clusMTD = new TH2D("","dr among clusters and Mtd hit for mct electrons",bin_dr,dr_min,dr_max,bin_pt,pt_min,pt_max);
	TH2D* photon_vert_dx = new TH2D("","dx of the real vertex wrt to the ECAL_MTD line",10,-15,15,50,0.0,0.4);
	TH2D* photon_vert_dy = new TH2D("","dy of the real vertex wrt to the ECAL_MTD line",35,-15,15,2,0,0.5);
	TH2D* photon_vert_dz = new TH2D("","dz of the real vertex wrt to the ECAL_MTD line",45,-15,15,bin_pt,pt_min,pt_max);




	//--loop on events
	for(i=0;i<cand.fChain->GetEntries();i++){
		cand.fChain->GetEntry(i);
		mct.fChain->GetEntry(i);
		//--loop on mct photons
		for(j=0;j<mct.pt->size();j++){
			dr =-1;
			//-- matching variables definitions
			std::pair <float,float>  matched_mct;
			std::pair <float,float>  pair_clus;
			struct pos matched_clus;
			struct pos matched_MTDhit;
					
	
			TVector3 *vertex=new TVector3();
			vertex->SetX(mct.primary_x->at(j));
			vertex->SetY(mct.primary_y->at(j));
			vertex->SetZ(mct.primary_z->at(j));
			if(abs(mct.eta->at(j))< 1.4){
		std::cout << "primaryvertexZ" << vertex->Z() << "___" << cand.ecalx->size() << std::endl;
//			//--matches mct phtons to clusters
			dr = matcher(mct.eta->at(j),mct.phi->at(j),cand.ecaleta, cand.ecalphi,cand.ecalx,cand.ecaly,cand.ecalz,cand.pt,&matched_clus,&matched_mct,vertex,true);
			if(dr <0.05 && dr != -1){
	std::cout << "quiiii______mct "<< matched_mct.first << std::endl;
				label = mct.PId->at(j);
				if (label == 4  && (mct.convRadius->at(j)>116 && mct.convRadius->at(j) < 119)) //--MTD converting photons
				{ 
				 dr_photon_clus->Fill(dr,mct.pt->at(j));
				AngularCoordinates(angular_dis,matched_mct,matched_clus,vertex);	//--angular coordinates histograms
				}else if(label == 2) dr_electron_clus->Fill(dr,mct.pt->at(j));
			
//		std::cout << "herewww"	 << std::endl; 	
			if(nonoise == 0){ NonoiseMTD(vertex,matched_clus,matched_mct,&matched_MTDhit);
			std::cout << "-----------------------------------------------clus " << matched_clus.x << "-------------------------" << matched_clus.y << "---------------------------" << matched_clus.z << std::endl;
			}else if (nonoise == 1)dr = matcher(matched_clus.eta,matched_clus.phi,cand.MTDeta, cand.MTDphi, cand.MTDx, cand.MTDy, cand.MTDz,cand.pt,&matched_MTDhit,&pair_clus, vertex, false);
			std::cout << "clus_mtd dr____" << dr << "label__" <<label <<std::endl;

			if(dr !=-1){
				if (label == 4 && (mct.convRadius->at(j)>116 && mct.convRadius->at(j) < 119)){ //--MTD converting photons
				 photon_dr_clusMTD->Fill(dr,matched_MTDhit.pt);
				std::cout << "---------------------------------------------------------------------------------clus rad  " << sqrt(pow(matched_clus.x,2)+pow(matched_clus.y,2)) << std::endl;
					if(dr < 0.05  && sqrt(pow(matched_clus.x,2)+pow(matched_clus.y,2))>129 && sqrt(pow(matched_MTDhit.x,2)+pow(matched_MTDhit.y,2))>116){
						//if(sqrt(pow(matched_clus.z,2)+pow(matched_clus.y,2))>129 && sqrt(pow(matched_MTDhit.z,2)+pow(matched_MTDhit.y,2))>116){
					/*	std::cout << "" << std::endl;
						std::cout << "before function" << std::endl;
						std::cout << "coord_MTDx" << matched_MTDhit.x <<"ECALX__" <<  matched_clus.x << std::endl;
						std::cout << "coord_MTDy" << matched_MTDhit.y << "ECALY__"<< matched_clus.y << std::endl;
						std::cout << "" << std::endl;
						std::cout << "in function" << std::endl;
						std::cout << "dx" << std::endl;*/
						float dx = VertexDistance(vertex,matched_MTDhit.x,matched_clus.x,matched_MTDhit.y,matched_clus.y,1,PLOTPATH);	
					//	std::cout << "dy" << std::endl;
						float dy = VertexDistance(vertex,matched_MTDhit.x,matched_clus.x,matched_MTDhit.y,matched_clus.y,2,PLOTPATH);	
						std::cout << "" << std::endl;
						std::cout << "dz" << std::endl;
						float dz = VertexDistance(vertex,matched_MTDhit.z,matched_clus.z,matched_MTDhit.y,matched_clus.y,3,PLOTPATH);	
						std::cout << "vert "<< vertex->z()  <<"____________________________________________________________________________________________________________dz=" << dz << std::endl;
						std::cout << "" << std::endl;
						if(dx != -1 && dy != -1 && dz != -1 ){	
					//	std::cout << "displacements" << dx << "____" << dz << std::endl;
					/*	if(sqrt(dxy*dxy + dzy*dzy) < 10)*/
						photon_vert_dx->Fill(dx,dr);
						photon_vert_dy->Fill(dy,dr);
						photon_vert_dz->Fill(dz,mct.pt->at(j));
				//		if(mct.pt->at(j) < 50) std::cout << "pt less thn expecteeed " << std::endl;
						std::cout << "holaaa" << std::endl;
						}
					
					}
				}
				else if(label == 2) ele_dr_clusMTD->Fill(dr,mct.pt->at(j));

			}



		}
		}

		}

	}
	system(("mkdir"+PLOTPATH+"/controlplots").c_str());
	Slicer(PLOTPATH,bin_pt,pt_min,pt_max,"dr",photon_dr_clusMTD,"Sphoton_dr_clusMTD");
	Slicer(PLOTPATH,bin_pt,pt_min,pt_max,"dr",ele_dr_clusMTD,"Sele_dr_clusMTD");
	Slicer(PLOTPATH,10,-1.4,pt_max,"dz(cm)",photon_vert_dz,"Sphoton_vert_dz");
	
	slices[0]=photon_vert_dx->ProjectionY("p",0,bin_dr-1);
	proj_full[0]=dr_photon_clus->ProjectionX("1",0,bin_pt-1);
	proj_full[1]=dr_electron_clus->ProjectionX("2",0,bin_pt-1);
        proj_full[2]=photon_dr_clusMTD->ProjectionX("3",0,bin_pt-1);
	proj_full[3]=ele_dr_clusMTD->ProjectionX("4",0,bin_pt-1);
	proj_full[4]=photon_vert_dx->ProjectionX("5",0,bin_dr-1);
	proj_full[5]=photon_vert_dy->ProjectionX("6",0,bin_dr-1);
	proj_full[6]=photon_vert_dz->ProjectionX("7",0,bin_dr-1);
	
	proj_full[0]->GetXaxis()->SetTitle("dr");
	proj_full[1]->GetXaxis()->SetTitle("dr");
	proj_full[2]->GetXaxis()->SetTitle("dr");
	proj_full[3]->GetXaxis()->SetTitle("dr");
	proj_full[4]->GetXaxis()->SetTitle("dx(cm)");
	proj_full[5]->GetXaxis()->SetTitle("dy(cm)");
	proj_full[6]->GetXaxis()->SetTitle("dz(cm)");


	for(i=0;i<7;i++){
	proj_full[i]->GetYaxis()->SetTitle("entries");
	}
/*	narrow_gaus->SetParLimits(0,proj_full[6]->GetMaximum()-50,1000000);
	//narrow_gaus->SetParameter(0,proj_full[6]->GetMaximum());
	narrow_gaus->SetParameter("Sigma",1);
	narrow_gaus->SetParLimits(2,1,6);
	wide_gaus->SetParameter("Sigma",3);
	wide_gaus->SetParLimits(2,2,1000000);
	wide_gaus->SetParLimits(0,0,60);
	
	proj_full[6]->Fit("narrowgaus","0R");
	proj_full[6]->Fit("widegaus","0R");
	narrow_gaus->GetParameters(&par[0]);
	wide_gaus->GetParameters(&par[3]);
	doublegaus_fit->SetParameters(par);
	doublegaus_fit->SetParLimits(0,proj_full[6]->GetMaximum()-100,proj_full[6]->GetMaximum()+100);
	proj_full[6]->Fit("doublegaus","0R");
*/	
	doublegaus_fit = N_gausFit(proj_full[6],2);
	for (i=0;i<4;i++){
	angular_fit[i]=N_gausFit(angular_dis[i],2);
	angular_dis[i]->GetXaxis()->SetTitle(angular_dis[i]->GetName());
	angular_dis[i]->GetYaxis()->SetTitle("entries");
	SavePlot((char*)angular_dis[i]->GetName(),angular_dis[i],(PLOTPATH+"/"+ang_filename[i]).c_str(),false,angular_fit[i]);
	


	}
//	SavePlot("Unconverted #gamma gen-cluster #Delta#eta",angular_dis[0],(PLOTPATH+"/gen_clus_pho_deta.pdf").c_str(),false,NULL);
//	SavePlot("Unconverted #gamma gen cluster  #Delta#phi ",angular_dis[1],(PLOTPATH+"/gen_clus_pho_dphi.pdf").c_str(),false,NULL);
//	SavePlot("Unconv #gamma gen-cluster geom #Delta#eta",angular_dis[2],(PLOTPATH+"/gen_clus_gpho_eta.pdf").c_str(),false,NULL);
//	SavePlot("Unconv #gamma gen-cluster geom #Delta#phi",angular_dis[3],(PLOTPATH+"/gen_clus_gpho_phi.pdf").c_str(),false,NULL);
	SavePlot("dr among matched mct photons and clusters",proj_full[0],(PLOTPATH+"/dr_photon_clus.pdf").c_str(),false,NULL);
	SavePlot("dr among matched mct electrons and clusters",proj_full[1],(PLOTPATH+"/dr_electron_clus.pdf").c_str(),false,NULL);
	SavePlot("dr among MTD hit and cluster for mct photons",proj_full[2],(PLOTPATH+"/photon_dr_clusMTD.pdf").c_str(),false,NULL);
	SavePlot("dr among MTD hit and cluster for mct electrons",proj_full[3],(PLOTPATH+"/ele_dr_clusMTD.pdf").c_str(),false,NULL);
	SavePlot("dx between real vertex and MTD_ECAL line",proj_full[4],(PLOTPATH+"/photon_vert_dx.pdf").c_str(),false,NULL);
	SavePlot("dx in 0, 0.02 dr interval",slices[0],(PLOTPATH+"/controlplots/photon_vert_proj.pdf").c_str(),false,NULL);
//	SavePlot("dx for dr >0.02",slices[1],(PLOTPATH+"/controlplots/photon_vert_dx0.4.pdf").c_str(),false);
	SavePlot("dy between real vertex and MTD_ECAL line",proj_full[5],(PLOTPATH+"/photon_vert_dy.pdf").c_str(),false,NULL);
	SavePlot("dz between real vertex and MTD_ECAL line",proj_full[6],(PLOTPATH+"/photon_vert_dz.pdf").c_str(),false,doublegaus_fit);
	SavePlot2D("dr vs dz",photon_vert_dz,(PLOTPATH+"/photon_vert_dz2.pdf").c_str(),false);
//	proj_full[2]->Scale(1/proj_full[2]->Integral());	
//	slices[0]->Scale(1/slices[0]->Integral());	
	slices[0]->SetLineColor(kRed);
	TCanvas* canvas = new TCanvas("","superposition",600,550);
	proj_full[2]->DrawNormalized("hist");
	slices[0]->DrawNormalized("histsame");
	canvas->SaveAs((PLOTPATH+"/superpos.pdf").c_str());
	canvas->SaveAs((PLOTPATH+"/superpos.root").c_str());
	




}


