#include "Tclass.h"


void Tclass::Init(TTree* tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   gen_DR = 0;
   Track_ECAL_x = 0;
   Track_ECAL_y = 0;
   Track_ECAL_z = 0;
   Track_ECAL_phi = 0;
   Track_ECAL_eta = 0;
   Track_In_eta = 0;
   Track_pt = 0;
   Track_In_phi = 0;
   Track_In_x = 0;
   Track_In_y = 0;
   Track_In_z = 0;
   Track_pt_Iphi = 0;
   Track_pt_Ieta = 0;
   Track_pt_Ephi = 0;
   Track_pt_Eeta = 0;
   Track_charge = 0;
   Track_nhits = 0;
   genpv_z = 0;
   dr = 0;
   pt = 0;
   label =0;
   eta = 0;
   phi = 0;
   ecalEnergy = 0;
   hcalEnergy = 0;
   particleId = 0;
   clus_n = 0;
   gen_pt = 0;
   gen_eta = 0;
   gen_phi = 0;
   gen_pdgId = 0;
   mct_pt = 0;
   mct_eta = 0;
   mct_phi = 0;
   mct_energy = 0;
   mct_convRadius = 0;
   mct_convZ = 0;
   mct_convPhi = 0;
   mct_ele1_pt = 0;
   mct_ele1_eta = 0;
   mct_nlegs = 0;
   mct_ele1_phi = 0;
   ele1_match = 0;
   mct_ele2_pt = 0;
   mct_ele2_eta = 0;
   mct_ele2_phi = 0;
   ele2_match = 0;
   mct_eles_dr = 0;
   mct_ele1_dr = 0;
   mct_ele2_dr = 0;
   genpv_t = 0;
   

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("gen_DR", &gen_DR, &b_gen_DR);
   fChain->SetBranchAddress("Track_ECAL_x", &Track_ECAL_x, &b_Track_ECAL_x);
   fChain->SetBranchAddress("Track_ECAL_y", &Track_ECAL_y, &b_Track_ECAL_y);
   fChain->SetBranchAddress("Track_ECAL_z", &Track_ECAL_z, &b_Track_ECAL_z);
   fChain->SetBranchAddress("Track_ECAL_phi", &Track_ECAL_phi, &b_Track_ECAL_phi);
   fChain->SetBranchAddress("Track_ECAL_eta", &Track_ECAL_eta, &b_Track_ECAL_eta);
   fChain->SetBranchAddress("Track_In_eta", &Track_In_eta, &b_Track_In_eta);
   fChain->SetBranchAddress("Track_pt", &Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_In_phi", &Track_In_phi, &b_Track_In_phi);
   fChain->SetBranchAddress("Track_In_x", &Track_In_x, &b_Track_In_x);
   fChain->SetBranchAddress("Track_In_y", &Track_In_y, &b_Track_In_y);
   fChain->SetBranchAddress("Track_In_z", &Track_In_z, &b_Track_In_z);
   fChain->SetBranchAddress("Track_pt_Iphi", &Track_pt_Iphi, &b_Track_pt_Iphi);
   fChain->SetBranchAddress("Track_pt_Ieta", &Track_pt_Ieta, &b_Track_pt_Ieta);
   fChain->SetBranchAddress("Track_pt_Ephi", &Track_pt_Ephi, &b_Track_pt_Ephi);
   fChain->SetBranchAddress("Track_pt_Eeta", &Track_pt_Eeta, &b_Track_pt_Eeta);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_nhits", &Track_nhits, &b_Track_nhits);
   fChain->SetMakeClass(0);
   fChain->SetBranchAddress("genpv_z", &genpv_z, &b_genpv_z);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("ecalEnergy", &ecalEnergy, &b_ecalEnergy);
   fChain->SetBranchAddress("hcalEnergy", &hcalEnergy, &b_hcalEnergy);
   fChain->SetBranchAddress("particleId", &particleId, &b_particleId);
   fChain->SetBranchAddress("clus_n", &clus_n, &b_clus_n);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("label", &label, &b_label);
   fChain->SetBranchAddress("dr", &dr, &b_dr);
   fChain->SetBranchAddress("gen_eta", &gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_pdgId", &gen_pdgId, &b_gen_pdgId);
   fChain->SetBranchAddress("mct_pt", &mct_pt, &b_mct_pt);
   fChain->SetBranchAddress("mct_eta", &mct_eta, &b_mct_eta);
   fChain->SetBranchAddress("mct_phi", &mct_phi, &b_mct_phi);
   fChain->SetBranchAddress("mct_energy", &mct_energy, &b_mct_energy);
   fChain->SetBranchAddress("mct_convRadius", &mct_convRadius, &b_mct_convRadius);
   fChain->SetBranchAddress("mct_convZ", &mct_convZ, &b_mct_convZ);
   fChain->SetBranchAddress("mct_convPhi", &mct_convPhi, &b_mct_convPhi);
   fChain->SetBranchAddress("mct_ele1_pt", &mct_ele1_pt, &b_mct_ele1_pt);
   fChain->SetBranchAddress("mct_ele1_eta", &mct_ele1_eta, &b_mct_ele1_eta);
   fChain->SetBranchAddress("mct_nlegs", &mct_nlegs, &b_mct_nlegs);
   fChain->SetBranchAddress("mct_ele1_phi", &mct_ele1_phi, &b_mct_ele1_phi);
   fChain->SetBranchAddress("ele1_match", &ele1_match, &b_ele1_match);
   fChain->SetBranchAddress("mct_ele2_pt", &mct_ele2_pt, &b_mct_ele2_pt);
   fChain->SetBranchAddress("mct_ele2_eta", &mct_ele2_eta, &b_mct_ele2_eta);
   fChain->SetBranchAddress("mct_ele2_phi", &mct_ele2_phi, &b_mct_ele2_phi);
   fChain->SetBranchAddress("ele2_match", &ele2_match, &b_ele2_match);
   fChain->SetBranchAddress("mct_eles_dr", &mct_eles_dr, &b_mct_eles_dr);
   fChain->SetBranchAddress("mct_ele1_dr", &mct_ele1_dr, &b_mct_ele1_dr);
   fChain->SetBranchAddress("mct_ele2_dr", &mct_ele2_dr, &b_mct_ele2_dr);
   fChain->SetBranchAddress("genpv_t", &genpv_t, &b_genpv_t);
   

   /*genpv_z->resize(2);
   pt->resize(2);
   eta->resize(2);
   phi->resize(2);
   ecalEnergy->resize(2);
   hcalEnergy->resize(2);
   particleId->resize(2);
   clus_n->resize(2);
   gen_pt->resize(2);
   dr->resize(2);
   label->resize(2);
   gen_eta->resize(2);
   gen_phi->resize(2);
   gen_pdgId->resize(2);
   mct_pt->resize(2);
   mct_eta->resize(2);
   mct_phi->resize(2);
   mct_energy->resize(2);
   mct_convRadius->resize(2);
   mct_convZ->resize(2);
   mct_convPhi->resize(2);
   mct_ele1_pt->resize(2);
   mct_ele1_eta->resize(2);
   mct_nlegs->resize(2);
   mct_ele1_phi->resize(2);
   ele1_match->resize(2);
   mct_ele2_pt->resize(2);
   mct_ele2_eta->resize(2);
   mct_ele2_phi->resize(2);
   ele2_match->resize(2);
   mct_eles_dr->resize(2);
   mct_ele1_dr->resize(2);
   mct_ele2_dr->resize(2);
   genpv_t ->resize(2);
for(int i=0;i<phi->size() ;i++){

   genpv_z->at(i) = 0;
   pt ->at(i) = 0;
   label ->at(i) =0;
   dr ->at(i) =0;
   eta ->at(i) = 0;
   phi ->at(i) = 0;
   ecalEnergy ->at(i) = 0;
   hcalEnergy ->at(i) = 0;
   particleId ->at(i) = 0;
   clus_n ->at(i) = 0;
   gen_pt ->at(i) = 0;
   gen_eta ->at(i) = 0;
   gen_phi ->at(i) = 0;
   gen_pdgId ->at(i) = 0;
   mct_pt ->at(i) = 0;
   mct_eta ->at(i) = 0;
   mct_phi ->at(i) = 0;
   mct_energy ->at(i) = 0;
   mct_convRadius ->at(i) = 0;
   mct_convZ ->at(i) = 0;
   mct_convPhi ->at(i) = 0;
   mct_ele1_pt->at(i)  = 0;
   mct_ele1_eta ->at(i) = 0;
   mct_nlegs ->at(i) = 0;
   mct_ele1_phi ->at(i) = 0;
   ele1_match ->at(i) = 0;
   mct_ele2_pt ->at(i) = 0;
   mct_ele2_eta ->at(i) = 0;
   mct_ele2_phi ->at(i) = 0;
   ele2_match ->at(i)= 0;
   mct_eles_dr ->at(i) = 0;
   mct_ele1_dr ->at(i) = 0;
   mct_ele2_dr ->at(i) = 0;
   genpv_t ->at(i) = 0;*/

//std::cout << "here " <<phi->size() <<  std::endl;
	
 //  Notify();
}
