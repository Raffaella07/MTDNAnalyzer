#include "NeutralsClass.h"

/*NeutralsClass::NeutralsClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("gen_tree",tree);

   }
   Init(tree);
}

NeutralsClass::~NeutralsClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NeutralsClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NeutralsClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
*/
void NeutralsClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Track_ECAL_x = 0;
   Track_ECAL_y = 0;
   Track_ECAL_z = 0;
   Track_ECAL_phi = 0;
   Track_ECAL_eta = 0;
   Track_In_eta = 0;
   Track_pt = 0;
   Track_mtdt = 0;
   Track_genVtx_t = 0;
   Track_In_phi = 0;
   Track_In_x = 0;
   Track_In_y = 0;
   Track_In_z = 0;
   Track_charge = 0;
   Track_nhits = 0;
   Track_pt_Iphi = 0;
   Track_pt_Ieta = 0;
   Track_pt_Ephi = 0;
   Track_pt_Eeta = 0;
   pt = 0;
   phi = 0;
   eta = 0;
   dr = 0;
   pdgId = 0;
   label = 0;
   PId = 0;
   convRadius = 0;
   convPhi = 0;
   convZ = 0;
   energy = 0;
   nlegs = 0;
   plabel = 0;
   ecalEnergy = 0;
   ecalx = 0;
   ecaly = 0;
   ecalz = 0;
   ecalphi = 0;
   ecaleta = 0;
   hcalEnergy = 0;
   clus_size =0;   
   clus_size_x=0; 
   clus_size_y=0;
   clus_energy= 0; 
   clus_time = 0;   
   clus_eta =0;    
   clus_phi =0;   
   clus_x=0;      
   clus_y=0;      
   clus_z=0;      
   MTDx=0;        
   MTDy=0;        
   MTDz=0;        
   MTDphi=0;      
   MTDeta=0;      
   genpv_t = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetMakeClass(kFALSE);
   fChain->SetBranchAddress("Track_ECAL_x", &Track_ECAL_x, &b_Track_ECAL_x);
   fChain->SetBranchAddress("Track_ECAL_y", &Track_ECAL_y, &b_Track_ECAL_y);
   fChain->SetBranchAddress("Track_ECAL_z", &Track_ECAL_z, &b_Track_ECAL_z);
   fChain->SetBranchAddress("Track_ECAL_phi", &Track_ECAL_phi, &b_Track_ECAL_phi);
   fChain->SetBranchAddress("Track_ECAL_eta", &Track_ECAL_eta, &b_Track_ECAL_eta);
   fChain->SetBranchAddress("Track_In_eta", &Track_In_eta, &b_Track_In_eta);
   fChain->SetBranchAddress("Track_pt", &Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_mtdt", &Track_mtdt, &b_Track_mtdt);
   fChain->SetBranchAddress("Track_genVtx_t", &Track_genVtx_t, &b_Track_genVtx_t);
   fChain->SetBranchAddress("Track_In_phi", &Track_In_phi, &b_Track_In_phi);
   fChain->SetBranchAddress("Track_In_x", &Track_In_x, &b_Track_In_x);
   fChain->SetBranchAddress("Track_In_y", &Track_In_y, &b_Track_In_y);
   fChain->SetBranchAddress("Track_In_z", &Track_In_z, &b_Track_In_z);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_nhits", &Track_nhits, &b_Track_nhits);
   fChain->SetBranchAddress("Track_pt_Iphi", &Track_pt_Iphi, &b_Track_pt_Iphi);
   fChain->SetBranchAddress("Track_pt_Ieta", &Track_pt_Ieta, &b_Track_pt_Ieta);
   fChain->SetBranchAddress("Track_pt_Ephi", &Track_pt_Ephi, &b_Track_pt_Ephi);
   fChain->SetBranchAddress("Track_pt_Eeta", &Track_pt_Eeta, &b_Track_pt_Eeta);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("dr", &dr, &b_dr);
   fChain->SetBranchAddress("pdgId", &pdgId, &b_pdgId);
   fChain->SetBranchAddress("label", &label, &b_label);
   fChain->SetBranchAddress("PId", &PId, &b_PId);
   fChain->SetBranchAddress("convRadius", &convRadius, &b_convRadius);
   fChain->SetBranchAddress("convPhi", &convPhi, &b_convPhi);
   fChain->SetBranchAddress("convZ", &convZ, &b_convZ);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("nlegs", &nlegs, &b_nlegs);
   fChain->SetBranchAddress("plabel", &plabel, &b_plabel);
   fChain->SetBranchAddress("ecalEnergy", &ecalEnergy, &b_ecalEnergy);
   fChain->SetBranchAddress("ecalx", &ecalx, &b_ecalx);
   fChain->SetBranchAddress("ecaly", &ecaly, &b_ecaly);
   fChain->SetBranchAddress("ecalz", &ecalz, &b_ecalz);
   fChain->SetBranchAddress("ecalphi", &ecalphi, &b_ecalphi);
   fChain->SetBranchAddress("ecaleta", &ecaleta, &b_ecaleta);
   fChain->SetBranchAddress("hcalEnergy", &hcalEnergy, &b_hcalEnergy);
   fChain->SetBranchAddress("clus_size", &clus_size, &b_clus_size);
   fChain->SetBranchAddress("clus_size_x", &clus_size_x, &b_clus_size_x);
   fChain->SetBranchAddress("clus_size_y", &clus_size_y, &b_clus_size_y);
   fChain->SetBranchAddress("clus_energy", &clus_energy, &b_clus_energy);
   fChain->SetBranchAddress("clus_time", &clus_time, &b_clus_time);
   fChain->SetBranchAddress("clus_x", &clus_x, &b_clus_x);
   fChain->SetBranchAddress("clus_y", &clus_y, &b_clus_y);
   fChain->SetBranchAddress("clus_z", &clus_z, &b_clus_z);
   fChain->SetBranchAddress("clus_eta", &clus_eta, &b_clus_eta);
   fChain->SetBranchAddress("clus_phi", &clus_phi, &b_clus_phi);
   fChain->SetBranchAddress("MTDx", &MTDx, &b_MTDx);
   fChain->SetBranchAddress("MTDy", &MTDy, &b_MTDy);
   fChain->SetBranchAddress("MTDz", &MTDz, &b_MTDz);
   fChain->SetBranchAddress("MTDphi", &MTDphi, &b_MTDphi);
   fChain->SetBranchAddress("MTDeta", &MTDeta, &b_MTDeta);
   fChain->SetBranchAddress("genpv_t", &genpv_t, &b_genpv_t);
 //  Notify();
}

/*Bool_t NeutralsClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NeutralsClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NeutralsClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NeutralsClass_cxx */
