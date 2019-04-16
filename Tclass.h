//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr  4 09:54:41 2019 by ROOT version 6.12/07
// from TTree tracks_tree/4D TOFPID studies
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef Tclass_h
#define Tclass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.
#include <vector>


class Tclass {
public :
   TTree* fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
  std::vector<float>    *gen_DR;
  std::vector<float>     *Track_ECAL_x;
  std::vector<float>     *Track_ECAL_y;
  std::vector<float>     *Track_ECAL_z;
  std::vector<float>     *Track_ECAL_phi;
  std::vector<float>     *Track_ECAL_eta;
  std::vector<float>     *Track_In_eta;
  std::vector<float>     *Track_pt;
  std::vector<float>     *Track_In_phi;
  std::vector<float>     *Track_In_x;
  std::vector<float>     *Track_In_y;
  std::vector<float>     *Track_In_z;
  std::vector<float>     *Track_pt_Iphi;
  std::vector<float>     *Track_pt_Ieta;
  std::vector<float>     *Track_pt_Ephi;
  std::vector<float>     *Track_pt_Eeta;
  std::vector<int>    *Track_charge;
  std::vector<short unsigned>    *Track_nhits;
   std::vector < std::vector<float> >  *genpv_z;
   std::vector <std::vector<float> >  *pt;
   std::vector <std::vector<float> >  *dr;
   std::vector <std::vector<float> >  *label;
   std::vector <std::vector<float> >  *eta;
   std::vector <std::vector<float> >  *phi;
   std::vector <std::vector<float> >  *ecalEnergy;
   std::vector <std::vector<float> >  *hcalEnergy;
   std::vector <std::vector<float> >  *particleId;
   std::vector <std::vector<float> >  *clus_n;
   std::vector <std::vector<float> >  *gen_pt;
   std::vector <std::vector<float> >  *gen_eta;
   std::vector <std::vector<float> >  *gen_phi;
   std::vector <std::vector<float> >  *gen_pdgId;
   std::vector <std::vector<float> >  *mct_pt;
   std::vector <std::vector<float> >  *mct_eta;
   std::vector <std::vector<float> >  *mct_phi;
   std::vector <std::vector<float> >  *mct_energy;
   std::vector <std::vector<float> >  *mct_convRadius;
   std::vector <std::vector<float> >  *mct_convZ;
   std::vector <std::vector<float> >  *mct_convPhi;
   std::vector <std::vector<float> >  *mct_ele1_pt;
   std::vector <std::vector<float> >  *mct_ele1_eta;
   std::vector <std::vector<float> >  *mct_nlegs;
   std::vector <std::vector<float> >  *mct_ele1_phi;
   std::vector <std::vector<float> >  *ele1_match;
   std::vector <std::vector<float> >  *mct_ele2_pt;
   std::vector <std::vector<float> >  *mct_ele2_eta;
   std::vector <std::vector<float> >  *mct_ele2_phi;
   std::vector <std::vector<float> >  *ele2_match;
   std::vector <std::vector<float> >  *mct_eles_dr;
   std::vector <std::vector<float> >  *mct_ele1_dr;
   std::vector <std::vector<float> >  *mct_ele2_dr;
   std::vector <std::vector<float> >  *genpv_t;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_gen_DR;   //!
   TBranch        *b_Track_ECAL_x;   //!
   TBranch        *b_Track_ECAL_y;   //!
   TBranch        *b_Track_ECAL_z;   //!
   TBranch        *b_Track_ECAL_phi;   //!
   TBranch        *b_Track_ECAL_eta;   //!
   TBranch        *b_Track_In_eta;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_In_phi;   //!
   TBranch        *b_Track_In_x;   //!
   TBranch        *b_Track_In_y;   //!
   TBranch        *b_Track_In_z;   //!
   TBranch        *b_Track_pt_Iphi;   //!
   TBranch        *b_Track_pt_Ieta;   //!
   TBranch        *b_Track_pt_Ephi;   //!
   TBranch        *b_Track_pt_Eeta;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_nhits;   //!
   TBranch        *b_genpv_z;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_label;   //!
   TBranch        *b_dr;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_ecalEnergy;   //!
   TBranch        *b_hcalEnergy;   //!
   TBranch        *b_particleId;   //!
   TBranch        *b_clus_n;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_pdgId;   //!
   TBranch        *b_mct_pt;   //!
   TBranch        *b_mct_eta;   //!
   TBranch        *b_mct_phi;   //!
   TBranch        *b_mct_energy;   //!
   TBranch        *b_mct_convRadius;   //!
   TBranch        *b_mct_convZ;   //!
   TBranch        *b_mct_convPhi;   //!
   TBranch        *b_mct_ele1_pt;   //!
   TBranch        *b_mct_ele1_eta;   //!
   TBranch        *b_mct_nlegs;   //!
   TBranch        *b_mct_ele1_phi;   //!
   TBranch        *b_ele1_match;   //!
   TBranch        *b_mct_ele2_pt;   //!
   TBranch        *b_mct_ele2_eta;   //!
   TBranch        *b_mct_ele2_phi;   //!
   TBranch        *b_ele2_match;   //!
   TBranch        *b_mct_eles_dr;   //!
   TBranch        *b_mct_ele1_dr;   //!
   TBranch        *b_mct_ele2_dr;   //!
   TBranch        *b_genpv_t;   //!

  //Tclass(TTree *tree=0);
    void Init(TTree*);
/*   virtual ~Tclass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);*/
};

#endif

/*#ifdef Tclass_cxx
Tclass::Tclass(TTree* tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("tracks_tree",tree);

   }
   Init(tree);


Tclass::~Tclass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tclass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tclass::LoadTree(Long64_t entry)
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
/*void Tclass::Init(TTree* tree)
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
   Track_charge = 0;
   genpv_z = 0;
   pt = 0;
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
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("genpv_z", &genpv_z, &b_genpv_z);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("ecalEnergy", &ecalEnergy, &b_ecalEnergy);
   fChain->SetBranchAddress("hcalEnergy", &hcalEnergy, &b_hcalEnergy);
   fChain->SetBranchAddress("particleId", &particleId, &b_particleId);
   fChain->SetBranchAddress("clus_n", &clus_n, &b_clus_n);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
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
 //  Notify();
}
Bool_t Tclass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tclass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tclass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}*/
