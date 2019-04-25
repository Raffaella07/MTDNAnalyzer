//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  9 10:13:01 2019 by ROOT version 6.12/07
// from TTree gen_tree/4D TOFPID studies
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef NeutralsClass_h
#define NeutralsClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "vector"
#include "vector"

class NeutralsClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   std::vector<float>   *Track_ECAL_x;
  std::vector<float>   *Track_ECAL_y;
  std::vector<float>   *Track_ECAL_z;
  std::vector<float>   *Track_ECAL_phi;
  std::vector<float>   *Track_ECAL_eta;
  std::vector<float>   *Track_In_eta;
  std::vector<float>   *Track_pt;
  std::vector<float>   *Track_mtdt;
  std::vector<float>   *Track_genVtx_t;
  std::vector<float>   *Track_In_phi;
  std::vector<float>   *Track_In_x;
  std::vector<float>   *Track_In_y;
  std::vector<float>   *Track_In_z;
  std::vector<int>     *Track_charge;
  std::vector<unsigned short> *Track_nhits;
  std::vector<float>   *Track_pt_Iphi;
  std::vector<float>   *Track_pt_Ieta;
  std::vector<float>   *Track_pt_Ephi;
  std::vector<float>   *Track_pt_Eeta;
  std::vector<float>   *pt;
  std::vector<float>   *phi;
  std::vector<float>   *eta;
  std::vector<float>   *dr;
  std::vector<float>   *pdgId;
  std::vector<float>   *label;
  std::vector<float>   *PId;
  std::vector<float>   *primary_x;
  std::vector<float>   *primary_y;
  std::vector<float>   *primary_z;
  std::vector<float>   *primary_phi;
  std::vector<float>   *primary_eta;
  std::vector<float>   *convRadius;
  std::vector<float>   *convPhi;
  std::vector<float>   *convZ;
  std::vector<float>   *conv_x;
  std::vector<float>   *conv_y;
  std::vector<float>   *conv_z;
  std::vector<float>   *conv_phi;
  std::vector<float>   *conv_eta;
  std::vector<float>   *energy;
  std::vector<float>   *nlegs;
  std::vector<float>   *plabel;
  std::vector<float>   *ecalEnergy;
  std::vector<float>   *ecalx;
  std::vector<float>   *ecaly;
  std::vector<float>   *ecalz;
  std::vector<float>   *ecalphi;
  std::vector<float>   *ecaleta;
  std::vector<float>   *hcalEnergy;
  std::vector<float>   *clus_size;
  std::vector<float>   *clus_energy;
  std::vector<float>   *clus_time;
  std::vector<float>   *clus_eta;
  std::vector<float>   *clus_phi;
  std::vector<float>   *clus_x;
  std::vector<float>   *clus_y;
  std::vector<float>   *clus_z;
  std::vector<float>   *MTDx;
  std::vector<float>   *MTDy;
  std::vector<float>   *MTDz;
  std::vector<float>   *MTDphi;
  std::vector<float>   *MTDeta;
  std::vector<float>   *clus_size_x;
  std::vector<float>   *clus_size_y;
  std::vector<float>   *genpv_t;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_Track_ECAL_x;   //!
   TBranch        *b_Track_ECAL_y;   //!
   TBranch        *b_Track_ECAL_z;   //!
   TBranch        *b_Track_ECAL_phi;   //!
   TBranch        *b_Track_ECAL_eta;   //!
   TBranch        *b_Track_In_eta;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_mtdt;   //!
   TBranch        *b_Track_genVtx_t;   //!
   TBranch        *b_Track_In_phi;   //!
   TBranch        *b_Track_In_x;   //!
   TBranch        *b_Track_In_y;   //!
   TBranch        *b_Track_In_z;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_nhits;   //!
   TBranch        *b_Track_pt_Iphi;   //!
   TBranch        *b_Track_pt_Ieta;   //!
   TBranch        *b_Track_pt_Ephi;   //!
   TBranch        *b_Track_pt_Eeta;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_dr;   //!
   TBranch        *b_pdgId;   //!
   TBranch        *b_label;   //!
   TBranch        *b_primary_x;
   TBranch        *b_primary_y;
   TBranch        *b_primary_z;
   TBranch        *b_primary_phi;
   TBranch        *b_primary_eta;
   TBranch        *b_PId;   //!
   TBranch        *b_convRadius;   //!
   TBranch        *b_convPhi;   //!
   TBranch        *b_convZ;   //!
   TBranch        *b_conv_x;
   TBranch        *b_conv_y;
   TBranch        *b_conv_z;
   TBranch        *b_conv_phi;
   TBranch        *b_conv_eta;
   TBranch        *b_energy;   //!
   TBranch        *b_nlegs;   //!
   TBranch        *b_plabel;   //!
   TBranch        *b_ecalEnergy;   //!
   TBranch        *b_ecalx;   //!
   TBranch        *b_ecaly;   //!
   TBranch        *b_ecalz;   //!
   TBranch        *b_ecalphi;   //!
   TBranch        *b_ecaleta;   //!
   TBranch        *b_hcalEnergy;   //!
   TBranch         *b_clus_size;     
   TBranch         *b_clus_energy;
   TBranch         *b_clus_time;
   TBranch         *b_clus_eta;
   TBranch         *b_clus_phi;
   TBranch         *b_clus_x;
   TBranch         *b_clus_y;
   TBranch         *b_clus_z;
   TBranch         *b_MTDx;
   TBranch         *b_MTDy;
   TBranch         *b_MTDz;
   TBranch         *b_MTDphi;
   TBranch         *b_MTDeta;
   TBranch         *b_clus_size_x;
   TBranch         *b_clus_size_y;
   TBranch        *b_genpv_t;   //!

  // NeutralsClass(TTree *tree=0);
 // ~NeutralsClass();
  //Int_t    Cut(Long64_t entry);
  //Int_t    GetEntry(Long64_t entry);
  //Long64_t LoadTree(Long64_t entry);
  void     Init(TTree *tree);
 // void     Loop();
 // Bool_t   Notify();
 // void     Show(Long64_t entry = -1);
};

#endif

