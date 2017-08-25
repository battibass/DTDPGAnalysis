//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 16 10:08:45 2017 by ROOT version 6.08/02
// from TTree DTTree/CMSSW DT tree
// found on file: /afs/cern.ch/work/f/ferrico/private/DT_update/CMSSW_9_0_0_pre4/src/UserCode/DTDPGAnalysis/test/DTNtuple.root
//////////////////////////////////////////////////////////

#ifndef combined_trigger_studies_h
#define combined_trigger_studies_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "TClonesArray.h"
#include "TVector.h"
#include "vector"
#include "TLorentzVector.h"

class combined_trigger_studies {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.
			
  Float_t pig = acos(-1.);
  
  //********************
  // Analysis parameters
  //********************

  // Files and number of events
  TString input_file_name  = "../DTNtuple.root";
  TString output_file_name = "analysis_results.root";

  Long64_t n_events = 500000;

  // RPC matching cuts
  Float_t cluster_size_cut = 3.;
  Float_t range_strips = 4.;

  // Trigger and RECO muon cuts   	

  Float_t min_TnP_mass = 50;
  Float_t max_TnP_mass = 200;

  TString trigger_filter_name = "hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09";
  Float_t muon_hlt_dR = 0.1;

  Float_t min_pt_tag    = 26.;
  Float_t min_pt_probe  = 20.;
  Float_t max_eta_probe = 1.2; 

  Float_t Dxy_cut = 0.2;
  Float_t Dz_cut =  0.5;

  Float_t N_hits_cut =   4; 
  Float_t muchi2_cut = 10.;
  Float_t npix_cut   =   1;
  Float_t ntkr_cut   =   5;

  Float_t tkr_iso_cut = 0.10;

  // Declaration of leaf types
  Int_t           runnumber;
  Int_t           lumiblock;
  Int_t           eventNumber;

  vector<TString> *hlt_path;
  vector<TString> *hlt_filter;
  vector<float> *hlt_filter_phi;
  vector<float> *hlt_filter_eta;
  vector<float> *hlt_filter_pt;

  vector<short>   *dtsegm4D_wheel;
  vector<short>   *dtsegm4D_sector;
  vector<short>   *dtsegm4D_station;
  vector<short>   *dtsegm4D_hasPhi;
  vector<short>   *dtsegm4D_hasZed;
  vector<float>   *dtsegm4D_x_pos_loc;
  vector<float>   *dtsegm4D_y_pos_loc;
  vector<float>   *dtsegm4D_z_pos_loc;
  vector<float>   *dtsegm4D_x_dir_loc;
  vector<float>   *dtsegm4D_y_dir_loc;
  vector<float>   *dtsegm4D_z_dir_loc;
  vector<float>   *dtsegm4D_cosx;
  vector<float>   *dtsegm4D_cosy;
  vector<float>   *dtsegm4D_cosz;
  vector<float>   *dtsegm4D_phi;
  vector<float>   *dtsegm4D_theta;
  vector<float>   *dtsegm4D_eta;
  vector<float>   *dtsegm4D_t0;
  vector<float>   *dtsegm4D_vdrift;
  vector<float>   *dtsegm4D_phinormchisq;
  vector<short>   *dtsegm4D_phinhits;
  vector<float>   *dtsegm4D_znormchisq;
  vector<short>   *dtsegm4D_znhits;
  TClonesArray    *dtsegm4D_hitsExpPos;
  TClonesArray    *dtsegm4D_hitsExpWire;
  TClonesArray    *dtsegm4D_phi_hitsPos;
  TClonesArray    *dtsegm4D_phi_hitsPosCh;
  TClonesArray    *dtsegm4D_phi_hitsPosErr;
  TClonesArray    *dtsegm4D_phi_hitsSide;
  TClonesArray    *dtsegm4D_phi_hitsWire;
  TClonesArray    *dtsegm4D_phi_hitsLayer;
  TClonesArray    *dtsegm4D_phi_hitsSuperLayer;
  TClonesArray    *dtsegm4D_phi_hitsTime;
  TClonesArray    *dtsegm4D_phi_hitsTimeCali;
  TClonesArray    *dtsegm4D_z_hitsPos;
  TClonesArray    *dtsegm4D_z_hitsPosCh;
  TClonesArray    *dtsegm4D_z_hitsPosErr;
  TClonesArray    *dtsegm4D_z_hitsSide;
  TClonesArray    *dtsegm4D_z_hitsWire;
  TClonesArray    *dtsegm4D_z_hitsLayer;
  TClonesArray    *dtsegm4D_z_hitsTime;
  TClonesArray    *dtsegm4D_z_hitsTimeCali;

  vector<short>   *ltTwinMuxIn_wheel;
  vector<short>   *ltTwinMuxIn_sector;
  vector<short>   *ltTwinMuxIn_station;
  vector<short>   *ltTwinMuxIn_quality;
  vector<short>   *ltTwinMuxIn_bx;
  vector<int>     *ltTwinMuxIn_phi;
  vector<int>     *ltTwinMuxIn_phiB;
  vector<short>   *ltTwinMuxIn_is2nd;
  vector<short>   *ltTwinMuxOut_wheel;
  vector<short>   *ltTwinMuxOut_sector;
  vector<short>   *ltTwinMuxOut_station;
  vector<short>   *ltTwinMuxOut_quality;
  vector<short>   *ltTwinMuxOut_rpcbit;
  vector<short>   *ltTwinMuxOut_bx;
  vector<int>     *ltTwinMuxOut_phi;
  vector<int>     *ltTwinMuxOut_phiB;
  vector<short>   *ltTwinMuxOut_is2nd;

  vector<short>   *Mu_isMuGlobal;
  vector<short>   *Mu_isMuTracker;
  vector<int>     *Mu_numberOfChambers_sta;
  vector<int>     *Mu_numberOfMatches_sta;
  vector<int>     *Mu_numberOfHits_sta;
  vector<int>     *Mu_segmentIndex_sta;
  vector<float>   *Mu_px;
  vector<float>   *Mu_py;
  vector<float>   *Mu_pz;
  vector<float>   *Mu_phi;
  vector<float>   *Mu_eta;
  vector<short>   *Mu_recHitsSize;
  vector<float>   *Mu_normchi2_sta;
  vector<short>   *Mu_charge;
  vector<float>   *Mu_dxy_sta;
  vector<float>   *Mu_dz_sta;
  vector<float>   *Mu_normchi2_glb;
  vector<float>   *Mu_dxy_glb;
  vector<float>   *Mu_dz_glb;
  vector<int>     *Mu_numberOfPixelHits_glb;
  vector<int>     *Mu_numberOfTrackerHits_glb;
  vector<float>   *Mu_tkIsoR03_glb;
  vector<float>   *Mu_ntkIsoR03_glb;
  vector<float>   *Mu_emIsoR03_glb;
  vector<float>   *Mu_hadIsoR03_glb;
  vector<float>   *STAMu_caloCompatibility;
  vector<float>   *Mu_z_mb2_mu;
  vector<float>   *Mu_phi_mb2_mu;
  vector<float>   *Mu_pseta_mb2_mu;

  vector<float>   *Mu_x_MB1;
  vector<float>   *Mu_y_MB1;
  vector<short>   *Mu_wheel_MB1;
  vector<short>   *Mu_sector_MB1;

  vector<float>   *Mu_x_MB2;
  vector<float>   *Mu_y_MB2;
  vector<short>   *Mu_wheel_MB2;
  vector<short>   *Mu_sector_MB2;

  vector<float>   *Mu_x_MB3;
  vector<float>   *Mu_y_MB3;
  vector<short>   *Mu_wheel_MB3;
  vector<short>   *Mu_sector_MB3;

  vector<float>   *Mu_x_MB4;
  vector<float>   *Mu_y_MB4;
  vector<short>   *Mu_wheel_MB4;
  vector<short>   *Mu_sector_MB4;

  vector<int>     *RpcRecHitTwinMuxRegion;
  vector<int>     *RpcRecHitTwinMuxClusterSize;
  vector<int>     *RpcRecHitTwinMuxStrip;
  vector<int>     *RpcRecHitTwinMuxBx;
  vector<int>     *RpcRecHitTwinMuxStation;
  vector<int>     *RpcRecHitTwinMuxSector;
  vector<int>     *RpcRecHitTwinMuxLayer;
  vector<int>     *RpcRecHitTwinMuxSubsector;
  vector<int>     *RpcRecHitTwinMuxRoll;
  vector<int>     *RpcRecHitTwinMuxRing;
  vector<float>   *RpcRechitTwinMuxLocX;
  vector<float>   *RpcRechitTwinMuxLocY;
  vector<float>   *RpcRechitTwinMuxLocZ;
  vector<float>   *RpcRechitTwinMuxGlobX;
  vector<float>   *RpcRechitTwinMuxGlobY;
  vector<float>   *RpcRechitTwinMuxGlobZ;
  vector<float>   *RpcRechitTwinMuxGlobPhi;

  vector<short>   *NDTsegmentonRPC;

  TClonesArray   *DTextrapolatedOnRPCBX;
  TClonesArray   *DTextrapolatedOnRPCLocX;
  TClonesArray   *DTextrapolatedOnRPCLocY;
  TClonesArray   *DTextrapolatedOnRPCLocZ;
  TClonesArray   *DTextrapolatedOnRPCGlobX;
  TClonesArray   *DTextrapolatedOnRPCGlobY;
  TClonesArray   *DTextrapolatedOnRPCGlobZ;
  TClonesArray   *DTextrapolatedOnRPCGlobEta;
  TClonesArray   *DTextrapolatedOnRPCGlobPhi;
  TClonesArray   *DTextrapolatedOnRPCRegion;
  TClonesArray   *DTextrapolatedOnRPCSector;
  TClonesArray   *DTextrapolatedOnRPCStation;
  TClonesArray   *DTextrapolatedOnRPCLayer;
  TClonesArray   *DTextrapolatedOnRPCRoll;
  TClonesArray   *DTextrapolatedOnRPCRing;
  TClonesArray   *DTextrapolatedOnRPCStripw;

  Short_t         Ndtsegments;
  Short_t         NdtltTwinMuxOut;
  Short_t         NdtltTwinMuxIn;
  Short_t         Nmuons;
  Short_t         NhltPaths;
  Short_t         NhltFilters;
  Short_t         NirpcrechitsTwinMux;
  
  // List of branches
  TBranch        *b_runnumber;   //!
  TBranch        *b_lumiblock;   //!
  TBranch        *b_eventNumber;   //!

  TBranch        *b_hlt_path;   //!
  TBranch        *b_hlt_filter;   //!
  TBranch        *b_hlt_filter_phi;   //!
  TBranch        *b_hlt_filter_eta;   //!
  TBranch        *b_hlt_filter_pt;    //!

  TBranch        *b_dtsegm4D_wheel;   //!
  TBranch        *b_dtsegm4D_sector;   //!
  TBranch        *b_dtsegm4D_station;   //!
  TBranch        *b_dtsegm4D_hasPhi;   //!
  TBranch        *b_dtsegm4D_hasZed;   //!
  TBranch        *b_dtsegm4D_x_pos_loc;   //!
  TBranch        *b_dtsegm4D_y_pos_loc;   //!
  TBranch        *b_dtsegm4D_z_pos_loc;   //!
  TBranch        *b_dtsegm4D_x_dir_loc;   //!
  TBranch        *b_dtsegm4D_y_dir_loc;   //!
  TBranch        *b_dtsegm4D_z_dir_loc;   //!
  TBranch        *b_dtsegm4D_cosx;   //!
  TBranch        *b_dtsegm4D_cosy;   //!
  TBranch        *b_dtsegm4D_cosz;   //!
  TBranch        *b_dtsegm4D_phi;   //!
  TBranch        *b_dtsegm4D_theta;   //!
  TBranch        *b_dtsegm4D_eta;   //!
  TBranch        *b_dtsegm4D_t0;   //!
  TBranch        *b_dtsegm4D_vdrift;   //!
  TBranch        *b_dtsegm4D_phinormchisq;   //!
  TBranch        *b_dtsegm4D_phinhits;   //!
  TBranch        *b_dtsegm4D_znormchisq;   //!
  TBranch        *b_dtsegm4D_znhits;   //!
  TBranch        *b_dtsegm4D_hitsExpPos;   //!
  TBranch        *b_dtsegm4D_hitsExpWire;   //!
  TBranch        *b_dtsegm4D_phi_hitsPos;   //!
  TBranch        *b_dtsegm4D_phi_hitsPosCh;   //!
  TBranch        *b_dtsegm4D_phi_hitsPosErr;   //!
  TBranch        *b_dtsegm4D_phi_hitsSide;   //!
  TBranch        *b_dtsegm4D_phi_hitsWire;   //!
  TBranch        *b_dtsegm4D_phi_hitsLayer;   //!
  TBranch        *b_dtsegm4D_phi_hitsSuperLayer;   //!
  TBranch        *b_dtsegm4D_phi_hitsTime;   //!
  TBranch        *b_dtsegm4D_phi_hitsTimeCali;   //!
  TBranch        *b_dtsegm4D_z_hitsPos;   //!
  TBranch        *b_dtsegm4D_z_hitsPosCh;   //!
  TBranch        *b_dtsegm4D_z_hitsPosErr;   //!
  TBranch        *b_dtsegm4D_z_hitsSide;   //!
  TBranch        *b_dtsegm4D_z_hitsWire;   //!
  TBranch        *b_dtsegm4D_z_hitsLayer;   //!
  TBranch        *b_dtsegm4D_z_hitsTime;   //!
  TBranch        *b_dtsegm4D_z_hitsTimeCali;   //!

  TBranch        *b_ltTwinMuxIn_wheel;   //!
  TBranch        *b_ltTwinMuxIn_sector;   //!
  TBranch        *b_ltTwinMuxIn_station;   //!
  TBranch        *b_ltTwinMuxIn_quality;   //!
  TBranch        *b_ltTwinMuxIn_bx;   //!
  TBranch        *b_ltTwinMuxIn_phi;   //!
  TBranch        *b_ltTwinMuxIn_phiB;   //!
  TBranch        *b_ltTwinMuxIn_is2nd;   //!
  TBranch        *b_ltTwinMuxOut_wheel;   //!
  TBranch        *b_ltTwinMuxOut_sector;   //!
  TBranch        *b_ltTwinMuxOut_station;   //!
  TBranch        *b_ltTwinMuxOut_quality;   //!
  TBranch        *b_ltTwinMuxOut_rpcbit;   //!
  TBranch        *b_ltTwinMuxOut_bx;   //!
  TBranch        *b_ltTwinMuxOut_phi;   //!
  TBranch        *b_ltTwinMuxOut_phiB;   //!
  TBranch        *b_ltTwinMuxOut_is2nd;   //!

  TBranch        *b_Mu_isMuGlobal;   //!
  TBranch        *b_Mu_isMuTracker;   //!
  TBranch        *b_Mu_numberOfChambers_sta;   //!
  TBranch        *b_Mu_numberOfMatches_sta;   //!
  TBranch        *b_Mu_numberOfHits_sta;   //!
  TBranch        *b_Mu_segmentIndex_sta;   //!
  TBranch        *b_Mu_px;   //!
  TBranch        *b_Mu_py;   //!
  TBranch        *b_Mu_pz;   //!
  TBranch        *b_Mu_phi;   //!
  TBranch        *b_Mu_eta;   //!
  TBranch        *b_Mu_recHitsSize;   //!
  TBranch        *b_Mu_normchi2_sta;   //!
  TBranch        *b_Mu_charge;   //!
  TBranch        *b_Mu_dxy_sta;   //!
  TBranch        *b_Mu_dz_sta;   //!
  TBranch        *b_Mu_normchi2_glb;   //!
  TBranch        *b_Mu_dxy_glb;   //!
  TBranch        *b_Mu_dz_glb;   //!
  TBranch        *b_Mu_numberOfPixelHits_glb;   //!
  TBranch        *b_Mu_numberOfTrackerHits_glb;   //!
  TBranch        *b_Mu_tkIsoR03_glb;   //!
  TBranch        *b_Mu_ntkIsoR03_glb;   //!
  TBranch        *b_Mu_emIsoR03_glb;   //!
  TBranch        *b_Mu_hadIsoR03_glb;   //!
  TBranch        *b_STAMu_caloCompatibility;   //!
  TBranch        *b_Mu_z_mb2_mu;   //!
  TBranch        *b_Mu_phi_mb2_mu;   //!
  TBranch        *b_Mu_pseta_mb2_mu;   //!
  
  TBranch   *b_Mu_x_MB1;
  TBranch   *b_Mu_y_MB1;
  TBranch   *b_Mu_wheel_MB1;
  TBranch   *b_Mu_sector_MB1;
  
  TBranch   *b_Mu_x_MB2;
  TBranch   *b_Mu_y_MB2;
  TBranch   *b_Mu_wheel_MB2;
  TBranch   *b_Mu_sector_MB2;
  
  TBranch   *b_Mu_x_MB3;
  TBranch   *b_Mu_y_MB3;
  TBranch   *b_Mu_wheel_MB3;
  TBranch   *b_Mu_sector_MB3;
  
  TBranch   *b_Mu_x_MB4;
  TBranch   *b_Mu_y_MB4;
  TBranch   *b_Mu_wheel_MB4;
  TBranch   *b_Mu_sector_MB4;

  TBranch        *b_RpcRecHitTwinMuxRegion;   //!
  TBranch        *b_RpcRecHitTwinMuxClusterSize;   //!
  TBranch        *b_RpcRecHitTwinMuxStrip;   //!
  TBranch        *b_RpcRecHitTwinMuxBx;   //!
  TBranch        *b_RpcRecHitTwinMuxStation;   //!
  TBranch        *b_RpcRecHitTwinMuxSector;   //!
  TBranch        *b_RpcRecHitTwinMuxLayer;   //!
  TBranch        *b_RpcRecHitTwinMuxSubsector;   //!
  TBranch        *b_RpcRecHitTwinMuxRoll;   //!
  TBranch        *b_RpcRecHitTwinMuxRing;   //!
  TBranch        *b_RpcRechitTwinMuxLocX;   //!
  TBranch        *b_RpcRechitTwinMuxLocY;   //!
  TBranch        *b_RpcRechitTwinMuxLocZ;   //!
  TBranch        *b_RpcRechitTwinMuxGlobX;   //!
  TBranch        *b_RpcRechitTwinMuxGlobY;   //!
  TBranch        *b_RpcRechitTwinMuxGlobZ;   //!
  TBranch        *b_RpcRechitTwinMuxGlobPhi;   //!

  TBranch        *b_DTextrapolatedOnRPCBX;   //!
  TBranch        *b_DTextrapolatedOnRPCLocX;   //!
  TBranch        *b_DTextrapolatedOnRPCLocY;   //!
  TBranch        *b_DTextrapolatedOnRPCLocZ;   //!
  TBranch        *b_DTextrapolatedOnRPCGlobX;   //!
  TBranch        *b_DTextrapolatedOnRPCGlobY;   //!
  TBranch        *b_DTextrapolatedOnRPCGlobZ;   //!
  TBranch        *b_DTextrapolatedOnRPCGlobEta;   //!
  TBranch        *b_DTextrapolatedOnRPCGlobPhi;   //!
  TBranch        *b_DTextrapolatedOnRPCRegion;   //!
  TBranch        *b_DTextrapolatedOnRPCSector;   //!
  TBranch        *b_DTextrapolatedOnRPCStation;   //!
  TBranch        *b_DTextrapolatedOnRPCLayer;   //!
  TBranch        *b_DTextrapolatedOnRPCRoll;   //!
  TBranch        *b_DTextrapolatedOnRPCRing;   //!
  TBranch        *b_DTextrapolatedOnRPCStripw;   //!

  TBranch        *b_Ndtsegments;   //!
  TBranch        *b_NdtltTwinMuxOut;   //!
  TBranch        *b_NdtltTwinMuxIn;   //!
  TBranch        *b_Nmuons;   //!

  TBranch        *b_NhltFilters;   //!
  TBranch        *b_NhltPaths;   //!   

  TBranch        *b_NirpcrechitsTwinMux;   //!
  TBranch        *b_NDTsegmentonRPC;   //!
  
   combined_trigger_studies(TTree *tree=0);
   virtual ~combined_trigger_studies();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(); 

   template<typename T> T getXY(TClonesArray * arr, int x, int y) { return static_cast<T>((*((TVectorT<float> *)(arr->At(x))))[y]); };
    
   float PhiConversion(int, int);
   std::vector<std::pair<int,int>> TnPSelection(Float_t, Float_t);
   bool HasTrigger(const TLorentzVector &, Float_t);

   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef combined_trigger_studies_cxx
combined_trigger_studies::combined_trigger_studies(TTree *tree) : fChain(0) 
{
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(input_file_name);
    
    if (!f || !f->IsOpen()) {
      f = new TFile(input_file_name);
    }
    f->GetObject("DTTree",tree);
    
  }
  Init(tree);
}

combined_trigger_studies::~combined_trigger_studies()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t combined_trigger_studies::GetEntry(Long64_t entry)
{
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t combined_trigger_studies::LoadTree(Long64_t entry)
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

void combined_trigger_studies::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set object pointer
  hlt_path = 0;
  hlt_filter = 0;
  hlt_filter_phi = 0;
  hlt_filter_eta = 0;
  hlt_filter_pt  = 0;

  dtsegm4D_wheel = 0;
  dtsegm4D_sector = 0;
  dtsegm4D_station = 0;
  dtsegm4D_hasPhi = 0;
  dtsegm4D_hasZed = 0;
  dtsegm4D_x_pos_loc = 0;
  dtsegm4D_y_pos_loc = 0;
  dtsegm4D_z_pos_loc = 0;
  dtsegm4D_x_dir_loc = 0;
  dtsegm4D_y_dir_loc = 0;
  dtsegm4D_z_dir_loc = 0;
  dtsegm4D_cosx = 0;
  dtsegm4D_cosy = 0;
  dtsegm4D_cosz = 0;
  dtsegm4D_phi = 0;
  dtsegm4D_theta = 0;
  dtsegm4D_eta = 0;
  dtsegm4D_t0 = 0;
  dtsegm4D_vdrift = 0;
  dtsegm4D_phinormchisq = 0;
  dtsegm4D_phinhits = 0;
  dtsegm4D_znormchisq = 0;
  dtsegm4D_znhits = 0;
  dtsegm4D_hitsExpPos = 0;
  dtsegm4D_hitsExpWire = 0;
  dtsegm4D_phi_hitsPos = 0;
  dtsegm4D_phi_hitsPosCh = 0;
  dtsegm4D_phi_hitsPosErr = 0;
  dtsegm4D_phi_hitsSide = 0;
  dtsegm4D_phi_hitsWire = 0;
  dtsegm4D_phi_hitsLayer = 0;
  dtsegm4D_phi_hitsSuperLayer = 0;
  dtsegm4D_phi_hitsTime = 0;
  dtsegm4D_phi_hitsTimeCali = 0;
  dtsegm4D_z_hitsPos = 0;
  dtsegm4D_z_hitsPosCh = 0;
  dtsegm4D_z_hitsPosErr = 0;
  dtsegm4D_z_hitsSide = 0;
  dtsegm4D_z_hitsWire = 0;
  dtsegm4D_z_hitsLayer = 0;
  dtsegm4D_z_hitsTime = 0;
  dtsegm4D_z_hitsTimeCali = 0;

  ltTwinMuxIn_wheel = 0;
  ltTwinMuxIn_sector = 0;
  ltTwinMuxIn_station = 0;
  ltTwinMuxIn_quality = 0;
  ltTwinMuxIn_bx = 0;
  ltTwinMuxIn_phi = 0;
  ltTwinMuxIn_phiB = 0;
  ltTwinMuxIn_is2nd = 0;
  ltTwinMuxOut_wheel = 0;
  ltTwinMuxOut_sector = 0;
  ltTwinMuxOut_station = 0;
  ltTwinMuxOut_quality = 0;
  ltTwinMuxOut_rpcbit = 0;
  ltTwinMuxOut_bx = 0;
  ltTwinMuxOut_phi = 0;
  ltTwinMuxOut_phiB = 0;
  ltTwinMuxOut_is2nd = 0;

  Mu_isMuGlobal = 0;
  Mu_isMuTracker = 0;
  Mu_numberOfChambers_sta = 0;
  Mu_numberOfMatches_sta = 0;
  Mu_numberOfHits_sta = 0;
  Mu_segmentIndex_sta = 0;
  Mu_px = 0;
  Mu_py = 0;
  Mu_pz = 0;
  Mu_phi = 0;
  Mu_eta = 0;
  Mu_recHitsSize = 0;
  Mu_normchi2_sta = 0;
  Mu_charge = 0;
  Mu_dxy_sta = 0;
  Mu_dz_sta = 0;
  Mu_normchi2_glb = 0;
  Mu_dxy_glb = 0;
  Mu_dz_glb = 0;
  Mu_numberOfPixelHits_glb = 0;
  Mu_numberOfTrackerHits_glb = 0;
  Mu_tkIsoR03_glb = 0;
  Mu_ntkIsoR03_glb = 0;
  Mu_emIsoR03_glb = 0;
  Mu_hadIsoR03_glb = 0;
  STAMu_caloCompatibility = 0;
  Mu_z_mb2_mu = 0;
  Mu_phi_mb2_mu = 0;
  Mu_pseta_mb2_mu = 0;
  Mu_x_MB1 = 0;
  Mu_y_MB1 = 0;
  Mu_wheel_MB1  = 0;
  Mu_sector_MB1 = 0;
  Mu_x_MB2 = 0;
  Mu_y_MB2 = 0;
  Mu_wheel_MB2  = 0;
  Mu_sector_MB2 = 0;
  Mu_x_MB3 = 0;
  Mu_y_MB3 = 0;
  Mu_wheel_MB3  = 0;
  Mu_sector_MB3 = 0;
  Mu_x_MB4 = 0;
  Mu_y_MB4 = 0;
  Mu_wheel_MB4  = 0;
  Mu_sector_MB4 = 0;
  

  RpcRecHitTwinMuxRegion = 0;
  RpcRecHitTwinMuxClusterSize = 0;
  RpcRecHitTwinMuxStrip = 0;
  RpcRecHitTwinMuxBx = 0;
  RpcRecHitTwinMuxStation = 0;
  RpcRecHitTwinMuxSector = 0;
  RpcRecHitTwinMuxLayer = 0;
  RpcRecHitTwinMuxSubsector = 0;
  RpcRecHitTwinMuxRoll = 0;
  RpcRecHitTwinMuxRing = 0;
  RpcRechitTwinMuxLocX = 0;
  RpcRechitTwinMuxLocY = 0;
  RpcRechitTwinMuxLocZ = 0;
  RpcRechitTwinMuxGlobX = 0;
  RpcRechitTwinMuxGlobY = 0;
  RpcRechitTwinMuxGlobZ = 0;
  RpcRechitTwinMuxGlobPhi = 0;

  DTextrapolatedOnRPCBX = 0;
  DTextrapolatedOnRPCLocX = 0;
  DTextrapolatedOnRPCLocY = 0;
  DTextrapolatedOnRPCLocZ = 0;
  DTextrapolatedOnRPCGlobX = 0;
  DTextrapolatedOnRPCGlobY = 0;
  DTextrapolatedOnRPCGlobZ = 0;
  DTextrapolatedOnRPCGlobEta = 0;
  DTextrapolatedOnRPCGlobPhi = 0;
  DTextrapolatedOnRPCRegion = 0;
  DTextrapolatedOnRPCSector = 0;
  DTextrapolatedOnRPCStation = 0;
  DTextrapolatedOnRPCLayer = 0;
  DTextrapolatedOnRPCRoll = 0;
  DTextrapolatedOnRPCRing = 0; 
  DTextrapolatedOnRPCStripw = 0;   

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
  fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
  fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);

  fChain->SetBranchAddress("hlt_path", &hlt_path, &b_hlt_path);
  fChain->SetBranchAddress("hlt_filter", &hlt_filter, &b_hlt_filter);
  fChain->SetBranchAddress("hlt_filter_phi", &hlt_filter_phi, &b_hlt_filter_phi);
  fChain->SetBranchAddress("hlt_filter_eta", &hlt_filter_eta, &b_hlt_filter_eta);
  fChain->SetBranchAddress("hlt_filter_pt",  &hlt_filter_pt, &b_hlt_filter_pt);

  fChain->SetBranchAddress("dtsegm4D_wheel", &dtsegm4D_wheel, &b_dtsegm4D_wheel);
  fChain->SetBranchAddress("dtsegm4D_sector", &dtsegm4D_sector, &b_dtsegm4D_sector);
  fChain->SetBranchAddress("dtsegm4D_station", &dtsegm4D_station, &b_dtsegm4D_station);
  fChain->SetBranchAddress("dtsegm4D_hasPhi", &dtsegm4D_hasPhi, &b_dtsegm4D_hasPhi);
  fChain->SetBranchAddress("dtsegm4D_hasZed", &dtsegm4D_hasZed, &b_dtsegm4D_hasZed);
  fChain->SetBranchAddress("dtsegm4D_x_pos_loc", &dtsegm4D_x_pos_loc, &b_dtsegm4D_x_pos_loc);
  fChain->SetBranchAddress("dtsegm4D_y_pos_loc", &dtsegm4D_y_pos_loc, &b_dtsegm4D_y_pos_loc);
  fChain->SetBranchAddress("dtsegm4D_z_pos_loc", &dtsegm4D_z_pos_loc, &b_dtsegm4D_z_pos_loc);
  fChain->SetBranchAddress("dtsegm4D_x_dir_loc", &dtsegm4D_x_dir_loc, &b_dtsegm4D_x_dir_loc);
  fChain->SetBranchAddress("dtsegm4D_y_dir_loc", &dtsegm4D_y_dir_loc, &b_dtsegm4D_y_dir_loc);
  fChain->SetBranchAddress("dtsegm4D_z_dir_loc", &dtsegm4D_z_dir_loc, &b_dtsegm4D_z_dir_loc);
  fChain->SetBranchAddress("dtsegm4D_cosx", &dtsegm4D_cosx, &b_dtsegm4D_cosx);
  fChain->SetBranchAddress("dtsegm4D_cosy", &dtsegm4D_cosy, &b_dtsegm4D_cosy);
  fChain->SetBranchAddress("dtsegm4D_cosz", &dtsegm4D_cosz, &b_dtsegm4D_cosz);
  fChain->SetBranchAddress("dtsegm4D_phi", &dtsegm4D_phi, &b_dtsegm4D_phi);
  fChain->SetBranchAddress("dtsegm4D_theta", &dtsegm4D_theta, &b_dtsegm4D_theta);
  fChain->SetBranchAddress("dtsegm4D_eta", &dtsegm4D_eta, &b_dtsegm4D_eta);
  fChain->SetBranchAddress("dtsegm4D_t0", &dtsegm4D_t0, &b_dtsegm4D_t0);
  fChain->SetBranchAddress("dtsegm4D_vdrift", &dtsegm4D_vdrift, &b_dtsegm4D_vdrift);
  fChain->SetBranchAddress("dtsegm4D_phinormchisq", &dtsegm4D_phinormchisq, &b_dtsegm4D_phinormchisq);
  fChain->SetBranchAddress("dtsegm4D_phinhits", &dtsegm4D_phinhits, &b_dtsegm4D_phinhits);
  fChain->SetBranchAddress("dtsegm4D_znormchisq", &dtsegm4D_znormchisq, &b_dtsegm4D_znormchisq);
  fChain->SetBranchAddress("dtsegm4D_znhits", &dtsegm4D_znhits, &b_dtsegm4D_znhits);
  fChain->SetBranchAddress("dtsegm4D_hitsExpPos", &dtsegm4D_hitsExpPos, &b_dtsegm4D_hitsExpPos);
  fChain->SetBranchAddress("dtsegm4D_hitsExpWire", &dtsegm4D_hitsExpWire, &b_dtsegm4D_hitsExpWire);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsPos", &dtsegm4D_phi_hitsPos, &b_dtsegm4D_phi_hitsPos);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsPosCh", &dtsegm4D_phi_hitsPosCh, &b_dtsegm4D_phi_hitsPosCh);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsPosErr", &dtsegm4D_phi_hitsPosErr, &b_dtsegm4D_phi_hitsPosErr);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsSide", &dtsegm4D_phi_hitsSide, &b_dtsegm4D_phi_hitsSide);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsWire", &dtsegm4D_phi_hitsWire, &b_dtsegm4D_phi_hitsWire);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsLayer", &dtsegm4D_phi_hitsLayer, &b_dtsegm4D_phi_hitsLayer);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsSuperLayer", &dtsegm4D_phi_hitsSuperLayer, &b_dtsegm4D_phi_hitsSuperLayer);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsTime", &dtsegm4D_phi_hitsTime, &b_dtsegm4D_phi_hitsTime);
  fChain->SetBranchAddress("dtsegm4D_phi_hitsTimeCali", &dtsegm4D_phi_hitsTimeCali, &b_dtsegm4D_phi_hitsTimeCali);
  fChain->SetBranchAddress("dtsegm4D_z_hitsPos", &dtsegm4D_z_hitsPos, &b_dtsegm4D_z_hitsPos);
  fChain->SetBranchAddress("dtsegm4D_z_hitsPosCh", &dtsegm4D_z_hitsPosCh, &b_dtsegm4D_z_hitsPosCh);
  fChain->SetBranchAddress("dtsegm4D_z_hitsPosErr", &dtsegm4D_z_hitsPosErr, &b_dtsegm4D_z_hitsPosErr);
  fChain->SetBranchAddress("dtsegm4D_z_hitsSide", &dtsegm4D_z_hitsSide, &b_dtsegm4D_z_hitsSide);
  fChain->SetBranchAddress("dtsegm4D_z_hitsWire", &dtsegm4D_z_hitsWire, &b_dtsegm4D_z_hitsWire);
  fChain->SetBranchAddress("dtsegm4D_z_hitsLayer", &dtsegm4D_z_hitsLayer, &b_dtsegm4D_z_hitsLayer);
  fChain->SetBranchAddress("dtsegm4D_z_hitsTime", &dtsegm4D_z_hitsTime, &b_dtsegm4D_z_hitsTime);
  fChain->SetBranchAddress("dtsegm4D_z_hitsTimeCali", &dtsegm4D_z_hitsTimeCali, &b_dtsegm4D_z_hitsTimeCali);
  
  fChain->SetBranchAddress("ltTwinMuxIn_wheel", &ltTwinMuxIn_wheel, &b_ltTwinMuxIn_wheel);
  fChain->SetBranchAddress("ltTwinMuxIn_sector", &ltTwinMuxIn_sector, &b_ltTwinMuxIn_sector);
  fChain->SetBranchAddress("ltTwinMuxIn_station", &ltTwinMuxIn_station, &b_ltTwinMuxIn_station);
  fChain->SetBranchAddress("ltTwinMuxIn_quality", &ltTwinMuxIn_quality, &b_ltTwinMuxIn_quality);
  fChain->SetBranchAddress("ltTwinMuxIn_bx", &ltTwinMuxIn_bx, &b_ltTwinMuxIn_bx);
  fChain->SetBranchAddress("ltTwinMuxIn_phi", &ltTwinMuxIn_phi, &b_ltTwinMuxIn_phi);
  fChain->SetBranchAddress("ltTwinMuxIn_phiB", &ltTwinMuxIn_phiB, &b_ltTwinMuxIn_phiB);
  fChain->SetBranchAddress("ltTwinMuxIn_is2nd", &ltTwinMuxIn_is2nd, &b_ltTwinMuxIn_is2nd);
  fChain->SetBranchAddress("ltTwinMuxOut_wheel", &ltTwinMuxOut_wheel, &b_ltTwinMuxOut_wheel);
  fChain->SetBranchAddress("ltTwinMuxOut_sector", &ltTwinMuxOut_sector, &b_ltTwinMuxOut_sector);
  fChain->SetBranchAddress("ltTwinMuxOut_station", &ltTwinMuxOut_station, &b_ltTwinMuxOut_station);
  fChain->SetBranchAddress("ltTwinMuxOut_quality", &ltTwinMuxOut_quality, &b_ltTwinMuxOut_quality);
  fChain->SetBranchAddress("ltTwinMuxOut_rpcbit", &ltTwinMuxOut_rpcbit, &b_ltTwinMuxOut_rpcbit);
  fChain->SetBranchAddress("ltTwinMuxOut_bx", &ltTwinMuxOut_bx, &b_ltTwinMuxOut_bx);
  fChain->SetBranchAddress("ltTwinMuxOut_phi", &ltTwinMuxOut_phi, &b_ltTwinMuxOut_phi);
  fChain->SetBranchAddress("ltTwinMuxOut_phiB", &ltTwinMuxOut_phiB, &b_ltTwinMuxOut_phiB);
  fChain->SetBranchAddress("ltTwinMuxOut_is2nd", &ltTwinMuxOut_is2nd, &b_ltTwinMuxOut_is2nd);
  
  fChain->SetBranchAddress("Mu_isMuGlobal", &Mu_isMuGlobal, &b_Mu_isMuGlobal);
  fChain->SetBranchAddress("Mu_isMuTracker", &Mu_isMuTracker, &b_Mu_isMuTracker);
  fChain->SetBranchAddress("Mu_numberOfChambers_sta", &Mu_numberOfChambers_sta, &b_Mu_numberOfChambers_sta);
  fChain->SetBranchAddress("Mu_numberOfMatches_sta", &Mu_numberOfMatches_sta, &b_Mu_numberOfMatches_sta);
  fChain->SetBranchAddress("Mu_numberOfHits_sta", &Mu_numberOfHits_sta, &b_Mu_numberOfHits_sta);
  fChain->SetBranchAddress("Mu_segmentIndex_sta", &Mu_segmentIndex_sta, &b_Mu_segmentIndex_sta);
  fChain->SetBranchAddress("Mu_px", &Mu_px, &b_Mu_px);
  fChain->SetBranchAddress("Mu_py", &Mu_py, &b_Mu_py);
  fChain->SetBranchAddress("Mu_pz", &Mu_pz, &b_Mu_pz);
  fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
  fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
  fChain->SetBranchAddress("Mu_recHitsSize", &Mu_recHitsSize, &b_Mu_recHitsSize);
  fChain->SetBranchAddress("Mu_normchi2_sta", &Mu_normchi2_sta, &b_Mu_normchi2_sta);
  fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
  fChain->SetBranchAddress("Mu_dxy_sta", &Mu_dxy_sta, &b_Mu_dxy_sta);
  fChain->SetBranchAddress("Mu_dz_sta", &Mu_dz_sta, &b_Mu_dz_sta);
  fChain->SetBranchAddress("Mu_normchi2_glb", &Mu_normchi2_glb, &b_Mu_normchi2_glb);
  fChain->SetBranchAddress("Mu_dxy_glb", &Mu_dxy_glb, &b_Mu_dxy_glb);
  fChain->SetBranchAddress("Mu_dz_glb", &Mu_dz_glb, &b_Mu_dz_glb);
  fChain->SetBranchAddress("Mu_numberOfPixelHits_glb", &Mu_numberOfPixelHits_glb, &b_Mu_numberOfPixelHits_glb);
  fChain->SetBranchAddress("Mu_numberOfTrackerHits_glb", &Mu_numberOfTrackerHits_glb, &b_Mu_numberOfTrackerHits_glb);
  fChain->SetBranchAddress("Mu_tkIsoR03_glb", &Mu_tkIsoR03_glb, &b_Mu_tkIsoR03_glb);
  fChain->SetBranchAddress("Mu_ntkIsoR03_glb", &Mu_ntkIsoR03_glb, &b_Mu_ntkIsoR03_glb);
  fChain->SetBranchAddress("Mu_emIsoR03_glb", &Mu_emIsoR03_glb, &b_Mu_emIsoR03_glb);
  fChain->SetBranchAddress("Mu_hadIsoR03_glb", &Mu_hadIsoR03_glb, &b_Mu_hadIsoR03_glb);
  fChain->SetBranchAddress("STAMu_caloCompatibility", &STAMu_caloCompatibility, &b_STAMu_caloCompatibility);
  fChain->SetBranchAddress("Mu_z_mb2_mu", &Mu_z_mb2_mu, &b_Mu_z_mb2_mu);
  fChain->SetBranchAddress("Mu_phi_mb2_mu", &Mu_phi_mb2_mu, &b_Mu_phi_mb2_mu);
  fChain->SetBranchAddress("Mu_pseta_mb2_mu", &Mu_pseta_mb2_mu, &b_Mu_pseta_mb2_mu);
  fChain->SetBranchAddress("TRKMu_x_MB1", &Mu_x_MB1, &b_Mu_x_MB1);
  fChain->SetBranchAddress("TRKMu_y_MB1", &Mu_y_MB1, &b_Mu_y_MB1);
  fChain->SetBranchAddress("TRKMu_wheel_MB1", &Mu_wheel_MB1, &b_Mu_wheel_MB1);
  fChain->SetBranchAddress("TRKMu_sector_MB1", &Mu_sector_MB1, &b_Mu_sector_MB1);
  fChain->SetBranchAddress("TRKMu_x_MB2", &Mu_x_MB2, &b_Mu_x_MB2);
  fChain->SetBranchAddress("TRKMu_y_MB2", &Mu_y_MB2, &b_Mu_y_MB2);
  fChain->SetBranchAddress("TRKMu_wheel_MB2", &Mu_wheel_MB2, &b_Mu_wheel_MB2);
  fChain->SetBranchAddress("TRKMu_sector_MB2", &Mu_sector_MB2, &b_Mu_sector_MB2);
  fChain->SetBranchAddress("TRKMu_x_MB3", &Mu_x_MB3, &b_Mu_x_MB3);
  fChain->SetBranchAddress("TRKMu_y_MB3", &Mu_y_MB3, &b_Mu_y_MB3);
  fChain->SetBranchAddress("TRKMu_wheel_MB3", &Mu_wheel_MB3, &b_Mu_wheel_MB3);
  fChain->SetBranchAddress("TRKMu_sector_MB3", &Mu_sector_MB3, &b_Mu_sector_MB3);
  fChain->SetBranchAddress("TRKMu_x_MB4", &Mu_x_MB4, &b_Mu_x_MB4);
  fChain->SetBranchAddress("TRKMu_y_MB4", &Mu_y_MB4, &b_Mu_y_MB4);
  fChain->SetBranchAddress("TRKMu_wheel_MB4", &Mu_wheel_MB4, &b_Mu_wheel_MB4);
  fChain->SetBranchAddress("TRKMu_sector_MB4", &Mu_sector_MB4, &b_Mu_sector_MB4);

  fChain->SetBranchAddress("RpcRecHitTwinMuxRegion", &RpcRecHitTwinMuxRegion, &b_RpcRecHitTwinMuxRegion);
  fChain->SetBranchAddress("RpcRecHitTwinMuxClusterSize", &RpcRecHitTwinMuxClusterSize, &b_RpcRecHitTwinMuxClusterSize);
  fChain->SetBranchAddress("RpcRecHitTwinMuxStrip", &RpcRecHitTwinMuxStrip, &b_RpcRecHitTwinMuxStrip);
  fChain->SetBranchAddress("RpcRecHitTwinMuxBx", &RpcRecHitTwinMuxBx, &b_RpcRecHitTwinMuxBx);
  fChain->SetBranchAddress("RpcRecHitTwinMuxStation", &RpcRecHitTwinMuxStation, &b_RpcRecHitTwinMuxStation);
  fChain->SetBranchAddress("RpcRecHitTwinMuxSector", &RpcRecHitTwinMuxSector, &b_RpcRecHitTwinMuxSector);
  fChain->SetBranchAddress("RpcRecHitTwinMuxLayer", &RpcRecHitTwinMuxLayer, &b_RpcRecHitTwinMuxLayer);
  fChain->SetBranchAddress("RpcRecHitTwinMuxSubsector", &RpcRecHitTwinMuxSubsector, &b_RpcRecHitTwinMuxSubsector);
  fChain->SetBranchAddress("RpcRecHitTwinMuxRoll", &RpcRecHitTwinMuxRoll, &b_RpcRecHitTwinMuxRoll);
  fChain->SetBranchAddress("RpcRecHitTwinMuxRing", &RpcRecHitTwinMuxRing, &b_RpcRecHitTwinMuxRing);
  fChain->SetBranchAddress("RpcRechitTwinMuxLocX", &RpcRechitTwinMuxLocX, &b_RpcRechitTwinMuxLocX);
  fChain->SetBranchAddress("RpcRechitTwinMuxLocY", &RpcRechitTwinMuxLocY, &b_RpcRechitTwinMuxLocY);
  fChain->SetBranchAddress("RpcRechitTwinMuxLocZ", &RpcRechitTwinMuxLocZ, &b_RpcRechitTwinMuxLocZ);
  fChain->SetBranchAddress("RpcRechitTwinMuxGlobX", &RpcRechitTwinMuxGlobX, &b_RpcRechitTwinMuxGlobX);
  fChain->SetBranchAddress("RpcRechitTwinMuxGlobY", &RpcRechitTwinMuxGlobY, &b_RpcRechitTwinMuxGlobY);
  fChain->SetBranchAddress("RpcRechitTwinMuxGlobZ", &RpcRechitTwinMuxGlobZ, &b_RpcRechitTwinMuxGlobZ);
  fChain->SetBranchAddress("RpcRechitTwinMuxGlobPhi", &RpcRechitTwinMuxGlobPhi, &b_RpcRechitTwinMuxGlobPhi);

  fChain->SetBranchAddress("DTextrapolatedOnRPCBX", &DTextrapolatedOnRPCBX, &b_DTextrapolatedOnRPCBX);
  fChain->SetBranchAddress("DTextrapolatedOnRPCLocX", &DTextrapolatedOnRPCLocX, &b_DTextrapolatedOnRPCLocX);
  fChain->SetBranchAddress("DTextrapolatedOnRPCLocY", &DTextrapolatedOnRPCLocY, &b_DTextrapolatedOnRPCLocY);
  fChain->SetBranchAddress("DTextrapolatedOnRPCLocZ", &DTextrapolatedOnRPCLocZ, &b_DTextrapolatedOnRPCLocZ);
  fChain->SetBranchAddress("DTextrapolatedOnRPCGlobX", &DTextrapolatedOnRPCGlobX, &b_DTextrapolatedOnRPCGlobX);
  fChain->SetBranchAddress("DTextrapolatedOnRPCGlobY", &DTextrapolatedOnRPCGlobY, &b_DTextrapolatedOnRPCGlobY);
  fChain->SetBranchAddress("DTextrapolatedOnRPCGlobZ", &DTextrapolatedOnRPCGlobZ, &b_DTextrapolatedOnRPCGlobZ);
  fChain->SetBranchAddress("DTextrapolatedOnRPCGlobEta", &DTextrapolatedOnRPCGlobEta, &b_DTextrapolatedOnRPCGlobEta);
  fChain->SetBranchAddress("DTextrapolatedOnRPCGlobPhi", &DTextrapolatedOnRPCGlobPhi, &b_DTextrapolatedOnRPCGlobPhi);
  fChain->SetBranchAddress("DTextrapolatedOnRPCRegion", &DTextrapolatedOnRPCRegion, &b_DTextrapolatedOnRPCRegion);
  fChain->SetBranchAddress("DTextrapolatedOnRPCSector", &DTextrapolatedOnRPCSector, &b_DTextrapolatedOnRPCSector);
  fChain->SetBranchAddress("DTextrapolatedOnRPCStation", &DTextrapolatedOnRPCStation, &b_DTextrapolatedOnRPCStation);
  fChain->SetBranchAddress("DTextrapolatedOnRPCLayer", &DTextrapolatedOnRPCLayer, &b_DTextrapolatedOnRPCLayer);
  fChain->SetBranchAddress("DTextrapolatedOnRPCRoll", &DTextrapolatedOnRPCRoll, &b_DTextrapolatedOnRPCRoll);
  fChain->SetBranchAddress("DTextrapolatedOnRPCRing", &DTextrapolatedOnRPCRing, &b_DTextrapolatedOnRPCRing);
  fChain->SetBranchAddress("DTextrapolatedOnRPCStripw", &DTextrapolatedOnRPCStripw, &b_DTextrapolatedOnRPCStripw);

  fChain->SetBranchAddress("Ndtsegments", &Ndtsegments, &b_Ndtsegments);
  fChain->SetBranchAddress("NdtltTwinMuxOut", &NdtltTwinMuxOut, &b_NdtltTwinMuxOut);
  fChain->SetBranchAddress("NdtltTwinMuxIn", &NdtltTwinMuxIn, &b_NdtltTwinMuxIn);
  fChain->SetBranchAddress("Nmuons", &Nmuons, &b_Nmuons);
  fChain->SetBranchAddress("NhltFilters", &NhltFilters, &b_NhltFilters);
  fChain->SetBranchAddress("NhltPaths",   &NhltPaths,   &b_NhltPaths);

  fChain->SetBranchAddress("NirpcrechitsTwinMux", &NirpcrechitsTwinMux, &b_NirpcrechitsTwinMux);
  fChain->SetBranchAddress("NDTsegmentonRPC", &NDTsegmentonRPC, &b_NDTsegmentonRPC);

   Notify();
}

Bool_t combined_trigger_studies::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void combined_trigger_studies::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

#endif // #ifdef combined_trigger_studies_cxx
