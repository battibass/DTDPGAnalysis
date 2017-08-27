#define combined_trigger_studies_cxx
#include "combined_trigger_studies.h"

#include "TMath.h"
#include "TFile.h"

#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"

#include <iostream>

void combined_trigger_studies::Loop()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //   In a ROOT session, you can do:
  //      root> .L combined_trigger_studies.C
  //      root> combined_trigger_studies t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries


  TFile *outputFile = new TFile(output_file_name,"recreate");
 
  if (!(outputFile->IsOpen())) {
    std::cout<<("Output File cannot be opened\n");
    return;
  }

  outputFile->mkdir("kinematics");
  outputFile->mkdir("efficiencies");
  outputFile->mkdir("distance");
  outputFile->mkdir("trigger");

  // muon kinematics

  outputFile->cd("/kinematics");

  TH1F *h_Zmumu_mass = new TH1F("h_Zmumu_mass", "Z boson candidate mass; m(tag_{#mu},probe_{#mu}) (GeV);# entries", 100, 50., 150.);
  TH1F *h_mu_eta = new TH1F("h_mu_eta", "muon #eta;muon #eta;# entries", 48, -2.4, 2.4);
  TH1F *h_mu_phi = new TH1F("h_mu_phi", "muon #phi;muon #phi (rad);# entries", 48, -pig,  pig);
  TH1F *h_mu_pt  = new TH1F("h_mu_pt" , "muon p_{T};muon p_{T} (GeV/c);# entries", 150, 0., 150.);
  
  // distance

  outputFile->cd("/distance");

  TH1F *h_dR_seg_muon_MB1 = new TH1F("h_dR_seg_muon_MB1", 
				     "Distance in MB1 between segment and muon;#DeltaR (cm);# entires", 
				     100, 0., 50.);

  TH1F *h_dR_seg_muon_MB2 = new TH1F("h_dR_seg_muon_MB2", 
				     "Distance in MB1 between segment and muon;#DeltaR (cm);# entires", 
				     100, 0., 50.);

  TH2F *h_dX_dY_seg_muon_MB1 = new TH2F("h_dX_dY_seg_muon_MB1", 
					"Distance in MB1 between segment and muon;#Deltax (cm);#Deltay (cm)", 
					100, -50., 50., 100, -50., 50.);

  TH2F *h_dX_dY_seg_muon_MB2 = new TH2F("h_dX_dY_seg_muon_MB2", 
					"Distance DT MB2 between segment and muon;#Deltax (cm);#Deltay (cm)", 
					100, -50., 50., 100, -50., 50.);

  TH1F *h_dPhi_seg_TrigIn_MB1 = new TH1F("h_dPhi_seg_TrigIn_MB1", 
					 "Distance in MB1 between segment and TwinMux In;#Delta#phi (rad);# entires", 
					 100, 0., .5);

  TH1F *h_dPhi_seg_TrigIn_MB2 = new TH1F("h_dPhi_seg_TrigIn_MB2", 
					 "Distance in MB2 between segment and TwinMux In;#Delta#phi (rad);# entires", 
					 100, 0., .5);  

  TH1F *h_dPhi_seg_TrigOut_MB1 = new TH1F("h_dPhi_seg_TrigOut_MB1", 
					 "Distance in MB1 between segment and TwinMux Out;#Delta#phi (rad);# entires", 
					 100, 0., .5);

  TH1F *h_dPhi_seg_TrigOut_MB2 = new TH1F("h_dPhi_seg_TrigOut_MB2", 
					 "Distance in MB2 between segment and TwinMux Out;#Delta#phi (rad);# entires", 
					 100, 0., .5);  

  TH1F *h_dX_MB1_layer_1 = new TH1F("h_dX_MB1_layer_1", "Distance on MB1 layer 1;#Deltax (cm); # entries", 100, 0, 10);
  TH1F *h_dX_MB1_layer_2 = new TH1F("h_dX_MB1_layer_2", "Distance on MB1 layer 2;#Deltax (cm); # entries", 100, 0, 10);
  
  TH1F *h_dX_MB2_layer_1 = new TH1F("h_dX_MB2_layer_1", "Distance on MB2 layer 1;#Deltax (cm); # entries", 100, 0, 10);
  TH1F *h_dX_MB2_layer_2 = new TH1F("h_dX_MB2_layer_2", "Distance on MB2 layer 2;#Deltax (cm); # entires", 100, 0, 10);

  

  
  // efficiency

  outputFile->cd("/efficiencies");

  TEfficiency* h_eff_twinmux_in_pt_MB1 = new TEfficiency("h_eff_twinmux_in_pt_MB1", 
							 "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 50,0.,100.);
  TEfficiency* h_eff_twinmux_in_pt_MB2 = new TEfficiency("h_eff_twinmux_in_pt_MB2", 
							 "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 50,0.,100.);    
  
  TEfficiency* h_eff_twinmux_in_eta_MB1 = new TEfficiency("h_eff_twinmux_in_eta_MB1",
							  "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  TEfficiency* h_eff_twinmux_in_eta_MB2 = new TEfficiency("h_eff_twinmux_in_eta_MB2",
							  "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  TEfficiency* h_eff_twinmux_in_phi_MB1 = new TEfficiency("h_eff_twinmux_in_phi_MB1",
							  "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  TEfficiency* h_eff_twinmux_in_phi_MB2 = new TEfficiency("h_eff_twinmux_in_phi_MB2",
							  "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  
  TEfficiency* h_eff_rpc_pt_MB1 = new TEfficiency("h_eff_rpc_pt_MB1", 
						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 50,0.,100.);
  TEfficiency* h_eff_rpc_pt_MB2 = new TEfficiency("h_eff_rpc_pt_MB2", 
						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 50,0.,100.);    
  
  TEfficiency* h_eff_rpc_eta_MB1 = new TEfficiency("h_eff_rpc_eta_MB1",
						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  TEfficiency* h_eff_rpc_eta_MB2 = new TEfficiency("h_eff_rpc_eta_MB2",
						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  TEfficiency* h_eff_rpc_phi_MB1 = new TEfficiency("h_eff_rpc_phi_MB1",
						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  TEfficiency* h_eff_rpc_phi_MB2 = new TEfficiency("h_eff_rpc_phi_MB2",
						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);

  TEfficiency* h_eff_combined_toy_pt_MB1 = new TEfficiency("h_eff_combined_toy_pt_MB1", 
						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 50,0.,100.);
  TEfficiency* h_eff_combined_toy_pt_MB2 = new TEfficiency("h_eff_combined_toy_pt_MB2", 
						  "Efficiency vs p_{T};p_{T} (GeV/c);#epsilon", 50,0.,100.);    
  
  TEfficiency* h_eff_combined_toy_eta_MB1 = new TEfficiency("h_eff_combined_toy_eta_MB1",
						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);
  TEfficiency* h_eff_combined_toy_eta_MB2 = new TEfficiency("h_eff_combined_toy_eta_MB2",
						   "Efficiency vs #eta;#eta;#epsilon", 24,-1.2,1.2);

  TEfficiency* h_eff_combined_toy_phi_MB1 = new TEfficiency("h_eff_combined_toy_phi_MB1",
						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);
  TEfficiency* h_eff_combined_toy_phi_MB2 = new TEfficiency("h_eff_combined_toy_phi_MB2",
						   "Efficiency vs #phi;#phi (rad);#epsilon", 48,-pig,pig);

  // Twin MuX quality and BX

  outputFile->cd("/trigger");

  TH1F *h_BX_twinmux_in_MB1 = new TH1F("h_BX_twinmux_in_MB1", "TwinMux In BX for MB1; BX; # entries", 7, -3.5, 3.5);
  TH1F *h_BX_twinmux_in_MB2 = new TH1F("h_BX_twinmux_in_MB2", "TwinMux In BX for MB2; BX; # entries", 7, -3.5, 3.5);

  TH1F *h_BX_twinmux_out_MB1 = new TH1F("h_BX_twinmux_out_MB1", "TwinMux Out BX for MB1; BX; # entries", 7, -3.5, 3.5);
  TH1F *h_BX_twinmux_out_MB2 = new TH1F("h_BX_twinmux_out_MB2", "TwinMux Out BX for MB2; BX; # entries", 7, -3.5, 3.5);

  TH1F *h_qual_twinmux_in_MB1 = new TH1F("h_qual_twinmux_in_MB1", "TwinMux In quality for MB1; quality; # entries", 7, -0.5, 6.5);
  TH1F *h_qual_twinmux_in_MB2 = new TH1F("h_qual_twinmux_in_MB2", "TwinMux In quality for MB2; quality; # entries", 7, -0.5, 6.5);
 


  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast() < n_events ?
                      fChain->GetEntriesFast() : n_events ; 

  Long64_t nbytes = 0, nb = 0;
   
  Int_t count_muon = 0;
  
  std::cout << "[combined_trigger_studies::Loop] processing : " 
	    << nentries << " entries " << std::endl;
 
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if(jentry % 10000 == 0) 
      std::cout << "[combined_trigger_studies::Loop] processed : " << jentry << " entries\r" << std::flush;
    
    auto tnpPairs = combined_trigger_studies::TnPSelection(min_TnP_mass,max_TnP_mass);

    for(const auto & pair : tnpPairs) {

      count_muon++;

      // the first  element of the pair is the index of the tag muon
      // the second element of the pair is the index of the probe muon
      // we want to study the probe

      Int_t iTag   = pair.first;
      Int_t iProbe = pair.second;
      
      TLorentzVector tag_vec;
      tag_vec.SetXYZM(Mu_px->at(iTag),
		      Mu_py->at(iTag),
		      Mu_pz->at(iTag),
		      0.106);
      
      TLorentzVector probe_vec;
      probe_vec.SetXYZM(Mu_px->at(iProbe),
			Mu_py->at(iProbe),
			Mu_pz->at(iProbe),
			0.106);
      
      h_mu_pt->Fill(probe_vec.Pt());
      h_mu_phi->Fill(probe_vec.Phi());
      h_mu_eta->Fill(probe_vec.Eta());
      
      Float_t mass = (probe_vec + tag_vec).M();
      h_Zmumu_mass->Fill(mass);
      
      
      // ************************************
      // We want to look for segmets closest 
      // to the muon in MB1 and MB2
      // ************************************

      // ************************************
      // 1st excercise step 1:
      // select the segment closest to the
      // propagate muon track and tune the
      // between segments and propagated RECO
      // muons to a reasonable value
      // ************************************

      bool has_match_DT_MB1_muon = false;
      bool has_match_DT_MB2_muon = false;
      Float_t minDrMB1 = 999.;
      Float_t minDrMB2 = 999.;
      
      Int_t dtsegment_index[2] = {-999, -999};

      /// Loop on DT 4D segments ///
      for(Int_t iDtSegm = 0; iDtSegm < Ndtsegments; ++iDtSegm) {
	
	if(dtsegm4D_station->at(iDtSegm) > 2)   continue; // MB3 and MB4 do not matter
	
	if(dtsegm4D_hasZed->at(iDtSegm)   <= 0) continue; // Miminal quality cuts:
	if(dtsegm4D_phinhits->at(iDtSegm) <  4) continue; // ask for z view and >= 4 hits in phi  			

	///// station 1 /////
	
	if(Mu_sector_MB1->at(iProbe)      > 0 && // non valid matches with MB1 are -999
	   dtsegm4D_station->at(iDtSegm) == 1 &&
	   dtsegm4D_sector->at(iDtSegm) == Mu_sector_MB1->at(iProbe) &&
	   dtsegm4D_wheel->at(iDtSegm)  == Mu_wheel_MB1->at(iProbe)  ) {

	  Float_t dXMB1 = (Mu_x_MB1->at(iProbe) - dtsegm4D_x_pos_loc->at(iDtSegm));
	  Float_t dYMB1 = (Mu_y_MB1->at(iProbe) - dtsegm4D_y_pos_loc->at(iDtSegm));
	  
	  Float_t dRMB1 = sqrt(dXMB1*dXMB1 + dYMB1*dYMB1);
	  
	  h_dX_dY_seg_muon_MB1->Fill(dXMB1,dYMB1);
	  
	  if(dRMB1 < minDrMB1) {
	    dtsegment_index[0] = iDtSegm;
	    minDrMB1 = dRMB1; 
	  }
	  
	}
	
	///// station 2 /////
	
 	if(Mu_sector_MB2->at(iProbe)      > 0  && // non valid matches with MB2 are -999
	   dtsegm4D_station->at(iDtSegm) == 2 &&
	   dtsegm4D_sector->at(iDtSegm) == Mu_sector_MB2->at(iProbe) &&
	   dtsegm4D_wheel->at(iDtSegm)  == Mu_wheel_MB2->at(iProbe)  ) {
	  
	  Float_t dXMB2 = (Mu_x_MB2->at(iProbe) - dtsegm4D_x_pos_loc->at(iDtSegm));
	  Float_t dYMB2 = (Mu_y_MB2->at(iProbe) - dtsegm4D_y_pos_loc->at(iDtSegm));
	    
	  Float_t dRMB2 = sqrt(dXMB2*dXMB2 + dYMB2*dYMB2);
	  
	  h_dX_dY_seg_muon_MB2->Fill(dXMB2,dYMB2);
	  
	  if(dRMB2 < minDrMB2) {  
	    dtsegment_index[1] = iDtSegm;
	    minDrMB2 = dRMB2; 
	  }
	  
	}
	
      }

      if(minDrMB1 < 900.) 
	h_dR_seg_muon_MB1->Fill(minDrMB1);

      if(minDrMB2 < 900.)
	h_dR_seg_muon_MB2->Fill(minDrMB2);

      if(minDrMB1 < 10.) // CB actually tune it!
	has_match_DT_MB1_muon = true;

      if(minDrMB2 < 10.) // CB actually tune it!
	has_match_DT_MB2_muon = true;

      // ************************************
      // We want to look for _all_ DT trigger
      // primitives (TwinMux input) close
      // to the segments found above
      // ************************************

      // ************************************
      // 1st excercise step 2:
      // tune the matching between segments
      // and trigger (can do the same for
      // TwinMux Out later on!)
      // ************************************
     
      std::vector<Int_t> twinmux_in_MB1;
      std::vector<Int_t> twinmux_in_MB2;

      for(Int_t iTrig = 0; iTrig < NdtltTwinMuxIn; ++iTrig) {
	
	if(has_match_DT_MB1_muon) {
	  
	  Int_t trigWh  = ltTwinMuxIn_wheel->at(iTrig);
	  Int_t trigSt  = ltTwinMuxIn_station->at(iTrig);
	  Int_t trigSec = ltTwinMuxIn_sector->at(iTrig);
	  
	  if(trigWh  == dtsegm4D_wheel->at(dtsegment_index[0])   &&
	     trigSt  == dtsegm4D_station->at(dtsegment_index[0]) &&
	     trigSec == dtsegm4D_sector->at(dtsegment_index[0]) ) {
	    
	    Float_t trigPhiGlb = PhiConversion(ltTwinMuxIn_phi->at(iTrig), trigSec);
	    
	    Float_t dPhi = acos(cos(dtsegm4D_phi->at(dtsegment_index[0]) - trigPhiGlb));
	    
	    h_dPhi_seg_TrigIn_MB1->Fill(dPhi);
	    
	    if(std::abs(dPhi) < 0.1) //CB tune it
	      twinmux_in_MB1.push_back(iTrig);
	  }
	}
	  
	if(has_match_DT_MB2_muon) {
	    
	  Int_t trigWh  = ltTwinMuxIn_wheel->at(iTrig);
	  Int_t trigSt  = ltTwinMuxIn_station->at(iTrig);
	  Int_t trigSec = ltTwinMuxIn_sector->at(iTrig);
	  
	  if(trigWh  == dtsegm4D_wheel->at(dtsegment_index[1])   &&
	     trigSt  == dtsegm4D_station->at(dtsegment_index[1]) &&
	     trigSec == dtsegm4D_sector->at(dtsegment_index[1]) ) {
	    
	    Float_t trigPhiGlb = PhiConversion(ltTwinMuxIn_phi->at(iTrig), trigSec);
	    
	    Float_t dPhi = acos(cos(dtsegm4D_phi->at(dtsegment_index[1]) - trigPhiGlb));
	    
	    h_dPhi_seg_TrigIn_MB2->Fill(dPhi);
	    
	    if(std::abs(dPhi) < 0.1) //CB tune it
	      twinmux_in_MB2.push_back(iTrig);
	  }
	    
	}

      } // loop on TwinMux In


      // ************************************
      // We want to look at the quality, BX
      // assignment and efficiency of trigger
      // primitives from "good" muons
      // ************************************

      // ************************************
      // 1st excercise step 3:
      // plot BX and quality for all the
      // matched trigger primitives and 
      // compute efficiency to trigger at 
      // BX == 0
      // ************************************

      bool has_twinmux_in_BX0_MB1 = false;
      bool has_twinmux_in_BX0_MB2 = false;

      for (const auto & iTrig : twinmux_in_MB1)
	{
	  
	  Int_t bx   = ltTwinMuxIn_bx->at(iTrig);
	  Int_t qual = ltTwinMuxIn_quality->at(iTrig);

	  if (bx == 0)
	    has_twinmux_in_BX0_MB1 = true;

	  h_BX_twinmux_in_MB1->Fill(bx);
	  h_qual_twinmux_in_MB1->Fill(qual);
	}	 

      for (const auto & iTrig : twinmux_in_MB2)
	{
	  
	  Int_t bx   = ltTwinMuxIn_bx->at(iTrig);
	  Int_t qual = ltTwinMuxIn_quality->at(iTrig);
	  
	  if (bx == 0)
	    has_twinmux_in_BX0_MB2 = true;

	  h_BX_twinmux_in_MB2->Fill(bx);
	  h_qual_twinmux_in_MB2->Fill(qual);
	}	 

      if (has_match_DT_MB1_muon)
	{
	  h_eff_twinmux_in_pt_MB1->Fill(has_twinmux_in_BX0_MB1,probe_vec.Pt());
	  h_eff_twinmux_in_eta_MB1->Fill(has_twinmux_in_BX0_MB1,probe_vec.Eta());
	  h_eff_twinmux_in_phi_MB1->Fill(has_twinmux_in_BX0_MB1,probe_vec.Phi());
	}

      if (has_match_DT_MB2_muon)
	{
	  h_eff_twinmux_in_pt_MB2->Fill(has_twinmux_in_BX0_MB2,probe_vec.Pt());
	  h_eff_twinmux_in_eta_MB2->Fill(has_twinmux_in_BX0_MB2,probe_vec.Eta());
	  h_eff_twinmux_in_phi_MB2->Fill(has_twinmux_in_BX0_MB2,probe_vec.Phi());
	}

      // ************************************
      // We want to study RPC trigger 
      // quantities
      // ************************************

      // ************************************
      // 2nd excercise step 1:
      // extrapolate DT segment on inner and
      // outer layer of the RPC and match
      // with the closest RPC recHit
      // ************************************

      bool has_dt_extrapolation_InOut_MB1 = false;
      bool has_dt_extrapolation_InOut_MB2 = false;

      bool has_rpc_match_MB1 = false;
      bool has_rpc_match_MB2 = false;

      // indexes in arrays are inner and outer layer
      bool has_dt_extrapolation_MB1[2] = {false, false};
      
      Float_t minDxMB1[2] = { 999., 999.};
	  
      Float_t rpcCluSizeMB1[2] = {0., 0.}; 
      Float_t stripWidthMB1[2] = {0., 0.};	    

      bool has_dt_extrapolation_MB2[2] = {false, false};
      
      Float_t minDxMB2[2] = { 999., 999.};
	  
      Float_t rpcCluSizeMB2[2] = {0., 0.}; 
      Float_t stripWidthMB2[2] = {0., 0.};	    

      // MB1

      if( has_match_DT_MB1_muon ) {
	for(Int_t iOnRpc = 0; iOnRpc < NDTsegmentonRPC->at(dtsegment_index[0]); ++iOnRpc) {	  

	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[0],iOnRpc) == 1) { // inner ring
	    
	    has_dt_extrapolation_MB1[0] = true;
	    
	    stripWidthMB1[0] = getXY<float>(DTextrapolatedOnRPCStripw,dtsegment_index[0],iOnRpc);
	    
	    for(int iRpc = 0; iRpc < NirpcrechitsTwinMux; ++iRpc) {
	      
	      if((RpcRecHitTwinMuxStation->at(iRpc) == getXY<int>(DTextrapolatedOnRPCStation,dtsegment_index[0],iOnRpc)) &&
		 (RpcRecHitTwinMuxSector->at(iRpc)  == getXY<int>(DTextrapolatedOnRPCSector,dtsegment_index[0],iOnRpc))  &&
		 (RpcRecHitTwinMuxRing->at(iRpc)    == getXY<int>(DTextrapolatedOnRPCRing,dtsegment_index[0],iOnRpc))    &&
		 (RpcRecHitTwinMuxLayer->at(iRpc)   == getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[0],iOnRpc)) ) {
		
		Float_t dXInner  = abs(getXY<float>(DTextrapolatedOnRPCLocX,dtsegment_index[0],iOnRpc) - RpcRechitTwinMuxLocX->at(iRpc));
		
		if(dXInner < minDxMB1[0]) {
		  minDxMB1[0] = dXInner;
		  rpcCluSizeMB1[0] = RpcRecHitTwinMuxClusterSize->at(iRpc);
		}
	      }
	    }
	    
	    h_dX_MB1_layer_1->Fill(minDxMB1[0]);
	    
	  }
	  
	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[0],iOnRpc) == 2) { // outer ring
	    
	    has_dt_extrapolation_MB1[1] = true;
	    
	    stripWidthMB1[1] = getXY<float>(DTextrapolatedOnRPCStripw,dtsegment_index[0],iOnRpc);
	    
	    for(int iRpc = 0; iRpc < NirpcrechitsTwinMux; ++iRpc) {
	      
	      if((RpcRecHitTwinMuxStation->at(iRpc) == getXY<int>(DTextrapolatedOnRPCStation,dtsegment_index[0],iOnRpc)) &&
	   	 (RpcRecHitTwinMuxSector->at(iRpc)  == getXY<int>(DTextrapolatedOnRPCSector,dtsegment_index[0],iOnRpc))  &&
	   	 (RpcRecHitTwinMuxRing->at(iRpc)    == getXY<int>(DTextrapolatedOnRPCRing,dtsegment_index[0],iOnRpc))    &&
	   	 (RpcRecHitTwinMuxLayer->at(iRpc)   == getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[0],iOnRpc)) ) {
		
	   	Float_t dXInner  = abs(getXY<float>(DTextrapolatedOnRPCLocX,dtsegment_index[0],iOnRpc) - RpcRechitTwinMuxLocX->at(iRpc));
		
	   	if(dXInner < minDxMB1[1]) {
	   	  minDxMB1[1] = dXInner;
	   	  rpcCluSizeMB1[1] = RpcRecHitTwinMuxClusterSize->at(iRpc);
	   	}
	      }
	    }
	    
	    h_dX_MB1_layer_2->Fill(minDxMB1[1]);
	    
	  }
	  
	}
      }
    	
      // MB2
      
      if( has_match_DT_MB2_muon ) {
	
	for(Int_t iOnRpc = 0; iOnRpc < NDTsegmentonRPC->at(dtsegment_index[1]); ++iOnRpc) {
	  
	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[1],iOnRpc) == 1) { // inner ring
	    
	    has_dt_extrapolation_MB2[0] = true;
	    
	    stripWidthMB2[0] = getXY<float>(DTextrapolatedOnRPCStripw,dtsegment_index[1],iOnRpc);
	    
	    for(int iRpc = 0; iRpc < NirpcrechitsTwinMux; ++iRpc) {
	      
	      if((RpcRecHitTwinMuxStation->at(iRpc) == getXY<int>(DTextrapolatedOnRPCStation,dtsegment_index[1],iOnRpc)) &&
		 (RpcRecHitTwinMuxSector->at(iRpc)  == getXY<int>(DTextrapolatedOnRPCSector,dtsegment_index[1],iOnRpc))  &&
		 (RpcRecHitTwinMuxRing->at(iRpc)    == getXY<int>(DTextrapolatedOnRPCRing,dtsegment_index[1],iOnRpc))    &&
		 (RpcRecHitTwinMuxLayer->at(iRpc)   == getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[1],iOnRpc)) ) {
		
		Float_t dXInner  = abs(getXY<float>(DTextrapolatedOnRPCLocX,dtsegment_index[1],iOnRpc) - RpcRechitTwinMuxLocX->at(iRpc));
		
		if(dXInner < minDxMB2[0]) {
		  minDxMB2[0] = dXInner;
		  rpcCluSizeMB2[0] = RpcRecHitTwinMuxClusterSize->at(iRpc);
		}
	      }
	    }
	    
	    h_dX_MB2_layer_1->Fill(minDxMB2[0]);
	  
	  }
	
	  if(getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[1],iOnRpc) == 2) { // outer ring
	    
	    has_dt_extrapolation_MB2[1] = true;
	    
	    stripWidthMB2[1] = getXY<float>(DTextrapolatedOnRPCStripw,dtsegment_index[1],iOnRpc);
	    
	    for(int iRpc = 0; iRpc < NirpcrechitsTwinMux; ++iRpc) {
	      
	      if((RpcRecHitTwinMuxStation->at(iRpc) == getXY<int>(DTextrapolatedOnRPCStation,dtsegment_index[1],iOnRpc)) &&
		 (RpcRecHitTwinMuxSector->at(iRpc)  == getXY<int>(DTextrapolatedOnRPCSector,dtsegment_index[1],iOnRpc))  &&
		 (RpcRecHitTwinMuxRing->at(iRpc)    == getXY<int>(DTextrapolatedOnRPCRing,dtsegment_index[1],iOnRpc))    &&
		 (RpcRecHitTwinMuxLayer->at(iRpc)   == getXY<int>(DTextrapolatedOnRPCLayer,dtsegment_index[1],iOnRpc)) ) {
		
		Float_t dXInner  = abs(getXY<float>(DTextrapolatedOnRPCLocX,dtsegment_index[1],iOnRpc) - RpcRechitTwinMuxLocX->at(iRpc));
		
		if(dXInner < minDxMB2[1]) {
		  minDxMB2[1] = dXInner;
		  rpcCluSizeMB2[1] = RpcRecHitTwinMuxClusterSize->at(iRpc);
		}
	      }
	    }
	    
	    h_dX_MB2_layer_2->Fill(minDxMB2[1]);
	    
	  }	  
	}
      }
    
      // ************************************
      // 2nd excercise step 2:
      // Tune matching between extrapolated
      // recHit and RPC recHit: our goal is 
      // to evaluate muon efficiency
      // ************************************      
      
      if (has_dt_extrapolation_MB1[0] && has_dt_extrapolation_MB1[1] ) {
	
	has_dt_extrapolation_InOut_MB1 = true;
	
	if( minDxMB1[0] < (range_strips + rpcCluSizeMB1[0]*0.5) * stripWidthMB1[0] && rpcCluSizeMB1[0] <= cluster_size_cut &&
	    minDxMB1[1] < (range_strips + rpcCluSizeMB1[1]*0.5) * stripWidthMB1[1] && rpcCluSizeMB1[1] <= cluster_size_cut ) {
	  has_rpc_match_MB1 = true;
	}
      }
      
      if (has_dt_extrapolation_MB2[0] && has_dt_extrapolation_MB2[1] ) {
	
	has_dt_extrapolation_InOut_MB2 = true;
	
	if( minDxMB2[0] < (range_strips + rpcCluSizeMB2[0]*0.5) * stripWidthMB2[0] && rpcCluSizeMB2[0] <= cluster_size_cut &&
	    minDxMB2[1] < (range_strips + rpcCluSizeMB2[1]*0.5) * stripWidthMB2[1] && rpcCluSizeMB2[1] <= cluster_size_cut ) {
	  has_rpc_match_MB2 = true;
	}
      }
      
      if (has_match_DT_MB1_muon && has_dt_extrapolation_InOut_MB1)
	{
	  h_eff_rpc_pt_MB1->Fill(has_rpc_match_MB1,probe_vec.Pt());
	  h_eff_rpc_eta_MB1->Fill(has_rpc_match_MB1,probe_vec.Eta());
	  h_eff_rpc_phi_MB1->Fill(has_rpc_match_MB1,probe_vec.Phi());
	}
      
      if (has_match_DT_MB2_muon && has_dt_extrapolation_InOut_MB2)
	{
	  h_eff_rpc_pt_MB2->Fill(has_rpc_match_MB2,probe_vec.Pt());
	  h_eff_rpc_eta_MB2->Fill(has_rpc_match_MB2,probe_vec.Eta());
	  h_eff_rpc_phi_MB2->Fill(has_rpc_match_MB2,probe_vec.Phi());
	}

      if (has_match_DT_MB1_muon)
	{
	  h_eff_combined_toy_pt_MB1->Fill(has_twinmux_in_BX0_MB1 || has_rpc_match_MB1,probe_vec.Pt());
	  h_eff_combined_toy_eta_MB1->Fill(has_twinmux_in_BX0_MB1 || has_rpc_match_MB1,probe_vec.Eta());
	  h_eff_combined_toy_phi_MB1->Fill(has_twinmux_in_BX0_MB1 || has_rpc_match_MB1,probe_vec.Phi());
	}

      if (has_match_DT_MB2_muon)
	{
	  h_eff_combined_toy_pt_MB2->Fill(has_twinmux_in_BX0_MB2 || has_rpc_match_MB2,probe_vec.Pt());
	  h_eff_combined_toy_eta_MB2->Fill(has_twinmux_in_BX0_MB2 || has_rpc_match_MB2,probe_vec.Eta());
	  h_eff_combined_toy_phi_MB2->Fill(has_twinmux_in_BX0_MB2 || has_rpc_match_MB2,probe_vec.Phi());
	}

      
      // ************************************
      // 2nd excercise step 3:
      // Evaluate how RPC improve primitive
      // BX looking at the TWinMux OUT
      // distribution 
      // ************************************	
	
      std::vector<Int_t> twinmux_out_MB1;
      std::vector<Int_t> twinmux_out_MB2;

      for(Int_t iTrig = 0; iTrig < NdtltTwinMuxOut; ++iTrig) {
	
	if(has_match_DT_MB1_muon) {
	  
	  Int_t trigWh  = ltTwinMuxOut_wheel->at(iTrig);
	  Int_t trigSt  = ltTwinMuxOut_station->at(iTrig);
	  Int_t trigSec = ltTwinMuxOut_sector->at(iTrig);
	  
	  if(trigWh  == dtsegm4D_wheel->at(dtsegment_index[0])   &&
	     trigSt  == dtsegm4D_station->at(dtsegment_index[0]) &&
	     trigSec == dtsegm4D_sector->at(dtsegment_index[0]) ) {
	    
	    Float_t trigPhiGlb = PhiConversion(ltTwinMuxOut_phi->at(iTrig), trigSec);
	    
	    Float_t dPhi = acos(cos(dtsegm4D_phi->at(dtsegment_index[0]) - trigPhiGlb));
	    
	    h_dPhi_seg_TrigOut_MB1->Fill(dPhi);
	    
	    if(std::abs(dPhi) < 0.1) //CB tune it
	      twinmux_out_MB1.push_back(iTrig);
	  }
	}
	  
	if(has_match_DT_MB2_muon) {
	    
	  Int_t trigWh  = ltTwinMuxOut_wheel->at(iTrig);
	  Int_t trigSt  = ltTwinMuxOut_station->at(iTrig);
	  Int_t trigSec = ltTwinMuxOut_sector->at(iTrig);
	  
	  if(trigWh  == dtsegm4D_wheel->at(dtsegment_index[1])   &&
	     trigSt  == dtsegm4D_station->at(dtsegment_index[1]) &&
	     trigSec == dtsegm4D_sector->at(dtsegment_index[1]) ) {
	    
	    Float_t trigPhiGlb = PhiConversion(ltTwinMuxOut_phi->at(iTrig), trigSec);
	    
	    Float_t dPhi = acos(cos(dtsegm4D_phi->at(dtsegment_index[1]) - trigPhiGlb));
	    
	    h_dPhi_seg_TrigOut_MB2->Fill(dPhi);
	    
	    if(std::abs(dPhi) < 0.1) //CB tune it
	      twinmux_out_MB2.push_back(iTrig);
	  }
	    
	}

      } // loop on TwinMux Out

      for (const auto & iTrig : twinmux_out_MB1)	  
	h_BX_twinmux_out_MB1->Fill(ltTwinMuxOut_bx->at(iTrig));

      for (const auto & iTrig : twinmux_out_MB2)
	h_BX_twinmux_out_MB2->Fill(ltTwinMuxOut_bx->at(iTrig));

    } //loop on muons			
    
  } // loop on event
  
  std::cout << std::endl;
  std::cout << "[combined_trigger_studies::Loop] analyzed " << count_muon << " muons" << std::endl;  

  outputFile->Write();
  outputFile->Close();

}


vector<std::pair<Int_t,Int_t>> combined_trigger_studies::TnPSelection(Float_t minMass,
							 Float_t maxMass)
{

  vector<std::pair<Int_t,Int_t>> pairs;

  for(Int_t iTag = 0; iTag < Nmuons; ++iTag) 
    {
    
      TLorentzVector tagVec;
      tagVec.SetXYZM(Mu_px->at(iTag),
		     Mu_py->at(iTag),
		     Mu_pz->at(iTag),
		     0.106);
      
      bool tagQuality = 
	Mu_isMuGlobal->at(iTag)  == 1 &&
	Mu_isMuTracker->at(iTag) == 1 &&
	std::abs(Mu_dxy_glb->at(iTag)) < Dxy_cut &&
	std::abs(Mu_dz_glb->at(iTag))  < Dz_cut  &&
	Mu_normchi2_glb->at(iTag)      < muchi2_cut &&
	Mu_numberOfHits_sta->at(iTag)        > N_hits_cut && 
	Mu_numberOfPixelHits_glb->at(iTag)   >= npix_cut   &&
	Mu_numberOfTrackerHits_glb->at(iTag) >= ntkr_cut   &&
 	Mu_tkIsoR03_glb->at(iTag) < tkr_iso_cut            &&
	tagVec.Pt() > min_pt_tag ;
      
      if(tagQuality && HasTrigger(tagVec,muon_hlt_dR)) 
	{
	  
	  for(Int_t iProbe = 0; iProbe < Nmuons; ++iProbe) 
	    {
	      
	      if (iTag == iProbe) 
		continue;

	      TLorentzVector probeVec;
	      probeVec.SetXYZM(Mu_px->at(iProbe),
			       Mu_py->at(iProbe),
			       Mu_pz->at(iProbe),
			       0.106);
	      
	      bool probeQuality =
		Mu_isMuTracker->at(iProbe) == 1 &&
		Mu_isMuGlobal->at(iTag)    == 1 &&
		fabs(Mu_dxy_glb->at(iProbe)) < Dxy_cut &&
		fabs(Mu_dz_glb->at(iProbe))  < Dz_cut  &&
		Mu_numberOfPixelHits_glb->at(iProbe)   >= npix_cut   &&
		Mu_numberOfTrackerHits_glb->at(iProbe) >= ntkr_cut   &&
		Mu_tkIsoR03_glb->at(iProbe) < tkr_iso_cut            &&
		std::abs(probeVec.Eta()) < max_eta_probe &&
		probeVec.Pt() > min_pt_probe;

	      if (probeQuality)
		{
		  Float_t mass = (tagVec + probeVec).M();

		  if (mass > minMass && mass < maxMass)
		    {
		      pairs.push_back(std::make_pair(iTag,iProbe));
		      break; // just one probe per tag
		    }
		} 
	    }
	}
    }
  
  return pairs;

}

bool combined_trigger_studies::HasTrigger(const TLorentzVector & muon,
			     Float_t deltaR) {

  for(Int_t iTrig = 0; iTrig < NhltFilters; ++iTrig)
    {

      if (hlt_filter->at(iTrig).Contains(trigger_filter_name))
	{
	  TLorentzVector trigVec;
	  trigVec.SetPtEtaPhiM(hlt_filter_pt->at(iTrig),
			       hlt_filter_eta->at(iTrig),
			       hlt_filter_phi->at(iTrig),
			       0.106);
	  
	  if (trigVec.DeltaR(muon) < deltaR)
	    return true;
  	}
    }
  
  return false;
  
}

Float_t combined_trigger_studies::PhiConversion(Int_t phi_In, Int_t sector)
{	
  Float_t locphi = phi_In / 4096.0;
  Float_t newphi = locphi + ((sector-1) * (pig/6.));
  if (newphi > pig) newphi -= 2 * pig;
  
  return newphi;
}
