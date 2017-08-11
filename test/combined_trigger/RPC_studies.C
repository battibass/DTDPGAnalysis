#define RPC_studies_cxx
#include "RPC_studies.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"

void RPC_studies::Loop()
{

  TString inputFileName  = "";
  TString outputFileName = "";
  Long64_t nEvents = 999999999;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //   In a ROOT session, you can do:
  //      root> .L RPC_studies.C
  //      root> RPC_studies t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries


  TFile *outputFile = new TFile("../DTNtuple.root","recreate");
 
  if (!(outputFile->IsOpen())) {
    std::cout<<("OutputFile cannot be opened\n");
    return;
  }

  if (fChain == 0) return;
   	
  Long64_t nentries = fChain->GetEntriesFast() < nEvents ?
                      fChain->GetEntriesFast() : nEvents ; 
  

  Long64_t nbytes = 0, nb = 0;
   
  std::cout << "[RPC_studies::Loop] processing : " 
	    << nentries << " entries " << std::endl;
 
  int i = 0;
   
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    if(jentry % 10000 == 0) 
      std::cout << "[RPC_studies::Loop] processed : " << jentry << " entries\r";
      
    muon_index = RPC_studies::Selection();
      
    for(int mu = 0; mu < muon_index.index.size(); mu++){
      
      MuPt = muon_index.pt.at(mu);
      MuEta = muon_index.eta.at(mu);
      MuPhi = muon_index.phi.at(mu);
      		
      probe_muon.SetPtEtaPhiM(MuPt, MuEta, MuPhi, 0.105);

// 			h_Zmumu_Mass->Fill(Mass(tag_muon, probe_muon));

      
//       }
//       for(int mu = 0; mu < Nmuons; mu++){

		//////// muon selection
//       	Mu_pt = sqrt(Mu_px->at(mu)*Mu_px->at(mu) + Mu_py->at(mu)*Mu_py->at(mu));
//       	
// 		if(Mu_isMuGlobal->at(mu) !=1 || Mu_isMuTracker->at(mu) !=1 ||
// 				fabs(Mu_dxy_glb->at(mu)) > D0_cut || fabs(Mu_dz_glb->at(mu)) > Dz_cut ||
// 				fabs(Mu_eta->at(mu)) > eta_cut || (Mu_normchi2_glb->at(mu)) > muchi2_cut ||
// 				Mu_pt < pt_cut || Mu_pt > pt_max || 
// 				Mu_numberOfHits_sta->at(mu) < N_hits_cut || 
// 				Mu_numberOfPixelHits_glb->at(mu) < npix_cut ||
// 				Mu_numberOfTrackerHits_glb->at(mu) < ntkr_cut) continue;
// 
// 		if(Mu_pt < pt_localtrig_cut) continue;
		
		//////// muon selection

      
      count_muon++;		
      
      min_DT_MB1_muon = 20;
      min_DT_MB2_muon = 20;
      match_DT_MB1_muon = false;
      match_DT_MB2_muon = false;
      dtsegment_index[0] = -999;
      dtsegment_index[1] = -999;
      
      rpc_extrapolation_mb1 = false;
      rpc_extrapolation_mb2 = false;
      
		      				
      /// for in 4dtSegment///
      for(int dtsegm = 0; dtsegm < Ndtsegments; dtsegm++){
	
	if(dtsegm4D_station->at(dtsegm) > 2) continue; // MB3 and MB4 do not matter
	
	if(dtsegm4D_hasZed->at(dtsegm) <= 0) continue;
	if(dtsegm4D_phinhits->at(dtsegm) <= 0) continue;  				
	
	delta_phi = abs(dtsegm4D_phi->at(dtsegm) - MuPhi);
	if(delta_phi > pig) delta_phi=2*pig-delta_phi;
	
	delta_eta = 0.;
	if(dtsegm4D_hasZed->at(dtsegm) != 0) 
	  delta_eta = abs(dtsegm4D_eta->at(dtsegm)-MuEta);
	
	///// station 1 /////
	
	if(dtsegm4D_station->at(dtsegm) == 1){
	  deltaR_MB1 = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
	  if(deltaR_MB1 < min_DT_MB1_muon){	
	    min_DT_MB1_muon = deltaR_MB1;
	  }
	  else 
	    continue;
	  
	  h_dist_DT_MB1_muon->Fill(min_DT_MB1_muon);
	  
	  if(min_DT_MB1_muon < 0.4){
	    dtsegment_index[0] = dtsegm;
	    match_DT_MB1_muon = true;
	    eta_mu[0] = MuEta;
	    phi_mu[0] = MuPhi;
	    pt_mu[0] = MuPt;
	    
	    min_layer[0] = 100;
	    min_layer[1] = 100;
	    
	    dt_extrapolation_mb1 = false;
	    for(int rpcExtra = 0; rpcExtra < NDTsegmentonRPC; rpcExtra++){
	      if((dtsegm4D_station->at(dtsegm) != DTextrapolatedOnRPCStation->at(rpcExtra)) || 
		 (dtsegm4D_sector->at(dtsegm) != DTextrapolatedOnRPCSector->at(rpcExtra)) || 
		 (dtsegm4D_wheel->at(dtsegm) != DTextrapolatedOnRPCRing->at(rpcExtra))) continue;
	      
	      if(DTextrapolatedOnRPCRegion->at(rpcExtra) != 0) std::cout<<" DT EXTRAPOLATION OUT OF BARREL !!! PLEASE CHECK "<<std::endl;
	      
	      if(DTextrapolatedOnRPCStation->at(rpcExtra) == 3){
		// 								std::cout<<"MB3: we don't care about it"<<std::endl;
		continue;
	      }
	      
	      dt_ring = DTextrapolatedOnRPCRing->at(rpcExtra)+2;
	      dt_station = DTextrapolatedOnRPCStation->at(rpcExtra)-1;
	      dt_sector = DTextrapolatedOnRPCSector->at(rpcExtra)-1;
	      
	      if(DTextrapolatedOnRPCLayer->at(rpcExtra) == 1){
		den[0][dt_ring][dt_station][dt_sector]++;
		stripw[0] = DTextrapolatedOnRPCStripw->at(rpcExtra);
		
		for(int rpc = 0; rpc < NirpcrechitsTwinMux; rpc++) {
		  if(	(RpcRecHitTwinMuxStation->at(rpc) == DTextrapolatedOnRPCStation->at(rpcExtra)) && 
			(RpcRecHitTwinMuxSector->at(rpc) == DTextrapolatedOnRPCSector->at(rpcExtra)) && 
			(RpcRecHitTwinMuxRing->at(rpc) == DTextrapolatedOnRPCRing->at(rpcExtra))) {
		    
		    if(RpcRecHitTwinMuxLayer->at(rpc) != DTextrapolatedOnRPCLayer->at(rpcExtra)) continue;
		    
		    dist_layer[0] = abs(DTextrapolatedOnRPCLocX->at(rpcExtra) - RpcRechitTwinMuxLocX->at(rpc));
		    
		    if(dist_layer[0] < min_layer[0]){
		      min_layer[0] = dist_layer[0];
		      RPC_cluSize[0] = RpcRecHitTwinMuxClusterSize->at(rpc);
		    }
		  }
		}
		h_dist_MB1_layer_1->Fill(min_layer[0]);
		
	      } // for on RPC rechit TwinMux  
	      
	      if(DTextrapolatedOnRPCLayer->at(rpcExtra) == 2){
		den[1][dt_ring][dt_station][dt_sector]++;
		stripw[1] = DTextrapolatedOnRPCStripw->at(rpcExtra);
		
		for(int rpc = 0; rpc < NirpcrechitsTwinMux; rpc++) {
		  if(	(RpcRecHitTwinMuxStation->at(rpc) == DTextrapolatedOnRPCStation->at(rpcExtra)) && 
			(RpcRecHitTwinMuxSector->at(rpc) == DTextrapolatedOnRPCSector->at(rpcExtra)) && 
			(RpcRecHitTwinMuxRing->at(rpc) == DTextrapolatedOnRPCRing->at(rpcExtra))) {
		    
		    if(RpcRecHitTwinMuxLayer->at(rpc) != DTextrapolatedOnRPCLayer->at(rpcExtra)) continue;
		    
		    dist_layer[1] = abs(DTextrapolatedOnRPCLocX->at(rpcExtra) - RpcRechitTwinMuxLocX->at(rpc));
		    
		    if(dist_layer[1] < min_layer[1]){
		      min_layer[1] = dist_layer[1];
		      RPC_cluSize[1] = RpcRecHitTwinMuxClusterSize->at(rpc);
		    }
		  }
		}
		h_dist_MB1_layer_2->Fill(min_layer[1]);
	      }
	      
	      
	      
	    }
	    if(		den[0][dt_ring][dt_station][dt_sector] > 0 && den[1][dt_ring][dt_station][dt_sector] > 0)	{
	      dt_extrapolation_mb1 = true;
	      den[0][dt_ring][dt_station][dt_sector] = 0;
	      den[1][dt_ring][dt_station][dt_sector] = 0;
	    }
	    
	    // 
	  } // if segment < 0.4 to muon
	  
	  if( min_layer[0] < (rangestrips + RPC_cluSize[0]*0.5)*stripw[0] && RPC_cluSize[0]<=clsCut &&
	      min_layer[1] < (rangestrips + RPC_cluSize[1]*0.5)*stripw[1] && RPC_cluSize[1]<=clsCut ){
	    rpc_extrapolation_mb1 = true;
	  }
	  
	}  // if station 1
	
	
	// if on station 2 //
	
	if(dtsegm4D_station->at(dtsegm) == 2){
	  deltaR_MB2 = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
	  h_dist_DT_MB2_muon->Fill(min_DT_MB2_muon);
	  if(deltaR_MB2 < min_DT_MB2_muon){	
	    min_DT_MB2_muon = deltaR_MB2;
	  }
	  else 
	    continue;
	  
	  if(min_DT_MB2_muon < 0.4){
	    dtsegment_index[1] = dtsegm;
	    match_DT_MB2_muon = true;
	    eta_mu[1] = MuEta;
	    phi_mu[1] = MuPhi;
	    pt_mu[1] = MuPt;
	    
	    min_layer[0] = 100;
	    min_layer[1] = 100;
	    
	    dt_extrapolation_mb2 = false;
	    for(int rpcExtra = 0; rpcExtra < NDTsegmentonRPC; rpcExtra++){
	      if((dtsegm4D_station->at(dtsegm) != DTextrapolatedOnRPCStation->at(rpcExtra)) || 
		 (dtsegm4D_sector->at(dtsegm) != DTextrapolatedOnRPCSector->at(rpcExtra)) || 
		 (dtsegm4D_wheel->at(dtsegm) != DTextrapolatedOnRPCRing->at(rpcExtra))) continue;
	      
	      if(DTextrapolatedOnRPCRegion->at(rpcExtra) != 0) std::cout<<" DT EXTRAPOLATION OUT OF BARREL !!! PLEASE CHECK "<<std::endl;
	      
	      if(DTextrapolatedOnRPCStation->at(rpcExtra) == 3){
		// 								std::cout<<"MB3: we don't care about it"<<std::endl;
		continue;
	      }
	      
	      dt_ring = DTextrapolatedOnRPCRing->at(rpcExtra)+2;
	      dt_station = DTextrapolatedOnRPCStation->at(rpcExtra)-1;
	      dt_sector = DTextrapolatedOnRPCSector->at(rpcExtra)-1;
	      
	      if(DTextrapolatedOnRPCLayer->at(rpcExtra) == 1){
		den[0][dt_ring][dt_station][dt_sector]++;
		stripw[0] = DTextrapolatedOnRPCStripw->at(rpcExtra);
		
		for(int rpc = 0; rpc < NirpcrechitsTwinMux; rpc++) {
		  if(	(RpcRecHitTwinMuxStation->at(rpc) == DTextrapolatedOnRPCStation->at(rpcExtra)) && 
			(RpcRecHitTwinMuxSector->at(rpc) == DTextrapolatedOnRPCSector->at(rpcExtra)) && 
			(RpcRecHitTwinMuxRing->at(rpc) == DTextrapolatedOnRPCRing->at(rpcExtra))) {
		    
		    if(RpcRecHitTwinMuxLayer->at(rpc) != DTextrapolatedOnRPCLayer->at(rpcExtra)) continue;
		    
		    dist_layer[0] = abs(DTextrapolatedOnRPCLocX->at(rpcExtra) - RpcRechitTwinMuxLocX->at(rpc));
		    
		    if(dist_layer[0] < min_layer[0]){
		      min_layer[0] = dist_layer[0];
		      RPC_cluSize[0] = RpcRecHitTwinMuxClusterSize->at(rpc);
		    }
		  }
		}
		h_dist_MB2_layer_1->Fill(min_layer[0]);
		
	      } // for on RPC rechit TwinMux  
	      
	      if(DTextrapolatedOnRPCLayer->at(rpcExtra) == 2){
		den[1][dt_ring][dt_station][dt_sector]++;
		stripw[1] = DTextrapolatedOnRPCStripw->at(rpcExtra);
		
		for(int rpc = 0; rpc < NirpcrechitsTwinMux; rpc++) {
		  if(	(RpcRecHitTwinMuxStation->at(rpc) == DTextrapolatedOnRPCStation->at(rpcExtra)) && 
			(RpcRecHitTwinMuxSector->at(rpc) == DTextrapolatedOnRPCSector->at(rpcExtra)) && 
			(RpcRecHitTwinMuxRing->at(rpc) == DTextrapolatedOnRPCRing->at(rpcExtra))) {
		    
		    if(RpcRecHitTwinMuxLayer->at(rpc) != DTextrapolatedOnRPCLayer->at(rpcExtra)) continue;
		    
		    dist_layer[1] = abs(DTextrapolatedOnRPCLocX->at(rpcExtra) - RpcRechitTwinMuxLocX->at(rpc));
		    
		    if(dist_layer[1] < min_layer[1]){
		      min_layer[1] = dist_layer[1];
		      RPC_cluSize[1] = RpcRecHitTwinMuxClusterSize->at(rpc);
		    }
		  }
		}
		h_dist_MB2_layer_2->Fill(min_layer[1]);
	      }
	      
	      
	      
	    }
	    if(den[0][dt_ring][dt_station][dt_sector] > 0 && den[1][dt_ring][dt_station][dt_sector] > 0) {
	      dt_extrapolation_mb2 = true;
	      den[0][dt_ring][dt_station][dt_sector] = 0;
	      den[1][dt_ring][dt_station][dt_sector] = 0;
	    }
	    
	    
	  } // if segment < 0.4 to muon
	  
	  if( min_layer[0] < (rangestrips + RPC_cluSize[0]*0.5)*stripw[0] && RPC_cluSize[0]<=clsCut &&
	      min_layer[1] < (rangestrips + RPC_cluSize[1]*0.5)*stripw[1] && RPC_cluSize[1]<=clsCut ){
	    rpc_extrapolation_mb2 = true;
	  }
	  
	} // station 2						
	
	
	
      } //for on 4dtSegment
      
      
      ////////////////////
      ///// TWIN MUX /////
      for(int in = 0; in < NdtltTwinMuxIn; ++in) {
	
	if(!match_DT_MB1_muon) continue;
	
	int iWh_In  = ltTwinMuxIn_wheel->at(in)+2;
	int iSt_In  = ltTwinMuxIn_station->at(in)-1;
	int iSec_In = ltTwinMuxIn_sector->at(in)-1;
	if(iWh_In != dtsegm4D_wheel->at(dtsegment_index[0]) || 
	   iSt_In != dtsegm4D_station->at(dtsegment_index[0]) ||
	   iSec_In != dtsegm4D_sector->at(dtsegment_index[0])) continue;
	min_dtSegm_In = 999;
	float phi_In = PhiConversion(ltTwinMuxIn_phi->at(in), ltTwinMuxIn_sector->at(in));
	h_DeltaPhi_InOut_MB1->Fill(phi_In);
	
	if(dtsegm4D_phi->at(dtsegment_index[0]) > 0)
	  delta_dtSegm_In = phi_In - dtsegm4D_phi->at(dtsegment_index[0]);
	else
	  delta_dtSegm_In = phi_In - (2*pig + dtsegm4D_phi->at(dtsegment_index[0]));
      	
	if(abs(delta_dtSegm_In) > pig) delta_dtSegm_In = 2 * pig - delta_dtSegm_In;
      	
	if(delta_dtSegm_In < min_dtSegm_In){
	  min_dtSegm_In = delta_dtSegm_In;
	  bx_In  = ltTwinMuxIn_bx->at(in);
	}
      } // for on TwinMux In
      if(abs(min_dtSegm_In) < 1){ // distanza da stabilire
	h_BX_IN_MB1->Fill(bx_In);
	for(int Out = 0; Out < NdtltTwinMuxOut; Out++){
	  int iWh_Out  = ltTwinMuxOut_wheel->at(Out)+2;
	  int iSt_Out  = ltTwinMuxOut_station->at(Out)-1;
	  int iSec_Out = ltTwinMuxOut_sector->at(Out)-1;
	  if(iWh_In != iWh_Out || iSt_In != iSt_Out || iSec_In != iSec_Out) continue;
	  bx_Out      = ltTwinMuxOut_bx->at(Out);
	  h_BX_OUT_MB1->Fill(bx_Out);
	} // for on TwinMux Out
      } // if on distance between dtSegment and TwinMux In
      
      
      
      for(int in = 0; in < NdtltTwinMuxIn; ++in) {
	
	if(!match_DT_MB2_muon) continue;
	
	int iWh_In  = ltTwinMuxIn_wheel->at(in)+2;
	int iSt_In  = ltTwinMuxIn_station->at(in)-1;
	int iSec_In = ltTwinMuxIn_sector->at(in)-1;
	if(iWh_In != dtsegm4D_wheel->at(dtsegment_index[1]) || 
	   iSt_In != dtsegm4D_station->at(dtsegment_index[1]) ||
	   iSec_In != dtsegm4D_sector->at(dtsegment_index[1])) continue;
	min_dtSegm_In = 999;
	float phi_In = PhiConversion(ltTwinMuxIn_phi->at(in), ltTwinMuxIn_sector->at(in));
	h_DeltaPhi_InOut_MB2->Fill(phi_In);
	if(dtsegm4D_phi->at(dtsegment_index[1]) > 0)
	  delta_dtSegm_In = phi_In - dtsegm4D_phi->at(dtsegment_index[1]);
	else
	  delta_dtSegm_In = phi_In - (2*pig + dtsegm4D_phi->at(dtsegment_index[1]));
      	
	if(abs(delta_dtSegm_In)> pig) delta_dtSegm_In = 2 * pig - delta_dtSegm_In;
	
	if(delta_dtSegm_In < -4)
	  std::cout<<phi_In<<"\t"<<dtsegm4D_phi->at(dtsegment_index[1])<<std::endl;
	
	if(delta_dtSegm_In < min_dtSegm_In){
	  min_dtSegm_In = delta_dtSegm_In;
	  bx_In  = ltTwinMuxIn_bx->at(in);
	}
      	
      } // for on TwinMux In
      if(abs(min_dtSegm_In) < 1){ // distanza da stabilire
	h_BX_IN_MB2->Fill(bx_In);
	for(int Out = 0; Out < NdtltTwinMuxOut; Out++){
	  iWh_Out  = ltTwinMuxOut_wheel->at(Out)+2;
	  iSt_Out  = ltTwinMuxOut_station->at(Out)-1;
	  iSec_Out = ltTwinMuxOut_sector->at(Out)-1;
	  if(iWh_In != iWh_Out || iSt_In != iSt_Out || iSec_In != iSec_Out) continue;
	  bx_Out      = ltTwinMuxOut_bx->at(Out);
	  h_BX_OUT_MB2->Fill(bx_Out);
	} // for on TwinMux Out
      } // if on distance between dtSegment and TwinMux In						
      ///// TWIN MUX /////
      ////////////////////
      
      
      
      
      if(dt_extrapolation_mb1)
	{
	  h_eff_eta_station_1->Fill(rpc_extrapolation_mb1, eta_mu[0]);
	  h_eff_phi_station_1->Fill(rpc_extrapolation_mb1, phi_mu[0]);
	  h_eff_pt_station_1->Fill(rpc_extrapolation_mb1,  pt_mu[0]);
	}

      if(dt_extrapolation_mb2)
	{
	  h_eff_eta_station_2->Fill(rpc_extrapolation_mb2, eta_mu[1]);
	  h_eff_phi_station_2->Fill(rpc_extrapolation_mb2, phi_mu[1]);
	  h_eff_pt_station_2->Fill(rpc_extrapolation_mb2,  pt_mu[1]);
	}
      
    } 
    //loop on muons			
    
  } // loop on event
  
  std::cout << std::endl;
  std::cout << "Muon = " << count_muon << std::endl;
  
  outputFile->cd();  
  
  h_Zmumu_Mass->Write("h_Zmumu_Mass");
  h_dist_DT_MB1_muon->Write("h_dist_DT_MB1_muon");
  h_dist_DT_MB2_muon->Write("h_dist_DT_MB2_muon");
  h_dist_MB1_layer_1->Write("h_dist_MB1_layer_1");
  h_dist_MB1_layer_2->Write("h_dist_MB1_layer_2");
  h_dist_MB2_layer_1->Write("h_dist_MB2_layer_1");
  h_dist_MB2_layer_2->Write("h_dist_MB2_layer_2");
  h_eff_eta_station_1->Write("h_eff_eta_station_1");
  h_eff_eta_station_2->Write("h_eff_eta_station_2");
  h_eff_pt_station_1->Write("h_eff_pt_station_1");
  h_eff_pt_station_2->Write("h_eff_pt_station_2");
  h_eff_phi_station_1->Write("h_eff_phi_station_1");
  h_eff_phi_station_2->Write("h_eff_phi_station_2");
  h_BX_IN_MB1->Write("h_BX_IN_MB1");
  h_BX_OUT_MB1->Write("h_BX_OUT_MB1");
  h_BX_IN_MB2->Write("h_BX_IN_MB2");
  h_BX_OUT_MB2->Write("h_BX_OUT_MB2");
  h_DeltaPhi_InOut_MB1->Write("h_DeltaPhi_InOut_MB1");
  h_DeltaPhi_InOut_MB2->Write("h_DeltaPhi_InOut_MB2");	
  
  outputFile->Close();

}


RPC_studies::muons_vector RPC_studies::Selection()
{

  vector<int> index;
  vector<float> pt;
  vector<float> eta;
  vector<float> phi;
  bool probe;
  index.clear();
  pt.clear();
  eta.clear();
  phi.clear();
  
  float Mu_pt = -999;
	
	   		
  for(int mu = 0; mu < Nmuons; mu++){
    
    Mu_pt = sqrt(Mu_px->at(mu)*Mu_px->at(mu) + Mu_py->at(mu)*Mu_py->at(mu));
    
    if(Mu_isMuGlobal->at(mu) !=1 || Mu_isMuTracker->at(mu) !=1 ||
       fabs(Mu_dxy_glb->at(mu)) > D0_cut || fabs(Mu_dz_glb->at(mu)) > Dz_cut ||
       fabs(Mu_eta->at(mu)) > eta_cut || (Mu_normchi2_glb->at(mu)) > muchi2_cut ||
       Mu_pt < pt_cut || Mu_pt > pt_max || 
       Mu_numberOfHits_sta->at(mu) < N_hits_cut || 
       Mu_numberOfPixelHits_glb->at(mu) < npix_cut ||
       Mu_numberOfTrackerHits_glb->at(mu) < ntkr_cut) continue;
    
    else{			
      index.push_back(mu);
      pt.push_back(Mu_pt);
      eta.push_back(Mu_eta->at(mu));
      phi.push_back(Mu_phi->at(mu));
      probe = true;
    }
  }
  
  muons_vector good_muons = {index, pt, eta, phi, probe};
					
  return good_muons;	

}

float RPC_studies::Mass(TLorentzVector v1, TLorentzVector v2){
	
  float massWW;
  float ener1, ener2;
  float px1, px2;
  float py1, py2;
  float pz1, pz2;
  
  ener1 = v1.E();
  px1 = v1.Px();
  py1 = v1.Py();
  pz1 = v1.Pz();
  
  ener2 = v2.E();
  px2 = v2.Px();
  py2 = v2.Py();
  pz2 = v2.Pz();
  
  massWW = sqrt((ener1+ener2)*(ener1+ener2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
  
  return massWW;
  
}

float RPC_studies::PhiConversion(float phi_In, int sector)
{	
  double locphi = (((double)phi_In)/4096.0);
  double newphi = locphi + ((sector-1) * (pig/6.));
  if (newphi > pig) newphi -= 2 * pig;
  
  return newphi;
}
