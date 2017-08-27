#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TFile.h"
#include "tdrstyle.C"

TCanvas* compareHistos(TString title, TH1F* histo1, TH1F* histo2, TString leg1, TString leg2)
{

  TCanvas *canvas = new TCanvas(title, title, 500, 500);
  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(1);
  histo1->SetMarkerColor(kRed);
  histo1->SetLineColor(kRed);
  histo1->Scale(1./histo1->Integral());
  histo2->SetMarkerStyle(21);
  histo2->SetMarkerSize(1);
  histo2->SetMarkerColor(kBlue);
  histo2->SetLineColor(kBlue);
  histo2->Scale(1./histo2->Integral());
  histo2->Draw();
  histo1->Draw("same");
  histo2->SetMaximum(1.01);
  TLegend *legend = new TLegend(0.65,0.75,0.9,0.9);   
  legend->SetBorderSize(1);
  legend->AddEntry(histo1, leg1, "lep");
  legend->AddEntry(histo2, leg2, "lep");    
  legend->Draw();

  return canvas;
  
}

TCanvas* compareEfficiencies(TString title, TEfficiency* eff1, TEfficiency* eff2, TString leg1, TString leg2)
{
  
  TCanvas *canvas = new TCanvas(title, title, 500, 500);
  eff1->SetMarkerStyle(20);
  eff2->SetMarkerStyle(21);
  eff1->SetMarkerSize(1);
  eff2->SetMarkerSize(1);
  eff1->SetMarkerColor(kRed);
  eff2->SetMarkerColor(kBlue);
  eff1->SetLineColor(kRed);
  eff2->SetLineColor(kBlue);
  eff1->Draw(); 
  canvas->Update(); 
  auto graph = eff1->GetPaintedGraph(); 
  graph->SetMinimum(0.);
  graph->SetMaximum(1.05); 
  canvas->Update();
  eff2->Draw("same");   
  TLegend *legend = new TLegend(0.75,0.5,0.9,0.6);
  legend->SetBorderSize(1);
  legend->AddEntry(eff1, leg1, "lep");
  legend->AddEntry(eff2, leg2, "lep");    
  legend->Draw();

  return canvas;
  
}


void plotter()
{

  bool save = false;
  TString inputFile = "analysis_results.root";
	
  TFile *file = new TFile(inputFile,"open");
  
  if (!(file->IsOpen())) {
    std::cout<<("File cannot be opened\n");
    return;
  }
	
  ///////////////////////////////////////////////
  // Recommended style macro used as strating /// 
  // point for many CMS plots                 /// 
  ///////////////////////////////////////////////

  setTDRStyle();
  

  //////////////////////////
  // Read the histograms ///
  //////////////////////////
  
  // From kinematics folder
  
  TH1F *h_Zmumu_mass = (TH1F*)gDirectory->Get("/kinematics/h_Zmumu_mass");
  TH1F *h_mu_eta = (TH1F*)gDirectory->Get("/kinematics/h_mu_eta");
  TH1F *h_mu_phi = (TH1F*)gDirectory->Get("/kinematics/h_mu_phi");
  TH1F *h_mu_pt = (TH1F*)gDirectory->Get("/kinematics/h_mu_pt");
  
  // From distance folder (matchings etc ...)
  
  TH1F *h_dR_seg_muon_MB1 = (TH1F*)gDirectory->Get("/distance/h_dR_seg_muon_MB1");
  TH1F *h_dR_seg_muon_MB2 = (TH1F*)gDirectory->Get("/distance/h_dR_seg_muon_MB2");
  
  TH2F *h_dX_dY_seg_muon_MB1 = (TH2F*)gDirectory->Get("/distance/h_dX_dY_seg_muon_MB1");
  TH2F *h_dX_dY_seg_muon_MB2 = (TH2F*)gDirectory->Get("/distance/h_dX_dY_seg_muon_MB2");
  
  TH1F *h_dPhi_seg_TrigIn_MB1 = (TH1F*)gDirectory->Get("/distance/h_dPhi_seg_TrigIn_MB1");
  TH1F *h_dPhi_seg_TrigIn_MB2 = (TH1F*)gDirectory->Get("/distance/h_dPhi_seg_TrigIn_MB2");
  
  TH1F *h_dX_MB1_layer_1 = (TH1F*)gDirectory->Get("/distance/h_dX_MB1_layer_1");
  TH1F *h_dX_MB1_layer_2 = (TH1F*)gDirectory->Get("/distance/h_dX_MB1_layer_2");
  TH1F *h_dX_MB2_layer_1 = (TH1F*)gDirectory->Get("/distance/h_dX_MB2_layer_1");
  TH1F *h_dX_MB2_layer_2 = (TH1F*)gDirectory->Get("/distance/h_dX_MB2_layer_2");
  
  TH1F *h_dPhi_TrigIn_TrigOut_MB1 = (TH1F*)gDirectory->Get("/distance/h_dPhi_TrigIn_TrigOut_MB1");
  TH1F *h_dPhi_TrigIn_TrigOut_MB2 = (TH1F*)gDirectory->Get("/distance/h_dPhi_TrigIn_TrigOut_MB2");
  
  // From efficiency folder
  
  TEfficiency *h_eff_twinmux_in_pt_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_twinmux_in_pt_MB1");
  TEfficiency *h_eff_twinmux_in_pt_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_twinmux_in_pt_MB2");
  
  TEfficiency *h_eff_twinmux_in_eta_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_twinmux_in_eta_MB1");
  TEfficiency *h_eff_twinmux_in_eta_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_twinmux_in_eta_MB2");
  
  TEfficiency *h_eff_twinmux_in_phi_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_twinmux_in_phi_MB1");
  TEfficiency *h_eff_twinmux_in_phi_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_twinmux_in_phi_MB2");
  
  TEfficiency *h_eff_rpc_pt_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_rpc_pt_MB1");
  TEfficiency *h_eff_rpc_pt_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_rpc_pt_MB2");
  
  TEfficiency *h_eff_rpc_eta_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_rpc_eta_MB1");
  TEfficiency *h_eff_rpc_eta_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_rpc_eta_MB2");
  
  TEfficiency *h_eff_rpc_phi_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_rpc_phi_MB1");
  TEfficiency *h_eff_rpc_phi_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_rpc_phi_MB2");

  TEfficiency *h_eff_combined_toy_pt_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_combined_toy_pt_MB1");
  TEfficiency *h_eff_combined_toy_pt_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_combined_toy_pt_MB2");
  
  TEfficiency *h_eff_combined_toy_eta_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_combined_toy_eta_MB1");
  TEfficiency *h_eff_combined_toy_eta_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_combined_toy_eta_MB2");
  
  TEfficiency *h_eff_combined_toy_phi_MB1 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_combined_toy_phi_MB1");
  TEfficiency *h_eff_combined_toy_phi_MB2 = (TEfficiency*)gDirectory->Get("/efficiencies/h_eff_combined_toy_phi_MB2");
  
  // From trigger folder (TwinMux quality and BX)
  
  TH1F *h_BX_twinmux_in_MB1 = (TH1F*)gDirectory->Get("/trigger/h_BX_twinmux_in_MB1");
  TH1F *h_BX_twinmux_in_MB2 = (TH1F*)gDirectory->Get("/trigger/h_BX_twinmux_in_MB2");
  
  TH1F *h_qual_twinmux_in_MB1 = (TH1F*)gDirectory->Get("/trigger/h_qual_twinmux_in_MB1");
  TH1F *h_qual_twinmux_in_MB2 = (TH1F*)gDirectory->Get("/trigger/h_qual_twinmux_in_MB2");
  
  TH1F *h_BX_twinmux_out_MB1 = (TH1F*)gDirectory->Get("/trigger/h_BX_twinmux_out_MB1");
  TH1F *h_BX_twinmux_out_MB2 = (TH1F*)gDirectory->Get("/trigger/h_BX_twinmux_out_MB2");
  
  
  /////////////////////////
  ////  Now plot the  /////
  //// the histograms /////
  /////////////////////////
  
  TCanvas *c_kinematcis = new TCanvas("Muon (pair) kinematcis", "Muon (pair) kinematcis",750,750);
  
  c_kinematcis->Divide(2,2);
  
  c_kinematcis->cd(1);
  h_Zmumu_mass->Draw();
  
  c_kinematcis->cd(2);
  h_mu_eta->Draw();	
  
  c_kinematcis->cd(3);
  h_mu_pt->Draw();
  
  c_kinematcis->cd(4);
  h_mu_phi->SetMinimum(0.);	
  h_mu_phi->Draw();	

  TCanvas *dX_dY_seg_muon_MB1 = new TCanvas("#DeltaX vs #DeltaY: MB1", "#DeltaX vs #DeltaY: MB1", 500 ,500);
  dX_dY_seg_muon_MB1->SetLogz();
  h_dX_dY_seg_muon_MB1->Draw("colz");
  TCanvas *dX_dY_seg_muon_MB2 = new TCanvas("#DeltaX vs #DeltaY: MB2", "#DeltaX vs #DeltaY: MB2", 500,500);
  dX_dY_seg_muon_MB2->SetLogz();
  h_dX_dY_seg_muon_MB2->Draw("colz");

  TCanvas * dR_seg_muon  = compareHistos("Distance muon - DT segment", h_dR_seg_muon_MB1, h_dR_seg_muon_MB2, "MB1", "MB2");
  dR_seg_muon->SetLogy();
  dR_seg_muon->Update();

  TCanvas * dPhi_seg_TrigIn = compareHistos("#Delta#phi DT segment - TwinMux In", h_dPhi_seg_TrigIn_MB1, h_dPhi_seg_TrigIn_MB2, "MB1", "MB2");
  dPhi_seg_TrigIn->SetLogy();
  dPhi_seg_TrigIn->Update();

  TCanvas *dX_MB1_layer_1 = new TCanvas("#DeltaX extrapolation - recHit: MB1, layer 1", "#DeltaX extrapolation - recHit: MB1, layer 1", 210,45,750,500);
  h_dX_MB1_layer_1->Draw();
  TCanvas *dX_MB1_layer_2 = new TCanvas("#DeltaX extrapolation - recHit: MB1, layer 2", "#DeltaX extrapolation - recHit: MB1, layer 2", 210,45,750,500);
  h_dX_MB1_layer_2->Draw();
  TCanvas *dX_MB2_layer_1 = new TCanvas("#DeltaX extrapolation - recHit: MB2, layer 1", "#DeltaX extrapolation - recHit: MB2, layer 1", 210,45,750,500);
  h_dX_MB2_layer_1->Draw();
  TCanvas *dX_MB2_layer_2 = new TCanvas("#DeltaX extrapolation - recHit: MB2, layer 2", "#DeltaX extrapolation - recHit: MB2, layer 2", 210,45,750,500);
  h_dX_MB2_layer_2->Draw();
    
  // TCanvas *dPhi_TrigIn_TrigOut_MB1 = new TCanvas("#Delta#phi TrigIn - TrigOut: MB1", "#Delta#phi TrigIn - TrigOut: MB1", 210,45,750,500);
  // h_dPhi_TrigIn_TrigOut_MB1->Draw();
    
  // TCanvas *dPhi_TrigIn_TrigOut_MB2 = new TCanvas("#Delta#phi TrigIn - TrigOut: MB2", "#Delta#phi TrigIn - TrigOut: MB2", 210,45,750,500);
  // h_dPhi_TrigIn_TrigOut_MB2->Draw();

  compareEfficiencies("RPC eff vs p_{T}", h_eff_rpc_pt_MB1, h_eff_rpc_pt_MB2, "MB1", "MB2");    
  compareEfficiencies("RPC eff vs #eta", h_eff_rpc_eta_MB1, h_eff_rpc_eta_MB2, "MB1", "MB2");    
  compareEfficiencies("RPC eff vs #phi", h_eff_rpc_phi_MB1, h_eff_rpc_phi_MB2, "MB1", "MB2");    

  compareEfficiencies("TwinMux In eff vs p_{T}", h_eff_twinmux_in_pt_MB1, h_eff_twinmux_in_pt_MB2, "MB1", "MB2");    
  compareEfficiencies("TwinMux In eff vs #eta", h_eff_twinmux_in_eta_MB1, h_eff_twinmux_in_eta_MB2, "MB1", "MB2");    
  compareEfficiencies("TwinMux In eff vs #phi", h_eff_twinmux_in_phi_MB1, h_eff_twinmux_in_phi_MB2, "MB1", "MB2");    

  compareEfficiencies("Eff vs p_{T} MB1", (TEfficiency*) h_eff_twinmux_in_pt_MB1->Clone(), h_eff_combined_toy_pt_MB1, "TwinMux", "TwinMux + RPC");    
  compareEfficiencies("Eff vs #eta MB1", (TEfficiency*) h_eff_twinmux_in_eta_MB1->Clone(), h_eff_combined_toy_eta_MB1, "TwinMux", "TwinMux + RPC");    
  compareEfficiencies("Eff vs #phi MB1", (TEfficiency*) h_eff_twinmux_in_phi_MB1->Clone(), h_eff_combined_toy_phi_MB1, "TwinMux", "TwinMux + RPC");    

  compareEfficiencies("Eff vs p_{T} MB2", (TEfficiency*) h_eff_twinmux_in_pt_MB2->Clone(), h_eff_combined_toy_pt_MB2, "TwinMux", "TwinMux + RPC");    
  compareEfficiencies("Eff vs #eta MB2", (TEfficiency*) h_eff_twinmux_in_eta_MB2->Clone(), h_eff_combined_toy_eta_MB2, "TwinMux", "TwinMux + RPC");    
  compareEfficiencies("Eff vs #phi MB2", (TEfficiency*) h_eff_twinmux_in_phi_MB2->Clone(), h_eff_combined_toy_phi_MB2, "TwinMux", "TwinMux + RPC");    
  
  compareHistos("TwinMux IN/OUT BX MB1", h_BX_twinmux_in_MB1, h_BX_twinmux_out_MB1, "TwinMux In", "TwinMux Out");
  compareHistos("TwinMux IN/OUT BX MB2", h_BX_twinmux_in_MB2, h_BX_twinmux_out_MB2, "TwinMux In", "TwinMux Out");
  
  compareHistos("TwinMux IN Quality", h_qual_twinmux_in_MB1, h_qual_twinmux_in_MB2, "MB1", "MB2");

}
