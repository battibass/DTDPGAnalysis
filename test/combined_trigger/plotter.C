#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TFile.h"
#include "tdrstyle.C"

TCanvas* compareHistos(TString title, TH1F* histo1, TH1F* histo2, TString leg1, TString leg2)
{

  TCanvas *canvas = new TCanvas(title, title, 210,45,750,500);
  //histo1->SetStats(0);
  //histo2->SetStats(0);
  //histo1->GetXaxis()->SetTitle("BX");
  //histo2->GetXaxis()->SetTitle("BX");
  //histo1->SetTitle("BX in MB2");
  //histo2->SetTitle("BX in MB2");
  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(1);
  histo1->SetMarkerColor(kRed);
  histo1->SetLineColor(kRed);
  histo1->Scale(1./histo1->Integral());
  histo2->SetMarkerStyle(20);
  histo2->SetMarkerSize(1);
  histo2->SetMarkerColor(kBlue);
  histo2->SetLineColor(kBlue);
  histo2->Scale(1./histo1->Integral());
  histo2->Draw();
  histo1->Draw("same");
  histo1->SetMaximum(1.01);
  histo2->SetMaximum(1.01);
  TLegend *legend = new TLegend(0.65,0.75,0.9,0.9);   
  legend->AddEntry(histo1, leg1, "lep");
  legend->AddEntry(histo2, leg2, "lep");    
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
  h_mu_phi->Draw();	

  /*
    TCanvas *dR_seg_muon_MB1 = new TCanvas("Distance muon - DT segment: MB1", "Distance muon - DT segment: MB1", 210,45,750,500);
    h_dR_seg_muon_MB1->Draw();
    TCanvas *dR_seg_muon_MB2 = new TCanvas("Distance muon - DT segment: MB2", "Distance muon - DT segment: MB2", 210,45,750,500);
    h_dR_seg_muon_MB2->Draw();
    
    TCanvas *dX_dY_seg_muon_MB1 = new TCanvas("#DeltaY vs #DeltaY: MB1", "#DeltaY vs #DeltaY: MB1", 210,45,750,500);
    h_dX_dY_seg_muon_MB1->Draw();
    TCanvas *dX_dY_seg_muon_MB2 = new TCanvas("#DeltaY vs #DeltaY: MB2", "#DeltaY vs #DeltaY: MB2", 210,45,750,500);
    h_dX_dY_seg_muon_MB2->Draw();
    
    TCanvas *dPhi_seg_TrigIn_MB1 = new TCanvas("#Delta#phi DT segment - TrigIn: MB1", "#Delta#phi DT segment - TrigIn: MB1", 210,45,750,500);
    h_dPhi_seg_TrigIn_MB1->Draw();
    TCanvas *dPhi_seg_TrigIn_MB2 = new TCanvas("#Delta#phi DT segment - TrigIn: MB2", "#Delta#phi DT segment - TrigIn: MB2", 210,45,750,500);
    h_dPhi_seg_TrigIn_MB2->Draw();
    
    TCanvas *dX_MB1_layer_1 = new TCanvas("#DeltaX extrapolation - recHit: MB1, layer 1", "#DeltaX extrapolation - recHit: MB1, layer 1", 210,45,750,500);
    h_dX_MB1_layer_1->Draw();
    TCanvas *dX_MB1_layer_2 = new TCanvas("#DeltaX extrapolation - recHit: MB1, layer 2", "#DeltaX extrapolation - recHit: MB1, layer 2", 210,45,750,500);
    h_dX_MB1_layer_2->Draw();
    TCanvas *dX_MB2_layer_1 = new TCanvas("#DeltaX extrapolation - recHit: MB2, layer 1", "#DeltaX extrapolation - recHit: MB2, layer 1", 210,45,750,500);
    h_dX_MB2_layer_1->Draw();
    TCanvas *dX_MB2_layer_2 = new TCanvas("#DeltaX extrapolation - recHit: MB2, layer 2", "#DeltaX extrapolation - recHit: MB2, layer 2", 210,45,750,500);
    h_dX_MB2_layer_2->Draw();
    
    TCanvas *dPhi_TrigIn_TrigOut_MB1 = new TCanvas("#Delta#phi TrigIn - TrigOut: MB1", "#Delta#phi TrigIn - TrigOut: MB1", 210,45,750,500);
    h_dPhi_TrigIn_TrigOut_MB1->Draw();
    
    TCanvas *dPhi_TrigIn_TrigOut_MB2 = new TCanvas("#Delta#phi TrigIn - TrigOut: MB2", "#Delta#phi TrigIn - TrigOut: MB2", 210,45,750,500);
    h_dPhi_TrigIn_TrigOut_MB2->Draw();
    
    TCanvas *eff_pt = new TCanvas("eff vs p_{T}", "eff vs p_{T}", 210,45,750,500);
    h_eff_rpc_pt_MB1->SetMarkerStyle(20);
    h_eff_rpc_pt_MB2->SetMarkerStyle(20);
    h_eff_rpc_pt_MB1->SetMarkerSize(1);
    h_eff_rpc_pt_MB2->SetMarkerSize(1);
    h_eff_rpc_pt_MB1->SetMarkerColor(kRed);
    h_eff_rpc_pt_MB2->SetMarkerColor(kBlue);
    h_eff_rpc_pt_MB1->SetLineColor(kRed);
    h_eff_rpc_pt_MB2->SetLineColor(kBlue);
    TLegend *legend_eff_pt = new TLegend(0.75,0.8,0.9,0.9);
    h_eff_rpc_pt_MB1->Draw(); 
    gPad->Update(); 
    auto graph_pt = h_eff_rpc_pt_MB1->GetPaintedGraph(); 
    graph_pt->SetMinimum(0.5);
    graph_pt->SetMaximum(1.03); 
    gPad->Update();
    h_eff_rpc_pt_MB2->Draw("same");   
    legend_eff_pt->AddEntry(h_eff_rpc_pt_MB1, "MB1", "lep");
    legend_eff_pt->AddEntry(h_eff_rpc_pt_MB2,"MB2", "lep");	
    legend_eff_pt->Draw();
    
    TCanvas *eff_eta = new TCanvas("eff vs #eta", "eff vs #eta", 210,45,750,500);
    h_eff_rpc_eta_MB1->SetMarkerStyle(20);
    h_eff_rpc_eta_MB2->SetMarkerStyle(20);
    h_eff_rpc_eta_MB1->SetMarkerSize(1);
    h_eff_rpc_eta_MB2->SetMarkerSize(1);
    h_eff_rpc_eta_MB1->SetMarkerColor(kRed);
    h_eff_rpc_eta_MB2->SetMarkerColor(kBlue);
    h_eff_rpc_eta_MB1->SetLineColor(kRed);
    h_eff_rpc_eta_MB2->SetLineColor(kBlue);
    TLegend *legend_eff_eta = new TLegend(0.75,0.8,0.9,0.9);
    h_eff_rpc_eta_MB1->Draw(); 
    gPad->Update(); 
    auto graph_eta = h_eff_rpc_eta_MB1->GetPaintedGraph(); 
    graph_eta->SetMinimum(0.5);
    graph_eta->SetMaximum(1.03); 
    gPad->Update();
    h_eff_rpc_eta_MB2->Draw("same");   
    legend_eff_eta->AddEntry(h_eff_rpc_eta_MB1, "MB1", "lep");
    legend_eff_eta->AddEntry(h_eff_rpc_eta_MB2,"MB2", "lep");	
    legend_eff_eta->Draw();
    
    TCanvas *eff_phi = new TCanvas("eff vs #phi", "eff vs #phi", 210,45,750,500);
    h_eff_rpc_phi_MB1->SetMarkerStyle(20);
    h_eff_rpc_phi_MB2->SetMarkerStyle(20);
    h_eff_rpc_phi_MB1->SetMarkerSize(1);
    h_eff_rpc_phi_MB2->SetMarkerSize(1);
    h_eff_rpc_phi_MB1->SetMarkerColor(kRed);
    h_eff_rpc_phi_MB2->SetMarkerColor(kBlue);
    h_eff_rpc_phi_MB1->SetLineColor(kRed);
    h_eff_rpc_phi_MB2->SetLineColor(kBlue);
    TLegend *legend_eff_phi = new TLegend(0.75,0.8,0.9,0.9);
    h_eff_rpc_phi_MB1->Draw(); 
    gPad->Update(); 
    auto graph_phi = h_eff_rpc_phi_MB1->GetPaintedGraph(); 
    graph_phi->SetMinimum(0.5);
    graph_phi->SetMaximum(1.03); 
    gPad->Update();
    h_eff_rpc_phi_MB2->Draw("same");   
    legend_eff_phi->AddEntry(h_eff_rpc_phi_MB1, "MB1", "lep");
    legend_eff_phi->AddEntry(h_eff_rpc_phi_MB2,"MB2", "lep");	
    legend_eff_phi->Draw();
    
    TCanvas *eff_pt_TwinMux = new TCanvas("TwinMux In eff vs p_{T}", "TwinMux In eff vs p_{T}", 210,45,750,500);
    h_eff_twinmux_in_pt_MB1->SetMarkerStyle(20);
    h_eff_twinmux_in_pt_MB2->SetMarkerStyle(20);
    h_eff_twinmux_in_pt_MB1->SetMarkerSize(1);
    h_eff_twinmux_in_pt_MB2->SetMarkerSize(1);
    h_eff_twinmux_in_pt_MB1->SetMarkerColor(kRed);
    h_eff_twinmux_in_pt_MB2->SetMarkerColor(kBlue);
    h_eff_twinmux_in_pt_MB1->SetLineColor(kRed);
    h_eff_twinmux_in_pt_MB2->SetLineColor(kBlue);
    h_eff_twinmux_in_pt_MB1->Draw(); 
    gPad->Update(); 
    auto graph_pt_TwinMux = h_eff_twinmux_in_pt_MB1->GetPaintedGraph(); 
    graph_pt_TwinMux->SetMinimum(0.5);
    graph_pt_TwinMux->SetMaximum(1.03); 
    gPad->Update();
    h_eff_twinmux_in_pt_MB2->Draw("same");   
    TLegend *legend_eff_pt_TwinMux = new TLegend(0.75,0.8,0.9,0.9);
    legend_eff_pt_TwinMux->AddEntry(h_eff_twinmux_in_pt_MB1, "MB1", "lep");
    legend_eff_pt_TwinMux->AddEntry(h_eff_twinmux_in_pt_MB2,"MB2", "lep");	
    legend_eff_pt_TwinMux->Draw();
    
    TCanvas *eff_eta_TwinMux = new TCanvas("TwinMux In eff vs #eta", "TwinMux In eff vs #eta", 210,45,750,500);
    h_eff_twinmux_in_eta_MB1->SetMarkerStyle(20);
    h_eff_twinmux_in_eta_MB2->SetMarkerStyle(20);
    h_eff_twinmux_in_eta_MB1->SetMarkerSize(1);
    h_eff_twinmux_in_eta_MB2->SetMarkerSize(1);
    h_eff_twinmux_in_eta_MB1->SetMarkerColor(kRed);
    h_eff_twinmux_in_eta_MB2->SetMarkerColor(kBlue);
    h_eff_twinmux_in_eta_MB1->SetLineColor(kRed);
    h_eff_twinmux_in_eta_MB2->SetLineColor(kBlue);
    h_eff_twinmux_in_eta_MB1->Draw(); 
    gPad->Update(); 
    auto graph_eta_TwinMux = h_eff_twinmux_in_eta_MB1->GetPaintedGraph(); 
    graph_eta_TwinMux->SetMinimum(0.5);
    graph_eta_TwinMux->SetMaximum(1.03); 
    gPad->Update();
    h_eff_twinmux_in_eta_MB2->Draw("same");   
    TLegend *legend_eff_eta_TwinMux = new TLegend(0.75,0.8,0.9,0.9);
    legend_eff_eta_TwinMux->AddEntry(h_eff_twinmux_in_eta_MB1, "MB1", "lep");
    legend_eff_eta_TwinMux->AddEntry(h_eff_twinmux_in_eta_MB2,"MB2", "lep");	
    legend_eff_eta_TwinMux->Draw();
   
    TCanvas *eff_phi_TwinMux = new TCanvas("TwinMux In eff vs #phi", "TwinMux In eff vs #phi", 210,45,750,500);
    h_eff_twinmux_in_phi_MB1->SetMarkerStyle(20);
    h_eff_twinmux_in_phi_MB2->SetMarkerStyle(20);
    h_eff_twinmux_in_phi_MB1->SetMarkerSize(1);
    h_eff_twinmux_in_phi_MB2->SetMarkerSize(1);
    h_eff_twinmux_in_phi_MB1->SetMarkerColor(kRed);
    h_eff_twinmux_in_phi_MB2->SetMarkerColor(kBlue);
    h_eff_twinmux_in_phi_MB1->SetLineColor(kRed);
    h_eff_twinmux_in_phi_MB2->SetLineColor(kBlue);
    h_eff_twinmux_in_phi_MB1->Draw(); 
    gPad->Update(); 
    auto graph_phi_TwinMux = h_eff_twinmux_in_phi_MB1->GetPaintedGraph(); 
    graph_phi_TwinMux->SetMinimum(0.5);
    graph_phi_TwinMux->SetMaximum(1.03); 
    gPad->Update();
    h_eff_twinmux_in_phi_MB2->Draw("same");   
    TLegend *legend_eff_phi_TwinMux = new TLegend(0.75,0.8,0.9,0.9);
    legend_eff_phi_TwinMux->AddEntry(h_eff_twinmux_in_phi_MB1, "MB1", "lep");
    legend_eff_phi_TwinMux->AddEntry(h_eff_twinmux_in_phi_MB2,"MB2", "lep");	
    legend_eff_phi_TwinMux->Draw();
  */
  
  
  TCanvas *BX_MB2 = new TCanvas("BX in MB2", "BX in MB2", 210,45,750,500);
  //    BX_MB2->SetLogy();
  h_BX_twinmux_in_MB2->SetStats(0);
  h_BX_twinmux_out_MB2->SetStats(0);
  h_BX_twinmux_in_MB2->GetXaxis()->SetTitle("BX");
  h_BX_twinmux_out_MB2->GetXaxis()->SetTitle("BX");
  h_BX_twinmux_in_MB2->SetTitle("BX in MB2");
  h_BX_twinmux_out_MB2->SetTitle("BX in MB2");
  h_BX_twinmux_in_MB2->SetMarkerStyle(20);
  h_BX_twinmux_in_MB2->SetMarkerSize(1);
  h_BX_twinmux_in_MB2->SetMarkerColor(kRed);
  h_BX_twinmux_in_MB2->SetLineColor(kRed);
  h_BX_twinmux_in_MB2->Scale(1/h_BX_twinmux_in_MB2->Integral());
  h_BX_twinmux_out_MB2->SetMarkerStyle(20);
  h_BX_twinmux_out_MB2->SetMarkerSize(1);
  h_BX_twinmux_out_MB2->SetMarkerColor(kBlue);
  h_BX_twinmux_out_MB2->SetLineColor(kBlue);
  h_BX_twinmux_out_MB2->Scale(1/h_BX_twinmux_out_MB2->Integral());
  h_BX_twinmux_out_MB2->Draw();
  h_BX_twinmux_in_MB2->Draw("same");
  h_BX_twinmux_out_MB2->SetMaximum(1.01);
  h_BX_twinmux_in_MB2->SetMaximum(1.01);
  TLegend *legend_BX_MB2 = new TLegend(0.65,0.75,0.9,0.9);   
  legend_BX_MB2->AddEntry(h_BX_twinmux_in_MB2, "TwinMux In", "lep");
  legend_BX_MB2->AddEntry(h_BX_twinmux_out_MB2,"TwinMux Out", "lep");    
  legend_BX_MB2->Draw();
  
  compareHistos(TString("TwinMux in Quality"),
		h_qual_twinmux_in_MB1, h_qual_twinmux_in_MB2,
		TString("MB1"), TString("MB2"))
  

}
