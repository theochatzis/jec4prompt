// Purpose: Read input JES data as TProfile2D and produce Closure Plots
#include "TFile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"

// Custom utils
#include "../../utils.C"

// JSON headers
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <vector>
#include <set>

#include "../../tdrstyle_mod22.C"

// --- Helper Struct for Eta-dependent Pt Binning ---
struct EtaBinConfig {
    double abs_eta_min;
    double abs_eta_max;
    std::vector<double> pt_bins;
};

void DrawClosures(
    TString channel="photonjet",
    TString outputDirectory="",
    TString inputPathMC="",
    TString inputPathData="",
    TString title = "Closure",
    TString luminosity = "1"
  ) 
  {
  // Don't display any graphics on the screen
  gROOT->SetBatch(kTRUE);
  
  // --- READ JSON CONFIGURATION ---
  boost::property_tree::ptree propertyTree;
  try {
      boost::property_tree::read_json("../../constants.json", propertyTree);
  } catch (const std::exception& e) {
      std::cerr << "Error reading constants.json: " << e.what() << std::endl;
      return; 
  }


  double jes_limitMin = propertyTree.get<double>("global.jes_limitMin", 0.82-0.20);
  double jes_limitMax = propertyTree.get<double>("global.jes_limitMax", 1.12+0.20);
  double jes_db_limitMin = propertyTree.get<double>("global.jes_db_limitMin", 0.82-0.20);
  double jes_db_limitMax = propertyTree.get<double>("global.jes_db_limitMax", 1.12+0.20);
  
  std::string ch = channel.Data();
  std::string chPath = "channels." + ch; // e.g., "channels.photonjet"
  
  // Verify the channel exists inside the JSON file
  if (!propertyTree.count("channels") || !propertyTree.get_child("channels").count(ch)) {
      std::cerr << "Error: Channel '" << ch << "' not found in constants.json!" << std::endl;
      return;
  }
  
  // Extract variables for the current channel
  std::string profileName = propertyTree.get<std::string>(chPath + ".profile_name_mpf");
  std::string profileNameDB = propertyTree.get<std::string>(chPath + ".profile_name_db");
  int colorData = propertyTree.get<int>(chPath + ".color_data");
  int colorMC = propertyTree.get<int>(chPath + ".color_mc");

  // Parse pt-eta dependent binning
  std::vector<EtaBinConfig> etaConfigs;
  for (auto& item : propertyTree.get_child(chPath + ".eta_binnings")) {
      EtaBinConfig config;
      config.abs_eta_min = item.second.get<double>("abs_eta_min");
      config.abs_eta_max = item.second.get<double>("abs_eta_max");
      for (auto& pt_item : item.second.get_child("pt_binning")) {
          config.pt_bins.push_back(pt_item.second.get<double>(""));
      }
      etaConfigs.push_back(config);
  }

  // Helper lambda to fetch the correct pt binning for a given eta
  auto getPtBinning = [&](double eta) -> std::vector<double> {
      double abs_eta = std::abs(eta);
      for (const auto& config : etaConfigs) {
          if (abs_eta >= config.abs_eta_min && abs_eta <= config.abs_eta_max + 1e-5) {
              return config.pt_bins;
          }
      }
      if (!etaConfigs.empty()) return etaConfigs.front().pt_bins;
      return {}; 
  };


  // --- BUILD CUSTOM ETA BINNING FROM JSON ---
  std::set<double> unique_abs_etas;
  for (const auto& config : etaConfigs) {
      unique_abs_etas.insert(config.abs_eta_min);
      unique_abs_etas.insert(config.abs_eta_max);
  }
  
  std::vector<double> custom_eta_edges;
  for (double val : unique_abs_etas) {
      if (val > 0) {
          custom_eta_edges.push_back(val);
          custom_eta_edges.push_back(-val);
      } else {
          custom_eta_edges.push_back(0.0);
      }
  }
  std::sort(custom_eta_edges.begin(), custom_eta_edges.end());
  int nCustomEtaBins = custom_eta_edges.size() - 1;
  
  // -------------------------------
  gROOT->ProcessLine(Form(".! mkdir -p %s/", outputDirectory.Data()));
  gROOT->ProcessLine(Form(".! touch %s/", outputDirectory.Data()));

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  // ----------------------------------------
  // INPUT FILES
  // ----------------------------------------  
  TFile *f = new TFile(inputPathData,"READ");
  if (!f || f->IsZombie()) {
      std::cerr << "Warning: Could not open data file " << inputPathData << std::endl;
      if (f) {
          f->Close(); 
          delete f;
      }
      return; 
  }



  TFile *fm = new TFile(inputPathMC, "READ");


  assert(fm && !fm->IsZombie());

  
  curdir->cd();
  
  // ----------------------------------------
  // LOAD MPF AND DB JES MEASURES
  // ----------------------------------------
  TProfile2D *p2_MPF = (TProfile2D*)f->Get("photonjet/MPF_2D"); assert(p2_MPF);
  p2_MPF->SetName("p2_MPF");
  TProfile2D *p2_MPF_MC = (TProfile2D*)fm->Get(profileName.c_str()); assert(p2_MPF_MC);
  p2_MPF_MC->SetName("p2_MPF_MC");


  TProfile2D *p2_DB = (TProfile2D*)f->Get("photonjet/DB_2D"); assert(p2_DB);
  p2_DB->SetName("p2_DB");
  TProfile2D *p2_DB_MC = (TProfile2D*)fm->Get(profileNameDB.c_str()); assert(p2_DB_MC);
  p2_DB_MC->SetName("p2_DB_MC");

  
  // ----------------------------------------
  // REFERENCE PLOTS
  // ----------------------------------------
  double refEtaLimit=0.0; 
  double refEtaLimitPrevious = 0.0;
  
  float fit_region_min = propertyTree.get<float>("global.ref_plot_pt_min", 40.0);
  float fit_region_max = propertyTree.get<float>("global.ref_plot_pt_max", 300.0);
  
  // --- Plots vs pt of the tag object for different eta regions ---
  for (auto& eta_lim : propertyTree.get_child("global.ref_plot_eta_slices")) {
    refEtaLimit = eta_lim.second.get_value<double>();
    
    std::vector<double> v_pt_bins_ref = getPtBinning(refEtaLimitPrevious + 1e-5);
    const Double_t* v_ref = v_pt_bins_ref.data();
    const int n_ref = v_pt_bins_ref.size() - 1;
    
    // --- MPF
    TProfile *p1_MPF_unbinned = GetFoldedPtProfile(p2_MPF, refEtaLimitPrevious, refEtaLimit, Form("p1_MPF_unbinned_%f", refEtaLimit));
    TProfile *p1_MPF_MC_unbinned = GetFoldedPtProfile(p2_MPF_MC, refEtaLimitPrevious, refEtaLimit, Form("p1_MPF_MC_unbinned_%f", refEtaLimit));


    TProfile *p1_MPF_rebin = (TProfile*)p1_MPF_unbinned->Rebin(n_ref, Form("p1_MPF_rebin_%f",refEtaLimit), v_ref);
    TProfile *p1_MPF_MC_rebin = (TProfile*)p1_MPF_MC_unbinned->Rebin(n_ref, Form("p1_MPF_MC_rebin_%f",refEtaLimit), v_ref);

    TH1D *h1_MPF = p1_MPF_rebin->ProjectionX(Form("h1_MPF_%f",refEtaLimit));
    TH1D *h1_MPF_MC = p1_MPF_MC_rebin->ProjectionX(Form("h1_MPF_MC_%f",refEtaLimit));

    h1_MPF->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    h1_MPF_MC->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 

    
    TH1D *h1_ratio = (TH1D*)h1_MPF->Clone(Form("h1_ratio_%f",refEtaLimit));
    TH1D *h1_mc_proj = (TH1D*)h1_MPF_MC->Clone(Form("h1_mc_proj_%f",refEtaLimit));
    h1_ratio->Divide(h1_mc_proj);

    // --- DB
    TProfile *p1_DB_unbinned = GetFoldedPtProfile(p2_DB, refEtaLimitPrevious, refEtaLimit, Form("p1_DB_unbinned_%f", refEtaLimit));
    TProfile *p1_DB_MC_unbinned = GetFoldedPtProfile(p2_DB_MC, refEtaLimitPrevious, refEtaLimit, Form("p1_DB_MC_unbinned_%f", refEtaLimit));


    TProfile *p1_DB_rebin = (TProfile*)p1_DB_unbinned->Rebin(n_ref, Form("p1_DB_rebin_%f",refEtaLimit), v_ref);
    TProfile *p1_DB_MC_rebin = (TProfile*)p1_DB_MC_unbinned->Rebin(n_ref, Form("p1_DB_MC_rebin_%f",refEtaLimit), v_ref);


    TH1D *h1_DB = p1_DB_rebin->ProjectionX(Form("h1_DB_%f",refEtaLimit));
    TH1D *h1_DB_MC = p1_DB_MC_rebin->ProjectionX(Form("h1_DB_MC_%f",refEtaLimit));


    h1_DB->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    h1_DB_MC->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 

    
    TH1D *h1_db_ratio = (TH1D*)h1_DB->Clone(Form("h1_db_ratio_%f",refEtaLimit));
    TH1D *h1_db_mc_proj = (TH1D*)h1_DB_MC->Clone(Form("h1_db_mc_proj_%f",refEtaLimit));
    h1_db_ratio->Divide(h1_db_mc_proj);
    
    // --- PLOTTING --- 
    auto [h1_ratio_min, h1_ratio_max] = GetHistMinMaxWithErrors(h1_ratio);
    TH1D *h_up = tdrHist("h_up", "JES (MPF)", jes_limitMin, jes_limitMax, "", fit_region_min, fit_region_max);
    TH1D *h_dw = tdrHist("h_dw", "Data / MC", h1_ratio_min*0.95, h1_ratio_max*1.05, "p^{tag}_{T} (GeV)", fit_region_min, fit_region_max);

    lumi_136TeV = Form(" %s fb^{-1}", luminosity.Data());
    TCanvas *cPt = tdrDiCanvas("cPt", h_up, h_dw, 8, 11);

    // --- TOP PAD ---
    cPt->cd(1);
    gPad->SetLogx();

    TLatex *tex_eta = new TLatex();
    tex_eta->SetNDC();
    tex_eta->SetTextFont(42);
    tex_eta->SetTextSize(0.045);
    tex_eta->SetTextAlign(11); 
    if (refEtaLimit > 1e-3){
        tex_eta->DrawLatex(0.16, 0.95, Form("%0.1f< |#eta| < %0.1f", refEtaLimitPrevious, refEtaLimit));
    }
    else{
        tex_eta->DrawLatex(0.16, 0.95, Form("|#eta| < %0.1f", refEtaLimit));
    }

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->SetLineColor(kBlack);
    l->DrawLine(fit_region_min,1,fit_region_max,1);


    tdrDraw(p1_MPF_MC_rebin,"HIST",kNone,colorMC,kSolid,-1,kNone,0);

    tdrDraw(p1_MPF_rebin,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0); 
    
    TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.90);
    leg->SetFillStyle(0);  
    leg->SetBorderSize(0); 
    leg->SetTextFont(42);  
    leg->SetTextSize(0.04);
    leg->AddEntry(p1_MPF_rebin, "Data", "pe");

    leg->AddEntry(p1_MPF_MC_rebin, "MC Online", "le");
    leg->Draw();
    
    // --- BOTTOM PAD (RATIO) ---
    cPt->cd(2);
    gPad->SetLogx();

    l->SetLineColor(kGray);
    l->DrawLine(fit_region_min,0.99,fit_region_max,0.99);
    l->SetLineColor(kBlack);
    l->DrawLine(fit_region_min,1,fit_region_max,1);
    l->SetLineColor(kGray);
    l->DrawLine(fit_region_min,1.01,fit_region_max,1.01);

    tdrDraw(h1_ratio, "Pz", kFullCircle, colorData, kSolid, -1, kNone, 0);

    cPt->SaveAs(Form("%s/L2L3Res_Ref_MPF_Et%dp%dto%dp%d.png", outputDirectory.Data(),
    int(floor(refEtaLimitPrevious)), int(round((refEtaLimitPrevious-floor(refEtaLimitPrevious))*1000.)), 
    int(floor(refEtaLimit)), int(round((refEtaLimit-floor(refEtaLimit))*1000.))));

    // Plots with DB
    auto [h1_db_ratio_min, h1_db_ratio_max] = GetHistMinMaxWithErrors(h1_db_ratio);
    TH1D *h_db_up = tdrHist("h_db_up", "JES (DB)", jes_db_limitMin, jes_db_limitMax, "", fit_region_min, fit_region_max);
    TH1D *h_db_dw = tdrHist("h_db_dw", "Data / MC", h1_ratio_min*0.95, h1_ratio_max*1.05, "p^{tag}_{T} (GeV)", fit_region_min, fit_region_max);

    TCanvas *cPt_db = tdrDiCanvas("cPt_db", h_db_up, h_db_dw, 8, 11);

    // --- TOP PAD ---
    cPt_db->cd(1);
    gPad->SetLogx();

    tex_eta->SetNDC();
    tex_eta->SetTextFont(42);
    tex_eta->SetTextSize(0.045);
    tex_eta->SetTextAlign(11); 
    if (refEtaLimit > 1e-3){ 
        tex_eta->DrawLatex(0.16, 0.95, Form("%0.1f< |#eta| < %0.1f", refEtaLimitPrevious, refEtaLimit));
    }
    else{
        tex_eta->DrawLatex(0.16, 0.95, Form("|#eta| < %0.1f", refEtaLimit));
    }

    l->SetLineStyle(kDashed);
    l->SetLineColor(kBlack);
    l->DrawLine(fit_region_min,1,fit_region_max,1);

    tdrDraw(p1_DB_MC_rebin,"HIST",kNone,colorMC,kSolid,-1,kNone,0);
    tdrDraw(p1_DB_rebin,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0); 
    
    TLegend *leg_db = new TLegend(0.55, 0.70, 0.88, 0.90);
    leg_db->SetFillStyle(0);  
    leg_db->SetBorderSize(0); 
    leg_db->SetTextFont(42);  
    leg_db->SetTextSize(0.04);
    leg_db->AddEntry(p1_MPF_rebin, "Data", "pe");

    leg_db->AddEntry(p1_MPF_MC_rebin, "MC Online", "le");
    leg_db->Draw();
    
    // --- BOTTOM PAD (RATIO) ---
    cPt_db->cd(2);
    gPad->SetLogx();

    l->SetLineColor(kGray);
    l->DrawLine(fit_region_min,0.99,fit_region_max,0.99);
    l->SetLineColor(kBlack);
    l->DrawLine(fit_region_min,1,fit_region_max,1);
    l->SetLineColor(kGray);
    l->DrawLine(fit_region_min,1.01,fit_region_max,1.01);

    tdrDraw(h1_db_ratio, "Pz", kFullCircle, colorData, kSolid, -1, kNone, 0);

    cPt_db->SaveAs(Form("%s/L2L3Res_Ref_DB_Et%dp%dto%dp%d.png", outputDirectory.Data(),
    int(floor(refEtaLimitPrevious)), int(round((refEtaLimitPrevious-floor(refEtaLimitPrevious))*1000.)), 
    int(floor(refEtaLimit)), int(round((refEtaLimit-floor(refEtaLimit))*1000.))));
    
    // Cleanup memory
    delete p1_MPF_unbinned; delete p1_MPF_MC_unbinned;
    delete p1_MPF_rebin; delete p1_MPF_MC_rebin;
    delete h1_MPF; delete h1_MPF_MC;
    delete h1_ratio; delete h1_mc_proj;

    delete p1_DB_unbinned; delete p1_DB_MC_unbinned;
    delete p1_DB_rebin; delete p1_DB_MC_rebin;
    delete h1_DB; delete h1_DB_MC;
    delete h1_db_ratio; delete h1_db_mc_proj;

    refEtaLimitPrevious = refEtaLimit;
  } 

  //--- JES vs eta for slices of pT ---
  std::vector<double> pt_cuts(0);
  for (auto& item : propertyTree.get_child("global.ref_plot_pt_slices")) {
      pt_cuts.push_back(item.second.get<double>(""));
  }
  std::vector<int> colors = {
    kBlack,
    kRed,
    kBlue,
    kGreen+2,
    kMagenta,
    kOrange+7,
    kCyan+2
  };

  std::vector<int> mc_colors = {
    kGray+1,
    kRed-9,
    kBlue-9,
    kGreen-9,
    kMagenta-9,
    kOrange-3,
    kCyan-9
  };
  
  // ---- MPF
  TH1D *h_eta_ref_up = tdrHist("h_eta_ref_up", "JES (MPF)", jes_limitMin, jes_limitMax, "#eta", -5.2, 5.2);
  double h_eta_ref_min_yaxis = 0.95;
  double h_eta_ref_max_yaxis = 1.05;
  TH1D *h_eta_ref_dw = tdrHist("h_eta_ref_dw", "Data / MC", h_eta_ref_min_yaxis, h_eta_ref_max_yaxis, "#eta", -5.2, 5.2);
  TCanvas *cEta = tdrDiCanvas("cEta", h_eta_ref_up, h_eta_ref_dw, 8, 11);
  
  cEta->cd(1);
  TLine *line_ref = new TLine();
  line_ref->SetLineStyle(kDashed); line_ref->SetLineColor(kBlack);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);

  TLegend *leg_eta = tdrLeg(0.15, 0.012, 0.90, 0.22);
  leg_eta->SetNColumns(2); 
  
  cEta->cd(2);
  line_ref->SetLineStyle(kDashed); 
  
  line_ref->SetLineColor(kGray);
  line_ref->DrawLine(-5.2, 0.99, 5.2, 0.99);
  line_ref->SetLineColor(kBlack);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);
  line_ref->SetLineColor(kGray);
  line_ref->DrawLine(-5.2, 1.01, 5.2, 1.01);

  for (size_t i = 0; i < pt_cuts.size(); ++i) {
      double pt_cut = pt_cuts[i];
      
      int y_bin_start = p2_MPF->GetYaxis()->FindBin(pt_cut);
      int y_bin_end   = p2_MPF->GetYaxis()->GetNbins() + 1;

      int y_bin_mc_start = p2_MPF_MC->GetYaxis()->FindBin(pt_cut);
      int y_bin_mc_end   = p2_MPF_MC->GetYaxis()->GetNbins() + 1;

      TProfile *p_eta_mc = p2_MPF_MC->ProfileX(Form("p_eta_mc_%d", (int)pt_cut), y_bin_mc_start, y_bin_mc_end);
      TProfile *p_eta_mc_rebin = (TProfile*)p_eta_mc->Rebin(nCustomEtaBins, Form("p_eta_mc_rebin_%d", (int)pt_cut), custom_eta_edges.data());

      TProfile *p_eta_data = p2_MPF->ProfileX(Form("p_eta_data_%d", (int)pt_cut), y_bin_start, y_bin_end);
      TProfile *p_eta_data_rebin = (TProfile*)p_eta_data->Rebin(nCustomEtaBins, Form("p_eta_data_rebin_%d", (int)pt_cut), custom_eta_edges.data());
      
      TH1D *h_eta_data = p_eta_data_rebin->ProjectionX(Form("p_eta_data_rebin_%d", (int)pt_cut));
      TH1D *h_eta_mc = p_eta_mc_rebin->ProjectionX(Form("p_eta_mc_rebin_%d", (int)pt_cut));
      TH1D *h_eta_ratio = (TH1D*) h_eta_data -> Clone(Form("h_eta_ratio_%d", (int)pt_cut));
      h_eta_ratio->Divide(h_eta_mc);
      
      auto [h1_eta_ref_ratio_min, h1_eta_ref_ratio_max] = GetHistMinMaxWithErrors(h_eta_ratio);
      h_eta_ref_min_yaxis = std::min(h_eta_ref_min_yaxis, h1_eta_ref_ratio_min*0.95);
      h_eta_ref_max_yaxis = std::max(h_eta_ref_max_yaxis, h1_eta_ref_ratio_max*1.05);
      h_eta_ref_dw -> GetYaxis() -> SetRangeUser(h_eta_ref_min_yaxis, h_eta_ref_max_yaxis);
      
      cEta->cd(1);
      tdrDraw(p_eta_data_rebin, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
      tdrDraw(p_eta_mc_rebin, "HIST", kNone, mc_colors[i], kSolid, mc_colors[i], 0, 0);
      
      leg_eta->AddEntry(p_eta_data_rebin, "Data", "pe");
      leg_eta->AddEntry(p_eta_mc_rebin, Form("MC   p^{tag}_{T} > %d GeV", (int)pt_cut), "le");
      
      cEta->cd(2);
      std::cout << Form("drawing ratio for pT cut: %d", (int)pt_cut) << std::endl;
      tdrDraw(h_eta_ratio, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
      
  }
  cEta->SaveAs(Form("%s/JES_MPF_vs_Eta_Slices.png", outputDirectory.Data()));

  // ---- DB
  TH1D *h_db_eta_ref_up = tdrHist("h_db_eta_ref_up", "JES (DB)", jes_db_limitMin, jes_db_limitMax, "#eta", -5.2, 5.2);
  double h_db_eta_ref_min_yaxis = 0.95;
  double h_db_eta_ref_max_yaxis = 1.05;
  TH1D *h_db_eta_ref_dw = tdrHist("h_db_eta_ref_dw", "Data / MC", h_db_eta_ref_min_yaxis, h_db_eta_ref_max_yaxis, "#eta", -5.2, 5.2);
  TCanvas *cEta_db = tdrDiCanvas("cEta_db", h_db_eta_ref_up, h_db_eta_ref_dw, 8, 11);
  
  cEta_db->cd(1);
  line_ref->SetLineStyle(kDashed); line_ref->SetLineColor(kBlack);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);

  TLegend *leg_db_eta = tdrLeg(0.15, 0.012, 0.90, 0.22);
  leg_db_eta->SetNColumns(2); 
  
  cEta_db->cd(2);

  line_ref->SetLineColor(kGray);
  line_ref->DrawLine(-5.2, 0.99, 5.2, 0.99);
  line_ref->SetLineColor(kBlack);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);
  line_ref->SetLineColor(kGray);
  line_ref->DrawLine(-5.2, 1.01, 5.2, 1.01);

  for (size_t i = 0; i < pt_cuts.size(); ++i) {
      double pt_cut = pt_cuts[i];
      
      int y_bin_start = p2_DB->GetYaxis()->FindBin(pt_cut);
      int y_bin_end   = p2_DB->GetYaxis()->GetNbins() + 1;

      int y_bin_mc_start = p2_DB_MC->GetYaxis()->FindBin(pt_cut);
      int y_bin_mc_end   = p2_DB_MC->GetYaxis()->GetNbins() + 1;

      TProfile *p_db_eta_mc = p2_DB_MC->ProfileX(Form("p_db_eta_mc_%d", (int)pt_cut), y_bin_mc_start, y_bin_mc_end);
      TProfile *p_db_eta_mc_rebin = (TProfile*)p_db_eta_mc->Rebin(nCustomEtaBins, Form("p_db_eta_mc_rebin_%d", (int)pt_cut), custom_eta_edges.data());

      TProfile *p_db_eta_data = p2_DB->ProfileX(Form("p_db_eta_data_%d", (int)pt_cut), y_bin_start, y_bin_end);
      TProfile *p_db_eta_data_rebin = (TProfile*)p_db_eta_data->Rebin(nCustomEtaBins, Form("p_db_eta_data_rebin_%d", (int)pt_cut), custom_eta_edges.data());
      
      TH1D *h_db_eta_data = p_db_eta_data_rebin->ProjectionX(Form("p_db_eta_data_rebin_%d", (int)pt_cut));
      TH1D *h_db_eta_mc = p_db_eta_mc_rebin->ProjectionX(Form("p_db_eta_mc_rebin_%d", (int)pt_cut));
      TH1D *h_db_eta_ratio = (TH1D*) h_db_eta_data -> Clone(Form("h_db_eta_ratio_%d", (int)pt_cut));
      h_db_eta_ratio->Divide(h_db_eta_mc);
      
      auto [h1_db_eta_ref_ratio_min, h1_db_eta_ref_ratio_max] = GetHistMinMaxWithErrors(h_db_eta_ratio);
      h_db_eta_ref_min_yaxis = std::min(h_db_eta_ref_min_yaxis, h1_db_eta_ref_ratio_min*0.95);
      h_db_eta_ref_max_yaxis = std::max(h_db_eta_ref_max_yaxis, h1_db_eta_ref_ratio_max*1.05);
      h_db_eta_ref_dw -> GetYaxis() -> SetRangeUser(h_db_eta_ref_min_yaxis, h_db_eta_ref_max_yaxis);
      
      cEta_db->cd(1);
      tdrDraw(p_db_eta_data_rebin, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
      tdrDraw(p_db_eta_mc_rebin, "HIST", kNone, mc_colors[i], kSolid, mc_colors[i], 0, 0);
      
      leg_db_eta->AddEntry(p_db_eta_data_rebin, "Data", "pe");
      leg_db_eta->AddEntry(p_db_eta_mc_rebin, Form("MC   p^{tag}_{T} > %d GeV", (int)pt_cut), "le");
      
      cEta_db->cd(2);
      tdrDraw(h_db_eta_ratio, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
  }
  cEta_db->SaveAs(Form("%s/JES_DB_vs_Eta_Slices.png", outputDirectory.Data()));
  
  // Close and delete files
  f->Close();    delete f;
  fm->Close();   delete fm;
  
}