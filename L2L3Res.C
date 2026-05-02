// Purpose: Read input JES data as TProfile2D and produce JEC L2L3Res fit
//          Focus on JEC stability and robustness
#include "TFile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

// Custom utils
#include "utils.C"

// JSON headers
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <vector>

#include "tdrstyle_mod22.C"


// --- Helper Struct for Eta-dependent Pt Binning ---
struct EtaBinConfig {
    double abs_eta_min;
    double abs_eta_max;
    std::vector<double> pt_bins;
};

void L2L3Res(int run = 398027, TString basePath="2025G", TString channel="photonjet") {
  // Don't display any graphics on the screen
  gROOT->SetBatch(kTRUE);
  // --- READ JSON CONFIGURATION ---
  boost::property_tree::ptree propertyTree;
  try {
      boost::property_tree::read_json("constants.json", propertyTree);
  } catch (const std::exception& e) {
      std::cerr << "Error reading constants.json: " << e.what() << std::endl;
      return; 
  }
  std::string outputBaseDirectory = propertyTree.get<std::string>("global.outputBaseDirectory", "");
  // Extract Global Variables
  // note : second argument acts as a fall-back in case it isn't found.
  std::string jsonWithLumis_path = propertyTree.get<std::string>("global.jsonWithLumis_path", "");
  std::string runsDirectoriesBase = propertyTree.get<std::string>("global.runsDirectoriesBase", "");
  //Get luminosity 
  float luminosity = 0.;
  try {
      boost::property_tree::ptree groupingJsonTree;
      boost::property_tree::read_json(jsonWithLumis_path, groupingJsonTree);
      luminosity = groupingJsonTree.get<float>(Form("runs.%d.recorded_lumi", run)); // it is pb^-1
      luminosity *= 0.001;
  } catch (const std::exception& e) {
      std::cerr << Form("Error reading %s: ", jsonWithLumis_path.c_str()) << e.what() << std::endl;
      return; 
  }


  double minimum_luminosity = propertyTree.get<double>("global.minimum_luminosity", 0.0);
  
  if (luminosity < minimum_luminosity){
    std::cout << Form("Won't derive JECs. Run has luminosity smaller than :%.3f" , minimum_luminosity) << std::endl;
    return;
  }

  double jes_limitMin = propertyTree.get<double>("global.jes_limitMin", 0.82-0.20);
  double jes_limitMax = propertyTree.get<double>("global.jes_limitMax", 1.12+0.20);
  double jes_db_limitMin = propertyTree.get<double>("global.jes_db_limitMin", 0.82-0.20);
  double jes_db_limitMax = propertyTree.get<double>("global.jes_db_limitMax", 1.12+0.20);

  double precision_tolerance = propertyTree.get<double>("global.precision_tolerance", 0.20);

  
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
  double xmin = propertyTree.get<double>(chPath + ".xmin");
  double xmax = propertyTree.get<double>(chPath + ".xmax");
  int colorData = propertyTree.get<int>(chPath + ".color_data");
  int colorMC = propertyTree.get<int>(chPath + ".color_mc");

  // Parse pt-eta dependent binning into the struct which has the eta bins (min,max) and pT as a vector 
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
      // Fallback to first configuration if no match is found
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
  gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%d", outputBaseDirectory.c_str(),basePath.Data(), run));
  gROOT->ProcessLine(Form(".! touch %s/%s/%d", outputBaseDirectory.c_str(),basePath.Data(), run));

  //gROOT->ProcessLine(Form(".! mkdir %s/L2L3Res", basePath.Data());
  //gROOT->ProcessLine(Form(".! touch %s/L2L3Res", basePath.Data()));

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  // ----------------------------------------
  // INPUT FILES
  // ----------------------------------------
  // Load input file, using the channel variable for the path as well
  TString dataPath = Form("%s/Run%s/run%d/%s/J4PHists_runs%dto%d_%s.root", 
                          runsDirectoriesBase.c_str(), basePath.Data(), run, channel.Data(), run, run, channel.Data());
  
  //TFile *f = new TFile("rootfiles/J4PHists_runs392175to392175_photonjet.root","READ");
  TFile *f = new TFile(dataPath,"READ");
  
  //assert(f && !f->IsZombie());
  if (!f || f->IsZombie()) {
      std::cerr << "Warning: Could not open data file for run " << run << ". Skipping..." << std::endl;
      if (f) {
          f->Close(); 
          delete f; // Prevent memory leaks
      }
      return; // Skips the rest of the code and goes to the next run in the loop
  }

  std::string fileMCOnline = propertyTree.get<std::string>(chPath + ".mc_online_file");
  std::string fileMCOffline = propertyTree.get<std::string>(chPath + ".mc_offline_file");
  std::string profileMCOffline = propertyTree.get<std::string>(chPath + ".profile_name_mpf_offline_mc");
  std::string profileMCOfflineDB = propertyTree.get<std::string>(chPath + ".profile_name_db_offline_mc");

  TFile *fm = new TFile(fileMCOnline.c_str(), "READ");
  TFile *fmOffline = new TFile(fileMCOffline.c_str(), "READ");

  assert(fm && !fm->IsZombie());
  assert(fmOffline && !fmOffline->IsZombie());
  
  curdir->cd();
  
  // ----------------------------------------
  // LOAD MPF AND DB JES MEASURES
  // ----------------------------------------
  // Load input data (MPF, DB) from file
  TProfile2D *p2_MPF = (TProfile2D*)f->Get(profileName.c_str()); assert(p2_MPF);
  p2_MPF->SetName("p2_MPF");
  TProfile2D *p2_MPF_MC = (TProfile2D*)fm->Get(profileName.c_str()); assert(p2_MPF_MC);
  p2_MPF_MC->SetName("p2_MPF_MC");
  TProfile2D *p2_MPF_OfflineMC = (TProfile2D*)fmOffline->Get(profileMCOffline.c_str()); assert(p2_MPF_OfflineMC);
  p2_MPF_OfflineMC->SetName("p2_MPF_OfflineMC");

  TProfile2D *p2_DB = (TProfile2D*)f->Get(profileNameDB.c_str()); assert(p2_DB);
  p2_DB->SetName("p2_DB");
  TProfile2D *p2_DB_MC = (TProfile2D*)fm->Get(profileNameDB.c_str()); assert(p2_DB_MC);
  p2_DB_MC->SetName("p2_DB_MC");
  TProfile2D *p2_DB_OfflineMC = (TProfile2D*)fmOffline->Get(profileMCOfflineDB.c_str()); assert(p2_DB_OfflineMC);
  p2_DB_OfflineMC->SetName("p2_DB_OfflineMC");
  
  // ----------------------------------------
  // REFERENCE PLOTS
  // ----------------------------------------
  // Make reference plots before the corrections.
  // --- Plots vs pt of the tag object for different eta regions ---------------------------------
  // Loop over reference eta regions
  double refEtaLimit=0.0; 
  double refEtaLimitPrevious = 0.0;
  
  float fit_region_min = propertyTree.get<float>("global.ref_plot_pt_min", 40.0);
  float fit_region_max = propertyTree.get<float>("global.ref_plot_pt_max", 300.0);
  
  for (auto& eta_lim : propertyTree.get_child("global.ref_plot_eta_slices")) {
    //std::cout << eta_lim.second.get_value<double>() << std::endl;
    // Find the indices of each eta limit
    refEtaLimit = eta_lim.second.get_value<double>();
    
    // Get pt binning
    std::vector<double> v_pt_bins_ref = getPtBinning(refEtaLimitPrevious + 1e-5);
    const Double_t* v_ref = v_pt_bins_ref.data();
    const int n_ref = v_pt_bins_ref.size() - 1;
    
    // --- MPF
    // Folded 1D slices of pT in abs eta bins
    TProfile *p1_MPF_unbinned = GetFoldedPtProfile(p2_MPF, refEtaLimitPrevious, refEtaLimit, Form("p1_MPF_unbinned_%f", refEtaLimit));
    TProfile *p1_MPF_MC_unbinned = GetFoldedPtProfile(p2_MPF_MC, refEtaLimitPrevious, refEtaLimit, Form("p1_MPF_MC_unbinned_%f", refEtaLimit));
    TProfile *p1_MPF_Offline_unbinned = GetFoldedPtProfile(p2_MPF_OfflineMC, refEtaLimitPrevious, refEtaLimit, Form("p1_MPF_Offline_unbinned_%f", refEtaLimit));

    // Rebin profiles
    TProfile *p1_MPF_rebin = (TProfile*)p1_MPF_unbinned->Rebin(n_ref, Form("p1_MPF_rebin_%f",refEtaLimit), v_ref);
    TProfile *p1_MPF_MC_rebin = (TProfile*)p1_MPF_MC_unbinned->Rebin(n_ref, Form("p1_MPF_MC_rebin_%f",refEtaLimit), v_ref);
    TProfile *p1_MPF_Offline_rebin = (TProfile*)p1_MPF_Offline_unbinned->Rebin(n_ref, Form("p1_MPF_Offline_rebin_%f",refEtaLimit), v_ref);

    // Project TProfiles to TH1D before dividing for correct error propagation
    TH1D *h1_MPF = p1_MPF_rebin->ProjectionX(Form("h1_MPF_%f",refEtaLimit));
    TH1D *h1_MPF_MC = p1_MPF_MC_rebin->ProjectionX(Form("h1_MPF_MC_%f",refEtaLimit));
    TH1D *h1_MPF_OfflineMC = p1_MPF_Offline_rebin->ProjectionX(Form("h1_MPF_MCOffline_%f",refEtaLimit));

    h1_MPF->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    h1_MPF_MC->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    h1_MPF_OfflineMC->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    
    // Calculate Ratio
    TH1D *h1_ratio = (TH1D*)h1_MPF->Clone(Form("h1_ratio_%f",refEtaLimit));
    TH1D *h1_mc_proj = (TH1D*)h1_MPF_MC->Clone(Form("h1_mc_proj_%f",refEtaLimit));
    h1_ratio->Divide(h1_mc_proj);

    // --- DB
    // Folded 1D slices of pT in abs eta bins
    TProfile *p1_DB_unbinned = GetFoldedPtProfile(p2_DB, refEtaLimitPrevious, refEtaLimit, Form("p1_DB_unbinned_%f", refEtaLimit));
    TProfile *p1_DB_MC_unbinned = GetFoldedPtProfile(p2_DB_MC, refEtaLimitPrevious, refEtaLimit, Form("p1_DB_MC_unbinned_%f", refEtaLimit));
    TProfile *p1_DB_Offline_unbinned = GetFoldedPtProfile(p2_DB_OfflineMC, refEtaLimitPrevious, refEtaLimit, Form("p1_DB_Offline_unbinned_%f", refEtaLimit));

    // Rebin profiles
    TProfile *p1_DB_rebin = (TProfile*)p1_DB_unbinned->Rebin(n_ref, Form("p1_DB_rebin_%f",refEtaLimit), v_ref);
    TProfile *p1_DB_MC_rebin = (TProfile*)p1_DB_MC_unbinned->Rebin(n_ref, Form("p1_DB_MC_rebin_%f",refEtaLimit), v_ref);
    TProfile *p1_DB_Offline_rebin = (TProfile*)p1_DB_Offline_unbinned->Rebin(n_ref, Form("p1_DB_Offline_rebin_%f",refEtaLimit), v_ref);

    // Project TProfiles to TH1D before dividing for correct error propagation
    TH1D *h1_DB = p1_DB_rebin->ProjectionX(Form("h1_DB_%f",refEtaLimit));
    TH1D *h1_DB_MC = p1_DB_MC_rebin->ProjectionX(Form("h1_DB_MC_%f",refEtaLimit));
    TH1D *h1_DB_OfflineMC = p1_DB_Offline_rebin->ProjectionX(Form("h1_DB_MCOffline_%f",refEtaLimit));

    h1_DB->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    h1_DB_MC->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max); 
    h1_DB_OfflineMC->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);
    
    // Calculate Ratio
    TH1D *h1_db_ratio = (TH1D*)h1_DB->Clone(Form("h1_db_ratio_%f",refEtaLimit));
    TH1D *h1_db_mc_proj = (TH1D*)h1_DB_MC->Clone(Form("h1_db_mc_proj_%f",refEtaLimit));
    h1_db_ratio->Divide(h1_db_mc_proj);
    
    // --- PLOTTING --- 
    auto [h1_ratio_min, h1_ratio_max] = GetHistMinMaxWithErrors(h1_ratio);
    // Create top and bottom dummy histograms for tdrDiCanvas
    TH1D *h_up = tdrHist("h_up", "JES (MPF)", jes_limitMin, jes_limitMax, "", fit_region_min, fit_region_max);
    TH1D *h_dw = tdrHist("h_dw", "Data / MC", h1_ratio_min*0.95, h1_ratio_max*1.05, "p^{tag}_{T} (GeV)", fit_region_min, fit_region_max);

    lumi_136TeV = Form("Run%d, %.2f fb^{-1}", run, luminosity);
    TCanvas *cPt = tdrDiCanvas("cPt", h_up, h_dw, 8, 11);

    // --- TOP PAD ---
    cPt->cd(1);
    gPad->SetLogx();

    TLatex *tex_eta = new TLatex();
    tex_eta->SetNDC();
    tex_eta->SetTextFont(42);
    tex_eta->SetTextSize(0.045);
    // Alignment 11: 1 = Left aligned horizontally, 1 = Bottom aligned vertically
    tex_eta->SetTextAlign(11); 
    // Draw at X=0.16 (aligned with the Y-axis) and Y=0.95 (in the top margin)
    if (refEtaLimit > 1e-3){ // handle the barrel case
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
    tdrDraw(p1_MPF_Offline_unbinned,"HIST",kNone,kGreen,kSolid,-1,kNone,0);
    //tdrDraw(h1_MPF_OfflineMC_corr,"HIST",kNone,kBlue,kSolid,-1,kNone,0);
    //tdrDraw(p1_MPF,"Pz",kOpenSquare,colorData,kSolid,-1,kNone,0);
    tdrDraw(p1_MPF_rebin,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0); //tdrDraw(p1_MPF_rebin,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
    //tdrDraw(h1_MPF_cut,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0);
    
    // Make legends
    TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.90);
    leg->SetFillStyle(0);   // Transparent background
    leg->SetBorderSize(0);  // No border
    leg->SetTextFont(42);   // Standard CMS TDR font
    leg->SetTextSize(0.04);

    //leg->AddEntry(h1_MPF_cut, "Data (fit region)", "pe");
    leg->AddEntry(p1_MPF_rebin, "Data", "pe");
    //leg->AddEntry(p1_MPF, "Data (unrebinned)", "pe");
    leg->AddEntry(p1_MPF_Offline_unbinned, "MC Offline", "le");
    leg->AddEntry(p1_MPF_MC_rebin, "MC Online", "le");
    
    leg->Draw();
    
    // --- BOTTOM PAD (RATIO) ---
    cPt->cd(2);
    gPad->SetLogx();

    l->DrawLine(fit_region_min,1,fit_region_max,1);
    tdrDraw(h1_ratio, "Pz", kFullCircle, colorData, kSolid, -1, kNone, 0);

    // Save to pdf
    cPt->SaveAs(Form("%s/%s/%d/L2L3Res_Ref_MPF_Et%dp%dto%dp%d_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run,
    int(floor(refEtaLimitPrevious)), int(round((refEtaLimitPrevious-floor(refEtaLimitPrevious))*1000.)), 
    int(floor(refEtaLimit)), int(round((refEtaLimit-floor(refEtaLimit))*1000.)), 
    run));

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
    // Alignment 11: 1 = Left aligned horizontally, 1 = Bottom aligned vertically
    tex_eta->SetTextAlign(11); 
    // Draw at X=0.16 (aligned with the Y-axis) and Y=0.95 (in the top margin)
    if (refEtaLimit > 1e-3){ // handle the barrel case
        tex_eta->DrawLatex(0.16, 0.95, Form("%0.1f< |#eta| < %0.1f", refEtaLimitPrevious, refEtaLimit));
    }
    else{
        tex_eta->DrawLatex(0.16, 0.95, Form("|#eta| < %0.1f", refEtaLimit));
    }

    l->SetLineStyle(kDashed);
    l->SetLineColor(kBlack);
    l->DrawLine(fit_region_min,1,fit_region_max,1);

    tdrDraw(p1_DB_MC_rebin,"HIST",kNone,colorMC,kSolid,-1,kNone,0);
    tdrDraw(p1_DB_Offline_unbinned,"HIST",kNone,kGreen,kSolid,-1,kNone,0);
    //tdrDraw(h1_DB_OfflineMC_corr,"HIST",kNone,kBlue,kSolid,-1,kNone,0);
    //tdrDraw(p1_DB,"Pz",kOpenSquare,colorData,kSolid,-1,kNone,0);
    tdrDraw(p1_DB_rebin,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0); //tdrDraw(p1_DB_rebin,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
    //tdrDraw(h1_DB_cut,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0);
    
    // Make legends
    TLegend *leg_db = new TLegend(0.55, 0.70, 0.88, 0.90);
    leg_db->SetFillStyle(0);   // Transparent background
    leg_db->SetBorderSize(0);  // No border
    leg_db->SetTextFont(42);   // Standard CMS TDR font
    leg_db->SetTextSize(0.04);

    //leg->AddEntry(h1_MPF_cut, "Data (fit region)", "pe");
    leg_db->AddEntry(p1_MPF_rebin, "Data", "pe");
    //leg->AddEntry(p1_MPF, "Data (unrebinned)", "pe");
    leg_db->AddEntry(p1_MPF_Offline_unbinned, "MC Offline", "le");
    leg_db->AddEntry(p1_MPF_MC_rebin, "MC Online", "le");
    
    leg_db->Draw();
    
    // --- BOTTOM PAD (RATIO) ---
    cPt_db->cd(2);
    gPad->SetLogx();

    l->DrawLine(fit_region_min,1,fit_region_max,1);
    tdrDraw(h1_db_ratio, "Pz", kFullCircle, colorData, kSolid, -1, kNone, 0);

    // Save to pdf
    cPt_db->SaveAs(Form("%s/%s/%d/L2L3Res_Ref_DB_Et%dp%dto%dp%d_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run,
    int(floor(refEtaLimitPrevious)), int(round((refEtaLimitPrevious-floor(refEtaLimitPrevious))*1000.)), 
    int(floor(refEtaLimit)), int(round((refEtaLimit-floor(refEtaLimit))*1000.)), 
    run));
    
    // Cleanup memory
    delete p1_MPF_unbinned; delete p1_MPF_MC_unbinned; delete p1_MPF_Offline_unbinned;
    delete p1_MPF_rebin; delete p1_MPF_MC_rebin; delete p1_MPF_Offline_rebin;
    delete h1_MPF; delete h1_MPF_MC; delete h1_MPF_OfflineMC; 
    delete h1_ratio; delete h1_mc_proj;

    delete p1_DB_unbinned; delete p1_DB_MC_unbinned; delete p1_DB_Offline_unbinned;
    delete p1_DB_rebin; delete p1_DB_MC_rebin; delete p1_DB_Offline_rebin;
    delete h1_DB; delete h1_DB_MC; delete h1_DB_OfflineMC; 
    delete h1_db_ratio; delete h1_db_mc_proj;

    // Set the previous eta limit to move on to next eta bin
    refEtaLimitPrevious = refEtaLimit;
  } // End of loop over reference eta bins


  //--- JES vs eta for slices of pT --------------------------------------------------
  std::vector<double> pt_cuts(0);
  for (auto& item : propertyTree.get_child("global.ref_plot_pt_slices")) {
      pt_cuts.push_back(item.second.get<double>(""));
  }
  std::vector<int> colors = {kBlack, kRed, kBlue};
  std::vector<int> mc_colors = {kGray+1, kRed-9, kBlue-9};
  // ---- MPF
  // Create dummy for the eta-response plot
  TH1D *h_eta_ref_up = tdrHist("h_eta_ref_up", "JES (MPF)", jes_limitMin, jes_limitMax, "#eta", -5.2, 5.2);
  double h_eta_ref_min_yaxis = 0.95;
  double h_eta_ref_max_yaxis = 1.05;
  TH1D *h_eta_ref_dw = tdrHist("h_eta_ref_dw", "Data / MC", h_eta_ref_min_yaxis, h_eta_ref_max_yaxis, "#eta", -5.2, 5.2);
  TCanvas *cEta = tdrDiCanvas("cEta", h_eta_ref_up, h_eta_ref_dw, 8, 11);
  
  // --- TOP PAD ---
  cEta->cd(1);
  TLine *line_ref = new TLine();
  line_ref->SetLineStyle(kDashed); line_ref->SetLineColor(kBlack);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);

  TLegend *leg_eta = tdrLeg(0.15, 0.012, 0.90, 0.22);
  leg_eta->SetNColumns(2); // Two columns to separate Data and MC nicely
  
  // --- BOTTOM PAD ---
  cEta->cd(2);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);

  for (size_t i = 0; i < pt_cuts.size(); ++i) {
      double pt_cut = pt_cuts[i];
      
      // Find Y-bins for pT integration (from cut bin to overflow)
      int y_bin_start = p2_MPF->GetYaxis()->FindBin(pt_cut);
      int y_bin_end   = p2_MPF->GetYaxis()->GetNbins() + 1;

      // Project and REBIN MC
      TProfile *p_eta_mc = p2_MPF_MC->ProfileX(Form("p_eta_mc_%d", (int)pt_cut), y_bin_start, y_bin_end);
      TProfile *p_eta_mc_rebin = (TProfile*)p_eta_mc->Rebin(nCustomEtaBins, Form("p_eta_mc_rebin_%d", (int)pt_cut), custom_eta_edges.data());

      // Project and REBIN Data
      TProfile *p_eta_data = p2_MPF->ProfileX(Form("p_eta_data_%d", (int)pt_cut), y_bin_start, y_bin_end);
      TProfile *p_eta_data_rebin = (TProfile*)p_eta_data->Rebin(nCustomEtaBins, Form("p_eta_data_rebin_%d", (int)pt_cut), custom_eta_edges.data());
      
      // Calculate Data/MC ratio
      TH1D *h_eta_data = p_eta_data_rebin->ProjectionX(Form("p_eta_data_rebin_%d", (int)pt_cut));
      TH1D *h_eta_mc = p_eta_mc_rebin->ProjectionX(Form("p_eta_mc_rebin_%d", (int)pt_cut));
      TH1D *h_eta_ratio = (TH1D*) h_eta_data -> Clone(Form("h_eta_ratio_%d", (int)pt_cut));
      h_eta_ratio->Divide(h_eta_mc);
      
      // automatic y-axis in ratio
      auto [h1_eta_ref_ratio_min, h1_eta_ref_ratio_max] = GetHistMinMaxWithErrors(h_eta_ratio);
      h_eta_ref_min_yaxis = std::min(h_eta_ref_min_yaxis, h1_eta_ref_ratio_min*0.95);
      h_eta_ref_max_yaxis = std::max(h_eta_ref_max_yaxis, h1_eta_ref_ratio_max*1.05);
      h_eta_ref_dw -> GetYaxis() -> SetRangeUser(h_eta_ref_min_yaxis, h_eta_ref_max_yaxis);
      //--- Top pad
      cEta->cd(1);

      tdrDraw(p_eta_data_rebin, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
      tdrDraw(p_eta_mc_rebin, "HIST", kNone, mc_colors[i], kSolid, mc_colors[i], 0, 0);
      
      // Pass the rebinned and styled objects to the legend, not the original ones!
      leg_eta->AddEntry(p_eta_data_rebin, "Data", "pe");
      leg_eta->AddEntry(p_eta_mc_rebin, Form("MC   p^{tag}_{T} > %d GeV", (int)pt_cut), "le");
      
      //--- Bottom pad
      cEta->cd(2);
      tdrDraw(h_eta_ratio, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);

  }

  cEta->SaveAs(Form("%s/%s/%d/JES_MPF_vs_Eta_Slices_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run, run));

  // ---- DB
  // Create dummy for the eta-response plot
  TH1D *h_db_eta_ref_up = tdrHist("h_db_eta_ref_up", "JES (DB)", jes_db_limitMin, jes_db_limitMax, "#eta", -5.2, 5.2);
  double h_db_eta_ref_min_yaxis = 0.95;
  double h_db_eta_ref_max_yaxis = 1.05;
  TH1D *h_db_eta_ref_dw = tdrHist("h_db_eta_ref_dw", "Data / MC", h_db_eta_ref_min_yaxis, h_db_eta_ref_max_yaxis, "#eta", -5.2, 5.2);
  TCanvas *cEta_db = tdrDiCanvas("cEta_db", h_db_eta_ref_up, h_db_eta_ref_dw, 8, 11);
  
  // --- TOP PAD ---
  cEta_db->cd(1);
  line_ref->SetLineStyle(kDashed); line_ref->SetLineColor(kBlack);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);

  TLegend *leg_db_eta = tdrLeg(0.15, 0.012, 0.90, 0.22);
  leg_db_eta->SetNColumns(2); // Two columns to separate Data and MC nicely
  
  // --- BOTTOM PAD ---
  cEta_db->cd(2);
  line_ref->DrawLine(-5.2, 1.0, 5.2, 1.0);

  for (size_t i = 0; i < pt_cuts.size(); ++i) {
      double pt_cut = pt_cuts[i];
      
      // Find Y-bins for pT integration (from cut bin to overflow)
      int y_bin_start = p2_DB->GetYaxis()->FindBin(pt_cut);
      int y_bin_end   = p2_DB->GetYaxis()->GetNbins() + 1;

      // Project and REBIN MC
      TProfile *p_db_eta_mc = p2_DB_MC->ProfileX(Form("p_db_eta_mc_%d", (int)pt_cut), y_bin_start, y_bin_end);
      TProfile *p_db_eta_mc_rebin = (TProfile*)p_db_eta_mc->Rebin(nCustomEtaBins, Form("p_db_eta_mc_rebin_%d", (int)pt_cut), custom_eta_edges.data());

      // Project and REBIN Data
      TProfile *p_db_eta_data = p2_DB->ProfileX(Form("p_db_eta_data_%d", (int)pt_cut), y_bin_start, y_bin_end);
      TProfile *p_db_eta_data_rebin = (TProfile*)p_db_eta_data->Rebin(nCustomEtaBins, Form("p_db_eta_data_rebin_%d", (int)pt_cut), custom_eta_edges.data());
      
      // Calculate Data/MC ratio
      TH1D *h_db_eta_data = p_db_eta_data_rebin->ProjectionX(Form("p_db_eta_data_rebin_%d", (int)pt_cut));
      TH1D *h_db_eta_mc = p_db_eta_mc_rebin->ProjectionX(Form("p_db_eta_mc_rebin_%d", (int)pt_cut));
      TH1D *h_db_eta_ratio = (TH1D*) h_db_eta_data -> Clone(Form("h_db_eta_ratio_%d", (int)pt_cut));
      h_db_eta_ratio->Divide(h_db_eta_mc);
      
      // automatic y-axis in ratio
      auto [h1_db_eta_ref_ratio_min, h1_db_eta_ref_ratio_max] = GetHistMinMaxWithErrors(h_db_eta_ratio);
      h_db_eta_ref_min_yaxis = std::min(h_db_eta_ref_min_yaxis, h1_db_eta_ref_ratio_min*0.95);
      h_db_eta_ref_max_yaxis = std::max(h_db_eta_ref_max_yaxis, h1_db_eta_ref_ratio_max*1.05);
      h_db_eta_ref_dw -> GetYaxis() -> SetRangeUser(h_db_eta_ref_min_yaxis, h_db_eta_ref_max_yaxis);
      //--- Top pad
      cEta_db->cd(1);

      tdrDraw(p_db_eta_data_rebin, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
      tdrDraw(p_db_eta_mc_rebin, "HIST", kNone, mc_colors[i], kSolid, mc_colors[i], 0, 0);
      
      // Pass the rebinned and styled objects to the legend, not the original ones!
      leg_db_eta->AddEntry(p_db_eta_data_rebin, "Data", "pe");
      leg_db_eta->AddEntry(p_db_eta_mc_rebin, Form("MC   p^{tag}_{T} > %d GeV", (int)pt_cut), "le");
      
      //--- Bottom pad
      cEta_db->cd(2);
      tdrDraw(h_db_eta_ratio, "Pz", kFullCircle, colors[i], kSolid, colors[i], 0, 0);
  }

  cEta_db->SaveAs(Form("%s/%s/%d/JES_DB_vs_Eta_Slices_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run, run));
  
  
  // ==========================================================
  // Production of L2L3Res text file & Fit Loop over |eta| bins
  // ==========================================================
  gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%d/fits", outputBaseDirectory.c_str(),basePath.Data(), run));
  gROOT->ProcessLine(Form(".! touch %s/%s/%d/fits", outputBaseDirectory.c_str(),basePath.Data(), run));
  // Define JEC fit formula string
  //TString fit_formula = "[2]*([3]*([4]+TMath::Log(max([0],min([1],x)))*([5]+TMath::Log(max([0],min([1],x)))*[6])+[7]/x))*1./([8]+[9]/x+[10]*log(x)/x+[11]*(pow(x/[12],[13])-1)/(pow(x/[12],[13])+1)+[14]*pow(x,-0.3051)+[15]*x)";
  // [0]=Low pT Clamp, [1]=High pT Clamp, [2]=Scale p0, [3]=Log Slope p1, [4]=Inverse pT offset p2
  //TString fit_formula = "[2] + [3]*log(max([0],min([1],x))) + )";
  // [0]=Low pT Clamp, [1]=High pT Clamp, [2]=Scale p0, [3]=Log p1, [4]=Log^2 p2, [5]=Inverse pT p3
  TString fit_formula = "[2] + [3]*log(max([0],min([1],x))) + [4]*pow(log(max([0],min([1],x))), 2) + [5]/max([0],min([1],x))";
  
  // using now the JERC database convention string (link for example: https://github.com/cms-jet/JECDatabase/tree/master/textFiles/Summer24Prompt24_V1_MC )
  std::string txtPrefix = propertyTree.get<std::string>("global.txt_output_prefix", "JEC_Campaign");
  std::string txtVersion = propertyTree.get<std::string>("global.txt_output_version", "V1");
  std::string txt_filename = Form("./txts/%sRun%d_%s_DATA_L2L3Residual_AK4PFPuppi.txt", txtPrefix.c_str(), run, txtVersion.c_str());
  std::ofstream out_file(txt_filename);
  out_file << "{ 1 JetEta 1 JetPt " << fit_formula << " Correction L2Relative}\n"; // header in .txt file with the parameters JetEta, JetPt and the fit function.
  
  // Dedicated TGraphErrors to track chi2/ndf performance
  TGraphErrors *g_chi2ndf = new TGraphErrors();
  g_chi2ndf->SetName("g_chi2ndf");
  g_chi2ndf->SetTitle(Form("#chi^{2}/ndf of L2L3Res Fit (Run%d);Probe #eta;#chi^{2}/ndf", run));
  g_chi2ndf->SetMarkerStyle(kFullCircle);
  g_chi2ndf->SetMarkerColor(kBlue);
  g_chi2ndf->SetLineColor(kBlue);
  int iPointChi2 = 0;
  

  
  curdir->cd(); // Ensure objects are written to the current active directory

  // Loop over |eta| bins, rebinning data to keep uncertainties controlled

  for (int ix = 0; ix < nCustomEtaBins; ++ix) { 
    // ----------------------------------------------------------
    // AGGREGATE INPUTS FOR CORRESPONDING ETA BIN
    // ----------------------------------------------------------
    double eta_min = custom_eta_edges[ix];
    double eta_max = custom_eta_edges[ix+1];
    double eta_center = (eta_min + eta_max) / 2.0;
    
    // --- FIND MATCHING JSON CONFIG ---
    double abs_eta_center = std::abs(eta_center);
    double json_abs_eta_min = -1.0;
    double json_abs_eta_max = -1.0;
    int json_bin_index = -1; // Useful if you want to know if it's bin_01, bin_02, etc.

    for (size_t i = 0; i < etaConfigs.size(); ++i) {
        // Adding 1e-5 handles floating-point precision at the boundary edges
        if (abs_eta_center >= etaConfigs[i].abs_eta_min && abs_eta_center <= etaConfigs[i].abs_eta_max + 1e-5) {
            json_abs_eta_min = etaConfigs[i].abs_eta_min;
            json_abs_eta_max = etaConfigs[i].abs_eta_max;
            json_bin_index = i;
            break; // We found the matching block, stop searching
        }
    }
    
    // Now you have the exact absolute boundaries from the JSON!
    // std::cout << "Plotting " << eta_min << " to " << eta_max 
    //           << " | Matches JSON Bin Index: " << json_bin_index 
    //           << " (" << json_abs_eta_min << " - " << json_abs_eta_max << ")\n";



    // Find corresponding underlying X-axis bins in the TProfile2D
    // Adding/subtracting 1e-5 ensures we don't accidentally grab neighboring edge bins
    int bin_start = p2_MPF->GetXaxis()->FindBin(eta_min + 1e-5);
    int bin_end   = p2_MPF->GetXaxis()->FindBin(eta_max - 1e-5);
    
    std::vector<double> current_pt_bins = getPtBinning(eta_center);
    const Double_t* v_pt = current_pt_bins.data();
    const int n_pt = current_pt_bins.size() - 1;
    
    // Use ProfileY over the custom range: (name, start_bin, end_bin)
    TProfile *p_data = p2_MPF->ProfileY(Form("p_data_%d", ix), bin_start, bin_end);
    TProfile *p_data_rebin = (TProfile*)p_data->Rebin(n_pt, Form("p_data_rebin_%d", ix), v_pt);
    TH1D *h_data = p_data_rebin->ProjectionX(Form("h_data_%d", ix));

    TProfile *p_mc = p2_MPF_MC->ProfileY(Form("p_mc_%d", ix), bin_start, bin_end);
    TProfile *p_mc_rebin = (TProfile*)p_mc->Rebin(n_pt, Form("p_mc_rebin_%d", ix), v_pt);
    TH1D *h_mc = p_mc_rebin->ProjectionX(Form("h_mc_%d", ix));

    TH1D *h_ratio = (TH1D*)h_mc->Clone(Form("h_ratio_%d", ix));
    h_ratio->Divide(h_data);

    TProfile *p_db_data = p2_DB->ProfileY(Form("p_db_data_%d", ix), bin_start, bin_end);
    TProfile *p_db_data_rebin = (TProfile*)p_db_data->Rebin(n_pt, Form("p_db_data_rebin_%d", ix), v_pt);
    TH1D *h_db_data = p_db_data_rebin->ProjectionX(Form("h_db_data_%d", ix));

    // ----------------------------------------------------------
    // COORDINATE TRANSFORMATION: pT_tag -> pT_jet , 
    // by doing pT_jet = pT_tag*<pT_jet/pT_tag> 
    // where <pT_jet/pT_tag> is the direct balance (DB)
    // ----------------------------------------------------------
    TGraphErrors *g_ratio_vs_jetpt = new TGraphErrors();
    TGraphErrors *g_fit_graph = new TGraphErrors(); // New graph for points that are precise enough to used in the fit
    
    int n_points = 0;
    int n_fit_points = 0;
    
    double min_ratio = 9999.0, max_ratio = -9999.0;
    double pt_min_fit = 999999.0, pt_max_fit = -1.0; // To track dynamic pT clamps
    
    for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
        double pt_tag = h_ratio->GetBinCenter(b);
        double ratio = h_ratio->GetBinContent(b);
        double ratio_err = h_ratio->GetBinError(b);
        double data_err = h_data-> GetBinError(b);
        double r_data = h_db_data->GetBinContent(b);
        
        if (ratio > 0 && r_data > 0) {
            double pt_jet = pt_tag * r_data; 
            g_ratio_vs_jetpt->SetPoint(n_points, pt_jet, ratio);
            g_ratio_vs_jetpt->SetPointError(n_points, 0.0, ratio_err); 
            n_points++;
            
            if (ratio < min_ratio) min_ratio = ratio;
            if (ratio > max_ratio) max_ratio = ratio;

            // --- Apply Precision Filter (< precision tolerance ) ---
            if ((ratio_err / ratio) < precision_tolerance && data_err > 1e-5) {
                g_fit_graph->SetPoint(n_fit_points, pt_jet, ratio);
                g_fit_graph->SetPointError(n_fit_points, 0.0, ratio_err);
                
                // Track dynamic pT bounds of accepted points for the clamps
                if (pt_jet < pt_min_fit) pt_min_fit = pt_jet;
                if (pt_jet > pt_max_fit) pt_max_fit = pt_jet;
                n_fit_points++;
            }
        }
    }
    
    // --- Set Titles and Dynamic Range ---
    h_ratio->SetTitle(";Reco Jet p_{T} (GeV);JES MC/Data");
    
    // Safety check: only apply dynamic range if valid points were found
    if (min_ratio < 9999.0 && max_ratio > -9999.0) {
        h_ratio->GetYaxis()->SetRangeUser(0.8 * min_ratio, 1.2 * max_ratio);
    }

    TF1* func = new TF1(Form("fit_eta_%d", ix), fit_formula.Data(), xmin, xmax);
    
    // Apply exact lowest and highest used points as stability clamps
    if (n_fit_points >= 4) {
        func->FixParameter(0, pt_min_fit);
        func->FixParameter(1, pt_max_fit);

    } else { 
        // Fallback if bin is empty or statistically starved
        func->FixParameter(0, propertyTree.get<double>(chPath + ".fallback_clamp_min", 30.0));
        func->FixParameter(1, propertyTree.get<double>(chPath + ".fallback_clamp_max", 140.0)); 
    }
    
    // Initialize Quadratic Log-Inverse seeds
    func->SetParameter(2, 1.0); // p0 (Scale)
    func->SetParameter(3, 0.0); // p1 (Log slope)
    func->SetParameter(4, 0.0); // p2 (Log^2 slope)
    func->SetParameter(5, 0.0); // p3 (Inverse pT offset)

   // Uncertainty band
    TH1D *h_band = nullptr;

    // Perform fit against mapped TGraphErrors. Q=Quiet, R=Range, S=Save result.
    if (n_fit_points >= 4) { 
        TFitResultPtr r = g_fit_graph->Fit(func, "RQSN");
        
        // Fill Chi2/ndf monitoring graph
        if (func->GetNDF() > 0) {
            double chi2ndf = func->GetChisquare() / func->GetNDF();
            g_chi2ndf->SetPoint(iPointChi2, eta_center, chi2ndf);
            g_chi2ndf->SetPointError(iPointChi2, (eta_max - eta_min)/2.0, 0.0);
            iPointChi2++;
        }

        // Generate the uncertainty band from TFitResult 
        h_band = new TH1D(Form("h_band_%d", ix), "", 500, xmin, xmax);
        
        // Check if the fit result actually exists in memory
        if (r.Get()) {
            int n_bins = h_band->GetNbinsX();
            std::vector<double> x_vals(n_bins);
            std::vector<double> ci_vals(n_bins);
            
            // Get the X coordinates from our empty band histogram
            for (int i = 0; i < n_bins; ++i) {
                x_vals[i] = h_band->GetBinCenter(i + 1);
            }
            
            // Extract the 68% (1-sigma) intervals directly from the fit result pointer
            r->GetConfidenceIntervals(n_bins, 1, 1, x_vals.data(), ci_vals.data(), 0.68, false);
            
            // Populate the band histogram
            for (int i = 0; i < n_bins; ++i) {
                h_band->SetBinContent(i + 1, func->Eval(x_vals[i])); // Center the band on the fit curve
                h_band->SetBinError(i + 1, ci_vals[i]);              // Expand by the uncertainty error
            }
        }
        
        // Style the band
        h_band->SetStats(kFALSE);
        h_band->SetFillColor(kOrange); 
        //h_band->SetFillStyle(3001);
        h_band->SetMarkerSize(0); 
        h_band->SetLineWidth(0);  
    } else {
        std::cerr << Form("Warning: Not enough points for fit in bin %d (%d points)\n", ix, n_points);
    }
    
    // --- Plotting ---
    // Determine dynamic Y-axis bounds
    double y_min = (min_ratio < 9999.0) ? 0.8 * min_ratio : 0.8;
    double y_max = (max_ratio > -9999.0) ? 1.4 * max_ratio : 1.2;

    // Create an empty dummy histogram to enforce TDR fonts, axes, and titles
    TH1D *h_dummy = tdrHist(Form("h_dummy_%d", ix), "JES MC/Data", y_min, y_max, "Reco Jet p_{T} (GeV)", xmin, xmax);
    
    // Create the TDR canvas (Period 8 = 13.6 TeV, Pos 11 = Top-Left CMS label)
    TCanvas *cFit = tdrCanvas(Form("cFit_%d", ix), h_dummy, 8, 33, kSquare);
    gPad->SetLogx();

    // Draw a dashed line at exactly 1.0 for reference
    TLine *ll = new TLine(); 
    ll->SetLineStyle(kDashed); ll->SetLineColor(kBlack);
    ll->DrawLine(xmin, 1.0, xmax, 1.0);

    // Draw uncertainty band
    // "E3" draws a filled error band
    if (h_band) {
       h_band->Draw("E3 SAME");
    }

    // Draw data points using tdrDraw for consistent CMS styling
    // Tag-based (black open squares)
    tdrDraw(h_ratio, "P", kOpenSquare, kBlack, kSolid, -1, kNone, 0); 
    // Mapped Jet-based (black full circles)
    // Draw all original mapped points in open red markers to show what was rejected
    tdrDraw(g_ratio_vs_jetpt, "P", kOpenCircle, kRed, kSolid, -1, kNone, 0); 
    // Draw accepted mapped points (<20% unc) in solid black markers
    tdrDraw(g_fit_graph, "P", kFullCircle, kBlack, kSolid, -1, kNone, 0);

    // Draw the fit
    func->SetLineColor(kBlue); 
    func->SetLineWidth(3);
    func->Draw("L SAME"); 

    // Draw the TDR-styled Legend
    TLegend *lt = tdrLeg(0.20, 0.70, 0.50, 0.88); 
    lt->AddEntry(h_ratio, "Tag-based", "pe");
    lt->AddEntry(g_ratio_vs_jetpt, "Jet-based", "pe");
    lt->AddEntry(g_fit_graph, "Jet-based (for fit)", "pe");
    lt->AddEntry(func, "L2L3 Fit", "l");
    if (h_band) {
        lt->AddEntry(h_band, "Fit Unc. (1#sigma)", "f");
    }
    // Add standard CMS TLatex for the specific eta bin info
    TLatex *tex_eta = new TLatex();
    tex_eta->SetNDC();
    tex_eta->SetTextFont(42);
    tex_eta->SetTextSize(0.045);
    tex_eta->DrawLatex(0.20, 0.65, Form("%1.3f < #eta < %1.3f", eta_min, eta_max));

    // Save individual canvas
    TString plot_name = Form("%s/%s/%d/fits/L2L3Res_Fit_eta_%s_bin%d_%dp%d_to_%dp%d.png", outputBaseDirectory.c_str(), basePath.Data(), run, 
    eta_min < 0. ? "neg" : "pos",
    json_bin_index,
    int(floor(eta_min)), int(round((eta_min-floor(eta_min))*1000.)), 
    int(floor(eta_max)), int(round((eta_max-floor(eta_max))*1000.))
    );
    cFit->SaveAs(plot_name);

    // --- Write Structured Text File Payload ---
    int n_params = func->GetNpar();
    int n_tokens = n_params + 2; 
    
    out_file << Form("%9.6f %9.6f %d %9.0f %9.0f", eta_min, eta_max, n_tokens, xmin, xmax);
    for (int p = 0; p < n_params; ++p) {
        out_file << Form(" %12.6f", func->GetParameter(p));
    }
    out_file << "\n";
    
    // Cleanup within loop to control memory usage
    delete tex_eta; delete h_dummy; delete h_band; 
    delete lt; delete ll; delete func; delete g_ratio_vs_jetpt; delete g_fit_graph; delete h_ratio;
    delete h_mc; delete p_mc_rebin; delete p_mc; delete h_data; delete p_data_rebin; delete p_data;
    delete h_db_data; delete p_db_data_rebin; delete p_db_data;
    delete cFit;
  } // End of loop over eta bins

  // Close the text file and save the fit grid canvas
  out_file.close();
  std::cout << "Successfully generated corrections txt: " << txt_filename << std::endl;


  // --- Plot of chi2/ndf vs Eta ---
  // --- Final Plot: Chi2/ndf vs Eta ---
  TH1D *hDummyChi2 = tdrHist("hDummyChi2","#chi^{2}/ndf", -0.1, 5.0, "Probe #eta", -5.2, 5.2);
  TCanvas *cChi2 = tdrCanvas("cChi2", hDummyChi2, 8, 0, kRectangular); // tdrCanvas creates the canvas for you!
  g_chi2ndf->Draw("PZ SAME");
  TLine *lChi2 = new TLine(); lChi2->SetLineStyle(kDashed); lChi2->DrawLine(-5.191, 1.0, 5.191, 1.0);
  cChi2->SaveAs(Form("%s/%s/%d/fits/L2L3Res_Chi2_vs_Eta_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run, run));
  
  // CLEANUP SECTION ( to avoid memory leaks)
  // Delete monitoring stuff
  delete cChi2;
  delete g_chi2ndf;
  delete hDummyChi2;

  // Close and delete files
  f->Close();    delete f;
  fm->Close();   delete fm;
  fmOffline->Close(); delete fmOffline;
  
  // combine eta bins plots in grids
  try {
      gROOT->ProcessLine(Form(".! python3 combine_plots.py %s/%s/%d/fits --pattern L2L3Res_Fit_eta_pos* --name Grid_L2L3Res_Fit_eta_pos", outputBaseDirectory.c_str(),basePath.Data(), run));
      gROOT->ProcessLine(Form(".! python3 combine_plots.py %s/%s/%d/fits --pattern L2L3Res_Fit_eta_neg* --name Grid_L2L3Res_Fit_eta_neg", outputBaseDirectory.c_str(),basePath.Data(), run));
  } catch (const std::exception& e) {
      std::cerr << "Error creating grid images: " << e.what() << std::endl;
      return; 
  }
} // L2L3Res
