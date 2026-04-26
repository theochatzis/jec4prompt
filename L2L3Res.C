// Purpose: Read input JES data as TProfile2D and produce JEC L2L3Res fit
//          Focus on JEC stability and robustness
#include "TFile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"

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

void L2L3Res(int run = 398600, TString basePath="2025G", TString channel="photonjet") {
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

  double jes_limitMin = propertyTree.get<double>("global.jes_limitMin", 0.82-0.20);
  double jes_limitMax = propertyTree.get<double>("global.jes_limitMax", 1.12+0.20);
  
  
  
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
  // -------------------------------
  
  gROOT->ProcessLine(Form(".! mkdir -p %s/%s", outputBaseDirectory.c_str(),basePath.Data()));
  gROOT->ProcessLine(Form(".! touch %s/%s", outputBaseDirectory.c_str(),basePath.Data()));
  //gROOT->ProcessLine(Form(".! mkdir %s/L2L3Res", basePath.Data());
  //gROOT->ProcessLine(Form(".! touch %s/L2L3Res", basePath.Data()));

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Load input file, using the channel variable for the path as well
  TString dataPath = Form("/eos/user/j/jecpcl/public/jec4prompt/runs/Run%s/run%d/%s/J4PHists_runs%dto%d_%s.root", 
                          basePath.Data(), run, channel.Data(), run, run, channel.Data());
  
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

  TFile *fm = new TFile("rootfiles/reweighted_J4PHists_photonjet_GJ-4Jets.root","READ");
  TFile *fmOffline = new TFile("rootfiles/GamHistosFill_mc_summer2024P8_no-pu_w73.root","READ");
  //TFile *fm = new TFile("../jecsys3/rootfiles/Prompt2024/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root","READ");
  assert(fm && !fm->IsZombie());
  assert(fmOffline && !fmOffline->IsZombie());
  
  curdir->cd();
  
  // Load input data (MPF, DB) from file
  TProfile2D *p2_MPF = (TProfile2D*)f->Get(profileName.c_str()); assert(p2_MPF);
  TProfile2D *p2_MPF_MC = (TProfile2D*)fm->Get(profileName.c_str()); assert(p2_MPF_MC);
  TProfile2D *p2_MPFOffline_MC = (TProfile2D*)fmOffline->Get("Gamjet2/p2m0"); assert(p2_MPFOffline_MC);
  //TProfile2D *p2corrm = (TProfile2D*)fm->Get("Gamjet2/p2corr"); assert(p2corrm);
  TProfile2D *p2_DB = (TProfile2D*)f->Get(profileNameDB.c_str()); assert(p2_DB);

  // Initial fit at |eta|<1.305 needed for scaling dijet data to L2L3Res level
  int i1 = p2_MPF->GetXaxis()->FindBin(-1.305);
  int i2 = p2_MPF->GetXaxis()->FindBin(+1.305)-1;
  double eta1 = p2_MPF->GetXaxis()->GetBinLowEdge(i1);
  double eta2 = p2_MPF->GetXaxis()->GetBinLowEdge(i2)+1;

  int i1m = p2_MPF_MC->GetXaxis()->FindBin(-1.305);
  int i2m = p2_MPF_MC->GetXaxis()->FindBin(+1.305)-1;
  double eta1m = p2_MPF_MC->GetXaxis()->GetBinLowEdge(i1m);
  double eta2m = p2_MPF_MC->GetXaxis()->GetBinLowEdge(i2m)+1;

  int i1mOffline = p2_MPFOffline_MC->GetXaxis()->FindBin(-1.305);
  int i2mOffline = p2_MPFOffline_MC->GetXaxis()->FindBin(+1.305)-1;
  double eta1mOffline = p2_MPFOffline_MC->GetXaxis()->GetBinLowEdge(i1m);
  double eta2mOffline = p2_MPFOffline_MC->GetXaxis()->GetBinLowEdge(i2m)+1;

  cout << Form("Fitting %1.3f #LT eta < %1.3f reference region\n",eta1,eta2);
  
//   std::string txt_filename = Form("%s/%s/Summer22EE-22Sep2023_Run%d_%s_DATA_L2L3Residual_AK4PFPuppi.txt", outputBaseDirectory.c_str(), basePath.Data(), run, basePath.Data());
//   std::ofstream out_file(txt_filename);
//   out_file << "{ 1 JetEta 1 JetPt " << fit_formula << " Correction L2Relative}\n";
  
  // Get reference binning (central region, eta ~ 0.0)
  std::vector<double> v_pt_bins_ref = getPtBinning(0.0);
  const Double_t* v13_ref = v_pt_bins_ref.data();
  const int n13_ref = v_pt_bins_ref.size() - 1;

  float fit_region_min  = 40;
  float fit_region_max = 300;

  TProfile *p1_MPF = p2_MPF->ProfileY("p1_MPF",i1,i2);
  TProfile *p1_MPF_rebin = (TProfile*)p1_MPF->Rebin(n13_ref,"p1_MPF_rebinned",v13_ref);
  TH1D *h1_MPF = p1_MPF_rebin->ProjectionX("h1_MPF");
  TH1D *h1_MPF_cut = (TH1D*)h1_MPF->Clone("h1_MPF_cut");
  h1_MPF_cut->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);

  TProfile *p1_MPF_MC = p2_MPF_MC->ProfileY("p1_MPF_MC",i1m,i2m);
  TProfile *p1_MPF_MC_rebin = (TProfile*)p1_MPF_MC->Rebin(n13_ref,"p1_MPF_MC_rebinned",v13_ref);
  //TProfile *p1_MPF_MC_rebin = (TProfile*)p1_MPF_MC->Clone("p1_MPF_MC_rebinned");
  TH1D *h1_MPF_MC = p1_MPF_MC_rebin->ProjectionX("h1_MPF_MC");
  TH1D *h1_MPF_MC_cut = (TH1D*)h1_MPF_MC->Clone("h1_MPF_MC_cut");
  h1_MPF_MC_cut->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);

  TProfile *p1_MPFOffline_MC = p2_MPFOffline_MC->ProfileY("p1_MPFOffline_MC",i1mOffline,i2mOffline);
  //TProfile *p1_MPFOffline_MC_rebin = (TProfile*)p1_MPFOffline_MC->Rebin(n13_ref,"p1_MPFOffline_MC_rebinned",v13_ref);
  //TProfile *p1_MPFOffline_MC_rebin = (TProfile*)p1_MPFOffline_MC->Clone("p1_MPFOffline_MC_rebinned");
  TH1D *h1_MPFOffline_MC = p1_MPFOffline_MC->ProjectionX("h1_MPFOffline_MC");
  TH1D *h1_MPFOffline_MC_cut = (TH1D*)h1_MPFOffline_MC->Clone("h1_MPFOffline_MC_cut");
  h1_MPFOffline_MC_cut->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);

  // Un-do the jec correction in offline
  //TProfile *p1corrm = p2corrm->ProfileY("pcorrm",i1m,i2m);
  //TH1D *h1corrm = p1corrm->ProjectionX("h1corrm");
  //TH1D *h1_MPFOffline_MC_corr = (TH1D*)h1_MPFOffline_MC->Clone("h1_MPFOffline_MC_corr");
  //h1_MPFOffline_MC_corr->Divide(h1corrm);
  
  TH1D *h = tdrHist("h","JES",jes_limitMin,jes_limitMax,"p_{T,#gamma} (GeV)",xmin,xmax);
  lumi_136TeV = Form("Run %d, %.3f fb^{-1}", run, luminosity);
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();
  //drawCustomLogXLabels(h);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(xmin,1,xmax,1);

  tdrDraw(p1_MPF_MC_rebin,"HIST",kNone,colorMC,kSolid,-1,kNone,0);
  tdrDraw(p1_MPFOffline_MC,"HIST",kNone,kGreen,kSolid,-1,kNone,0);
  //tdrDraw(h1_MPFOffline_MC_corr,"HIST",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(p1_MPF,"Pz",kOpenSquare,kRed-9,kSolid,-1,kNone,0);
  tdrDraw(p1_MPF_rebin,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
  tdrDraw(h1_MPF_cut,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0);
  
  // Make legends
  TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.90);
  leg->SetFillStyle(0);   // Transparent background
  leg->SetBorderSize(0);  // No border
  leg->SetTextFont(42);   // Standard CMS TDR font
  leg->SetTextSize(0.04);

  leg->AddEntry(h1_MPF_cut, "Data (fit region)", "pe");
  leg->AddEntry(p1_MPF_rebin, "Data (rebinned)", "pe");
  leg->AddEntry(p1_MPF, "Data (unrebinned)", "pe");
  leg->AddEntry(p1_MPFOffline_MC, "MC Offline", "le");
  leg->AddEntry(p1_MPF_MC_rebin, "MC Online", "le");
  
  leg->Draw();

  // Save to pdf
  c1->SaveAs(Form("%s/%s/L2L3Res_Eta13_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run));
  
  // ==========================================================
  // Production of L2L3Res text file & Fit Loop over |eta| bins
  // ==========================================================
  // Define JEC fit formula string
  
  TString fit_formula = "[2]*([3]*([4]+TMath::Log(max([0],min([1],x)))*([5]+TMath::Log(max([0],min([1],x)))*[6])+[7]/x))*1./([8]+[9]/x+[10]*log(x)/x+[11]*(pow(x/[12],[13])-1)/(pow(x/[12],[13])+1)+[14]*pow(x,-0.3051)+[15]*x)";
  
  // using now the JERC database convention string (link for example: https://github.com/cms-jet/JECDatabase/tree/master/textFiles/Summer24Prompt24_V1_MC )
  std::string txt_filename = Form("./txts/Summer24Prompt24JEC4PromptRun%d_V1_DATA_L2L3Residual_AK4PFPuppi.txt", run);//outputBaseDirectory.c_str(), basePath.Data(), run);
  std::ofstream out_file(txt_filename);
  out_file << "{ 1 JetEta 1 JetPt " << fit_formula << " Correction L2Relative}\n"; // header in .txt file with the parameters JetEta, JetPt and the fit function.
  
  int nEtaBins = p2_MPF->GetXaxis()->GetNbins();

//   // --- Monitoring Plot Setup ---
//   // Create multi-canvas grid for monitoring plots (square, fits roughly 4x4 or 5x5 grid based on typical nEtaBins)
//   TCanvas *cGrid = new TCanvas("cGrid", "Fits Grid", 2000, 2000);
//   int gridDim = std::ceil(std::sqrt(nEtaBins));
//   cGrid->Divide(gridDim, gridDim);
  // --- PAGINATED PLOT SETUP ---
  int maxPadsPerCanvas = 16; // 4x4 grid
  int gridX = 4;
  int gridY = 4;
  int iCanvas = 1;
  int iPad = 1;
  
  TCanvas *cGrid = nullptr;
  std::vector<TObject*> objectsToClean; // To hold drawn objects safely until canvas saves
  
  // Dedicated TGraphErrors to track chi2/ndf performance
  TGraphErrors *g_chi2ndf = new TGraphErrors();
  g_chi2ndf->SetName("g_chi2ndf");
  g_chi2ndf->SetTitle("#chi^{2}/ndf of L2L3Res Fit;Probe #eta;#chi^{2}/ndf");
  g_chi2ndf->SetMarkerStyle(kFullCircle);
  g_chi2ndf->SetMarkerColor(kBlue);
  g_chi2ndf->SetLineColor(kBlue);
  int iPointChi2 = 0;

  curdir->cd(); // Ensure objects are written to the current active directory

  // Loop over |eta| bins, rebinning data to keep uncertainties controlled
  for (int ix = 1; ix <= nEtaBins; ++ix) { 
    // // Select correct multicanvas grid panel
    // cGrid->cd(ix);
    // gPad->SetLogx();
    // gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
    // --- CREATE NEW CANVAS IF NEEDED ---
    if ((ix - 1) % maxPadsPerCanvas == 0) {
        cGrid = new TCanvas(Form("cGrid_part%d", iCanvas), Form("Fits Grid Part %d", iCanvas), 2000, 2000);
        cGrid->Divide(gridX, gridY);
        iPad = 1; // Reset pad counter for the new canvas
    }

    cGrid->cd(iPad);
    gPad->SetLogx();
    gPad->SetRightMargin(0.05); gPad->SetTopMargin(0.10);
    // ----------------------------------------------------------
    // CREATION OF INPUTS PER ETA BIN
    // ----------------------------------------------------------
    double eta_min = p2_MPF->GetXaxis()->GetBinLowEdge(ix);
    double eta_max = p2_MPF->GetXaxis()->GetBinUpEdge(ix);
    double eta_center = p2_MPF->GetXaxis()->GetBinCenter(ix);
    
    // Dynamically fetch the correct pt binning for this specific eta slice
    std::vector<double> current_pt_bins = getPtBinning(eta_center);
    const Double_t* v_pt = current_pt_bins.data();
    const int n_pt = current_pt_bins.size() - 1;
    
    // Extract Profiles, rebin using the dynamically fetched array and project
    TProfile *p_data = p2_MPF->ProfileY(Form("p_data_%d", ix), ix, ix);
    TProfile *p_data_rebin = (TProfile*)p_data->Rebin(n_pt, Form("p_data_rebin_%d", ix), v_pt);
    TH1D *h_data = p_data_rebin->ProjectionX(Form("h_data_%d", ix));

    TProfile *p_mc = p2_MPF_MC->ProfileY(Form("p_mc_%d", ix), ix, ix);
    TProfile *p_mc_rebin = (TProfile*)p_mc->Rebin(n_pt, Form("p_mc_rebin_%d", ix), v_pt);
    TH1D *h_mc = p_mc_rebin->ProjectionX(Form("h_mc_%d", ix));

    // TH1D *h_ratio = (TH1D*)h_data->Clone(Form("h_ratio_%d", ix));
    // h_ratio->Divide(h_mc);
    TH1D *h_ratio = (TH1D*)h_mc->Clone(Form("h_ratio_%d", ix));
    h_ratio->Divide(h_data);

    // Take also Direct Balance in order to make the re-binning to reco Jet pT on x-axis
    TProfile *p_db_data = p2_DB->ProfileY(Form("p_db_data_%d", ix), ix, ix);
    TProfile *p_db_data_rebin = (TProfile*)p_db_data->Rebin(n_pt, Form("p_db_data_rebin_%d", ix), v_pt);
    TH1D *h_db_data = p_db_data_rebin->ProjectionX(Form("h_db_data_%d", ix));

    // ----------------------------------------------------------
    // COORDINATE TRANSFORMATION: pT_tag -> pT_jet , 
    // by doing pT_jet = pT_tag*<pT_jet/pT_tag> 
    // where <pT_jet/pT_tag> is the direct balance (DB)
    // ----------------------------------------------------------
    TGraphErrors *g_ratio_vs_jetpt = new TGraphErrors();
    int n_points = 0;
    
    for (int b = 1; b <= h_ratio->GetNbinsX(); ++b) {
        double pt_tag = h_ratio->GetBinCenter(b);
        double ratio = h_ratio->GetBinContent(b);
        double ratio_err = h_ratio->GetBinError(b);
        double r_data = h_db_data->GetBinContent(b);
        
        if (ratio > 0 && r_data > 0) {
            double pt_jet = pt_tag * r_data; 
            g_ratio_vs_jetpt->SetPoint(n_points, pt_jet, ratio);
            g_ratio_vs_jetpt->SetPointError(n_points, 0.0, ratio_err); 
            n_points++;
        }
    }

      TF1* func = new TF1(Form("fit_eta_%d", ix), fit_formula.Data(), xmin, xmax);
      
      // Fix stability clamps [0] and [1] globally
      func->FixParameter(0, 30.0);
      
      // Dynamic clamp for p1 (upper logarithmic stability anchor)
      double p1_clamp = 140.0;
      double abs_eta = std::abs(eta_center);
      if (abs_eta <= 1.305) p1_clamp = 1290.0;
      else if (abs_eta <= 2.5) p1_clamp = 700.0;
      else if (abs_eta <= 3.0) p1_clamp = 330.0;
      else p1_clamp = 140.0;
      func->FixParameter(1, p1_clamp);
      
      // Initial seeds for L2L3 residual parameters
      func->SetParameter(2, 1.0); // Global scale
      func->SetParameter(3, 1.0); // Secondary scale 
      func->SetParameter(4, 1.0); // Log offset
      func->SetParameter(5, 0.0);
      func->SetParameter(6, 0.0);
      func->SetParameter(7, 0.0);
      
      // Denominator anchored to standard Run 3 MC parameterization (from standard constants)
      func->SetParameter(8, 0.9461);
      func->SetParameter(9, 0.6498);
      func->SetParameter(10, 0.06043);
      func->SetParameter(11, 0.05852);
      func->SetParameter(12, 221.98);
      func->SetParameter(13, 0.8940);
      func->SetParameter(14, -0.17906);
      func->SetParameter(15, -0.00002410);

      // Perform fit against mapped TGraphErrors. Q=Quiet, R=Range.
      if (n_points >= 4) { // Require minimum points for a robust multiparameter fit
          g_ratio_vs_jetpt->Fit(func, "RQ"); 
          
          // Fill Chi2/ndf monitoring graph
          if (func->GetNDF() > 0) {
              double chi2ndf = func->GetChisquare() / func->GetNDF();
              g_chi2ndf->SetPoint(iPointChi2, eta_center, chi2ndf);
              g_chi2ndf->SetPointError(iPointChi2, (eta_max - eta_min)/2.0, 0.0);
              iPointChi2++;
          }
      } else {
          std::cerr << Form("Warning: Not enough points for fit in bin %d (%d points)\n", ix, n_points);
      }
      
      // --- Formatting and Plotting ---
      h_ratio->Draw("AXIS"); // Draw axes first
      //h_ratio->Draw("P SAME"); // Draw original tag-based ratio
      g_ratio_vs_jetpt->Draw("P SAME"); // Draw mapped jet-based ratio
      func->SetLineColor(kBlue); func->SetLineWidth(3);
      func->Draw("L SAME"); // Draw fit curve
      
    //   TLegend *lt = new TLegend(0.15, 0.70, 0.45, 0.88); lt->SetFillStyle(0); lt->SetBorderSize(0); lt->SetTextSize(0.04);
    //   lt->AddEntry(h_ratio, "Tag-based", "pe");
    //   lt->AddEntry(g_ratio_vs_jetpt, "Mapped Jet-based", "pe");
    //   lt->AddEntry(func, "L2L3 Fit", "l");
    //   lt->Draw();
      TLegend *lt = new TLegend(0.15, 0.70, 0.45, 0.88); lt->SetFillStyle(0); lt->SetBorderSize(0); lt->SetTextSize(0.04);
      lt->AddEntry(h_ratio, Form("%1.3f < #eta < %1.3f", eta_min, eta_max), "");
      //lt->AddEntry(h_ratio, "Tag-based", "pe");
      lt->AddEntry(g_ratio_vs_jetpt, "MPF", "pe");
      lt->AddEntry(func, "L2L3 Fit", "l");
      lt->Draw();

      TLine *ll = new TLine(); ll->SetLineStyle(kDashed); ll->DrawLine(gPad->GetUxmin(), 1.0, gPad->GetUxmax(), 1.0);

      // --- Write Structured Text File Payload ---
      int n_params = func->GetNpar();
      int n_tokens = n_params + 2; 
      
      out_file << Form("%9.6f %9.6f %d %9.0f %9.0f", eta_min, eta_max, n_tokens, xmin, xmax);
      for (int p = 0; p < n_params; ++p) {
          out_file << Form(" %12.6f", func->GetParameter(p));
      }
      out_file << "\n";
      
      // Cleanup within loop to control memory usage
      // CRITICAL: Do NOT delete drawn objects (lt, ll, func, g_ratio_vs_jetpt, h_ratio)
      // here because cGrid needs them to be alive when SaveAs() is called!
      //   delete lt;
      //   delete ll;
      //   delete func;
      //   delete g_ratio_vs_jetpt;
      //   delete h_ratio;
      delete h_mc;
      delete p_mc_rebin;
      delete p_mc;
      delete h_data;
      delete p_data_rebin;
      delete p_data;
      delete h_db_data;
      delete p_db_data_rebin;
      delete p_db_data;

      // 3. IF canvas is full OR it's the last bin: Save, clean drawn objects, and increment canvas counter
      if (ix % maxPadsPerCanvas == 0 || ix == nEtaBins) {
          cGrid->SaveAs(Form("%s/%s/L2L3Res_Fits_Grid_run%d_part%d.png", outputBaseDirectory.c_str(), basePath.Data(), run, iCanvas));          
          delete cGrid;
          cGrid = nullptr;
          iCanvas++;
      }
      
      iPad++; // Increment pad for next iteration
  } // End of loop over eta bins

  // Close the text file and save the fit grid canvas
  out_file.close();
  std::cout << "Successfully generated corrections txt: " << txt_filename << std::endl;
//   cGrid->SaveAs(Form("%s/%s/L2L3Res_Fits_Grid_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run));

  // --- Final Plot: Chi2/ndf vs Eta ---
//   TCanvas *cChi2 = new TCanvas("cChi2", "Chi2 Monitoring", 1200, 800);
//   TH1D *hDummyChi2 = tdrHist("hDummyChi2","#chi^{2}/ndf", -0.1, 5.0, "Probe #eta", -5.2, 5.2);
//   tdrCanvas("cChi2", hDummyChi2, 8, 11, kSquare);
//   g_chi2ndf->Draw("PZ SAME");
//   TLine *lChi2 = new TLine(); lChi2->SetLineStyle(kDashed); lChi2->DrawLine(-5.191, 1.0, 5.191, 1.0);
//   cChi2->SaveAs(Form("%s/%s/L2L3Res_Chi2_vs_Eta_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run));


  
  // CLEANUP SECTION ( to avoid memory leaks)
  // Delete monitoring stuff
  //delete cGrid;
  //delete cChi2;
  //delete g_chi2ndf;
  //delete hDummyChi2;

  // Delete the canvas
  delete c1;
  delete h; 
  delete l;

  // elete the histograms/profiles created
  delete p1_MPF;
  delete p1_MPF_rebin;
  // delete other projections/clones...

  // Close and delete files
  f->Close();    delete f;
  fm->Close();   delete fm;
  fmOffline->Close(); delete fmOffline;

} // L2L3Res