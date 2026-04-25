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
  std::string profileName = propertyTree.get<std::string>(chPath + ".profile_name");
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
  TProfile2D *p2m0 = (TProfile2D*)f->Get(profileName.c_str()); assert(p2m0);
  TProfile2D *p2m0m = (TProfile2D*)fm->Get(profileName.c_str()); assert(p2m0m);
  TProfile2D *p2m0mOffline = (TProfile2D*)fmOffline->Get("Gamjet2/p2m0"); assert(p2m0mOffline);
  //TProfile2D *p2corrm = (TProfile2D*)fm->Get("Gamjet2/p2corr"); assert(p2corrm);

  // Initial fit at |eta|<1.305 needed for scaling dijet data to L2L3Res level
  int i1 = p2m0->GetXaxis()->FindBin(-1.305);
  int i2 = p2m0->GetXaxis()->FindBin(+1.305)-1;
  double eta1 = p2m0->GetXaxis()->GetBinLowEdge(i1);
  double eta2 = p2m0->GetXaxis()->GetBinLowEdge(i2)+1;

  int i1m = p2m0m->GetXaxis()->FindBin(-1.305);
  int i2m = p2m0m->GetXaxis()->FindBin(+1.305)-1;
  double eta1m = p2m0m->GetXaxis()->GetBinLowEdge(i1m);
  double eta2m = p2m0m->GetXaxis()->GetBinLowEdge(i2m)+1;

  int i1mOffline = p2m0mOffline->GetXaxis()->FindBin(-1.305);
  int i2mOffline = p2m0mOffline->GetXaxis()->FindBin(+1.305)-1;
  double eta1mOffline = p2m0mOffline->GetXaxis()->GetBinLowEdge(i1m);
  double eta2mOffline = p2m0mOffline->GetXaxis()->GetBinLowEdge(i2m)+1;

  cout << Form("Fitting %1.3f #LT eta < %1.3f reference region\n",eta1,eta2);

  // Get reference binning (central region, eta ~ 0.0)
  std::vector<double> v_pt_bins_ref = getPtBinning(0.0);
  const Double_t* v13_ref = v_pt_bins_ref.data();
  const int n13_ref = v_pt_bins_ref.size() - 1;

  float fit_region_min  = 40;
  float fit_region_max = 300;

  TProfile *p1m0 = p2m0->ProfileY("p1m0",i1,i2);
  TProfile *p1m0_rebin = (TProfile*)p1m0->Rebin(n13_ref,"p1m0_rebinned",v13_ref);
  TH1D *h1m0 = p1m0_rebin->ProjectionX("h1m0");
  TH1D *h1m0_cut = (TH1D*)h1m0->Clone("h1m0_cut");
  h1m0_cut->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);

  TProfile *p1m0m = p2m0m->ProfileY("p1m0m",i1m,i2m);
  TProfile *p1m0m_rebin = (TProfile*)p1m0m->Rebin(n13_ref,"p1m0m_rebinned",v13_ref);
  //TProfile *p1m0m_rebin = (TProfile*)p1m0m->Clone("p1m0m_rebinned");
  TH1D *h1m0m = p1m0m_rebin->ProjectionX("h1m0m");
  TH1D *h1m0m_cut = (TH1D*)h1m0m->Clone("h1m0m_cut");
  h1m0m_cut->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);

  TProfile *p1m0mOffline = p2m0mOffline->ProfileY("p1m0mOffline",i1mOffline,i2mOffline);
  //TProfile *p1m0mOffline_rebin = (TProfile*)p1m0mOffline->Rebin(n13_ref,"p1m0mOffline_rebinned",v13_ref);
  //TProfile *p1m0mOffline_rebin = (TProfile*)p1m0mOffline->Clone("p1m0mOffline_rebinned");
  TH1D *h1m0mOffline = p1m0mOffline->ProjectionX("h1m0mOffline");
  TH1D *h1m0mOffline_cut = (TH1D*)h1m0mOffline->Clone("h1m0mOffline_cut");
  h1m0mOffline_cut->GetXaxis()->SetRangeUser(fit_region_min, fit_region_max);

  //TProfile *p1corrm = p2corrm->ProfileY("pcorrm",i1m,i2m);
  //TH1D *h1corrm = p1corrm->ProjectionX("h1corrm");
  //TH1D *h1m0m_corr = (TH1D*)h1m0m->Clone("h1m0m_corr");
  //h1m0m_corr->Divide(h1corrm);
  
  TH1D *h = tdrHist("h","JES",jes_limitMin,jes_limitMax,"p_{T,#gamma} (GeV)",xmin,xmax);
  lumi_136TeV = Form("Run %d, %.3f fb^{-1}", run, luminosity);
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();
  //drawCustomLogXLabels(h);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(xmin,1,xmax,1);

  tdrDraw(p1m0m_rebin,"HIST",kNone,colorMC,kSolid,-1,kNone,0);
  tdrDraw(p1m0mOffline,"HIST",kNone,kGreen,kSolid,-1,kNone,0);
  //tdrDraw(h1m0m_corr,"HIST",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(p1m0,"Pz",kOpenSquare,kRed-9,kSolid,-1,kNone,0);
  tdrDraw(p1m0_rebin,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
  tdrDraw(h1m0_cut,"Pz",kFullCircle,colorData,kSolid,-1,kNone,0);
  
  // Make legends
  TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.90);
  leg->SetFillStyle(0);   // Transparent background
  leg->SetBorderSize(0);  // No border
  leg->SetTextFont(42);   // Standard CMS TDR font
  leg->SetTextSize(0.04);

  leg->AddEntry(h1m0_cut, "Data (fit region)", "pe");
  leg->AddEntry(p1m0_rebin, "Data (rebinned)", "pe");
  leg->AddEntry(p1m0, "Data (unrebinned)", "pe");
  leg->AddEntry(p1m0mOffline, "MC Offline", "le");
  leg->AddEntry(p1m0m_rebin, "MC Online", "le");
  
  leg->Draw();

  // Save to pdf
  c1->SaveAs(Form("%s/%s/L2L3Res_Eta13_run%d.png", outputBaseDirectory.c_str(), basePath.Data(), run));
  
  // ==========================================================
  // Production of L2L3Res text file & Fit Loop over |eta| bins
  // ==========================================================
  // Define JEC fit formula string
  TString fit_formula = "[2]*([3]*([4]+TMath::Log(max([0],min([1],x)))*([5]+TMath::Log(max([0],min([1],x)))*[6])+[7]/x))*1./([8]+[9]/x+[10]*log(x)/x+[11]*(pow(x/[12],[13])-1)/(pow(x/[12],[13])+1)+[14]*pow(x,-0.3051)+[15]*x)";
  
  // Loop over |eta| bins, rebinning data to keep uncertainties controlled
  
  
  // Monitoring plots for data-fit


  // Monitoring plots for fit chi2 vs |eta|


  // Production of L2L3Res text file
  // (Later, add also mapping from pTtag to <pTprobe,raw> before this)
  
  // CLEANUP SECTION ( to avoid memory leaks)
  // 1. Delete the canvas to free up graphics memory
  delete c1;
  delete h; 
  delete l;

  // 2. Delete the histograms/profiles created
  delete p1m0;
  delete p1m0_rebin;
  // delete other projections/clones...

  // 3. Close and delete files
  f->Close();    delete f;
  fm->Close();   delete fm;
  fmOffline->Close(); delete fmOffline;

} // L2L3Res