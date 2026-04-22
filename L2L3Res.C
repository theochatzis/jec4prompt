// Purpose: Read input JES data as TProfile2D and produce JEC L2L3Res fit
//          Focus on JEC stability and robustness
#include "TFile.h"
#include "TProfile2D.h"

#include "tdrstyle_mod22.C"

void L2L3Res() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/L2L3Res");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/L2L3Res");

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  // Load input file
  int run = 398801;
  //TFile *f = new TFile("rootfiles/J4PHists_runs392175to392175_photonjet.root","READ");
  TFile *f = new TFile("/eos/user/j/jecpcl/public/jec4prompt/runs/Run2025C/run392175/photonjet/J4PHists_runs392175to392175_photonjet.root","READ");
  
  assert(f && !f->IsZombie());
  TFile *fm = new TFile("rootfiles/reweighted_J4PHists_photonjet_GJ-4Jets.root","READ");
  TFile *fmOffline = new TFile("rootfiles/GamHistosFill_mc_summer2024P8_no-pu_w73.root","READ");
  //TFile *fm = new TFile("../jecsys3/rootfiles/Prompt2024/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root","READ");
  assert(fm && !fm->IsZombie());
  assert(fmOffline && !fmOffline->IsZombie());
  
  curdir->cd();

  // Define rebinning to ensure sufficient statistics
  // Original binning
  // {5, 7, 9, 11, 13, 15, 17, 20, 24, 28, 32, 36, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000}
  const Double_t v13[] =
    {20, 28, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153,
     174, 220, 300, 430, 638, 1032, 2000, 7000};
  const int n13 = sizeof(v13)/sizeof(v13[0])-1;
  
  // Load input data (MPF, DB) from file
  TProfile2D *p2m0 = (TProfile2D*)f->Get("DB_2D"); assert(p2m0);
  TProfile2D *p2m0m = (TProfile2D*)fm->Get("DB_2D"); assert(p2m0m);
  TProfile2D *p2m0mOffline = (TProfile2D*)fmOffline->Get("Gamjet2/p2m2"); assert(p2m0mOffline);
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
  TProfile *p1m0 = p2m0->ProfileY("p1m0",i1,i2);
  TProfile *p1m0_rebin = (TProfile*)p1m0->Rebin(n13,"p1m0_rebinned",v13);
  TH1D *h1m0 = p1m0_rebin->ProjectionX("h1m0");
  TH1D *h1m0_cut = (TH1D*)h1m0->Clone("h1m0_cut");
  h1m0_cut->GetXaxis()->SetRangeUser(40,300);

  TProfile *p1m0m = p2m0m->ProfileY("p1m0m",i1m,i2m);
  TProfile *p1m0m_rebin = (TProfile*)p1m0m->Rebin(n13,"p1m0m_rebinned",v13);
  //TProfile *p1m0m_rebin = (TProfile*)p1m0m->Clone("p1m0m_rebinned");
  TH1D *h1m0m = p1m0m_rebin->ProjectionX("h1m0m");
  TH1D *h1m0m_cut = (TH1D*)h1m0m->Clone("h1m0m_cut");
  h1m0m_cut->GetXaxis()->SetRangeUser(40,300);

  TProfile *p1m0mOffline = p2m0mOffline->ProfileY("p1m0mOffline",i1mOffline,i2mOffline);
  //TProfile *p1m0mOffline_rebin = (TProfile*)p1m0mOffline->Rebin(n13,"p1m0mOffline_rebinned",v13);
  //TProfile *p1m0mOffline_rebin = (TProfile*)p1m0mOffline->Clone("p1m0mOffline_rebinned");
  TH1D *h1m0mOffline = p1m0mOffline->ProjectionX("h1m0mOffline");
  TH1D *h1m0mOffline_cut = (TH1D*)h1m0mOffline->Clone("h1m0mOffline_cut");
  h1m0mOffline_cut->GetXaxis()->SetRangeUser(40,300);

  //TProfile *p1corrm = p2corrm->ProfileY("pcorrm",i1m,i2m);
  //TH1D *h1corrm = p1corrm->ProjectionX("h1corrm");
  //TH1D *h1m0m_corr = (TH1D*)h1m0m->Clone("h1m0m_corr");
  //h1m0m_corr->Divide(h1corrm);
  
  double xmin(15), xmax(4500);
  TH1D *h = tdrHist("h","JES",0.82-0.20,1.12+0.20,"p_{T,#gamma} (GeV)",15,4500);
  lumi_136TeV = Form("Run %d, 1 fb^{-1}",run);
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();
  //drawCustomLogXLabels(h);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(xmin,1,xmax,1);

  tdrDraw(p1m0m_rebin,"HIST",kNone,kBlue-9,kSolid,-1,kNone,0);
  tdrDraw(p1m0mOffline,"HIST",kNone,kGreen,kSolid,-1,kNone,0);
  //tdrDraw(h1m0m_corr,"HIST",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(p1m0,"Pz",kOpenSquare,kRed-9,kSolid,-1,kNone,0);
  tdrDraw(p1m0_rebin,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
  tdrDraw(h1m0_cut,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
  
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
  c1->SaveAs("pdf/L2L3Res_c1_Eta13.pdf");
  
  // Loop over |eta| bins, rebinning data to keep uncertainties controlled
  
  
  // Monitoring plots for data-fit


  // Monitoring plots for fit chi2 vs |eta|


  // Production of L2L3Res text file
  // (Later, add also mapping from pTtag to <pTprobe,raw> before this)
  
} // L2L3Res
