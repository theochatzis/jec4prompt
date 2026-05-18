import ROOT
import os
import json
import math
import array

ROOT.TGaxis.SetMaxDigits(3)
ROOT.gErrorIgnoreLevel = ROOT.kError # Suppress Sumw2 warnings

def generate_db_dummy_data(is_mc=False, mc_index=1):
    """Generates dummy DB data to test against."""
    name = f"DB_2D_Dummy_{mc_index}" if is_mc else "DB_2D_Data_Dummy"
    h2 = ROOT.TProfile2D(name, f"{name} vs Probe #eta and Tag p_{{T}};Probe #eta;Tag p_{{T}} (GeV);DB Response", 
                         100, -5.2, 5.2, 2000, 10.0, 2000.0)  
    h2.SetDirectory(0)
    
    true_N = 2.5 if is_mc else 3.0
    true_S = 0.8 if is_mc else 1.0
    true_C = 0.03 + (mc_index * 0.01) if is_mc else 0.04

    for _ in range(5000000):
        eta = ROOT.gRandom.Uniform(-5.2, 5.2)
        pt = 55.0 + ROOT.gRandom.Exp(50.0) if ROOT.gRandom.Uniform() < 0.2 else 55.0 / math.pow(1.0 - ROOT.gRandom.Uniform(), 1.0/2.5)
        if pt > 3000: continue
        
        res = math.sqrt((true_N/pt)**2 + (true_S/math.sqrt(pt))**2 + true_C**2)
        h2.Fill(eta, pt, ROOT.gRandom.Gaus(0.95, 0.95 * res))
    return h2

def extract_jer_profile(h2_db, name_suffix, eta_min, eta_max, pt_bins):
    """
    Pure mathematical extraction function. 
    Takes a TProfile2D, folds negative eta, rebins, and computes sigma/mean.
    """
    b_low = h2_db.GetXaxis().FindBin(eta_min + 1e-5)
    b_high = h2_db.GetXaxis().FindBin(eta_max - 1e-5)
    
    h_prof_raw = h2_db.ProfileY(f"h_prof_raw_{name_suffix}", b_low, b_high)
    h_prof_raw.SetDirectory(0)
    
    # Fold negative eta
    if h2_db.GetXaxis().GetXmin() < 0:
        b_neg_low = h2_db.GetXaxis().FindBin(-eta_max + 1e-5)
        b_neg_high = h2_db.GetXaxis().FindBin(-eta_min - 1e-5)
        h_neg = h2_db.ProfileY(f"h_neg_{name_suffix}", b_neg_low, b_neg_high)
        h_neg.SetDirectory(0)
        h_prof_raw.Add(h_neg)
        h_neg.Delete()

    # Rebin
    if pt_bins:
        bin_array = array.array('d', pt_bins)
        h_prof = h_prof_raw.Rebin(len(pt_bins) - 1, f"h_prof_{name_suffix}", bin_array)
        h_prof.SetDirectory(0)
        h_prof_raw.Delete()
    else:
        h_prof = h_prof_raw
        h_prof.SetName(f"h_prof_{name_suffix}")

    # Create empty resolution histogram matching the profile's binning
    xbins = h_prof.GetXaxis().GetXbins()
    title = f"Extracted JER |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Tag p_{{T}} (GeV); JER(MPFx) "
    if xbins.GetSize() > 0: 
        h_res = ROOT.TH1D(f"h_res_{name_suffix}", title, h_prof.GetNbinsX(), xbins.GetArray())
    else: 
        h_res = ROOT.TH1D(f"h_res_{name_suffix}", title, h_prof.GetNbinsX(), h_prof.GetXaxis().GetXmin(), h_prof.GetXaxis().GetXmax())
    h_res.SetDirectory(0)

    # Compute sigma / mean
    for iy in range(1, h_prof.GetNbinsX() + 1):
        entries = h_prof.GetBinEntries(iy)
        err = h_prof.GetBinError(iy)
        mean = h_prof.GetBinContent(iy)
        
        if entries > 10 and mean > 0.0: 
            sigma = err * math.sqrt(entries)
            rel_res = sigma / mean
            err_rel_res = (sigma / math.sqrt(2 * entries)) / mean
            
            h_res.SetBinContent(iy, rel_res)
            h_res.SetBinError(iy, err_rel_res) 
            
    h_prof.Delete()
    return h_res

def process_and_plot_jers(input_sources, config_list, out_dir="./jer_fits"):
    """
    Loops over eta bins, extracts JERs for all inputs, fits the first one (baseline), 
    and plots them all together.
    """
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    jer_results = []

    for i, bin_config in enumerate(config_list):
        eta_min, eta_max = bin_config["eta_edges"][0], bin_config["eta_edges"][1]
        min_pt = bin_config["min_pt"]
        fit_min_pt = bin_config["fit_min_pt"]
        pt_bins = bin_config.get("pt_bins", None)
        fit_max = pt_bins[-1] if pt_bins else 1000.0 
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")
        
        # Initial fit seeds for the baseline
        init_N = bin_config.get("N", 3.0)
        init_S = bin_config.get("S", 1.0)
        init_C = bin_config.get("C", 0.05)

        # -------------------------------------------------------------
        # 1. Extraction Loop
        # -------------------------------------------------------------
        extracted_hists = []
        for src_idx, source in enumerate(input_sources):
            # Fetch the 2D profile from the file
            h2 = source["file"].Get(source["key"])
            if not h2:
                print(f"ERROR: Could not find '{source['key']}' in '{source['label']}' file. Skipping this source.")
                continue
            
            # Use our mathematical helper
            h_res = extract_jer_profile(h2, f"src{src_idx}_eta{i}", eta_min, eta_max, pt_bins)
            
            # Apply Visual Styling
            h_res.SetMarkerStyle(source.get("marker", ROOT.kFullCircle))
            h_res.SetMarkerSize(0.7)
            h_res.SetMarkerColor(source.get("color", ROOT.kBlack))
            h_res.SetLineColor(source.get("color", ROOT.kBlack))
            h_res.SetStats(0)
            
            extracted_hists.append((source["label"], h_res))

        if not extracted_hists:
            continue # Skip if no data could be loaded

        # -------------------------------------------------------------
        # 2. Baseline Fit (Only on the FIRST input source)
        # -------------------------------------------------------------
        baseline_label, h_baseline = extracted_hists[0]
        
        f_nsc = ROOT.TF1(f"f_nsc_{i}", "sqrt(([0]/x)**2 + ([1]/sqrt(x))**2 + [2]**2)", fit_min_pt, fit_max)
        f_nsc.SetParameters(init_N, init_S, init_C)
        f_nsc.SetParLimits(0, 0.5, 100.0) 
        f_nsc.SetParLimits(1, 0.0, 10.0)  
        f_nsc.SetParLimits(2, 0.01, 1.0)  

        fit_res = h_baseline.Fit(f_nsc, "Q R S 0")
        fit_status = int(fit_res) if fit_res.Get() else -1
        N, S, C = f_nsc.GetParameter(0), f_nsc.GetParameter(1), f_nsc.GetParameter(2)

        # Generate the 1-Sigma Band for Baseline
        h_band = ROOT.TH1D(f"h_band_{i}", "", 1000, min_pt, fit_max * 1.1)
        h_band.SetDirectory(0)
        if fit_status == 0:
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_band)

        # Save baseline parameters for subsequent Adaptive Binning steps
        jer_results.append({
            "eta_min": eta_min,
            "eta_max": eta_max,
            "N": N, "S": S, "C": C
        })

        # -------------------------------------------------------------
        # 3. Plotting
        # -------------------------------------------------------------
        c = ROOT.TCanvas(f"c_jer_{i}", "JER Overlay", 800, 600)
        c.SetLogx(); c.SetGrid()
        
        # Setup the frame based on the baseline histogram
        h_baseline.GetXaxis().SetRangeUser(min_pt, fit_max * 1.1)
        h_baseline.GetXaxis().SetNoExponent(True) 
        h_baseline.GetXaxis().SetMoreLogLabels(True)
        h_baseline.GetYaxis().SetRangeUser(0.0, 0.5) 
        h_baseline.Draw("AXIS")
        
        # Draw 1-sigma band
        if fit_status == 0:
            h_band.SetFillColor(ROOT.kOrange)
            h_band.SetLineColor(ROOT.kOrange)
            h_band.SetMarkerSize(0)
            h_band.Draw("E3 SAME")

        # Draw fit line
        f_nsc.SetLineColor(extracted_hists[0][1].GetLineColor())
        f_nsc.SetLineWidth(2)
        f_nsc.Draw("L SAME")
        
        # Draw all histograms (baseline and alternatives)
        for label, h in reversed(extracted_hists): 
            h.Draw("PE SAME")
            
        # Build Legend Dynamically
        leg = ROOT.TLegend(0.50, 0.68, 0.88, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.035)
        for label, h in extracted_hists:
            leg.AddEntry(h, label, "pe")
        leg.AddEntry(f_nsc, f"{baseline_label} Fit", "l")
        if fit_status == 0:
            leg.AddEntry(h_band, "Fit #pm 1#sigma", "f")
        leg.Draw()
        
        # Parameter text box for the baseline
        pave = ROOT.TPaveText(0.50, 0.40, 0.88, 0.65, "NDC")
        pave.SetFillStyle(0); pave.SetBorderSize(0); pave.SetTextFont(42); pave.SetTextSize(0.035)
        pave.AddText(f"Baseline: {baseline_label}")
        pave.AddText(r'#frac{#sigma}{#mu} = #sqrt{(N/p_{T})^{2} + (S/#sqrt{p_{T}})^{2} + C^{2}}')
        pave.AddText(f'N = {N:.2f} #pm {f_nsc.GetParError(0):.2f}')
        pave.AddText(f'S = {S:.3f} #pm {f_nsc.GetParError(1):.3f}')
        pave.AddText(f'C = {C:.4f} #pm {f_nsc.GetParError(2):.4f}')
        pave.Draw()
        
        c.SaveAs(os.path.join(out_dir, f"JER_Overlay_{eta_str}.png"))
        
        # Memory Cleanup
        c.Close()
        f_nsc.Delete()
        h_band.Delete()
        for _, h in extracted_hists:
            h.Delete()

    return jer_results

if __name__ == "__main__":
    
    # =======================================================================
    # 1. SETUP INPUT SOURCES
    # Define your files, the profile key, and visual style here.
    # IMPORTANT: The FIRST item in this list is treated as the Baseline/Data.
    # The fit and 1-sigma band will ONLY be applied to the first item!
    # =======================================================================
    
    input_sources = []
    
    # 1st Source (Baseline / Data)
    f_data = ROOT.TFile("../HistogramsMaker/testHistos/test_closure_data.root")
    if not f_data or f_data.IsZombie():
        print("WARNING: Could not open Data file. Generating dummy Data.")
        f_data = ROOT.TFile("dummy_data.root", "RECREATE")
        generate_db_dummy_data(is_mc=False).Write("photonjet/MPFx_2D")
    
    input_sources.append({
        "file": f_data,
        "key": "photonjet/MPFx_2D",
        "label": "Data",
        "color": ROOT.kBlack,
        "marker": ROOT.kFullCircle
    })
    
    # 2nd Source (MC)
    f_mc = ROOT.TFile("../../rootfiles/reweighted_J4PHists_photonjet_GJ-4Jets.root")
    if not f_mc or f_mc.IsZombie():
        print("WARNING: Could not open MC file. Generating dummy MC.")
        f_mc = ROOT.TFile("dummy_mc.root", "RECREATE")
        generate_db_dummy_data(is_mc=True, mc_index=1).Write("photonjet/MPFx_2D")
        
    input_sources.append({
        "file": f_mc,
        "key": "MPFx_2D",
        "label": "MC (MG5)",
        "color": ROOT.kRed,
        "marker": ROOT.kOpenSquare
    })
    
    # You could theoretically add a 3rd or 4th source here!
    # input_sources.append({"file": f_alt, "key": "photonjet/MPFx_2D", "label": "Alt MC", "color": ROOT.kBlue, "marker": ROOT.kOpenTriangleUp})

    # =======================================================================
    # 2. DEFINE KINEMATIC BINS
    # =======================================================================
    custom_pt_bins = []
    regions = [
        (30, 300, 5),      
        (300, 500, 10),   
        (500, 700, 20), 
        (700, 1000, 50)  
    ]
    for start, end, step in regions:
        for val in range(start, end, step):
            custom_pt_bins.append(val)

    custom_pt_bins = [ 30, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 220, 300, 400, 500, 600, 800, 1000]
    adaptive_config = [
        {"eta_edges": [0.0, 0.261], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [0.261, 0.522], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [0.522, 0.783], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [0.783, 1.044], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [1.044, 1.305], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [1.305, 1.479], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [1.479, 1.653], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [1.653, 1.930], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [1.930, 2.322], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [2.322, 2.650], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [2.650, 2.964], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "fit_min_pt":60.0},
        {"eta_edges": [2.964, 5.191], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt":60.0}
    ]
    
    # =======================================================================
    # 3. EXECUTE
    # =======================================================================
    print("Extracting and Fitting JERs...")
    jer_data = process_and_plot_jers(input_sources, adaptive_config)
    
    with open("jer_parameters.json", "w") as outfile:
        json.dump(jer_data, outfile, indent=4)
        
    # Cleanup files
    f_data.Close()
    f_mc.Close()
        
    print("Done! Check 'jer_parameters.json' and the './jer_fits' folder.")
