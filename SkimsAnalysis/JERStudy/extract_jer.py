import ROOT
import os
import json
import math
import array

ROOT.TGaxis.SetMaxDigits(3)

import subprocess
import sys

def extract_and_fit_jer(h2_mpfx, config_list, out_dir="./jer_fits"):
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    
    jer_results = []

    for i, bin_config in enumerate(config_list):
        eta_min, eta_max = bin_config["eta_edges"][0], bin_config["eta_edges"][1]
        
        min_pt = bin_config["min_pt"]
        fit_min_pt = bin_config["fit_min_pt"]

        init_N = bin_config.get("N", 30.0)
        init_S = bin_config.get("S", 1.0)
        init_C = bin_config.get("C", 0.1)
        
        pt_bins = bin_config.get("pt_bins", None)
        
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        b_low = h2_mpfx.GetXaxis().FindBin(eta_min + 1e-5)
        b_high = h2_mpfx.GetXaxis().FindBin(eta_max - 1e-5)
        
        h_prof_raw = h2_mpfx.ProfileY(f"h_prof_raw_{i}", b_low, b_high)
        h_prof_raw.SetDirectory(0)
        
        if h2_mpfx.GetXaxis().GetXmin() < 0:
            b_neg_low = h2_mpfx.GetXaxis().FindBin(-eta_max + 1e-5)
            b_neg_high = h2_mpfx.GetXaxis().FindBin(-eta_min - 1e-5)
            h_neg = h2_mpfx.ProfileY(f"h_neg_{i}", b_neg_low, b_neg_high)
            h_neg.SetDirectory(0)
            h_prof_raw.Add(h_neg)

        # --- THE FIX: Apply Custom Binning ---
        if pt_bins:
            bin_array = array.array('d', pt_bins)
            h_prof = h_prof_raw.Rebin(len(pt_bins) - 1, f"h_prof_{i}", bin_array)
            h_prof.SetDirectory(0)
        else:
            h_prof = h_prof_raw
            h_prof.SetName(f"h_prof_{i}")

        xbins = h_prof.GetXaxis().GetXbins()
        if xbins.GetSize() > 0: 
            h_res = ROOT.TH1D(f"h_res_{i}", f"Extracted JER |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Tag p_{{T}} (GeV); JER(MPFx) ",
                              h_prof.GetNbinsX(), xbins.GetArray())
        else: 
            h_res = ROOT.TH1D(f"h_res_{i}", f"Extracted JER |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Tag p_{{T}} (GeV); JER(MPFx) ",
                              h_prof.GetNbinsX(), h_prof.GetXaxis().GetXmin(), h_prof.GetXaxis().GetXmax())
        
        for iy in range(1, h_prof.GetNbinsX() + 1):
            entries = h_prof.GetBinEntries(iy)
            err = h_prof.GetBinError(iy)
            mean = h_prof.GetBinContent(iy)
            
            if entries > 10 and mean > 0.0: 
                sigma = err * math.sqrt(entries)
                rel_res = sigma / mean
                
                # Propagate the statistical error on sigma (sigma / sqrt(2N)) to the relative resolution
                err_rel_res = (sigma / math.sqrt(2 * entries)) / mean
                
                h_res.SetBinContent(iy, rel_res)
                h_res.SetBinError(iy, err_rel_res) 

        # Fit with NSC Formula
        # Extended fit maximum automatically if custom bins go higher
        fit_max = pt_bins[-1] if pt_bins else 1000.0 
        f_nsc = ROOT.TF1(f"f_nsc_{i}", "sqrt(([0]/x)**2 + ([1]/sqrt(x))**2 + [2]**2)", fit_min_pt, fit_max)
        
        f_nsc.SetParameters(init_N, init_S, init_C)
        f_nsc.SetParLimits(0, 1.0, 100.0) 
        f_nsc.SetParLimits(1, 0.0, 10.0)  
        f_nsc.SetParLimits(2, 0.001, 0.3)  
        
        fit_res = h_res.Fit(f_nsc, "Q R S 0")
        fit_status = int(fit_res) if fit_res.Get() else -1

        # --- Generate the 1-Sigma Band ---
        # Using a highly binned (1000 bins) empty histogram to make the curve look smooth
        h_band = ROOT.TH1D(f"h_band_{i}", "", 1000, min_pt, fit_max * 1.1)
        h_band.SetDirectory(0)
        
        if fit_status == 0:
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_band)

        N, S, C = f_nsc.GetParameter(0), f_nsc.GetParameter(1), f_nsc.GetParameter(2)
        jer_results.append({
            "eta_min": eta_min,
            "eta_max": eta_max,
            "N": N, "S": S, "C": C
        })

        c = ROOT.TCanvas(f"c_jer_{i}", "JER Fit", 800, 600)
        c.SetLogx()
        c.SetGrid()
        
        h_res.SetMarkerStyle(ROOT.kFullCircle)
        h_res.SetMarkerSize(0.7)
        h_res.SetLineColor(ROOT.kBlack)
        h_res.GetXaxis().SetRangeUser(min_pt, fit_max * 1.1)
        h_res.GetXaxis().SetNoExponent(True) 
        h_res.GetXaxis().SetMoreLogLabels(True)
        h_res.GetYaxis().SetRangeUser(0.0, 0.5) 
        h_res.SetStats(0)
        # Draw frame/axes first
        h_res.Draw("AXIS")
        
        if fit_status == 0:
            h_band.SetFillColor(ROOT.kOrange)
            h_band.SetLineColor(ROOT.kOrange) # Hides the outline of the histogram
            h_band.SetMarkerSize(0)
            h_band.Draw("E3 SAME")
        
        f_nsc.SetLineColor(ROOT.kBlue)
        f_nsc.SetLineWidth(3)
        f_nsc.Draw("L SAME")
        
        # Draw the data points back on top
        h_res.Draw("PE SAME")
        
        # Add legend
        leg = ROOT.TLegend(0.50, 0.72, 0.88, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.035)
        leg.AddEntry(h_res, "Data", "pe")
        leg.AddEntry(f_nsc, "Fit curve", "l")
        if fit_status == 0:
            leg.AddEntry(h_band, "1#sigma Uncertainty Band", "f")
        leg.Draw()
        
        # Shifted the parameter block slightly down to accommodate the legend
        pave = ROOT.TPaveText(0.50, 0.45, 0.88, 0.70, "NDC")
        pave.SetFillStyle(0); pave.SetBorderSize(0); pave.SetTextFont(42); pave.SetTextSize(0.035)
        pave.AddText(r'JER(p_{T}) = #sqrt{(N/p_{T})^{2} + (S/#sqrt{p_{T}})^{2} + C^{2}}')
        pave.AddText(f'N = {N:.2f}')
        pave.AddText(f'S = {S:.3f}')
        pave.AddText(f'C = {C:.4f}')
        pave.Draw()
        
        c.SaveAs(os.path.join(out_dir, f"JER_Fit_{eta_str}.png"))

    return jer_results

if __name__ == "__main__":
    print("Loading inputs JER Extraction...")
    
    f = ROOT.TFile("../HistogramsMaker/testHistos/test_closure_data.root")
    
    h2_mpfx = f.Get("photonjet/MPFx_2D") 
    if not h2_mpfx:
        print("ERROR: 2D profile not found. Exiting.")
        exit(1)
    h2_mpfx.SetDirectory(0)
    
    #f = ROOT.TFile("../../rootfiles/reweighted_J4PHists_photonjet_GJ-4Jets.root")
    
        
    # --- Custom Binning ---

        
    # custom_pt_bins = []
    
    # regions = [
    #     (30, 300, 5),      
    #     (300, 500, 10),   
    #     (500, 700, 20), 
    #     (700, 1000, 50)  
    # ]
    
    # for start, end, step in regions:
    #     for val in range(start, end, step):
    #         custom_pt_bins.append(val)
    
    #custom_pt_bins = [20, 28, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 220, 300, 430, 638, 1032, 2000, 7000]
    custom_pt_bins = [ 30, 40, 44, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 220, 300, 400, 500, 600, 800, 1000]

    adaptive_config = [
        {"eta_edges": [0.0, 0.261], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt":30.0},
        {"eta_edges": [0.261, 0.522], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt":30.0},
        {"eta_edges": [0.522, 0.783], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt":30.0},
        {"eta_edges": [0.783, 1.044], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt":30.0},
        {"eta_edges": [1.044, 1.305], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt":30.0},
        {"eta_edges": [1.305, 1.479], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt":30.0},
        {"eta_edges": [1.479, 1.653], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt":30.0},
        {"eta_edges": [1.653, 1.930], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt":30.0},
        {"eta_edges": [1.930, 2.322], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt":30.0},
        {"eta_edges": [2.322, 2.650], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt":30.0},
        {"eta_edges": [2.650, 2.964], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt":30.0},
        {"eta_edges": [2.964, 5.191], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt":30.0}
    ]
    
    # adaptive_config = [
    #     {"eta_edges": [0.0, 0.261], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [0.261, 0.522], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [0.522, 0.783], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [0.783, 1.044], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [1.044, 1.305], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "fit_min_pt": 60.0},

    #     {"eta_edges": [1.305, 1.479], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [1.479, 1.653], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [1.653, 1.930], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [1.930, 2.172], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [2.172, 2.322], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [2.322, 2.500], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "fit_min_pt": 60.0},

    #     {"eta_edges": [2.500, 2.650], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [2.650, 2.853], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [2.853, 2.964], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [2.964, 3.139], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "fit_min_pt": 60.0},

    #     {"eta_edges": [3.139, 3.489], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [3.489, 3.839], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt": 60.0},
    #     {"eta_edges": [3.839, 5.191], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt": 60.0}
    # ]
    

    # adaptive_config = [
    #     {"eta_edges": [0.0, 0.261], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [0.261, 0.522], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [0.522, 0.783], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [0.783, 1.044], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [1.044, 1.305], "pt_bins": custom_pt_bins, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 80.0},

    #     {"eta_edges": [1.305, 1.479], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [1.479, 1.653], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [1.653, 1.930], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [1.930, 2.172], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [2.172, 2.322], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [2.322, 2.500], "pt_bins": custom_pt_bins, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 80.0},

    #     {"eta_edges": [2.500, 2.650], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [2.650, 2.853], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [2.853, 2.964], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [2.964, 3.139], "pt_bins": custom_pt_bins, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 80.0},

    #     {"eta_edges": [3.139, 3.489], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [3.489, 3.839], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt": 80.0},
    #     {"eta_edges": [3.839, 5.191], "pt_bins": custom_pt_bins, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 30.0, "fit_min_pt": 80.0}
    # ]
    
    jer_data = extract_and_fit_jer(h2_mpfx, adaptive_config)
    
    with open("jer_parameters.json", "w") as outfile:
        json.dump(jer_data, outfile, indent=4)

    grid_plot_command = [
    sys.executable,      # Uses the current Python interpreter
    "../../combine_plots.py",  # The script to run
    "./jer_fits", # Parsed arguments
    "--pattern", "JER_Fit*",   
    "--name", "JER_Fit_Grid"       
    ]

    subprocess.run(grid_plot_command)
        
    print("Done! Check 'jer_parameters.json' and the './jer_fits' folder.")