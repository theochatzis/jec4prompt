import ROOT
import os
import json
import math
import array # Needed to pass double arrays to ROOT in Python

ROOT.TGaxis.SetMaxDigits(3)

def plot_jer_profiles(config_list, out_dir="./adaptive_bins"):
    """
    Reads the configuration list and generates plots comparing the 
    physical Jet Energy Resolution to the target statistical uncertainty profile.
    """
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)

    for i, bin_config in enumerate(config_list):
        eta_min = bin_config["eta_edges"][0]
        eta_max = bin_config["eta_edges"][1]
        N, S, C = bin_config["N"], bin_config["S"], bin_config["C"]
        unc_profile = bin_config["uncertainty_profile"]
        
        # We need a finite upper limit just for drawing the plots
        draw_min_pt = bin_config.get("min_pt", 30.0)
        draw_max_pt = 3000.0 

        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        c = ROOT.TCanvas(f"c_prof_{i}", "JER Profile", 800, 600)
        c.SetLogx()
        c.SetGrid()

        # 1. Draw the Physical JER Curve using a TF1
        f_jer = ROOT.TF1(f"jer_{i}", "sqrt(([0]/x)**2 + ([1]/sqrt(x))**2 + [2]**2)", draw_min_pt, draw_max_pt)
        f_jer.SetParameters(N, S, C)
        f_jer.SetLineColor(ROOT.kBlue + 1)
        f_jer.SetLineWidth(3)
        f_jer.SetTitle(f"Resolution & Constraints |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Jet p_{{T}} (GeV);Resolution / Target Unc")
        
        # Set Y-axis scale to max 30% so we can see the HF region comfortably
        f_jer.SetMinimum(0.0)
        f_jer.SetMaximum(0.30)
        f_jer.GetXaxis().SetNoExponent(True)
        f_jer.GetXaxis().SetMoreLogLabels(True)
        
        f_jer.Draw("L")

        # 2. Draw the Target Uncertainty Step Profile using a TGraph
        g_unc = ROOT.TGraph()
        pt_cursor = draw_min_pt
        
        point_idx = 0
        for step in unc_profile:
            step_max = step["pt_max"] if step["pt_max"] != float('inf') else draw_max_pt
            step_unc = step["unc"]
            
            # Left point of the step (e.g., 30 GeV, 1.5%)
            g_unc.SetPoint(point_idx, pt_cursor, step_unc)
            point_idx += 1
            
            # Right point of the step (e.g., 300 GeV, 1.5%)
            g_unc.SetPoint(point_idx, step_max, step_unc)
            point_idx += 1
            
            # Move cursor for the next step
            pt_cursor = step_max
            if pt_cursor >= draw_max_pt: 
                break

        g_unc.SetLineColor(ROOT.kRed)
        g_unc.SetLineWidth(3)
        g_unc.SetLineStyle(2) # Dashed line for target
        g_unc.Draw("L SAME")

        # 3. Add a standard TLegend
        leg = ROOT.TLegend(0.40, 0.75, 0.88, 0.88)
        leg.SetFillStyle(0)   # Transparent background
        leg.SetBorderSize(0)
        leg.SetTextFont(42)   # TDR font
        leg.AddEntry(f_jer, f"Physical JER (N={N:.1f}, S={S:.1f}, C={C:.2f})", "l")
        leg.AddEntry(g_unc, "Target Stat. Uncertainty", "l")
        leg.Draw()

        # Save the plot
        c.SaveAs(os.path.join(out_dir, f"ResolutionProfile_{eta_str}.png"))

def generate_dummy_mpf_data():
    """Generates a dummy TProfile2D of MPF Response to test the logic."""
    h2 = ROOT.TProfile2D("h2_dummy", "MPF vs Probe #eta and Tag p_{T};Probe #eta;Tag p_{T} (GeV);MPF Response", 
                         100, -5.2, 5.2,      
                         2000, 10.0, 2000.0)  
    
    for _ in range(1000000):
        eta = ROOT.gRandom.Uniform(-5.2, 5.2)
        pt = ROOT.gRandom.Landau(100, 30) 
        mpf_response = ROOT.gRandom.Gaus(1.0, 0.15)
        h2.Fill(eta, pt, mpf_response)
        
    return h2

def format_custom_json(data_list):
    """Custom JSON formatter to keep the pt_binning array on a single line."""
    lines = ["[\n"]
    for i, d in enumerate(data_list):
        lines.append('  {\n')
        lines.append(f'    "abs_eta_min": {d["abs_eta_min"]},\n')
        lines.append(f'    "abs_eta_max": {d["abs_eta_max"]},\n')
        lines.append(f'    "pt_binning": {d["pt_binning"]}\n')
        if i < len(data_list) - 1:
            lines.append('  },\n')
        else:
            lines.append('  }\n')
    lines.append("]")
    return "".join(lines)

def find_adaptive_uncertainty_binning(h2_prof, config_list, out_dir="./adaptive_bins"):
    
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    json_output_data = []

    for i, bin_config in enumerate(config_list):
        eta_min = bin_config["eta_edges"][0]
        eta_max = bin_config["eta_edges"][1]
        unc_profile = bin_config["uncertainty_profile"] # <--- NEW: Grab the profile list
        
        N = bin_config["N"]
        S = bin_config["S"]
        C = bin_config["C"]
        
        min_pt_limit = bin_config.get("min_pt", None)
        max_pt_limit = bin_config.get("max_pt", None) 
        
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        # 1. Project the Positive Eta Side
        bin_pos_low = h2_prof.GetXaxis().FindBin(eta_min + 1e-5)
        bin_pos_high = h2_prof.GetXaxis().FindBin(eta_max - 1e-5)
        h_pt = h2_prof.ProfileY(f"h_pt_{i}", bin_pos_low, bin_pos_high)

        # 2. Fold in the Negative Eta Side
        if h2_prof.GetXaxis().GetXmin() < 0:
            bin_neg_low = h2_prof.GetXaxis().FindBin(-eta_max + 1e-5)
            bin_neg_high = h2_prof.GetXaxis().FindBin(-eta_min - 1e-5)
            h_neg = h2_prof.ProfileY(f"h_neg_{i}", bin_neg_low, bin_neg_high)
            h_pt.Add(h_neg) 

        # 3. Handle Min/Max Limits
        if min_pt_limit is not None:
            first_bin_to_process = h_pt.GetXaxis().FindBin(min_pt_limit + 1e-5)
        else:
            first_bin_to_process = 1

        if max_pt_limit is not None:
            last_bin_to_process = h_pt.GetXaxis().FindBin(max_pt_limit - 1e-5)
        else:
            last_bin_to_process = h_pt.GetNbinsX()

        edges = [h_pt.GetXaxis().GetBinLowEdge(first_bin_to_process)]
        W_M, SumWY_M, SumWY2_M = 0.0, 0.0, 0.0

        for iy in range(first_bin_to_process, last_bin_to_process + 1):
            W_i = h_pt.GetBinEntries(iy)
            if W_i == 0: 
                continue
                
            Y_i = h_pt.GetBinContent(iy)
            E_i = h_pt.GetBinError(iy)
            
            if W_i > 1:
                s_i = (E_i**2) * W_i * (W_i - 1)
            else:
                s_i = 0.0 
                
            SumWY2_i = s_i + (Y_i**2) * W_i
            SumWY_i = Y_i * W_i
            
            W_M += W_i
            SumWY_M += SumWY_i
            SumWY2_M += SumWY2_i
            
            if W_M > 1:
                s_M = SumWY2_M - (SumWY_M**2) / W_M
                if s_M < 0: s_M = 0 
                
                E_M = math.sqrt(s_M / (W_M * (W_M - 1)))
                
                current_edge_start = edges[-1]
                current_edge_end = h_pt.GetXaxis().GetBinUpEdge(iy)
                bin_width = current_edge_end - current_edge_start
                
                # --- NEW: DYNAMIC TARGET UNCERTAINTY LOOKUP ---
                # Check the current starting pT against our profile ranges
                current_target_unc = unc_profile[-1]["unc"] # Fallback to the last defined tolerance
                for prof in unc_profile:
                    if current_edge_start < prof["pt_max"]:
                        current_target_unc = prof["unc"]
                        break
                
                # Dynamic NSC JER Resolution Constraint
                pt_val = max(current_edge_start, 1e-3) 
                dynamic_jer_fraction = math.sqrt((N / pt_val)**2 + (S / math.sqrt(pt_val))**2 + C**2)
                min_required_width = dynamic_jer_fraction * pt_val 
                
                if E_M <= current_target_unc and bin_width >= min_required_width:
                    edges.append(current_edge_end)
                    W_M, SumWY_M, SumWY2_M = 0.0, 0.0, 0.0 
        
        # (Tail Protection safely removed per your previous request!)

        clean_edges = [int(round(e)) for e in edges]

        json_output_data.append({
            "abs_eta_min": eta_min,
            "abs_eta_max": eta_max,
            "pt_binning": clean_edges
        })

        # --- RE-BIN AND DRAW THE CLOSURE ---
        c = ROOT.TCanvas(f"c_{i}", "Adaptive Binning", 800, 600)
        c.SetLogx()
        c.SetGrid()
        
        edges_array = array.array('d', clean_edges)
        h_rebin = h_pt.Rebin(len(clean_edges) - 1, f"h_rebinned_{i}", edges_array)
        
        # Updated Title to show it's a variable profile
        h_rebin.SetTitle(f"MPF |#eta| #in [{eta_min:.3f}, {eta_max:.3f}] (Variable Target Unc, JER: N={N:.1f}, S={S:.2f}, C={C:.3f});Tag p_{{T}} (GeV);MPF Response")
        h_rebin.SetMarkerStyle(ROOT.kFullCircle)
        h_rebin.SetMarkerSize(0.8)
        h_rebin.SetMarkerColor(ROOT.kBlack)
        h_rebin.SetLineColor(ROOT.kBlack)
        h_rebin.GetXaxis().SetNoExponent(True)
        h_rebin.GetXaxis().SetMoreLogLabels(True)
        h_rebin.GetYaxis().SetRangeUser(0.8, 1.2) 
        h_rebin.SetStats(0)
        
        h_rebin.Draw("PE")

        drawn_lines = [] 
        for e in clean_edges[1:-1]:
            l_vert = ROOT.TLine(e, 0.8, e, 1.2)
            l_vert.SetLineColor(ROOT.kRed)
            l_vert.SetLineStyle(2)
            l_vert.SetLineWidth(2)
            l_vert.Draw("SAME")
            drawn_lines.append(l_vert)

        c.SaveAs(os.path.join(out_dir, f"AdaptiveClosure_{eta_str}.png"))

    return json_output_data

if __name__ == "__main__":
    print("Loading 2D MPF TProfile2D...\n")
    
    # 1. Open the ROOT file
    f = ROOT.TFile("../HistogramsMaker/testHistos/test_closure_data.root")
    
    # 2. Take the TProfile2D with MPF on z axis and probe eta on x-axis , tag pt on y-axis
    h2_prof = f.Get("photonjet/MPF_2D") 
    
    # ---------------------------------------------------------
    # CENTRAL CONFIGURATION LIST
    # ---------------------------------------------------------
    # ---------------------------------------------------------
    # DYNAMIC UNCERTAINTY PROFILES
    # ---------------------------------------------------------
    
    # BARREL (|eta| < 1.3)
    barrel_unc_profile = [
        {"pt_max": 300.0,  "unc": 0.010},  # 1.0% up to 300 GeV
        {"pt_max": 800.0,  "unc": 0.015},  # 1.5% up to 800 GeV
        {"pt_max": 2000.0, "unc": 0.025},  # 2.5% up to 2 TeV
        {"pt_max": 3000.0, "unc": 0.040},  # 4.0% up to 3 TeV
        {"pt_max": float('inf'), "unc": 0.060} # 6.0% for the extreme tail
    ]

    # INNER ENDCAP (1.3 <= |eta| < 2.5)
    inner_endcap_unc_profile = [
        {"pt_max": 300.0,  "unc": 0.015},  # 1.5% up to 300 GeV
        {"pt_max": 800.0,  "unc": 0.020},  # 2.0% up to 800 GeV
        {"pt_max": 1500.0, "unc": 0.030},  # 3.0% up to 1.5 TeV
        {"pt_max": float('inf'), "unc": 0.050} # 5.0% for the tail
    ]

    # OUTER ENDCAP (2.5 <= |eta| < 3.0)
    outer_endcap_unc_profile = [
        {"pt_max": 200.0,  "unc": 0.020},  # 2.0% up to 200 GeV
        {"pt_max": 500.0,  "unc": 0.030},  # 3.0% up to 500 GeV
        {"pt_max": 1000.0, "unc": 0.050},  # 5.0% up to 1 TeV
        {"pt_max": float('inf'), "unc": 0.060} # 6.0% for the tail
    ]

    # HADRONIC FORWARD / HF (3.0 <= |eta| < 5.2)
    hf_unc_profile = [
        {"pt_max": 200.0,  "unc": 0.030},  # 3.0% up to 200 GeV
        {"pt_max": 500.0,  "unc": 0.050},  # 5.0% up to 500 GeV
        {"pt_max": float('inf'), "unc": 0.080} # 8.0% for the tail
    ]

    adaptive_config = [
        # --- BARREL ---
        # Target Unc: 1%, N=3.0, S=0.9, C=0.04
        {"eta_edges": [0.0, 0.261], "uncertainty_profile": barrel_unc_profile, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0},
        {"eta_edges": [0.261, 0.522], "uncertainty_profile": barrel_unc_profile, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0},
        {"eta_edges": [0.522, 0.783], "uncertainty_profile": barrel_unc_profile, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0},
        {"eta_edges": [0.783, 1.044], "uncertainty_profile": barrel_unc_profile, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0},
        {"eta_edges": [1.044, 1.305], "uncertainty_profile": barrel_unc_profile, "N": 3.0, "S": 1.0, "C": 0.04, "min_pt": 30.0, "max_pt": 2000.0},

        # --- INNER ENDCAP ---
        # Target Unc: 2%, N=4.0, S=1.0, C=0.05
        {"eta_edges": [1.305, 1.479], "uncertainty_profile": inner_endcap_unc_profile, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0},
        {"eta_edges": [1.479, 1.653], "uncertainty_profile": inner_endcap_unc_profile, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0},
        {"eta_edges": [1.653, 1.930], "uncertainty_profile": inner_endcap_unc_profile, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0},
        {"eta_edges": [1.930, 2.172], "uncertainty_profile": inner_endcap_unc_profile, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0},
        {"eta_edges": [2.172, 2.322], "uncertainty_profile": inner_endcap_unc_profile, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0},
        {"eta_edges": [2.322, 2.500], "uncertainty_profile": inner_endcap_unc_profile, "N": 4.0, "S": 1.0, "C": 0.05, "min_pt": 30.0, "max_pt": 1500.0},

        # --- OUTER ENDCAP  ---
        # Loosening Target Unc to 3%, N=6.0, S=1.2, C=0.06
        {"eta_edges": [2.500, 2.650], "uncertainty_profile": outer_endcap_unc_profile, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 40.0, "max_pt": 1000.0},
        {"eta_edges": [2.650, 2.853], "uncertainty_profile": outer_endcap_unc_profile, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 40.0, "max_pt": 1000.0},
        {"eta_edges": [2.853, 2.964], "uncertainty_profile": outer_endcap_unc_profile, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 40.0, "max_pt": 1000.0},
        {"eta_edges": [2.964, 3.139], "uncertainty_profile": outer_endcap_unc_profile, "N": 6.0, "S": 1.2, "C": 0.06, "min_pt": 40.0, "max_pt": 1000.0},

        # --- HADRONIC FORWARD / HF---
        # Loosening Target Unc to 5%, N=8.0, S=1.5, C=0.08
        {"eta_edges": [3.139, 3.489], "uncertainty_profile": hf_unc_profile, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 50.0},
        {"eta_edges": [3.489, 3.839], "uncertainty_profile": hf_unc_profile, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 50.0},
        {"eta_edges": [3.839, 5.191], "uncertainty_profile": hf_unc_profile, "N": 8.0, "S": 1.5, "C": 0.08, "min_pt": 50.0}
    ]
    # ---------------------------------------------------------
    
    print("\nDrawing JER vs Uncertainty profiles...")
    plot_jer_profiles(adaptive_config)

    # Run the algorithm using our config dictionary and your actual data
    json_results = find_adaptive_uncertainty_binning(h2_prof, config_list=adaptive_config)
    
    # Make better the json file
    formatted_json_str = format_custom_json(json_results)

    # print("================ ADAPTIVE JSON OUTPUT ================\n")
    # print(json.dumps(json_results, indent=2))
    # formatted_json_str = format_custom_json(json_results)
    
    with open("adaptive_eta_pt_bins.json", "w") as outfile:
        outfile.write(formatted_json_str)
        
    # Keep good memory hygiene 
    f.Close()
