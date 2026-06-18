import ROOT

# Force scientific notation globally for numbers > 3 digits
ROOT.TGaxis.SetMaxDigits(3)

import os
import json

def generate_dummy_data():
    """Generates a dummy 2D histogram with signed Eta to test the logic."""
    h2 = ROOT.TH2D("h2_dummy", "Tag p_{T} vs Probe #eta;Probe #eta;Tag p_{T} (GeV)", 
                   100, -5.2, 5.2,     # Signed Eta bins
                   2000, 10.0, 2000.0)  # Very fine pT bins (0.5 GeV per bin)
    
    for _ in range(500000):
        eta = ROOT.gRandom.Uniform(-5.2, 5.2)
        pt = ROOT.gRandom.Landau(100, 30) 
        h2.Fill(eta, pt)
        
    return h2

def find_dynamic_binning_json(h2, events_per_bin, 
                eta_edges = [0.0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 
                 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191], 
                 out_dir="./cdf_plots"):
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)

    json_output_data = []

    for i in range(len(eta_edges) - 1):
        eta_min = eta_edges[i]
        eta_max = eta_edges[i+1]
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        # 1. Project the Positive Eta Side
        bin_pos_low = h2.GetXaxis().FindBin(eta_min + 1e-5)
        bin_pos_high = h2.GetXaxis().FindBin(eta_max - 1e-5)
        h_pt = h2.ProjectionY(f"h_pt_{i}", bin_pos_low, bin_pos_high)

        # 2. Fold in the Negative Eta Side (if the histogram is signed)
        if h2.GetXaxis().GetXmin() < 0:
            bin_neg_low = h2.GetXaxis().FindBin(-eta_max + 1e-5)
            bin_neg_high = h2.GetXaxis().FindBin(-eta_min - 1e-5)
            h_neg = h2.ProjectionY(f"h_neg_{i}", bin_neg_low, bin_neg_high)
            h_pt.Add(h_neg) # Combine them!

        total_events = h_pt.Integral()
        if total_events == 0:
            continue
            
        # 3. Get the Cumulative Distribution
        h_cdf = h_pt.GetCumulative()
        h_cdf.SetTitle(f"CDF of Tag p_{{T}} for |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Tag p_{{T}} (GeV);Cumulative Events")
        h_cdf.SetLineColor(ROOT.kBlue + 1)
        h_cdf.SetLineWidth(2)
        h_cdf.SetStats(0)
        h_cdf.GetXaxis().SetNoExponent(True)      # Writes '100' instead of '10^2'
        h_cdf.GetXaxis().SetMoreLogLabels(True)   # Shows 20, 30, 40... 200, 300...

        # 4. Find the crossings
        edges = [h_pt.GetXaxis().GetBinLowEdge(1)] # Start with absolute minimum edge
        thresholds = []

        # Grab the target number of events specifically for THIS eta slice
        target_step = events_per_bin[i]
        current_target = target_step

        for iy in range(1, h_cdf.GetNbinsX() + 1):
            if h_cdf.GetBinContent(iy) >= current_target:
                edges.append(h_cdf.GetXaxis().GetBinUpEdge(iy))
                thresholds.append(current_target)
                current_target += target_step

        # # Check the "tail" (The leftover events after the last mark)
        # max_edge = h_pt.GetXaxis().GetBinUpEdge(h_pt.GetNbinsX())
        # leftover_events = total_events - (current_target - events_per_bin)
        
        # # If the high-pT tail has fewer than 50% of our target events, merge it 
        # # with the previous bin so our fits don't fail due to low statistics
        # if leftover_events < (0.5 * events_per_bin) and len(edges) > 2:
        #     edges[-1] = max_edge 
        # else:
        #     if edges[-1] != max_edge:
        #         edges.append(max_edge)

        # 5. Append to our JSON structure (Casting to int for clean outputs like [20, 28, 40])
        json_output_data.append({
            "abs_eta_min": eta_min,
            "abs_eta_max": eta_max,
            "pt_binning": [int(round(e)) for e in edges]
        })

        # 6. Plot the CDF and TLines for visual debugging
        c = ROOT.TCanvas(f"c_{i}", "CDF", 800, 600)
        c.SetGrid()
        c.SetLogx()
        h_cdf.Draw("HIST")

        drawn_lines = [] 
        for t, e in zip(thresholds, edges[1:-1]):
            # Horizontal line
            l_horiz = ROOT.TLine(h_cdf.GetXaxis().GetXmin(), t, e, t)
            l_horiz.SetLineColor(ROOT.kRed)
            l_horiz.SetLineStyle(2)
            l_horiz.SetLineWidth(3)
            l_horiz.Draw("SAME")
            drawn_lines.append(l_horiz)

            # Vertical line
            l_vert = ROOT.TLine(e, 0, e, t)
            l_vert.SetLineColor(ROOT.kGreen + 2)
            l_vert.SetLineStyle(2)
            l_vert.SetLineWidth(3)
            l_vert.Draw("SAME")
            drawn_lines.append(l_vert)

        c.SaveAs(os.path.join(out_dir, f"CDF_{eta_str}.png"))

    return json_output_data

if __name__ == "__main__":
    print("Loading 2D histogram...\n")
    
    # Loading your the TH2D with probe eta on x-axis and tag pt on y-axis:
    f = ROOT.TFile("../HistogramsMaker/testHistos/test_closure_data.root")
    h2_data = f.Get("photonjet/Probe_eta_Tag_pt")
    
    # To play with dummy data
    # h2_data = generate_dummy_data() 
    
    #eta_edges = [0.0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 
    #            2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191]
    eta_edges = [0.0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 
                2.172, 2.650, 2.964, 3.489, 5.191]
    
    events_per_bin = [
        10000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 
                5000, 5000, 5000, 5000, 5000
    ]
    
    json_results = find_dynamic_binning_json(h2_data, events_per_bin=events_per_bin, eta_edges=eta_edges)
    
    
    # This prints the exact formatted JSON to terminal
    #print("================ JSON OUTPUT ================\n")
    #print(json.dumps(json_results, indent=2))
    
    # Saving the eta bins to a json file
    with open("dynamic_eta_pt_bins.json", "w") as outfile:
        json.dump(json_results, outfile, indent=2)
