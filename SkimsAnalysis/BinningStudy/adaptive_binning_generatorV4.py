import ROOT
import os
import json
import math
import array

ROOT.TGaxis.SetMaxDigits(3)

def generate_complex_dummy_data():
    """Generates dummy data using Exponential (QCD) + Power Law (Signal)."""
    h2 = ROOT.TProfile2D("h2_dummy", "MPF vs Probe #eta and Tag p_{T};Probe #eta;Tag p_{T} (GeV);MPF Response", 
                         100, -5.2, 5.2, 2000, 10.0, 2000.0)  
    h2.SetDirectory(0)
    for _ in range(5000000):
        eta = ROOT.gRandom.Uniform(-5.2, 5.2)
        u = ROOT.gRandom.Uniform(0, 1)
        
        # Simulate QCD (Exponential) and Signal (Power Law)
        if u < 0.2:
            pt = 55.0 + ROOT.gRandom.Exp(40.0) # b = 1/40 = 0.025
        else:      
            pt = 55.0 / math.pow(1.0 - ROOT.gRandom.Uniform(0, 1), 1.0/(3.5 - 1.0))
        
        if pt > 3000: continue
        h2.Fill(eta, pt, ROOT.gRandom.Gaus(1.0, 0.15))
    return h2

def format_custom_json(data_list):
    lines = ["[\n"]
    for i, d in enumerate(data_list):
        lines.append('  {\n')
        lines.append(f'    "abs_eta_min": {d["abs_eta_min"]},\n')
        lines.append(f'    "abs_eta_max": {d["abs_eta_max"]},\n')
        lines.append(f'    "pt_binning": {d["pt_binning"]}\n')
        if i < len(data_list) - 1: lines.append('  },\n')
        else: lines.append('  }\n')
    lines.append("]")
    return "".join(lines)

def derive_complex_fit(h_pt, min_pt, max_pt, eta_str):
    """
    Fits N_tot * (B * PDF_Exp + A * PDF_PL).
    Returns: (f_fit, h_band, true_min_pt, fit_status, safe_N)
    """
    xbins = h_pt.GetXaxis().GetXbins()
    if xbins.GetSize() > 0: h_fit = ROOT.TH1D(f"h_fit_{eta_str}", "", h_pt.GetNbinsX(), xbins.GetArray())
    else: h_fit = ROOT.TH1D(f"h_fit_{eta_str}", "", h_pt.GetNbinsX(), h_pt.GetXaxis().GetXmin(), h_pt.GetXaxis().GetXmax())
    h_fit.SetDirectory(0)

    b1, b2 = h_pt.GetXaxis().FindBin(min_pt + 1e-5), h_pt.GetXaxis().FindBin(max_pt - 1e-5)
    
    first_data_bin = b1
    while first_data_bin <= b2 and h_pt.GetBinEntries(first_data_bin) == 0:
        first_data_bin += 1

    exp_pdf = "([0] * ([1]*exp(-[1]*x)) / (exp(-[1]*[5]) - exp(-[1]*[6])))"
    pl_pdf = "([2] * ((1.0-[3])*pow(x, -[3])) / (pow([6], 1.0-[3]) - pow([5], 1.0-[3])))"
    full_pdf_model = f"[4] * ( {exp_pdf} + {pl_pdf} )"

    if first_data_bin > b2: 
        f_fail = ROOT.TF1("f_fail", full_pdf_model, min_pt, max_pt)
        f_fail.SetParameters(0.0, 0.02, 1.0, 4.0, 1.0, min_pt, max_pt)
        return f_fail, None, min_pt, -1, 1.0

    true_min_pt = h_pt.GetXaxis().GetBinLowEdge(first_data_bin)
    
    N_tot = 0.0
    for iy in range(first_data_bin, b2 + 1):
        entries, width = h_pt.GetBinEntries(iy), h_pt.GetXaxis().GetBinWidth(iy)
        N_tot += entries
        h_fit.SetBinContent(iy, entries / width)
        h_fit.SetBinError(iy, math.sqrt(entries) / width if entries > 0 else 0)

    safe_N = max(N_tot, 1.0)
    
    f_fit = ROOT.TF1(f"f_fit_{eta_str}", full_pdf_model, true_min_pt, max_pt)
    
    # Initialize parameters
    init_b = 0.025
    init_alpha = 4.0
    init_B, init_A = 0.2, 0.8 # Exp fraction vs PL fraction

    f_fit.SetParameters(init_B, init_b, init_A, init_alpha, safe_N, true_min_pt, max_pt)
    
    f_fit.FixParameter(4, safe_N)
    f_fit.FixParameter(5, true_min_pt)
    f_fit.FixParameter(6, max_pt)

    # Hard bound the fractions between 0 and 1
    f_fit.SetParLimits(0, 0.0, 1.0)      # B (Exp Fraction)
    f_fit.SetParLimits(1, 0.001, 0.2)    # b (Exp slope)
    f_fit.SetParLimits(2, 0.0, 1.0)      # A (PL Fraction)
    f_fit.SetParLimits(3, 2.5, 6.0)      # alpha (PL slope)

    fit_res = h_fit.Fit(f_fit, "Q R S 0") 
    fit_status = int(fit_res) if fit_res.Get() else -1

    h_band = h_fit.Clone(f"h_band_{eta_str}")
    h_band.SetDirectory(0)
    h_band.Reset()

    if fit_status == 0:
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_band)
        
    return f_fit, h_band, true_min_pt, fit_status, safe_N

def plot_analytic_profiles(h2_prof, config_list, out_dir="./adaptive_bins"):
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)

    for i, bin_config in enumerate(config_list):
        eta_min, eta_max = bin_config["eta_edges"][0], bin_config["eta_edges"][1]
        N, S, C = bin_config["N"], bin_config["S"], bin_config["C"]
        draw_min_pt, max_pt_bound = bin_config.get("min_pt", 30.0), 3000.0 
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        b_low, b_high = h2_prof.GetXaxis().FindBin(eta_min + 1e-5), h2_prof.GetXaxis().FindBin(eta_max - 1e-5)
        h_pt_base = h2_prof.ProfileY(f"h_pt_prof_{i}", b_low, b_high)
        h_pt_base.SetDirectory(0)
        
        if h2_prof.GetXaxis().GetXmin() < 0:
            b_neg_low, b_neg_high = h2_prof.GetXaxis().FindBin(-eta_max + 1e-5), h2_prof.GetXaxis().FindBin(-eta_min - 1e-5)
            h_neg = h2_prof.ProfileY(f"h_neg_prof_{i}", b_neg_low, b_neg_high)
            h_neg.SetDirectory(0)
            h_pt_base.Add(h_neg)

        fit_min, fit_max = bin_config.get("fit_min_pt", draw_min_pt), bin_config.get("fit_max_pt", max_pt_bound)
        f_fit, _, true_fit_min, fit_status, safe_N = derive_complex_fit(h_pt_base, fit_min, fit_max, f"prof_{i}")
        
        B, b, A, alpha = f_fit.GetParameter(0), f_fit.GetParameter(1), f_fit.GetParameter(2), f_fit.GetParameter(3)
        p_min_fit, p_max_fit = f_fit.GetParameter(5), f_fit.GetParameter(6)
        pow_alpha = 1.0 - alpha

        c = ROOT.TCanvas(f"c_prof_{i}", "JER Profile", 800, 600)
        c.SetLogx(); c.SetGrid()

        f_jer = ROOT.TF1(f"jer_{i}", "sqrt(([0]/x)**2 + ([1]/sqrt(x))**2 + [2]**2)", draw_min_pt, max_pt_bound)
        f_jer.SetParameters(N, S, C)
        f_jer.SetLineColor(ROOT.kBlue + 1); f_jer.SetLineWidth(3)
        f_jer.SetTitle(f"Analytic Profile |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Jet p_{{T}} (GeV);Resolution / Target Unc")
        f_jer.SetMinimum(0.0); f_jer.SetMaximum(0.50)
        f_jer.GetXaxis().SetNoExponent(True); f_jer.GetXaxis().SetMoreLogLabels(True)
        f_jer.Draw("L")

        g_unc = ROOT.TGraph()
        for pt_step in range(int(draw_min_pt), int(max_pt_bound), 10):
            jer_val = math.sqrt((N/pt_step)**2 + (S/math.sqrt(pt_step))**2 + C**2)
            width = jer_val * pt_step
            if pt_step + width > max_pt_bound or width <= 0: continue

            if fit_status == 0:
                denom_exp = math.exp(-b * p_min_fit) - math.exp(-b * p_max_fit)
                num_exp = math.exp(-b * pt_step) - math.exp(-b * (pt_step + width))
                n_exp_comp = safe_N * B * (num_exp / denom_exp) if denom_exp > 0 else 0.0
                
                denom_pl = math.pow(p_max_fit, pow_alpha) - math.pow(p_min_fit, pow_alpha)
                num_pl = math.pow(pt_step + width, pow_alpha) - math.pow(pt_step, pow_alpha)
                n_pl_comp = safe_N * A * (num_pl / denom_pl) if denom_pl != 0 else 0.0
                
                n_exp = n_exp_comp + n_pl_comp
            else:
                denom_pl = math.pow(p_max_fit, pow_alpha) - math.pow(p_min_fit, pow_alpha)
                num_pl = math.pow(pt_step + width, pow_alpha) - math.pow(pt_step, pow_alpha)
                n_exp = safe_N * 1.0 * (num_pl / denom_pl) if denom_pl != 0 else 0.0
            
            stat_unc = 1.0 / math.sqrt(n_exp) if n_exp > 0 else 1.0
            
            # Cap the Target Uncertainty so it never exceeds the Resolution!
            target = min(stat_unc, jer_val)
            g_unc.SetPoint(g_unc.GetN(), pt_step, target)

        g_unc.SetLineColor(ROOT.kRed); g_unc.SetLineWidth(3); g_unc.SetLineStyle(2) 
        g_unc.Draw("L SAME")

        leg = ROOT.TLegend(0.35, 0.75, 0.88, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.03)
        leg.AddEntry(f_jer, f"Physical JER (N={N:.1f}, S={S:.1f}, C={C:.2f})", "l")
        leg.AddEntry(g_unc, f"Target (min(Stat, JER))" if fit_status == 0 else "Target (Fallback)", "l")
        leg.Draw()

        c.SaveAs(os.path.join(out_dir, f"ResolutionProfile_{eta_str}.png"))

def plot_complex_spectrum(h2_prof, config_list, out_dir="./adaptive_bins"):
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)

    for i, bin_config in enumerate(config_list):
        eta_min, eta_max = bin_config["eta_edges"][0], bin_config["eta_edges"][1]
        draw_min_pt, max_pt_limit = bin_config.get("min_pt", 30.0), bin_config.get("max_pt", 3000.0)
        fit_min, fit_max = bin_config.get("fit_min_pt", draw_min_pt), bin_config.get("fit_max_pt", max_pt_limit)
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        b_low, b_high = h2_prof.GetXaxis().FindBin(eta_min + 1e-5), h2_prof.GetXaxis().FindBin(eta_max - 1e-5)
        h_pt = h2_prof.ProfileY(f"h_pt_spec_{i}", b_low, b_high)
        h_pt.SetDirectory(0)
        
        if h2_prof.GetXaxis().GetXmin() < 0:
            b_neg_low, b_neg_high = h2_prof.GetXaxis().FindBin(-eta_max + 1e-5), h2_prof.GetXaxis().FindBin(-eta_min - 1e-5)
            h_neg = h2_prof.ProfileY(f"h_neg_spec_{i}", b_neg_low, b_neg_high)
            h_neg.SetDirectory(0)
            h_pt.Add(h_neg)

        f_fit, h_band, true_fit_min, fit_status, safe_N = derive_complex_fit(h_pt, fit_min, fit_max, f"spec_{i}")
        
        xbins = h_pt.GetXaxis().GetXbins()
        if xbins.GetSize() > 0: h_spectrum = ROOT.TH1D(f"h_spec_{i}", "", h_pt.GetNbinsX(), xbins.GetArray())
        else: h_spectrum = ROOT.TH1D(f"h_spec_{i}", "", h_pt.GetNbinsX(), h_pt.GetXaxis().GetXmin(), h_pt.GetXaxis().GetXmax())
        h_spectrum.SetDirectory(0)

        b1, b2 = h_pt.GetXaxis().FindBin(draw_min_pt + 1e-5), h_pt.GetXaxis().FindBin(max_pt_limit - 1e-5)
        for iy in range(b1, b2 + 1):
            entries, width = h_pt.GetBinEntries(iy), h_pt.GetXaxis().GetBinWidth(iy)
            h_spectrum.SetBinContent(iy, entries / width)
            h_spectrum.SetBinError(iy, math.sqrt(entries) / width if entries > 0 else 0)

        c = ROOT.TCanvas(f"c_spec_{i}", "pT Spectrum", 800, 600)
        c.SetLogx(); c.SetLogy(); c.SetGrid()

        h_spectrum.SetTitle(f"Tag Jet Spectrum |#eta| #in [{eta_min:.3f}, {eta_max:.3f}];Jet p_{{T}} (GeV);dN/dp_{{T}} (Events / GeV)")
        h_spectrum.SetMarkerStyle(ROOT.kFullCircle); h_spectrum.SetMarkerSize(0.6)
        h_spectrum.SetLineColor(ROOT.kBlack); h_spectrum.SetMarkerColor(ROOT.kBlack)
        h_spectrum.GetXaxis().SetRangeUser(draw_min_pt, max_pt_limit)
        h_spectrum.GetXaxis().SetNoExponent(True); h_spectrum.GetXaxis().SetMoreLogLabels(True)
        h_spectrum.SetStats(0)
        
        spec_max = h_spectrum.GetMaximum() * 10
        spec_min = 0.1
        h_spectrum.SetMaximum(spec_max); h_spectrum.SetMinimum(spec_min)
        h_spectrum.Draw("PE")

        if fit_status == 0:
            h_band.SetFillColor(ROOT.kOrange)
            h_band.SetMarkerSize(0)
            h_band.Draw("E3 SAME")
            
            f_fit.SetLineColor(ROOT.kBlue)
            f_fit.SetLineWidth(3)
            f_fit.Draw("L SAME")
            
            exp_str = f"[4] * ( ([0] * [1]*exp(-[1]*x)) / (exp(-[1]*[5]) - exp(-[1]*[6])) )"
            pl_str  = f"[4] * ( ([2] * (1.0-[3])*pow(x, -[3])) / (pow([6], 1.0-[3]) - pow([5], 1.0-[3])) )"
            
            f_qcd = ROOT.TF1(f"f_qcd_{i}", exp_str, true_fit_min, fit_max)
            f_qcd.SetParameters(f_fit.GetParameter(0), f_fit.GetParameter(1), f_fit.GetParameter(2), f_fit.GetParameter(3), safe_N, f_fit.GetParameter(5), f_fit.GetParameter(6))
            f_qcd.SetLineColor(ROOT.kGreen + 2); f_qcd.SetLineStyle(2); f_qcd.SetLineWidth(2)
            f_qcd.Draw("L SAME")
            
            f_pl = ROOT.TF1(f"f_pl_{i}", pl_str, true_fit_min, fit_max)
            f_pl.SetParameters(f_fit.GetParameter(0), f_fit.GetParameter(1), f_fit.GetParameter(2), f_fit.GetParameter(3), safe_N, f_fit.GetParameter(5), f_fit.GetParameter(6))
            f_pl.SetLineColor(ROOT.kRed); f_pl.SetLineStyle(2); f_pl.SetLineWidth(2)
            f_pl.Draw("L SAME")
            
            #h_spectrum.Draw("PE SAME")

            B, b, A, alpha = f_fit.GetParameter(0), f_fit.GetParameter(1), f_fit.GetParameter(2), f_fit.GetParameter(3)
            
            pave = ROOT.TPaveText(0.40, 0.42, 0.88, 0.70, "NDC")
            pave.SetFillStyle(0); pave.SetBorderSize(0); pave.SetTextFont(42); pave.SetTextSize(0.03)
            pave.AddText(r'Model: N_{tot} #times (B exp(-b p_{T}) + A p^{- #alpha}_{T})')
            pave.AddText(f'N_{{tot}} = {safe_N:.2e}')
            pave.AddText(f'B = {B:.3f} #pm {f_fit.GetParError(0):.3f}')
            pave.AddText(f'b = {b:.4f} #pm {f_fit.GetParError(1):.4f}')
            pave.AddText(f'A = {A:.3f} #pm {f_fit.GetParError(2):.3f}')
            pave.AddText(f'#alpha = {alpha:.2f} #pm {f_fit.GetParError(3):.2f}')
            pave.Draw()

        leg = ROOT.TLegend(0.40, 0.72, 0.88, 0.88)
        leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42); leg.SetTextSize(0.03)
        leg.AddEntry(h_spectrum, "Data Spectrum (dN/dp_{T})", "pe")
        if fit_status == 0:
            leg.AddEntry(h_band, "Fit Uncertainty (#pm 1#sigma)", "f")
            leg.AddEntry(f_fit, "Normalized PDF Fit", "l")
            leg.AddEntry(f_qcd, "QCD (Exp PDF)", "l")
            leg.AddEntry(f_pl, "Signal (PL PDF)", "l")
        leg.Draw()

        c.SaveAs(os.path.join(out_dir, f"SpectrumValidation_{eta_str}.png"))

def find_adaptive_uncertainty_binning(h2_prof, config_list, out_dir="./adaptive_bins"):
    os.makedirs(out_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    json_output_data = []

    for i, bin_config in enumerate(config_list):
        eta_min, eta_max = bin_config["eta_edges"][0], bin_config["eta_edges"][1]
        N, S, C = bin_config["N"], bin_config["S"], bin_config["C"]
        eta_str = f"abs_eta_{eta_min:.3f}_to_{eta_max:.3f}".replace(".", "p")

        b_low, b_high = h2_prof.GetXaxis().FindBin(eta_min + 1e-5), h2_prof.GetXaxis().FindBin(eta_max - 1e-5)
        h_pt = h2_prof.ProfileY(f"h_pt_{i}", b_low, b_high)
        h_pt.SetDirectory(0)
        
        if h2_prof.GetXaxis().GetXmin() < 0:
            b_neg_low, b_neg_high = h2_prof.GetXaxis().FindBin(-eta_max + 1e-5), h2_prof.GetXaxis().FindBin(-eta_min - 1e-5)
            h_neg = h2_prof.ProfileY(f"h_neg_{i}", b_neg_low, b_neg_high)
            h_neg.SetDirectory(0)
            h_pt.Add(h_neg)

        min_pt_limit, max_pt_limit = bin_config.get("min_pt", 30.0), bin_config.get("max_pt", 3000.0)
        fit_min, fit_max = bin_config.get("fit_min_pt", min_pt_limit), bin_config.get("fit_max_pt", max_pt_limit)
        
        f_fit, _, true_min_pt, fit_status, safe_N = derive_complex_fit(h_pt, fit_min, fit_max, f"alg_{i}")
        B, b, A, alpha = f_fit.GetParameter(0), f_fit.GetParameter(1), f_fit.GetParameter(2), f_fit.GetParameter(3)
        p_min_fit, p_max_fit = f_fit.GetParameter(5), f_fit.GetParameter(6)
        pow_alpha = 1.0 - alpha

        first_bin_to_process = h_pt.GetXaxis().FindBin(true_min_pt + 1e-5)
        last_bin_to_process = h_pt.GetXaxis().FindBin(max_pt_limit - 1e-5)

        edges = [h_pt.GetXaxis().GetBinLowEdge(first_bin_to_process)]
        W_M, SumWY_M, SumWY2_M = 0.0, 0.0, 0.0

        for iy in range(first_bin_to_process, last_bin_to_process + 1):
            W_i = h_pt.GetBinEntries(iy)
            if W_i == 0: continue
                
            Y_i, E_i = h_pt.GetBinContent(iy), h_pt.GetBinError(iy)
            s_i = (E_i**2) * W_i * (W_i - 1) if W_i > 1 else 0.0 
            SumWY2_i, SumWY_i = s_i + (Y_i**2) * W_i, Y_i * W_i
            W_M, SumWY_M, SumWY2_M = W_M + W_i, SumWY_M + SumWY_i, SumWY2_M + SumWY2_i
            
            if W_M > 1:
                s_M = max(0, SumWY2_M - (SumWY_M**2) / W_M)
                E_M = math.sqrt(s_M / (W_M * (W_M - 1)))
                
                pt_val, pt_end = edges[-1], h_pt.GetXaxis().GetBinUpEdge(iy)
                bin_width = pt_end - pt_val
                
                pt_val_safe = max(pt_val, 1e-3) 
                dynamic_jer_fraction = math.sqrt((N / pt_val_safe)**2 + (S / math.sqrt(pt_val_safe))**2 + C**2)
                min_required_width = dynamic_jer_fraction * pt_val_safe 
                
                if fit_status == 0:
                    denom_exp = math.exp(-b * p_min_fit) - math.exp(-b * p_max_fit)
                    num_exp = math.exp(-b * pt_val) - math.exp(-b * pt_end)
                    n_exp_comp = safe_N * B * (num_exp / denom_exp) if denom_exp > 0 else 0.0
                    
                    denom_pl = math.pow(p_max_fit, pow_alpha) - math.pow(p_min_fit, pow_alpha)
                    num_pl = math.pow(pt_end, pow_alpha) - math.pow(pt_val, pow_alpha)
                    n_pl_comp = safe_N * A * (num_pl / denom_pl) if denom_pl != 0 else 0.0
                    
                    n_exp = n_exp_comp + n_pl_comp
                else:
                    denom_pl = math.pow(p_max_fit, pow_alpha) - math.pow(p_min_fit, pow_alpha)
                    num_pl = math.pow(pt_end, pow_alpha) - math.pow(pt_val, pow_alpha)
                    n_exp = safe_N * 1.0 * (num_pl / denom_pl) if denom_pl != 0 else 0.0
                
                stat_unc = 1.0 / math.sqrt(n_exp) if n_exp > 0 else 1.0
                
                # The target uncertainty is capped by the detector resolution (JER)
                current_target_unc = min(stat_unc, dynamic_jer_fraction)
                
                if E_M <= current_target_unc and bin_width >= min_required_width:
                    edges.append(pt_end)
                    W_M, SumWY_M, SumWY2_M = 0.0, 0.0, 0.0 

        clean_edges = [int(round(e)) for e in edges]
        json_output_data.append({"abs_eta_min": eta_min, "abs_eta_max": eta_max, "pt_binning": clean_edges})

        c = ROOT.TCanvas(f"c_{i}", "Adaptive Binning", 800, 600)
        c.SetLogx(); c.SetGrid()
        
        if len(clean_edges) > 1:
            h_rebin = h_pt.Rebin(len(clean_edges) - 1, f"h_rebinned_{i}", array.array('d', clean_edges))
            label_complex = f"Fit (#alpha={alpha:.2f}, b={b:.3f})" if fit_status == 0 else "PL Fallback Model (#alpha=4.0)"
            h_rebin.SetTitle(f"MPF |#eta| #in [{eta_min:.3f}, {eta_max:.3f}] ({label_complex}, N={N:.1f}, S={S:.2f}, C={C:.3f});Tag p_{{T}} (GeV);MPF Response")
            h_rebin.SetMarkerStyle(ROOT.kFullCircle); h_rebin.SetMarkerSize(0.8)
            h_rebin.SetMarkerColor(ROOT.kBlack); h_rebin.SetLineColor(ROOT.kBlack)
            h_rebin.GetXaxis().SetNoExponent(True); h_rebin.GetXaxis().SetMoreLogLabels(True)
            h_rebin.GetYaxis().SetRangeUser(0.8, 1.2); h_rebin.SetStats(0)
            h_rebin.Draw("PE")

            c.SaveAs(os.path.join(out_dir, f"AdaptiveClosure_{eta_str}.png"))

    return json_output_data

if __name__ == "__main__":
    print("Loading 2D MPF TProfile2D...\n")
    
    f = ROOT.TFile("../HistogramsMaker/testHistos/test_closure_data.root")
    if not f or f.IsZombie():
        print("WARNING: Could not open file. Falling back to dummy data.")
        h2_prof = generate_complex_dummy_data()
    else:
        h2_prof = f.Get("photonjet/MPF_2D") 
        if not h2_prof:
            print("ERROR: Histogram not found. Exiting.")
            exit(1)
        h2_prof.SetDirectory(0)
    
    # -------------------------------------------------------------
    # 1. Base Configuration (No hardcoded N, S, C)
    # -------------------------------------------------------------
    base_config = [
        {"eta_edges": [0.0, 0.261],   "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 30.0},
        {"eta_edges": [0.261, 0.522], "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 30.0},
        {"eta_edges": [0.522, 0.783], "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 30.0},
        {"eta_edges": [0.783, 1.044], "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 30.0},
        {"eta_edges": [1.044, 1.305], "min_pt": 30.0, "max_pt": 2000.0, "fit_min_pt": 30.0},

        {"eta_edges": [1.305, 1.479], "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 30.0},
        {"eta_edges": [1.479, 1.653], "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 30.0},
        {"eta_edges": [1.653, 1.930], "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 30.0},
        {"eta_edges": [1.930, 2.172], "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 30.0},
        {"eta_edges": [2.172, 2.322], "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 30.0},
        {"eta_edges": [2.322, 2.500], "min_pt": 30.0, "max_pt": 1500.0, "fit_min_pt": 30.0},

        {"eta_edges": [2.500, 2.650], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0},
        {"eta_edges": [2.650, 2.853], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0},
        {"eta_edges": [2.853, 2.964], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0},
        {"eta_edges": [2.964, 3.139], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0},

        {"eta_edges": [3.139, 3.489], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0},
        {"eta_edges": [3.489, 3.839], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0},
        {"eta_edges": [3.839, 5.191], "min_pt": 30.0, "max_pt": 1000.0, "fit_min_pt": 30.0}
    ]

    # -------------------------------------------------------------
    # 2. Inject N, S, C from JSON
    # -------------------------------------------------------------
    if not os.path.exists("jer_parameters.json"):
        print("ERROR: jer_parameters.json not found! Run the extraction script first.")
        exit(1)
        
    with open("jer_parameters.json", "r") as jfile:
        jer_data = json.load(jfile)
        
    for conf in base_config:
        # Match by checking if the min and max eta are practically identical
        match = next((item for item in jer_data if abs(item["eta_min"] - conf["eta_edges"][0]) < 0.001 
                                              and abs(item["eta_max"] - conf["eta_edges"][1]) < 0.001), None)
        if match:
            conf["N"] = match["N"]
            conf["S"] = match["S"]
            conf["C"] = match["C"]
        else:
            print(f"WARNING: No JER parameters found for {conf['eta_edges']}. Using some defaults.")
            conf["N"], conf["S"], conf["C"] = 2.0, 1.0, 0.05
    
    # -------------------------------------------------------------
    # 3. Execute Algorithms
    # -------------------------------------------------------------
    json_results = find_adaptive_uncertainty_binning(h2_prof, config_list=base_config)
    with open("adaptive_eta_pt_bins.json", "w") as outfile: 
        outfile.write(format_custom_json(json_results))

    print("Drawing Analytic Resolution vs Dynamic Target Uncertainty profiles...")
    plot_analytic_profiles(h2_prof, config_list=base_config)
    
    print("Drawing pT Spectra vs Complex Fit (Exp + PL) profiles...")
    plot_complex_spectrum(h2_prof, config_list=base_config)
    
    print("Done! Check the './adaptive_bins' folder.")