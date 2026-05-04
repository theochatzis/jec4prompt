import correctionlib
import numpy as np
import matplotlib.pyplot as plt
import os  

# =========================================================
# Inputs and parameters
# =========================================================

json_file = "./jsons/Summer24Prompt24JEC4PromptRun398027_JECs.json" # Change to your actual JSON name

campaign = "Summer24Prompt24JEC4PromptRun398027"
algo = "AK4PFPuppi"

# The exact names of the corrections inside your JSON
name_L1 = f"{campaign}_V1_DATA_L1FastJet_{algo}"
name_L2 = f"{campaign}_V1_DATA_L2Relative_{algo}"
name_L3Abs = f"{campaign}_V1_DATA_L3Absolute_{algo}"
name_L2L3Res = f"{campaign}_V1_DATA_L2L3Residual_{algo}"
name_Compound = f"{campaign}_V1_DATA_L1L2L3Res_{algo}"

# Fixed constants
rho_val = 30.0
area_val = 0.5
fixed_phi = 0.0  

# Load the evaluator
evaluator = correctionlib.CorrectionSet.from_file(json_file)

# The different cases you want to plot
eta_cases = [0.0, 2.0, 2.6, 3.2]
pt_cases = [30.0, 50.0, 100.0, 500.0, 1000.0]

# =========================================================
# Output Directory
# =========================================================
output_dir = f"/eos/user/t/tchatzis/php-plots/JEC4Prompt/jecs_values_checks/{campaign}"

# Create the folder if it doesn't already exist
os.makedirs(output_dir, exist_ok=True)
print(f"Saving all plots to: ./{output_dir}/")

# Standard colors for consistency across plots
colors = {
    "L1FastJet": "tab:blue",
    "L2Relative": "tab:orange",
    "L3Absolute": "tab:green",
    "L2L3Residual": "tab:red",
    "Total_Compound": "black"
}

# =========================================================
# Corrections vs pT (Loop over Etas)
# =========================================================
print("Generating vs pT Plots...")
pt_range = np.logspace(np.log10(15), np.log10(3000), 200)

for fixed_eta in eta_cases:
    l1_vs_pt, l2_vs_pt, l3abs_vs_pt, l2l3res_vs_pt, comp_vs_pt = [], [], [], [], []

    for raw_pt in pt_range:
        l1 = evaluator[name_L1].evaluate(fixed_eta, raw_pt, rho_val, area_val)
        pt1 = raw_pt * l1
        
        l2 = evaluator[name_L2].evaluate(fixed_eta, fixed_phi, pt1)
        pt2 = pt1 * l2
        
        l3abs = evaluator[name_L3Abs].evaluate(fixed_eta, pt2)
        pt3 = pt2 * l3abs
        
        l2l3res = evaluator[name_L2L3Res].evaluate(fixed_eta, pt3)
        
        comp = evaluator.compound[name_Compound].evaluate(fixed_eta, raw_pt, fixed_phi, rho_val, area_val)
        
        l1_vs_pt.append(l1)
        l2_vs_pt.append(l2)
        l3abs_vs_pt.append(l3abs)
        l2l3res_vs_pt.append(l2l3res)
        comp_vs_pt.append(comp)

    # --- 1. Combined Plot ---
    plt.figure(figsize=(10, 6))
    plt.plot(pt_range, l1_vs_pt, label="L1FastJet", linestyle="--", color=colors["L1FastJet"])
    plt.plot(pt_range, l2_vs_pt, label="L2Relative", linestyle="--", color=colors["L2Relative"])
    plt.plot(pt_range, l3abs_vs_pt, label="L3Absolute", linestyle="--", color=colors["L3Absolute"])
    plt.plot(pt_range, l2l3res_vs_pt, label="L2L3Residual", linestyle="--", color=colors["L2L3Residual"])
    plt.plot(pt_range, comp_vs_pt, label="Total Compound", color=colors["Total_Compound"], linewidth=2)

    plt.xscale("log")
    plt.xlabel("Raw Jet $p_{T}$ [GeV]")
    plt.ylabel("JEC Correction Factor")
    plt.title(f"JEC Breakdown vs pT ($\eta$ = {fixed_eta}, $\phi$ = {fixed_phi}, $\\rho$ = {rho_val}, Area = {area_val})")
    plt.grid(True, which="both", linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    filename = f"JEC_Combined_vs_pT_eta{fixed_eta}.png".replace(".", "p").replace("ppng", ".png")
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close() 

    # --- 2. Individual Plots ---
    components = {
        "L1FastJet": l1_vs_pt, "L2Relative": l2_vs_pt, 
        "L3Absolute": l3abs_vs_pt, "L2L3Residual": l2l3res_vs_pt, 
        "Total_Compound": comp_vs_pt
    }
    
    for label, data in components.items():
        plt.figure(figsize=(10, 6))
        linestyle = "-" if label == "Total_Compound" else "--"
        linewidth = 2 if label == "Total_Compound" else 1.5
        
        plt.plot(pt_range, data, label=label.replace("_", " "), linestyle=linestyle, color=colors[label], linewidth=linewidth)
        plt.xscale("log")
        plt.xlabel("Raw Jet $p_{T}$ [GeV]")
        plt.ylabel(f"{label.replace('_', ' ')} Factor")
        plt.title(f"{label.replace('_', ' ')} vs pT ($\eta$ = {fixed_eta}, $\phi$ = {fixed_phi})")
        plt.grid(True, which="both", linestyle=":", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        
        filename = f"JEC_{label}_vs_pT_eta{fixed_eta}.png".replace(".", "p").replace("ppng", ".png")
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()

# =========================================================
# Corrections vs Eta (Loop over pTs)
# =========================================================
print("Generating vs Eta Plots...")
eta_range = np.linspace(-5.0, 5.0, 5000)

for fixed_pt in pt_cases:
    l1_vs_eta, l2_vs_eta, l3abs_vs_eta, l2l3res_vs_eta, comp_vs_eta = [], [], [], [], []

    for eta in eta_range:
        l1 = evaluator[name_L1].evaluate(eta, fixed_pt, rho_val, area_val)
        pt1 = fixed_pt * l1
        l2 = evaluator[name_L2].evaluate(eta, fixed_phi, pt1)
        pt2 = pt1 * l2
        l3abs = evaluator[name_L3Abs].evaluate(eta, pt2)
        pt3 = pt2 * l3abs
        l2l3res = evaluator[name_L2L3Res].evaluate(eta, pt3)
        comp = evaluator.compound[name_Compound].evaluate(eta, fixed_pt, fixed_phi, rho_val, area_val)
        
        l1_vs_eta.append(l1)
        l2_vs_eta.append(l2)
        l3abs_vs_eta.append(l3abs)
        l2l3res_vs_eta.append(l2l3res)
        comp_vs_eta.append(comp)

    # --- 1. Combined Plot ---
    plt.figure(figsize=(10, 6))
    plt.plot(eta_range, l1_vs_eta, label="L1FastJet", linestyle="--", color=colors["L1FastJet"])
    plt.plot(eta_range, l2_vs_eta, label="L2Relative", linestyle="--", color=colors["L2Relative"])
    plt.plot(eta_range, l3abs_vs_eta, label="L3Absolute", linestyle="--", color=colors["L3Absolute"])
    plt.plot(eta_range, l2l3res_vs_eta, label="L2L3Residual", linestyle="--", color=colors["L2L3Residual"])
    plt.plot(eta_range, comp_vs_eta, label="Total Compound", color=colors["Total_Compound"], linewidth=2)

    plt.xlabel("Jet $\eta$")
    plt.ylabel("Correction Factor")
    plt.title(f"JEC Breakdown vs Eta (Raw pT = {fixed_pt} GeV, $\phi$ = {fixed_phi}, $\\rho$ = {rho_val}, Area = {area_val})")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    filename = f"JEC_Combined_vs_Eta_pT{fixed_pt}.png".replace(".", "p").replace("ppng", ".png")
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()

    # --- 2. Individual Plots ---
    components = {
        "L1FastJet": l1_vs_eta, "L2Relative": l2_vs_eta, 
        "L3Absolute": l3abs_vs_eta, "L2L3Residual": l2l3res_vs_eta, 
        "Total_Compound": comp_vs_eta
    }
    
    for label, data in components.items():
        plt.figure(figsize=(10, 6))
        linestyle = "-" if label == "Total_Compound" else "--"
        linewidth = 2 if label == "Total_Compound" else 1.5
        
        plt.plot(eta_range, data, label=label.replace("_", " "), linestyle=linestyle, color=colors[label], linewidth=linewidth)
        plt.xlabel("Jet $\eta$")
        plt.ylabel(f"{label.replace('_', ' ')} Factor")
        plt.title(f"{label.replace('_', ' ')} vs Eta (Raw pT = {fixed_pt} GeV, $\phi$ = {fixed_phi})")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        
        filename = f"JEC_{label}_vs_Eta_pT{fixed_pt}.png".replace(".", "p").replace("ppng", ".png")
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()

# =========================================================
# Eta Asymmetry (Loop over pTs)
# =========================================================
print("Generating Asymmetry Plots...")

abs_eta_range = np.linspace(0.0, 5.0, 2500)

for fixed_pt in pt_cases:
    l1_asym, l2_asym, l3abs_asym, l2l3res_asym, comp_asym = [], [], [], [], []

    for abs_eta in abs_eta_range:
        # Positive Eta
        l1_pos = evaluator[name_L1].evaluate(abs_eta, fixed_pt, rho_val, area_val)
        pt1_pos = fixed_pt * l1_pos
        l2_pos = evaluator[name_L2].evaluate(abs_eta, fixed_phi, pt1_pos)
        pt2_pos = pt1_pos * l2_pos
        l3abs_pos = evaluator[name_L3Abs].evaluate(abs_eta, pt2_pos)
        pt3_pos = pt2_pos * l3abs_pos
        l2l3res_pos = evaluator[name_L2L3Res].evaluate(abs_eta, pt3_pos)
        comp_pos = evaluator.compound[name_Compound].evaluate(abs_eta, fixed_pt, fixed_phi, rho_val, area_val)
        
        # Negative Eta
        l1_neg = evaluator[name_L1].evaluate(-abs_eta, fixed_pt, rho_val, area_val)
        pt1_neg = fixed_pt * l1_neg
        l2_neg = evaluator[name_L2].evaluate(-abs_eta, fixed_phi, pt1_neg)
        pt2_neg = pt1_neg * l2_neg
        l3abs_neg = evaluator[name_L3Abs].evaluate(-abs_eta, pt2_neg)
        pt3_neg = pt2_neg * l3abs_neg
        l2l3res_neg = evaluator[name_L2L3Res].evaluate(-abs_eta, pt3_neg)
        comp_neg = evaluator.compound[name_Compound].evaluate(-abs_eta, fixed_pt, fixed_phi, rho_val, area_val)
        
        # Ratios
        l1_asym.append(l1_pos / l1_neg if l1_neg != 0 else 1.0)
        l2_asym.append(l2_pos / l2_neg if l2_neg != 0 else 1.0)
        l3abs_asym.append(l3abs_pos / l3abs_neg if l3abs_neg != 0 else 1.0)
        l2l3res_asym.append(l2l3res_pos / l2l3res_neg if l2l3res_neg != 0 else 1.0)
        comp_asym.append(comp_pos / comp_neg if comp_neg != 0 else 1.0)

    # --- 1. Combined Plot ---
    plt.figure(figsize=(10, 6))
    plt.plot(abs_eta_range, l1_asym, label="L1FastJet", linestyle="--", color=colors["L1FastJet"])
    plt.plot(abs_eta_range, l2_asym, label="L2Relative", linestyle="--", color=colors["L2Relative"])
    plt.plot(abs_eta_range, l3abs_asym, label="L3Absolute", linestyle="--", color=colors["L3Absolute"])
    plt.plot(abs_eta_range, l2l3res_asym, label="L2L3Residual", linestyle="--", color=colors["L2L3Residual"])
    plt.plot(abs_eta_range, comp_asym, label="Total Compound", color=colors["Total_Compound"], linewidth=2)

    plt.axhline(1.0, color='gray', linestyle='-', alpha=0.5)
    plt.xlabel("Jet |$\eta$|")
    plt.ylabel("Asymmetry (+$\eta$ / -$\eta$)")
    plt.title(f"JEC Asymmetry vs |$\eta$| (Raw pT = {fixed_pt} GeV, $\phi$ = {fixed_phi})")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    filename = f"JEC_Combined_Asym_pT{fixed_pt}.png".replace(".", "p").replace("ppng", ".png")
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()

    # --- 2. Individual Plots ---
    components = {
        "L1FastJet": l1_asym, "L2Relative": l2_asym, 
        "L3Absolute": l3abs_asym, "L2L3Residual": l2l3res_asym, 
        "Total_Compound": comp_asym
    }
    
    for label, data in components.items():
        plt.figure(figsize=(10, 6))
        linestyle = "-" if label == "Total_Compound" else "--"
        linewidth = 2 if label == "Total_Compound" else 1.5
        
        plt.plot(abs_eta_range, data, label=label.replace("_", " "), linestyle=linestyle, color=colors[label], linewidth=linewidth)
        plt.axhline(1.0, color='gray', linestyle='-', alpha=0.5)
        
        plt.xlabel("Jet |$\eta$|")
        plt.ylabel(f"Asymmetry (+$\eta$ / -$\eta$)")
        plt.title(f"{label.replace('_', ' ')} Asymmetry vs |$\eta$| (Raw pT = {fixed_pt} GeV, $\phi$ = {fixed_phi})")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        
        filename = f"JEC_{label}_Asym_pT{fixed_pt}.png".replace(".", "p").replace("ppng", ".png")
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()

print("All plots generated successfully!")