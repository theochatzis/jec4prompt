import os
import re
import argparse
import yaml
import ROOT
import fnmatch
import numpy as np
import correctionlib
from tqdm import tqdm  # Progress bar

import time

t0 = time.time()

ROOT.ROOT.EnableImplicitMT()
print("Threads enabled:", ROOT.ROOT.GetThreadPoolSize())

def book_histograms(dataframe, config):
    """Books all histograms from the YAML for a given RDataFrame node."""
    pointers = []
    for hist_name, hist_info in config.items():
        title = hist_info["title"]
        hist_type = hist_info.get("type", "TH1D")

        if hist_type == "TH1D":
            variable = hist_info["variable"]
            if "edges" in hist_info:
                if hist_info["edges"][0]=="bins": edges = np.linspace(hist_info["edges"][2], hist_info["edges"][3], hist_info["edges"][1])
                else: edges = np.array(hist_info["edges"], dtype=np.float64)
                model = (hist_name, title, len(edges)-1, edges)
            else:
                bins = hist_info["bins"]
                model = (hist_name, title, bins[0], bins[1], bins[2])
            pointers.append(dataframe.Histo1D(model, variable))

        elif hist_type in ["TProfile", "Profile1D"]:
            var_x, var_y = hist_info["variable_x"], hist_info["variable_y"]
            if "edges" in hist_info:
                if hist_info["edges"][0]=="bins": edges = np.linspace(hist_info["edges"][2], hist_info["edges"][3], hist_info["edges"][1])
                else: edges = np.array(hist_info["edges"], dtype=np.float64)
                model = (hist_name, title, len(edges)-1, edges)
            else:
                bins = hist_info["bins"]
                model = (hist_name, title, bins[0], bins[1], bins[2])
            pointers.append(dataframe.Profile1D(model, var_x, var_y))

        elif hist_type == "TH2D":
            var_x, var_y = hist_info["variable_x"], hist_info["variable_y"]
            if "edges_x" in hist_info and "edges_y" in hist_info:
                if hist_info["edges_x"][0]=="bins": edges_x = np.linspace(hist_info["edges_x"][2], hist_info["edges_x"][3], hist_info["edges_x"][1])
                else: edges_x = np.array(hist_info["edges_x"], dtype=np.float64)
                if hist_info["edges_y"][0]=="bins": edges_y = np.linspace(hist_info["edges_y"][2], hist_info["edges_y"][3], hist_info["edges_y"][1])
                else: edges_y = np.array(hist_info["edges_y"], dtype=np.float64)
                model = (hist_name, title, len(edges_x)-1, edges_x, len(edges_y)-1, edges_y)
            else:
                bins = hist_info["bins"]
                model = (hist_name, title, bins[0], bins[1], bins[2], bins[3], bins[4], bins[5])
            pointers.append(dataframe.Histo2D(model, var_x, var_y))

        elif hist_type in ["TProfile2D", "Profile2D"]:
            var_x, var_y, var_z = hist_info["variable_x"], hist_info["variable_y"], hist_info["variable_z"]
            if "edges_x" in hist_info and "edges_y" in hist_info:
                if hist_info["edges_x"][0]=="bins": edges_x = np.linspace(hist_info["edges_x"][2], hist_info["edges_x"][3], hist_info["edges_x"][1])
                else: edges_x = np.array(hist_info["edges_x"], dtype=np.float64)
                if hist_info["edges_y"][0]=="bins": edges_y = np.linspace(hist_info["edges_y"][2], hist_info["edges_y"][3], hist_info["edges_y"][1])
                else: edges_y = np.array(hist_info["edges_y"], dtype=np.float64)
                model = (hist_name, title, len(edges_x)-1, edges_x, len(edges_y)-1, edges_y)
            else:
                bins = hist_info["bins"]
                model = (hist_name, title, bins[0], bins[1], bins[2], bins[3], bins[4], bins[5])
            pointers.append(dataframe.Profile2D(model, var_x, var_y, var_z))
        else:
            print(f"WARNING: Unknown histogram type '{hist_type}'. Skipping.")
            
    return pointers

# ---------- Argument Parsing ----------
parser = argparse.ArgumentParser(description="RDataFrame Analysis with YAML configs")
parser.add_argument("--input-files-dir", required=True, help="Directory with subdirectories of input ROOT files")
parser.add_argument("--file-pattern", required=False, default="*Skim*.root", help="Pattern to match ROOT file names (e.g., '*Skim*.root')")
parser.add_argument("--output-dir", required=True, help="Directory to save output ROOT files")
parser.add_argument("--output-name", required=False, default="", help="Output ROOT files name")
parser.add_argument("--skip", default="", help="Comma-separated regex to skip subdirectories")
parser.add_argument("--include-only", default="", help="Comma-separated regex to include only specific subdirectories")
parser.add_argument("--histograms-defs", required=True, help="YAML file defining histograms")
parser.add_argument("--regions-defs", required=True, help="YAML file defining selection regions")
parser.add_argument("--tree-name", required=False, default="Events", help="Name of TTree to get from files")
parser.add_argument("--skip-first-nevents", required=False, type=int, default=0, help="Skip first N events from TTree")
parser.add_argument("--max-events", required=False, type=int, default=-1, help="Process only max events entries from TTree")
parser.add_argument("--input-files-depth", required=False, type=int, default=0, help="Subfolder depth to process from --input-files-dir, default is 0 i.e. no subdirectory process")
parser.add_argument("--add-no-selection", required=False, type=bool, default=False, help="Add histograms in the file without any selection")

args = parser.parse_args()
# === C++ HELPERS ===
# ---------- Compile C++ helper functions ----------
# Note : In this way the code will compile it here. You can use pre-compiled functions external .so from cc files or header (.h) files.
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

lib_path = os.path.join(parent_dir, "Common")
header_path = os.path.join(parent_dir, "Common/interface")

# CORRECTIONLIB
# Dynamically locate correctionlib
corr_path = correctionlib.__path__[0]

# Load the C++ shared libraries into ROOT
ROOT.gSystem.Load(os.path.join(corr_path, "lib", "libcorrectionlib.so"))

# Load the compiled library
ROOT.gSystem.Load(os.path.join(lib_path,"libJECUtils.so"))

# Load the header into the ROOT interpreter
ROOT.gInterpreter.ProcessLine(f'#include "{os.path.join(header_path, "JECUtils.h")}"')
ROOT.gInterpreter.ProcessLine(f'#include "{os.path.join(script_dir, "../Common/interface/utils.h")}"')

# --------------------------------------------------
# === CONFIGURATIONS ===
# Initialize the JECs json payload
json_path = "../../jsons/Summer24Prompt24JEC4PromptRun398027_JECs.json" 
payload = "Summer24Prompt24JEC4PromptRun398027_V1_DATA_L2L3Residual_AK4PFPuppi" # Make automatic that it understands which is the L1L2L3Res

cset = correctionlib.CorrectionSet.from_file(json_path)
if 'L1L2L3Res' in payload:
    jec = cset.compound[payload]
else:
    jec = cset[payload]

print(f"\nThe payload {payload} expects these inputs:")
for inp in jec.inputs:
    print(f" - {inp.name} ({inp.type})")

ROOT.initJEC(json_path, payload)

# ---------- Load YAML configuration ----------
with open(args.histograms_defs, "r") as f:
    hist_config = yaml.safe_load(f)
print('loaded histograms definitions YAML...')
with open(args.regions_defs, "r") as f:
    region_config = yaml.safe_load(f)
print('loaded regions YAML...')
# ---------- Compile regex filters ----------
def compile_patterns(raw):
    return [re.compile(p.strip()) for p in raw.split(",") if p.strip()]

skip_patterns = compile_patterns(args.skip)
include_patterns = compile_patterns(args.include_only)

def should_process(subdir):
    if include_patterns and not any(p.search(subdir) for p in include_patterns):
        return False
    if skip_patterns and any(p.search(subdir) for p in skip_patterns):
        return False
    return True

def get_subdirs_at_depth(base_dir, target_depth):
    """
    Return a sorted list of relative subdirectory paths
    that are exactly `target_depth` levels below base_dir.
    Example:
        depth=1 -> Base directory
        depth=1 -> direct children
        depth=2 -> grandchildren, etc.
    """
    if target_depth == 0:
        return ['']
    base_dir = os.path.abspath(base_dir)
    subdirs = []
    base_depth = base_dir.rstrip(os.sep).count(os.sep)

    for root, dirs, _ in os.walk(base_dir):
        current_depth = root.count(os.sep) - base_depth

        # Only collect subdirectories that will be exactly at target_depth
        if current_depth + 1 == target_depth:
            for d in dirs:
                rel_path = os.path.relpath(os.path.join(root, d), base_dir)
                subdirs.append(rel_path)

        # Stop walking deeper once beyond target depth
        if current_depth >= target_depth:
            dirs[:] = []

    return sorted(subdirs)
            
# ---------- Walk subdirectories ----------
# subdirs = [d for d in sorted(os.listdir(args.input_files_dir))
#            if os.path.isdir(os.path.join(args.input_files_dir, d)) and should_process(d)]
subdirs = [d for d in get_subdirs_at_depth(args.input_files_dir, args.input_files_depth)
           if should_process(d)]
print(subdirs)
for subdir in tqdm(subdirs, desc="Processing samples"):
    print("================================================================")
    full_subdir_path = os.path.join(args.input_files_dir, subdir)
    
    if subdir == '':
        subdir = args.input_files_dir.split(os.sep)[-1]
    print(f"Processing {subdir}")
    print(f'full_subdir_path {full_subdir_path}')
    input_files = []
    for root, _, files in os.walk(full_subdir_path):
        for f in files:
            # Use fnmatch to check if the file matches the requested pattern
            if fnmatch.fnmatch(f, args.file_pattern):
                input_files.append(os.path.join(root, f))

    if not input_files:
        print(f"Skipping {subdir}: No ROOT files found.")
        continue

    
    df = ROOT.RDataFrame(args.tree_name, input_files)
    
    # If you want to process a fraction of events
    if args.max_events > 0:
        df = df.Range(args.skip_first_nevents, args.max_events)  # skip first [skip_first_nevents], take next [max_events]

    # ---------- Define derived variables ----------    
    # Tag
    df = df.Define("Tag_PolarVec" , "ROOT::Math::Polar2DVector(Tag_pt, Tag_phi)")

    # Probe 
    df = df.Define("ProbeMC_PolarVec", "ROOT::Math::Polar2DVector(Probe_mcPt, Probe_phi)")
    df = df.Define("Probe_PolarVec", "ROOT::Math::Polar2DVector(Probe_pt, Probe_phi)")
    df = df.Define("Probe_jec", "getJEC(Probe_area, Probe_eta, Probe_phi, Probe_pt, Rho_fixedGridRhoFastjetAll)")
    df = df.Define("Probe_corPt" , "Probe_jec*Probe_pt")
    df = df.Define("Probe_corDB", "Probe_corPt/Tag_pt")
    df = df.Define("Probe_corPolarVec", "ROOT::Math::Polar2DVector(Probe_corPt, Probe_phi)")

    # Re-Apply correction of Probe Jet to MET
    df = df.Define("MET_polarVec", "ROOT::Math::Polar2DVector(T1MET_mc_pt, T1MET_mc_phi)")
    df = df.Define("corMETvec", "getCorrectedMET(T1MET_pt, PuppiMET_phi, Probe_pt, Probe_corPt, Probe_phi)")
    df = df.Define("T1MET_corPt", "corMETvec.Pt()")
    df = df.Define("T1MET_corPhi", "corMETvec.Phi()")
    df = df.Define("T1MET_corPolarVec", "ROOT::Math::Polar2DVector(T1MET_corPt, T1MET_corPhi)")
    

    # Make corrected MPF
    df = df.Define("corMPF", "1 + T1MET_corPolarVec.Dot(Tag_PolarVec)/Tag_PolarVec.Mag2()")
    
    # MC based definitions
    df = df.Define("MPF_mc", "1 + MET_polarVec.Dot(Tag_PolarVec)/Tag_PolarVec.Mag2()")
    df = df.Define("DB_mc", "ProbeMC_PolarVec.R()/Tag_PolarVec.R()")
    
    # Define HDM inputs

    #### Temporary fix for JetActivity: ## Note: Here JetActivity also has residuals
    df = df.Define("JetActivity_PolarVec", "ROOT::Math::Polar2DVector(JetActivity_pt, JetActivity_phi)")
    # Define unclustered component -> add to MET all the jets
    df = df.Define("Unclustered_PolarVec", "MET_polarVec + JetActivity_PolarVec") 
    df = df.Redefine("JetActivity_PolarVec", "JetActivity_PolarVec + Probe_PolarVec")
    ####

    df = (
    df
    .Define(
        "HDM_r0",
        "hdm_r0(Tag_pt, Tag_phi, MET_polarVec.R(), MET_polarVec.Phi())"
    )
    .Define(
        "HDM_r1",
        "hdm_r1(Tag_pt, Tag_phi, ProbeMC_PolarVec.R(), ProbeMC_PolarVec.Phi())"
    )
    .Define(
        "HDM_rn",
        "hdm_rn_from_scalar(Tag_pt, Tag_phi, JetActivity_PolarVec.R(), JetActivity_PolarVec.Phi())"
        #"hdm_rn_from_scalar(Tag_pt, Tag_phi, JetActivity_pt, JetActivity_phi)"
    )
    .Define(
        "HDM_ru",
        "hdm_rn_from_scalar(Tag_pt, Tag_phi, Unclustered_PolarVec.R(), Unclustered_PolarVec.Phi())"
    )
    .Define(
        "HDM_MPD_diff",
        "hdm_closure(HDM_r0, HDM_r1, HDM_rn, HDM_ru)"
    )
    )

    # ---------- Create output file ----------
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Get the output path
    if args.input_files_depth == 0:
        if args.output_name != "":
            output_path = os.path.join(args.output_dir, f"{args.output_name}")
        else:
            output_path = os.path.join(args.output_dir, f"{subdir}.root")
    else: # Get the first-level directory relative to input_files_dir which is the process.
        rel_path = os.path.relpath(full_subdir_path, args.input_files_dir)
        first_level_dir = rel_path.split(os.sep)[0]  # Get first subdirectory in the relative path ,os.sep is separator
        if args.output_name != "":
            output_path = os.path.join(args.output_dir, f"{args.output_name}")
        else:
            output_path = os.path.join(args.output_dir, f"{first_level_dir}.root")

    output = ROOT.TFile(output_path, "RECREATE")
    
    # Master list to collect all execution nodes
    all_hist_pointers = []
    reports = {}
    # --- Optional :No selection histograms
    if args.add_no_selection:
        output.cd() # Ensure we are at the top level of the ROOT file
        print("\nBooking no selection histograms (No Cuts, Top Directory)")
        
        # Book histograms directly on the base 'df' with no filters
        baseline_pointers = book_histograms(df, hist_config)
        # Store as a tuple: (Target Directory, Histogram Pointer)
        all_hist_pointers.extend([("", ptr) for ptr in baseline_pointers])
    
    # ---------- Region loop ----------
    for region_name, region_info in region_config.items():
        # Initialize a region_df having the initial df
        region_df = df 
        print(f"\nBooking histograms for region: {region_name}\n")
        cuts = region_info.get("cuts", []) # Doing it like this so if we want we can completely skip the cuts
        # Apply cuts from regions definitions sequentially 
        for selection in cuts:
            region_df = region_df.Filter(selection, selection)
        
        # Save the report pointer to print later
        reports[region_name] = region_df.Report()

        output.mkdir(region_name)
        output.cd(region_name)
        
        # Book histograms on the filtered region_df
        region_pointers = book_histograms(region_df, hist_config)
        # Store as a tuple: (Target Directory, Histogram Pointer)
        all_hist_pointers.extend([(region_name, ptr) for ptr in region_pointers])
        
    print("\nExecuting RDataFrame Graph...")
    
    for target_dir, hist in all_hist_pointers:
        # Change directory before writing
        if target_dir == "":
            output.cd()             # Go to the top level for baseline
        else:
            output.cd(target_dir)   # Go to the specific region folder
            
        hist.Write() # Evaluates the whole graph at once on the first call
        
    for reg, rep in reports.items():
        print(f"\n--- Cuts report for {reg} ---")
        rep.Print()
        
    output.Close()
    print(f"\nOutput written: {output_path}")

print(f"\nTotal runtime: {time.time() - t0:.2f} seconds")