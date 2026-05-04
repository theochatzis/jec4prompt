import os
import re
import argparse
import yaml
import ROOT
import fnmatch
import correctionlib
from tqdm import tqdm  # Progress bar

import time

t0 = time.time()

ROOT.ROOT.EnableImplicitMT()
print("Threads enabled:", ROOT.ROOT.GetThreadPoolSize())

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
args = parser.parse_args()

# ---------- Compile C++ helper functions ----------
# Note : In this way the code will compile it here. You can use pre-compiled functions external .so from cc files or header (.h) files.
script_dir = os.path.dirname(os.path.abspath(__file__))
# Load the header into the ROOT interpreter
ROOT.gInterpreter.ProcessLine(f'#include "{os.path.join(script_dir, "../Common/interface/utils.h")}"')

# An example:
# ROOT.gInterpreter.Declare("""
# float compute_mt(float lep_pt, float lep_phi, float met_pt, float met_phi) {
#     float dphi=TVector2::Phi_mpi_pi(lep_phi - met_phi);
#     return std::sqrt(2.*lep_pt*met_pt*(1. - std::cos(dphi)));
# }                       
# """)

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
    # Probe pt
    df = df.Define("Probe4Vec","GetObject4Vec(Probe_pt, Probe_eta, Probe_phi, Probe_mass)")

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
    
    # ---------- Region loop ----------
    for region_name, region_info in region_config.items():
        # Initialize a region_df having the initial df
        region_df = df 
        print(f"\nApplying cuts for region: {region_name}\n")
        cuts = region_info["cuts"]
        # Apply cuts from regions definitions sequentially 
        for selection in cuts:
            region_df = region_df.Filter(selection, selection)

        output.mkdir(region_name)
        output.cd(region_name)

        hist_pointers = []
        for hist_name, hist_info in hist_config.items():
            title = hist_info["title"]
            bins = hist_info["bins"]
            variable = hist_info["variable"]

            hist = region_df.Histo1D(
                (hist_name, title, bins[0], bins[1], bins[2]),
                variable
            )
            hist_pointers.append(hist)
            
        # Trigger the event loop once and write all histograms
        for hist in hist_pointers:
            hist.Write()
        print('Cuts report:')
        # Print cuts report
        region_df.Report().Print()
        
        output.cd()
    
    output.Close()
    print(f"Output written: {output_path}")

print(f"\nTotal runtime: {time.time() - t0:.2f} seconds")