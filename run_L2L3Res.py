#!/usr/bin/env python3
import argparse
import sys
import os

def main():
    # Set up the argument parser to match the arguments of C++ function
    parser = argparse.ArgumentParser(description="Run L2L3Res JEC derivation")
    
    # Matching the arguments from void L2L3Res(...)
    parser.add_argument("--run", type=int, default=398027, help="Run number")
    parser.add_argument("--basePath", type=str, default="2025G", help="Base path/era")
    parser.add_argument("--channel", type=str, default="photonjet", help="Analysis channel")
    parser.add_argument("--outputBaseDir", type=str, default="", help="Output directory")
    parser.add_argument("--jsonWithLumis", type=str, default="", help="Path to lumis JSON")
    parser.add_argument("--runsDirBase", type=str, default="", help="Base directory for runs")

    # Booleans in argparse are best handled with action='store_true' or 'store_false'
    parser.add_argument("--ignore_min_lumi", action="store_false", dest="use_min_lumi", 
                        help="Pass this flag to set use_minimum_luminosity_=false")
    
    parser.add_argument("--min_lumi", type=float, default=-1.0, help="Minimum luminosity")
    parser.add_argument("--l1_txt", type=str, default="", help="L1 JEC txt path")
    parser.add_argument("--l2_txt", type=str, default="", help="L2 JEC txt path")
    parser.add_argument("--l3abs_txt", type=str, default="", help="L3 JEC txt path")
    parser.add_argument("--outputJson", type=str, default="", help="Output JSON path")

    args = parser.parse_args()

    print("Executing via PyROOT...")
    import ROOT
    
    # Run in batch mode (no graphics)
    ROOT.gROOT.SetBatch(True)
    
    # Load the C++ macro
    # Ensure "L2L3Res.C" is in the same directory, or provide the full path
    ROOT.gSystem.Load("libTree") # Load standard libraries if needed
    ROOT.gROOT.LoadMacro("L2L3Res.C") 

    # Call the C++ function exactly as if it were a Python function
    # PyROOT automatically handles the conversion of Python strings to C++ TString/std::string
    ROOT.L2L3Res(
        args.run,
        args.basePath,
        args.channel,
        args.outputBaseDir,
        args.jsonWithLumis,
        args.runsDirBase,
        args.use_min_lumi,
        args.min_lumi,
        args.l1_txt,
        args.l2_txt,
        args.l3abs_txt,
        args.outputJson
    )
if __name__ == "__main__":
    main()