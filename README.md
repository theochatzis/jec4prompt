# jec4prompt
Jet energy corrections for Prompt reconstruction. This setup is made to be compatible with the outputs produced by base repositories
of `JEC4Prompt`:
- [Analysis](https://gitlab.cern.ch/cms-analysis/jme/jec4prompt/jec4prompt-analysis)
- [Automation](https://gitlab.cern.ch/cms-analysis/jme/jec4prompt/jec4prompt-automation)

All the constants and parameters to be changed manually are set in the JSON file `constants.json`.

# Contents
- [Setup](https://github.com/theochatzis/jec4prompt/tree/main#setup)
- [Deriving the JECs](https://github.com/theochatzis/jec4prompt/tree/main#deriving-the-jecs)
- [Useful links](https://github.com/theochatzis/jec4prompt/tree/main#useful-links)
- [Plotting scripts](https://github.com/theochatzis/jec4prompt/tree/main#plotting-scripts)

# Setup
Clone the repository
```bash
git clone git@github.com:theochatzis/jec4prompt.git 
```
then to run the correction you can do :
```bash
root -l -b -q 'L2L3Res.C(<run(integer)>, <era(string)>, <channel(string)>)'
```

For example for photonjet in 2025G run 398600:
```bash
root -l -b -q 'L2L3Res.C(398600, "2025G", "photonjet")'
```

You can find the constants like input paths, binnings for the channels etc in a single json file `constants.json`.

# Deriving the JECs
For now the $`\gamma`$ + jet channel is only used. Corrections are derived in L2L3Residuals for all etas of probe jet.

By running the `L2L3Res.C` file automatically a `txt` with the `L2L3Residual` corrections is produced. 
```bash
root -l -b -q 'L2L3Res.C(398600, "2025G", "photonjet")'
```
ToDos : 
- Add also to make the JSON automatically.
- Add checks of corrections and closures.

## Making correction `.json` file from `.txt`
To derive the JSON file with the L2L3Resdiual from JEC4Prompt one needs to combine them with L1FastJet, L2Relative and L3Absolute and make a JSON file.
Those files can be found in `offline_txts`. 
For the L2Relative use the same as the MC. For the L1FastJet and L3Absolute dummys can be used (found in `offline_txts/dummy_txts`).

For the files to be compatbile they all need to have the same structure in the end. That means:
```
[Campaign]_[Version]_[DATA/MC]_[CorrectionType]_[Algorithm]
```
So only the `[CorrectionType]` is different per file. 
For example:
```
Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi.txt
Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi.txt
Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi.txt
Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi.txt
```
here:
- Campaign : `Summer24Prompt24`
- Version : `V1` ( always written as V and a number)
- MC/DATA : Here used `MC` for data is `DATA`
- Algorithm : `AK4PFPuppi` is used for PUPPI based AK4 jets.

In order to make a `.json` file that one can use to apply the corrections you can use the `makeJECsJSON.py` script.
The script will automatically take the campaign and version etc from residuals. Also if `MC` is in the name of txts will convert to `DATA`, as it is made for data corrections.

For example can be used as follows:
```bash
python3 makeJECsJSON.py \
    --l1 ./offline_txts/dummy_txts/JEC_L1FastJet_Dummy.txt \
    --l2 ./offline_txts/Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi.txt \
    --l3abs ./offline_txts/dummy_txts/JEC_L3Absolute_Dummy.txt \
    --l2l3res ./txts/Summer24Prompt24JEC4PromptRun398027_V1_DATA_L2L3Residual_AK4PFPuppi.txt \
    --output jsons/Summer24Prompt24JEC4PromptRun398027_JECs.json
```
which takes as inputs the txts and outputs the json file stated as output. You can find an example to run in `makeJSON.sh`.

### Plotting the JSON corrections
If you want to plot the corrections from your JSON file to check everything is correct and make sense you can use the `jec_values_plotter.py` script. This script can plot the correction factors from all levels for the various pt and etas that are given, and asymmetries in eta as well.

Open the script and change the configurable parameters to set the JSON name, values of pT to check, etas etc. :
```python
# =========================================================
# Inputs and parameters
# =========================================================
json_file = "./jsons/Summer24Prompt24JEC4PromptRun398027_JECs.json" # Change to your actual JSON name

# The exact names of the corrections inside your JSON
name_L1 = "Summer24Prompt24JEC4PromptRun398027_V1_DATA_L1FastJet_AK4PFPuppi"
name_L2 = "Summer24Prompt24JEC4PromptRun398027_V1_DATA_L2Relative_AK4PFPuppi"
name_L3Abs = "Summer24Prompt24JEC4PromptRun398027_V1_DATA_L3Absolute_AK4PFPuppi"
name_L2L3Res = "Summer24Prompt24JEC4PromptRun398027_V1_DATA_L2L3Residual_AK4PFPuppi"
name_Compound = "Summer24Prompt24JEC4PromptRun398027_V1_DATA_L1L2L3Res_AK4PFPuppi"

# Fixed constants
rho_val = 30.0
area_val = 0.5
fixed_phi = 0.0  

# Load the evaluator
evaluator = correctionlib.CorrectionSet.from_file(json_file)

# The different cases you want to plot
eta_cases = [0.0, 2.0, 2.6, 3.2]
pt_cases = [30.0, 50.0, 100.0, 500.0, 1000.0]
```

The standard JERC group [framework](https://gitlab.cern.ch/cms-jetmet/jerc2json).

# Descriptions of scripts
- `L2L3Res.C`: Derivation of L2L3Res fit for calibrations. 
- `L2L3ResTestRuns.C`: Makes usage of `L2L3Res.C` and loops the procedure over all runs of one era stored in `jecpcl`.

# Useful links
- [JERC application tutorial](https://gitlab.cern.ch/cms-analysis/jme/jerc-application-tutorial)

# Plotting scripts
## `combine_plots.py`( Combine plots in one page )
Python script which gathers plots from a directory and merged them into arrays of plots. Useful for results presentation.

```bash
python3 combine_plots.py /path/to/input/images/ [options]
```

*Arguments*
- `dir` (Required): The input directory containing the .png plots to combine.
- `-n`, `--name` : The filename prefix for the generated grids. Default is "Grid_Summary_Part".
- `-c`, `--cols` : Number of columns in the output grid. Default is 5.
- `-r`, `--rows` : Number of rows in the output grid. Default is 5.
- `-p`, `--pattern` : File search pattern. Default is "L2L3Res_Fit_eta_*.png".

## CMS TDR Style Macro (`tdrstyle_mod22.C`)

This file is a unified ROOT macro for generating publication-ready plots that comply with the CMS Technical Design Report (TDR) style guidelines.

### Quick Start

Include the macro at the top of your ROOT script:

```cpp
#include "tdrstyle_mod22.C"
```

### Core functions
1. `tdrHist` (Axis Formatting)
Creates an empty dummy histogram to enforce strict axis limits, titles, and TDR fonts before drawing your data.

```cpp
TH1D* h_dummy = tdrHist("name", "Y-Axis Label", y_min, y_max, "X-Axis Label", x_min, x_max);
```

2. `tdrCanvas` (Standard Canvas)
Creates a TDR-formatted canvas, draws the dummy axes, and automatically calls CMS_lumi to place the "CMS Preliminary" and energy labels.
```cpp
TCanvas* c2 = tdrDiCanvas("c2", h_up, h_down, iPeriod, iPos);
```
Note: `square` means set to `kSquare` (true) for a 600x600 canvas, or `kRectangular` (false) for 800x600.

3. `tdrDiCanvas` (Ratio/Sub-pad Canvas)
Creates a vertically split canvas (typically used for Data/MC ratio plots) with properly scaled text and margins for both pads.
```cpp
TCanvas* c2 = tdrDiCanvas("c2", h_up, h_down, iPeriod, iPos);
```

4. `tdrDraw` (Data Drawing)
A wrapper for TGraph and TH1 objects that sets markers, line styles, colors, and draws the object in one line.
```cpp
// Signature: obj, drawOpt, markerStyle, markerColor, lineStyle, lineColor, fillStyle, fillColor
tdrDraw(myGraph, "P", kFullCircle, kRed, kSolid, -1, kNone, 0);
```

5. `tdrLeg` (Legend Creation)
Creates a transparent, borderless legend using the standard CMS font.
```cpp
// Coordinates are Normalized Device Coordinates (NDC)
TLegend* leg = tdrLeg(0.20, 0.70, 0.50, 0.88);
leg->AddEntry(myGraph, "Data", "pe");
```

### Configuration Variables (`iPeriod` and `iPos`)

iPeriod (Center of Mass Energy & Luminosity)
Value,Output
1 -> 7 TeV
2 -> 8 TeV
3 -> 7 TeV + 8 TeV
4 -> 13 TeV
7 -> 7 TeV + 8 TeV + 13 TeV
8 -> 13.6 TeV (Run 3)
12 -> 8 TeV (No lumi text)

iPos (CMS Label Positioning)
Defines where the "CMS" and "Preliminary" (or extraText) labels are placed.

0 -> Out of Frame,"Top-left, physically above the plotting box. Leaves the plot area completely clean."

11 -> Top Left,"Inside the frame, aligned left."

22 -> Top Center,"Inside the frame, centered."

33 -> Top Right,"Inside the frame, aligned right."

### Minimal Example:
```cpp
#include "tdrstyle_mod22.C"

void simple_plot() {
    // 1. Create dummy axes
    TH1D *h_axes = tdrHist("axes", "Response", 0.8, 1.2, "Jet p_{T} (GeV)", 30, 3000);
    
    // 2. Create Canvas (Run 3 Energy, CMS Label out-of-frame)
    TCanvas *c = tdrCanvas("c", h_axes, 8, 0, kSquare);
    gPad->SetLogx();

    // 3. Assume `g_data` is a populated TGraphErrors
    tdrDraw(g_data, "P", kFullCircle, kBlack, kSolid, -1, kNone, 0);

    // 4. Draw Legend
    TLegend *leg = tdrLeg(0.60, 0.75, 0.90, 0.88);
    leg->AddEntry(g_data, "Data 2025G", "pe");
    
    c->SaveAs("my_plot.png");
}
```


