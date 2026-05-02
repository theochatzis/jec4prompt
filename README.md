# jec4prompt
Jet energy corrections for Prompt reconstruction. This setup is made to be compatible with the outputs produced by base repositories
of `JEC4Prompt`:
- [Analysis](https://gitlab.cern.ch/cms-analysis/jme/jec4prompt/jec4prompt-analysis)
- [Automation](https://gitlab.cern.ch/cms-analysis/jme/jec4prompt/jec4prompt-automation)

All the constants and parameters to be changed manually are set in the JSON file `constants.json`.

# Setup
Clone the repository
```
git clone git@github.com:theochatzis/jec4prompt.git 
```
then to run the correction you can do :
```
root -l -b -q 'L2L3Res.C(<run(integer)>, <era(string)>, <channel(string)>)'
```

For example for photonjet in 2025G run 398600:
```
root -l -b -q 'L2L3Res.C(398600, "2025G", "photonjet")'
```

You can find the constants like input paths, binnings for the channels etc in a single json file `constants.json`.

# Deriving the JECs
```
# Placeholder for instructions
```
# Making correction JSON files from `.txt`
Using the standard JERC group [framework](https://gitlab.cern.ch/cms-jetmet/jerc2json).

```
# Placeholder for instructions
```

# Descriptions of scripts
- `L2L3Res.C`: Derivation of L2L3Res fit for calibrations. 
- `L2L3ResTestRuns.C`: Makes usage of `L2L3Res.C` and loops the procedure over all runs of one era stored in `jecpcl`.


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