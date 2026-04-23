# jec4prompt
Jet energy corrections for Prompt reconstruction

All the constants and parameters to be changed manually are set in the JSON file `constants.json`.

# Descriptions of scripts
- `L2L3Res.C`: Derivation of L2L3Res fit for calibrations. 
- `L2L3ResTestRuns.C`: Makes usage of `L2L3Res.C` and loops the procedure over all runs of one era stored in `jecpcl`.
