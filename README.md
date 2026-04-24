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

# Descriptions of scripts
- `L2L3Res.C`: Derivation of L2L3Res fit for calibrations. 
- `L2L3ResTestRuns.C`: Makes usage of `L2L3Res.C` and loops the procedure over all runs of one era stored in `jecpcl`.
