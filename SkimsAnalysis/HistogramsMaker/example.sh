python3 run_analysis.py --input-files-dir /eos/user/j/jecpcl/public/jec4prompt/runs/Run2025G/run398027/photonjet \
--output-dir ./testHistos \
--histograms-defs hist_defs.yaml \
--regions-defs regions.yaml \
--file-pattern "*Skim*.root" \
--output-name "test.root"