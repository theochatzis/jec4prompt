// Purpose: For testing. Usage of L2L3Res.C which runs on single runs of jec4prompt output and loop over all runs an era.

#include "TFile.h"
#include "TProfile2D.h"

#include "tdrstyle_mod22.C"

// Include function for calibrations in one single run.
#include "L2L3Res.C"

void L2L3ResTestRuns() {
    // SILENT MODE: Prevent GUI popups
    gROOT->SetBatch(kTRUE);
    
    // MEMORY MANAGEMENT: Don't let ROOT globally track every histogram
    TH1::AddDirectory(kFALSE);
    
    // The base directory containing the "runXXXXXX" folders
    TString era = "2025G";
    const char *basePath = Form("/eos/user/j/jecpcl/public/jec4prompt/runs/Run%s/", era.Data());
    
    TSystemDirectory dir("runDir", basePath);
    TList *files = dir.GetListOfFiles();

    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        
        cout << "--- Starting Batch Processing ---" << endl;

        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            
            // Check if the entry is a directory and matches the pattern "runXXXXXX"
            if (file->IsDirectory() && fname.BeginsWith("run")) {
                
                // Extract the numeric part of the run (e.g., "run398801" -> 398801)
                TString runStr = fname;
                runStr.ReplaceAll("run", "");
                
                if (runStr.IsDigit()) {
                    int runNumber = runStr.Atoi();
                    cout << "\n>> Processing Run: " << runNumber << endl;
                    
                    // Call the single-run function
                    L2L3Res(runNumber, era);
                }
            }
        }
    } else {
        cout << "Error: Could not open directory " << basePath << endl;
    }
    
    cout << "\n--- All runs processed ---" << endl;
} // L2L3ResTestRuns