#include "JECUtils.h"
#include "correction.h"
#include <cmath>
#include <memory>
#include <string>
#include <iostream>

std::unique_ptr<correction::CorrectionSet> global_cset;

// Pointers for standard and compound correction cpp types.
correction::Correction::Ref global_jec = nullptr;
correction::CompoundCorrection::Ref global_compound_jec = nullptr;

void initJEC(const char* filepath, const char* correction_name) {
    std::cout << "JEC initialized from: " << filepath << std::endl;
    std::cout << "key used: " << correction_name << std::endl;
    
    global_cset = correction::CorrectionSet::from_file(std::string(filepath));
    
    // Route the payload to the correct pointer type
    if (std::string(correction_name).find("L1L2L3Res") != std::string::npos) {
        global_compound_jec = global_cset->compound().at(std::string(correction_name));
    } else {
        global_jec = global_cset->at(std::string(correction_name));
    }
}

float getJEC(float area, float eta, float phi, float pt , float rho) {
    // Create an empty vector with correctionlib inputs type to hold the dynamically ordered arguments.
    std::vector<correction::Variable::Type> args;
    
    // Ask the correct pointer for its list of required inputs
    auto required_inputs = (global_compound_jec != nullptr) ? 
                            global_compound_jec->inputs() : 
                            global_jec->inputs();
    

    // Loop through the required inputs and map them to your C++ variables
    for (const auto& input : required_inputs) {
        std::string name = input.name();

        // Note: You must ensure these string names exactly match what your JSON uses!
        // Standard CMS conventions are usually "JetEta", "JetPt", "Rho", "JetA", etc.
        if (name == "JetEta" || name == "eta") {
            args.push_back(static_cast<float>(eta));
        } 
        else if (name == "JetPt" || name == "pt") {
            args.push_back(static_cast<float>(pt));
        } 
        else if (name == "JetA" || name == "area") {
            args.push_back(static_cast<float>(area));
        } 
        else if (name == "Rho" || name == "rho") {
            args.push_back(static_cast<float>(rho));
        } 
        else if (name == "JetPhi" || name == "phi") {
            args.push_back(static_cast<float>(phi));
        } 
        else {
            std::cerr << "[CRITICAL ERROR]: Unknown JEC input parameter requested by JSON: " << name << std::endl;
            return 1.0; // Return uncorrected scale factor to avoid crashing, or throw an error.
        }
    }

    // Evaluate
    if (global_compound_jec != nullptr) {
        return global_compound_jec->evaluate(args);
    } else {
        return global_jec->evaluate(args); 
    }
}

ROOT::Math::PtEtaPhiMVector getCorrectedMET(float met_pt, float met_phi, 
                                            float old_jet_pt, float new_jet_pt, float jet_phi) {
    float met_px = met_pt * std::cos(met_phi);
    float met_py = met_pt * std::sin(met_phi);
    
    float dpt = new_jet_pt - old_jet_pt;
    float dpx = dpt * std::cos(jet_phi);
    float dpy = dpt * std::sin(jet_phi);
    
    float new_met_px = met_px - dpx;
    float new_met_py = met_py - dpy;
    
    float new_met_pt = std::hypot(new_met_px, new_met_py);
    float new_met_phi = std::atan2(new_met_py, new_met_px);
    
    return ROOT::Math::PtEtaPhiMVector(new_met_pt, 0.0, new_met_phi, 0.0);
}