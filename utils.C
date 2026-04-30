// --------------------------------
// Set of useful utils functions   
// --------------------------------

#include <utility> // For std::pair
#include <iostream>
#include <map>
#include <string>

// Helper function to find min max (including the error bar) of TH1(includes TH1D, TH1F, TProfile)and TGraphErrors  
std::pair<double, double> GetHistMinMaxWithErrors(TObject* obj) {
    // Safety check
    if (!obj) return {0.0, 0.0}; 

    double min_val = 999999.0;
    double max_val = -999999.0;

    // --- CHECK 1: Is it a Histogram or Profile? ---
    // (TH1D, TH1F, and TProfile all inherit from TH1)
    if (TH1* hist = dynamic_cast<TH1*>(obj)) {
        if (hist->GetNbinsX() == 0) return {0.0, 0.0};
        
        for (int i = 1; i <= hist->GetNbinsX(); ++i) {
            if (hist->GetBinContent(i) == 0 && hist->GetBinError(i) == 0) continue;
            double high_edge = hist->GetBinContent(i) + hist->GetBinError(i);
            double low_edge  = hist->GetBinContent(i) - hist->GetBinError(i);
            
            if (high_edge > max_val) max_val = high_edge;
            if (low_edge < min_val)  min_val = low_edge;
        }
    } 
    // --- CHECK 2: Is it a Graph with Errors? ---
    else if (TGraphErrors* graph = dynamic_cast<TGraphErrors*>(obj)) {
        if (graph->GetN() == 0) return {0.0, 0.0};
        
        for (int i = 0; i < graph->GetN(); ++i) {
            double y = graph->GetPointY(i); 
            double ey = graph->GetErrorY(i);
            
            double high_edge = y + ey;
            double low_edge  = y - ey;
            
            if (high_edge > max_val) max_val = high_edge;
            if (low_edge < min_val)  min_val = low_edge;
        }
    } 
    // --- CHECK 3: Unsupported object ---
    else {
        std::cerr << "Warning: Object passed to GetHistMinMaxWithErrors is neither a TH1 nor a TGraphErrors!" << std::endl;
        return {0.0, 0.0};
    }

    return {min_val, max_val};
}


// Helper function to dynamically evaluate JECs from a standard txt payload
double getJEC(const std::string& txt_filename, double jetEta, double jetPt) {
    // Internal struct to hold the parsed bin data
    struct JECLine {
        double eta_min, eta_max;
        std::vector<double> params;
    };
    
    // Internal struct to hold the entire file's payload in memory
    struct JECPayload {
        std::string formula;
        std::vector<JECLine> lines;
        TF1* func = nullptr;
    };
    
    // STATIC cache: Persists in memory across multiple calls to this function
    static std::map<std::string, JECPayload> jec_cache;

    // --- STEP 1: Parse file if it is NOT in the cache yet ---
    if (jec_cache.find(txt_filename) == jec_cache.end()) {
        std::ifstream in_file(txt_filename);
        if (!in_file.is_open()) {
            std::cerr << "Error: Could not open JEC file " << txt_filename << std::endl;
            return 1.0; 
        }

        JECPayload payload;
        std::string line;
        
        // Read Header and extract the formula
        if (std::getline(in_file, line)) {
            // Find the boundaries of the formula in the string
            size_t pos_start = line.find("JetPt ");
            size_t pos_end   = line.find(" Correction");
            
            if (pos_start != std::string::npos && pos_end != std::string::npos) {
                // Extract everything between "JetPt " (offset by 6 chars) and " Correction"
                payload.formula = line.substr(pos_start + 6, pos_end - (pos_start + 6));
            } else {
                std::cerr << "Error: Malformed header. Could not parse formula!" << std::endl;
                return 1.0;
            }
        }

        // Initialize the TF1 with the dynamically read formula
        payload.func = new TF1(Form("func_%s", txt_filename.c_str()), payload.formula.c_str(), 1.0, 15000.0);

        // Read the actual parameters
        while (std::getline(in_file, line)) {
            std::istringstream iss(line);
            double e_min, e_max, x_min, x_max;
            int n_tok;
            
            if (!(iss >> e_min >> e_max >> n_tok >> x_min >> x_max)) continue;
            
            JECLine jl;
            jl.eta_min = e_min;
            jl.eta_max = e_max;
            
            double p_val;
            while (iss >> p_val) {
                jl.params.push_back(p_val);
            }
            payload.lines.push_back(jl);
        }
        in_file.close();
        
        // Save to cache so we never have to read this file again
        jec_cache[txt_filename] = payload;
    }

    // --- STEP 2: Evaluate the JEC from memory ---
    JECPayload& payload = jec_cache[txt_filename];
    
    for (const auto& jl : payload.lines) {
        // Find the matching eta bin
        if (jetEta > jl.eta_min && jetEta < jl.eta_max) {
            if (!jl.params.empty()) {
                // Load parameters into the cached TF1
                for (size_t i = 0; i < jl.params.size(); ++i) {
                    payload.func->SetParameter(i, jl.params[i]);
                }
                // Return the evaluated correction!
                return payload.func->Eval(jetPt);
            }
        }
    }
    
    // Fallback if the requested eta is completely outside the payload bounds
    return 1.0; 
}
