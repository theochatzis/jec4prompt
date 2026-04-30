// --------------------------------
// Set of useful utils functions   
// --------------------------------

#include <utility> // For std::pair
#include <iostream>

// function to find min max (including the error bar) of TH1(includes TH1D, TH1F, TProfile)and TGraphErrors  
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