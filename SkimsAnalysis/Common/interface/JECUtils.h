#ifndef JECUTILS_H
#define JECUTILS_H

#include "Math/Vector4D.h"

// Changed std::string to const char* to bypass the ABI border!
void initJEC(const char* filepath, const char* correction_name);
float getJEC(float area, float eta,float phi, float pt, float rho);
ROOT::Math::PtEtaPhiMVector getCorrectedMET(float met_pt, float met_phi, float old_jet_pt, float new_jet_pt, float jet_phi);

#endif
