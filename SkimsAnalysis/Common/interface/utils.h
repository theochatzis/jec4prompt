// SkimsAnalysis/Common/interface/utils.h

#ifndef UTILS_H
#define UTILS_H

#include "TLorentzVector.h"
#include <cmath>
#include <vector>
#include "ROOT/RVec.hxx"  // Needed for RVec

using namespace ROOT::VecOps;


//  ============ OBJECTS FUNCTIONS ============

// Generic
/*
Functions to create TLorentzVectors from kinematic variables.
*/
TLorentzVector GetObject4Vec(float pt, float eta, float phi, float mass){
    TLorentzVector p4_object;
    p4_object.SetPtEtaPhiM(pt,eta,phi,mass);
    return p4_object;
}

TLorentzVector GetObject4VecNoMass(float pt, float eta, float phi){ // application for leptons
    TLorentzVector p4_object;
    p4_object.SetPtEtaPhiM(pt,eta,phi,0.);
    return p4_object;
}

TLorentzVector GetObject4VecTransverse(float pt, float phi){ // application for MET
    TLorentzVector p4_object;
    p4_object.SetPtEtaPhiM(pt,0.,phi,0.);
    return p4_object;
}

// HDM Helpers
/*
Functions to define components needed for  Hadron Decomposition Method (HDM)
*/

// Project a pT,phi vector onto the tag axis.
// vec (dot) n_tag = pT * cos(phi - tag_phi)
double project_on_tag_axis(double pt, double phi, double tag_phi) {
    return pt * std::cos(phi - tag_phi);
}

// Sum a collection of jets vectorially, then project onto tag axis.
double project_jet_activity_on_tag_axis(
    const RVec<float>& pts,
    const RVec<float>& phis,
    double tag_phi
) {
    double projection = 0.0;

    for (std::size_t i = 0; i < pts.size(); ++i) {
        projection += pts[i] * std::cos(phis[i] - tag_phi);
    }

    return projection;
}

// r0 = MPF-like response
double hdm_r0(
    double tag_pt,
    double tag_phi,
    double met_pt,
    double met_phi
) {
    if (tag_pt <= 0.0) return -999.0;

    const double met_parallel = project_on_tag_axis(met_pt, met_phi, tag_phi);

    return 1.0 + met_parallel / tag_pt;
}

// r1 = leading tag-probe / DB-like term
double hdm_r1(
    double tag_pt,
    double tag_phi,
    double probe_pt,
    double probe_phi
) {
    if (tag_pt <= 0.0) return -999.0;

    const double probe_parallel = project_on_tag_axis(probe_pt, probe_phi, tag_phi);

    return -probe_parallel / tag_pt;
}

// rn = additional-jet / jet-activity term
double hdm_rn_from_collection(
    double tag_pt,
    double tag_phi,
    const RVec<float>& jet_activity_pt,
    const RVec<float>& jet_activity_phi
) {
    if (tag_pt <= 0.0) return -999.0;

    const double jet_activity_parallel =
        project_jet_activity_on_tag_axis(jet_activity_pt, jet_activity_phi, tag_phi);

    return -jet_activity_parallel / tag_pt;
}

// If JetActivity is already stored as a single vector-summed object
double hdm_rn_from_scalar(
    double tag_pt,
    double tag_phi,
    double jet_activity_pt,
    double jet_activity_phi
) {
    if (tag_pt <= 0.0) return -999.0;

    const double jet_activity_parallel =
        project_on_tag_axis(jet_activity_pt, jet_activity_phi, tag_phi);

    return -jet_activity_parallel / tag_pt;
}

// Closure definition:
// r0 = r1 + rn + ru
double hdm_ru_closure(double r0, double r1, double rn) {
    return r0 - r1 - rn;
}

// Closure diagnostic, should be zero by construction
double hdm_closure(double r0, double r1, double rn, double ru) {
    return r0 - r1 - rn - ru;
}

#endif
