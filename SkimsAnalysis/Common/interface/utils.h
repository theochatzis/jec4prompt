// SkimsAnalysis/Common/interface/utils.h

#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

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

#endif