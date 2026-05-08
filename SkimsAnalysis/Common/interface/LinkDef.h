// ROOT Dictionary Binder
#ifdef __CINT__

// By default, ignore absolutely everything in the C++ header file. Keep it hidden.
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Expose only these specific functions to the Python/ROOT interpreter 
// i.e. create a Python/Cling binding only for these functions bellow:
#pragma link C++ function initJEC;
#pragma link C++ function getJEC;
#pragma link C++ function getCorrectedMET;

#endif
