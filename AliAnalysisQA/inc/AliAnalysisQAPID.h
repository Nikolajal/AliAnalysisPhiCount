# ifndef    ALIANALYSISQAPID_h
# define    ALIANALYSISQAPID_h
//
#include "AliAnalysisQAVariables.h"
//
//>->->->->->->->->->->->->->->
//
//  PID Variables
//
//>->->->->->->->->->->->->->->
//
// >-> Detectors
TString     fDetectors[]        =   {"TOF","TPC"};
Double_t    fTOFCut1[]          =   {-3,+3,1,100};
//
// >-> Particles
TString     fParticles[]        =   {"TOF","TPC"};
//
#endif
