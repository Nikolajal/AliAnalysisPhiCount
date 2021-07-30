#include "../inc/AliAnalysisPhiPair.h"
#include "./Systematics/SE_Analysis.C"
#include "./Systematics/SE_Production.C"
#include "./Systematics/PID_Analysis.C"
#include "./Systematics/PID_Production.C"
// !TODO: All Set!

void runSystematics ( string fFileNameDT = "", string fFileNameMC = "", Int_t nEventsCut = -1., TString fOption = "", TString fFolder = "" )    {
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    gROOT->ProcessLine(".! mkdir -p ./result");
    gROOT->ProcessLine(".! mkdir -p ./result/Syst_SE");
    gROOT->ProcessLine(".! mkdir -p ./result/Syst_PID");
    
    runAnalysis ( fFileNameDT, fFileNameMC, nEventsCut, fOption = "", fFolder = "" );
    SE_Production();
    SE_Analysis();
    PID_Production();
    PID_Analysis();
}
