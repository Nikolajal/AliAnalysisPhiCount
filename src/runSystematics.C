#include "../inc/AliAnalysisPhiPair.h"
#include "./Systematics/SE_Analysis.C"
#include "./Systematics/SE_Production.C"
// !TODO: All Set!

void runSystematics ( )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    gROOT->ProcessLine(".! mkdir -p ./result");
    gROOT->ProcessLine(".! mkdir -p ./result/Syst_SE");
    gROOT->ProcessLine(".! mkdir -p ./result/Syst_PID");
    
    SE_Production();
    SE_Analysis();
    
    //PID_Production();
    //PID_Analysis();
}
