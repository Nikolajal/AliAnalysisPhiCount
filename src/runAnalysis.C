#include "../inc/AliAnalysisPhiPair.h"
#include "./Analysis/Analysis_SignalExtraction.C"
#include "./Analysis/Analysis_SignalCorrections.C"
#include "./Analysis/Analysis_FinalPlots.C"
#include "./PreProcessing.C"
// !TODO: All Set!

void runAnalysis ( string fFileNameDT = "", string fFileNameMC = "", Int_t nEventsCut = -1., TString fOption = "", TString fFolder = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    gROOT->ProcessLine(".! mkdir -p ./result");
    gROOT->ProcessLine(".! mkdir -p ./result/trigger");
    gROOT->ProcessLine(".! mkdir -p ./result/yield");
    gROOT->ProcessLine(".! mkdir -p ./result/multiplicity");
    gROOT->ProcessLine(".! mkdir -p ./result/rapidity");
    gROOT->ProcessLine(".! mkdir -p ./result/tmp");
    gROOT->ProcessLine(".! mkdir -p ./result/SEFitCheck");
    gROOT->ProcessLine(".! mkdir -p ./result/check");
    
    if ( fFileNameMC == "" )
    {
        cout << "[ERROR] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing(\"MonteCarloFile.root\",\"DataFile.root\",\"AnalysisOption\")" << endl;
        return;
    }
    if ( fFileNameMC != "" && fFileNameDT == "" )
    {
        cout << "[WARNING] Data File not specified, will try to use the MC file provided" << endl;
        fFileNameDT = fFileNameMC;
    }

    PreProcessing(fFileNameDT,fFileNameMC,nEventsCut,fOption);
    Analysis_SignalExtraction(true,fOption);
    Analysis_SignalCorrections(true,fOption);
}
