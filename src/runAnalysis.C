#include "../inc/AliAnalysisPhiPair.h"
#include "./Analysis/SignalExtraction.C"
#include "./Analysis/SignalCorrections.C"
#include "./Analysis/MassResolution.C"
#include "./PreProcessing.C"
// !TODO: All Set!

void runAnalysis ( string fFileNameDT = "", string fFileNameMC = "", Int_t nEventsCut = -1., TString fOption = "", TString fFolder = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    if ( fFileNameMC == "" )
    {
        cout << "[ERROR] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing(\"MonteCarloFile.root\",\"DataFile.root\",\"AnalysisOption\",\"TargetFolder\",nEvents)" << endl;
        return;
    }
    if ( fFileNameMC != "" && fFileNameDT == "" )
    {
        cout << "[WARNING] Data File not specified, will try to use the MC file provided" << endl;
        fFileNameDT = fFileNameMC;
    }
    //  Pre-Processing the Trees into histograms
    PreProcessing(fFileNameDT,fFileNameMC,fOption,nEventsCut);
    //  Evaluating the resolution
    MassResolution();
    //  Signal Extraction
    SignalExtraction(true,fOption);
    //  Sepctra correction and integration
    SignalCorrections(true,fOption);
}
