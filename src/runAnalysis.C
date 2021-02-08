#include "../inc/AliAnalysisPhiPair.h"
#include "./Analysis/Analysis_SignalExtraction.C"
#include "./Analysis/Analysis_SignalCorrections.C"
#include "./Analysis/Analysis_FinalPlots.C"
#include "./PreProcessing.C"
// !TODO: All Set!

void runAnalysis ( string fFileNameDT = "", string fFileNameMC = "", Int_t nEventsCut = -1., string fOption = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
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

    PreProcessing(fFileNameDT,fFileNameMC,nEventsCut);
    Analysis_SignalExtraction();
}
