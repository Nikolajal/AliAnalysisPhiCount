#include "../inc/AliAnalysisPhiPair.h"
#include "./PreProcessing/PreProcessing_MC.C"
#include "./PreProcessing/PreProcessing_Data.C"
// !TODO: All Set!

void PreProcessing ( string fFileNameDT = "", string fFileNameMC = "", TString fOption = "", Int_t nEventsCut = -1.)
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    if ( fFileNameDT == "" )
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
    cout << "[INFO] Starting the Data PreProcessing" << endl;
    PreProcessing_Data(fFileNameDT,fOption,nEventsCut);
    cout << "[INFO] Starting the Monte Carlo PreProcessing" << endl;
    PreProcessing_MC(fFileNameMC,fOption,nEventsCut);
    cout << "[INFO] Finished Pre-Processing" << endl;
    return;
}
