#include "../inc/AliAnalysisPhiPair.h"
#include "./PreProcessing/PP_Data.cxx"
#include "./PreProcessing/PP_MC.cxx"
//#include "./PreProcessing/PP_Resl.cxx"
#include "./Analysis/AN_SigExtraction.cxx"
#include "./Analysis/AN_SigCorrections.cxx"
/*
 .x src/runAnalysis.C("/Users/nikolajal/alice/MC/","/Users/nikolajal/alice/DATA/", "LHC10_FL/LHC10_FL_STD.root",{{"LHC10_FL","LHC14j4_FL"}}, "yield", -1, "_FL")
 
 */
// !TODO: All Set!

void runAnalysis ( TString fFolderMC = "/Volumes/\[HD\]\[Nikolajal\]_Toshiba\ 2/Dataset/_Sim/20210805_LHC14j4X/Partial/", TString fFolderDT = "/Volumes/\[HD\]\[Nikolajal\]_Toshiba\ 2/Dataset/_Data/20210805_LHC10X/Partial/", TString kDataFile = "", std::vector<pair<TString,TString>> kDataSet = kpp7TeVDataset, TString fOption = "yield", Int_t nEventsCut = -1., TString kFolder = "" ) {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  Data PreProcessing
    PP_Data             ( fFolderDT + kDataFile, fOption, nEventsCut, kFolder );
    PP_MC               ( fFolderDT, fFolderMC, kDataSet, fOption, nEventsCut, kFolder );
    PP_Resl             ( fOption, kFolder, true );
    AN_SigExtraction    ( fOption, kFolder, true );
    //AN_SigCorrections   ( fOption, kFolder, true );
}

/*
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
    MassResolution(fOption);
    //  Signal Extraction
    SignalExtraction(fOption,true);
    //  Sepctra correction and integration
    SignalCorrections(fOption,true);
    //FinalPlots();
}
*/
