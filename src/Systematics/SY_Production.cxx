// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "../runAnalysis.C"
#include "SY_Production/SY_PR_SigExtraction.cxx"
#include "RooMsgService.h"

void
SY_Production
 ( TString kType = "SEX", TString fFolderMC = "/Volumes/NRUBINI_DATASTASH/Dataset/_Sim/ARCHIVE/20210805_LHC14j4X/Partial/", TString fFolderDT = "/Volumes/NRUBINI_DATASTASH/Dataset/_Data/ARCHIVE/20210805_LHC10X/Partial/", TString kPrefixDataFile = "/../20210805_LHC10X_STD.root", std::vector<pair<TString,TString>> kDataSet = kpp7TeVDataset, TString fOption = "all", Int_t nEventsCut = -1., TString kFolder = "_p_p__5TeV" )    {
    //
    //  --- Signal Extraction Production
    if ( kType.Contains("SEX") )    { SY_PR_SigExtraction( fOption, kFolder, true ); return; }
    //
    //  --- Else Production
    auto iTer = 0;
    while ( true )  { iTer++; runAnalysis( fFolderMC, fFolderDT, kPrefixDataFile+TString("_")+kType+TString(Form("%i.root",iTer)), kDataSet, fOption, nEventsCut, kFolder + TString("/") + kType + TString(Form("%i",iTer)), kType+TString(Form("%i",iTer)) ); }
}

