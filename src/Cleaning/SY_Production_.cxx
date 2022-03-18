// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "../runAnalysis.C"
#include "RooMsgService.h"

void SY_Production ( TString fFolderMC = "/Volumes/NRUBINI_DATASTASH/Dataset/_Sim/ARCHIVE/20210805_LHC14j4X/Partial/", TString fFolderDT = "/Volumes/NRUBINI_DATASTASH/Dataset/_Data/ARCHIVE/20210805_LHC10X/Partial/", TString kPrefixDataFile = "/../20210805_LHC10X_STD.root", std::vector<pair<TString,TString>> kDataSet = kpp7TeVDataset, TString fOption = "yield", Int_t nEventsCut = -1., TString kFolder = "_pp_7TeV_old", TString kType = "STD" )
{
    for ( Int_t iTYP = 1; iTYP < 100; iTYP++ ) {
        if ( kType.Contains("PID") && ( iTYP > 6  ) ) break;
        if ( kType.Contains("TRK") && ( iTYP > 13 ) ) break;
        runAnalysis( fFolderMC, fFolderDT, kPrefixDataFile + TString("_") + kType + TString(Form("%i.root",iTYP)), kDataSet, fOption, nEventsCut, kFolder + kType + TString(Form("%i",iTYP)) , kType + TString(Form("%i",iTYP)) );
    }
}

