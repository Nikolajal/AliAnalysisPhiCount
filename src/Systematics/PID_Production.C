// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "../runAnalysis.C"
#include "RooMsgService.h"

void PID_Production ( TString fPrefixDT = "", TString fPrefixMC = "" )
{
    // Retrieving PreProcessed data histograms
    gROOT->ProcessLine(Form(".! mkdir -p ./result/Syst_PID_%i",0));
    gROOT->ProcessLine(Form(".! mv ./result/yield/* ./result/Syst_PID_%i",0));
    gROOT->ProcessLine(Form(".! mv ./result/SEFitCheck ./result/Syst_PID_%i",0));
    for ( Int_t iPID = 1; iPID <= nPIDFiles; iPID++ ) {
        gROOT->ProcessLine(Form(".! mkdir -p ./result/Syst_PID_%i",iPID));
        runAnalysis(Form((fPrefixDT+sPID_DT_Name).Data(),iPID),Form((fPrefixMC+sPID_MC_Name).Data(),iPID),-1,"Yield");
        gROOT->ProcessLine(Form(".! mv ./result/yield/* ./result/Syst_PID_%i",iPID));
        gROOT->ProcessLine(Form(".! mv ./result/SEFitCheck ./result/Syst_PID_%i",iPID));
    }
}
