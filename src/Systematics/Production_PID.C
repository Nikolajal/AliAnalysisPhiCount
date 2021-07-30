// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "../runAnalysis.C"
#include "RooMsgService.h"

void Production_PID ( TString fPrefixDT = "", TString fPrefixMC = "" )
{
    // Retrieving PreProcessed data histograms
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/PIDSystematics/Standard/ExtractionCheck/"));
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/PIDSystematics/Standard/ExtrapolateCheck/"));
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/PIDSystematics/Standard/Util_MassResolution/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/*.root                      ./result/yield/PIDSystematics/Standard/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/ExtractionCheck/*           ./result/yield/PIDSystematics/Standard/ExtractionCheck/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/ExtrapolateCheck/*          ./result/yield/PIDSystematics/Standard/ExtrapolateCheck/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/ExtrapolateCheck/*          ./result/yield/PIDSystematics/Standard/Util_MassResolution/"));
    for ( Int_t iPID = 1; iPID <= nPIDFiles; iPID++ ) {
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/PIDSystematics/PID_%i/ExtractionCheck/",iPID));
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/PIDSystematics/PID_%i/ExtrapolateCheck/",iPID));
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/PIDSystematics/PID_%i/Util_MassResolution/",iPID));
        runAnalysis(Form((fPrefixDT+sPID_DT_Name).Data(),iPID),Form((fPrefixMC+sPID_MC_Name).Data(),iPID),-1,"Yield");
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/*.root                  ./result/yield/PIDSystematics/PID_%i",iPID));
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/ExtractionCheck/*       ./result/yield/PIDSystematics/PID_%i/ExtractionCheck/",iPID));
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/ExtrapolateCheck/*      ./result/yield/PIDSystematics/PID_%i/ExtrapolateCheck/",iPID));
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/ExtrapolateCheck/*      ./result/yield/PIDSystematics/PID_%i/Util_MassResolution/",iPID));
    }
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/PIDSystematics/Standard/    ./result/yield/"));
}
