// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "../runAnalysis.C"
#include "RooMsgService.h"

void Production_Else ( TString fPrefixDT = "", TString fPrefixMC = "", TString fType = "PID" )
{
    // Retrieving PreProcessed data histograms
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/Standard/SignalExtraction/"));
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/Standard/SignalExtrapolation/"));
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/Standard/MassResolution/"));
    gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/Standard/PreProcessing/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/SignalExtraction/*              ./result/yield/Systematics/Standard/SignalExtraction/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/SignalExtrapolation/*           ./result/yield/Systematics/Standard/SignalExtrapolation/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/MassResolution/*                ./result/yield/Systematics/Standard/MassResolution/"));
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/PreProcessing/*                 ./result/yield/Systematics/Standard/PreProcessing/"));
    gROOT->ProcessLine(Form(".! rm -r ./result/yield/PreProcessing/*                 "));
    gROOT->ProcessLine(Form(".! rm -r ./result/yield/SignalExtraction/*              "));
    gROOT->ProcessLine(Form(".! rm -r ./result/yield/SignalExtrapolation/*           "));
    gROOT->ProcessLine(Form(".! rm -r ./result/yield/MassResolution/*                "));
    for ( Int_t iTYP = 0; iTYP < 100; iTYP++ ) {
        if ( fType.Contains("PID") && iTYP > 6 ) break;
        if ( fType.Contains("TRK") && iTYP > 12 ) break;
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/%s/%s_%i/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/%s/%s_%i/SignalExtraction/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/%s/%s_%i/SignalExtrapolation/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/%s/%s_%i/MassResolution/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! mkdir -p ./result/yield/Systematics/%s/%s_%i/PreProcessing/",fType.Data(),fType.Data(),iTYP));
        runAnalysis(Form((fPrefixDT+sPID_DT_Name).Data(),Form("_%s%i",fType.Data(),iTYP)),Form((fPrefixMC+sPID_MC_Name).Data(),Form("_%s%i",fType.Data(),iTYP)),-1,"Yield");
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/PreProcessing/*             ./result/yield/Systematics/%s/%s_%i/PreProcessing/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/SignalExtraction/*          ./result/yield/Systematics/%s/%s_%i/SignalExtraction/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/SignalExtrapolation/*       ./result/yield/Systematics/%s/%s_%i/SignalExtrapolation/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/MassResolution/*            ./result/yield/Systematics/%s/%s_%i/MassResolution/",fType.Data(),fType.Data(),iTYP));
        gROOT->ProcessLine(Form(".! rm -r ./result/yield/PreProcessing/*             "));
        gROOT->ProcessLine(Form(".! rm -r ./result/yield/SignalExtraction/*          "));
        gROOT->ProcessLine(Form(".! rm -r ./result/yield/SignalExtrapolation/*       "));
        gROOT->ProcessLine(Form(".! rm -r ./result/yield/MassResolution/*            "));
    }
    gROOT->ProcessLine(Form(".! cp -p -r ./result/yield/Systematics/Standard/*       ./result/yield/"));
}
