#include "../inc/AliAnalysisPhiPair.h"
#include "./Analysis/SignalExtraction.C"
#include "./PreProcessing.C"
// !TODO: Can there really be a TODO in a Test Macro?

void
TestMacro
()  {
    TFile*  insFile_DT_Multi            =   new TFile   (Form(kASigExtr_FitCheckRst,"Multiplicity"));
    TFile*  insFile_EF_Multi            =   new TFile   (Form(kAnalysis_MCTruthHist,"Multiplicity"));
    //
    TH1D*   hUtilEventMultiplicity;
    //
    hName                   =   "fQC_Event_Enum_V0M";
    hUtilEventMultiplicity  =   (TH1D*)(insFile_DT_Multi->Get(hName));
    //
    TH1D       *hEvntEff;
    //
    hName       =   "fQC_Event_Enum_FLL";
    hEvntEff    =   (TH1D*)(insFile_DT_Multi->Get(hName));
    //
    auto        kN_Trg          =   (hEvntEff->GetBinContent(kEventCount::kTrigger));
    auto        kN_Vtx          =   (hEvntEff->GetBinContent(kEventCount::kVertex));
    auto        kN_MB           =   (hEvntEff->GetBinContent(kEventCount::kVertex10));
    Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(1./kSignalMiss1D) *(kN_Vtx/kN_Trg);
    //
    fSetAllBins();
    //
    auto kTOTScale = 0.;
    cout << "FULL ENTRIES: " << hUtilEventMultiplicity->GetEntries() << endl;
    cout << "FULL INTEGRAL: " << fEvaluateINELgt0(-1,hUtilEventMultiplicity) << endl;
    cout << "START CYCLE: " << endl;
    for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
        //
        auto kBINMULT = fEvaluateINELgt0(iMlt-1,hUtilEventMultiplicity);
        if ( iMlt > 0 ) kTOTScale += kBINMULT;
        cout << "BIN " << iMlt << " : " << kBINMULT << endl;
        //
    }
    cout << "END CYCLE" << endl;
    cout << "SUM:" << kTOTScale << endl;
    cout << (1.)/(fEvaluateINELgt0(-1,hUtilEventMultiplicity) * kBR ) << endl;
    cout << f1DCorrection << endl;
    
    fStartTimer("test");
    for ( int i = 0; i < 10; i++ )  {
        fPrintLoopTimer("test",i,10000,1);
        sleep(2);
    }
    fStopTimer("test");
}
