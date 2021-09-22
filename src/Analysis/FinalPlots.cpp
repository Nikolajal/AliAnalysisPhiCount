// File for final combinations of results:
// !TODO: N/A
#include "../../inc/AliAnalysisPhiPair.h"

void
FinalPlots
() {
    TFile      *fIn_Data    =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield"));
    //
    TGraphMultiErrors      *gYields =   (TGraphMultiErrors*)(fIn_Data->Get("gYieldResult"));
    //
    auto    Y_1_yield   =   gYields->GetPointY(0);
    auto    Y_1_stat_   =   gYields->GetErrorY(0,0);
    auto    Y_1_syst_   =   gYields->GetErrorY(0,1);
    auto    Y_2_yield   =   gYields->GetPointY(1);
    auto    Y_2_stat_   =   gYields->GetErrorY(1,0);
    auto    Y_2_syst_   =   gYields->GetErrorY(1,1);
    //
    cout << Y_1_yield << " +" << Y_1_stat_ << " - "<< Y_1_syst_ << endl;
    cout << Y_2_yield << " +" << Y_2_stat_ << " - "<< Y_2_syst_ << endl;
    //
    cout << " - - - - - - - - " << endl;
    cout << "- gamma_{phi} -" << endl;
    cout << " - - - - - - - - " << endl;
    cout << " - " << endl;
    cout << " - " << fGammaPhiValue(Y_1_yield,Y_2_yield) << " +-" << fGammaPhiError(Y_1_yield,Y_1_stat_,Y_2_yield,Y_2_stat_) << " +-" << fGammaPhiError(Y_1_yield,Y_1_syst_,Y_2_yield,Y_2_syst_) << endl;
    cout << " - " << endl;
    cout << " - - - - - - - - " << endl;
    cout << "- sigma_{phi} -" << endl;
    cout << " - - - - - - - - " << endl;
    cout << " - " << endl;
    cout << " - " << fSigmaPhiValue(Y_1_yield,Y_2_yield) << " +-" << fSigmaPhiError(Y_1_yield,Y_1_stat_,Y_2_yield,Y_2_stat_) << " +-" << fSigmaPhiError(Y_1_yield,Y_1_syst_,Y_2_yield,Y_2_syst_) << endl;
    cout << " - " << endl;
    cout << " - - - - - - - - " << endl;
    cout << "- boh -" << endl;
    cout << " - - - - - - - - " << endl;
    cout << " - " << endl;
    cout << " - " << Y_2_yield/(Y_1_yield*Y_1_yield) << " +-" << fGammaPhiError(Y_1_yield,Y_1_stat_,Y_2_yield,Y_2_stat_) << " +-" << fGammaPhiError(Y_1_yield,Y_1_syst_,Y_2_yield,Y_2_syst_) << endl;
    cout << " - " << endl;
}
