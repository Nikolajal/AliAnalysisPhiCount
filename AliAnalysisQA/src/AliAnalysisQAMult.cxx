#include "../inc/AliAnalysisQAMult.h"
//
void
uMultiplicityEventDistribution  (  )  {
    //
    //  --- Batch Mode
    gROOT->SetBatch(kTRUE);
    //
    //  --- Output Directory
    gROOT->ProcessLine(Form(".! mkdir -p %s/%s/%s",fOutput_Directory.Data(),fOutput_Directory_Mult.Data(),"EVT"));
    //
    //  --- Recover Histogram
    TH1F   *fQC_Event_Mult_Dist_DT  = (TH1F*)fInputList_DT->FindObject(fQC_Event_Enum_Mult);
    TH1F   *fQC_Event_Mult_Dist_MC  = (TH1F*)fInputList_MC->FindObject(fQC_Event_Enum_Mult);
    //
    //  --- Set Histograms
    fSetColorScheme(fQC_Event_Mult_Dist_DT,kTRUE);
    fSetColorScheme(fQC_Event_Mult_Dist_MC,kFALSE);
    fQC_Event_Mult_Dist_DT  ->  Rebin(100);
    fQC_Event_Mult_Dist_MC  ->  Rebin(100);
    fQC_Event_Mult_Dist_DT  ->  Scale( 1./fQC_Event_Mult_Dist_DT->GetEntries() );
    fQC_Event_Mult_Dist_MC  ->  Scale( 1./fQC_Event_Mult_Dist_MC->GetEntries() );
    fQC_Event_Mult_Dist_DT  ->  GetXaxis()->SetRangeUser( 0, 100 );
    fQC_Event_Mult_Dist_DT  ->  SetMinimum  ( 0 );
    fQC_Event_Mult_Dist_DT  ->  SetMaximum  ( 1.25*fQC_Event_Mult_Dist_DT->GetMaximum() );
    //
    //  --- Legend
    TLegend*    lLegend =   new TLegend(0.75,0.25,0.9,0.1);
    lLegend             ->  AddEntry(fQC_Event_Mult_Dist_DT, "Data", "L");
    lLegend             ->  AddEntry(fQC_Event_Mult_Dist_MC, "MC",   "L");
    //
    //  --- Plots
    TCanvas*    cDrawMultDist   =   new TCanvas();
    gStyle                      ->  SetOptStat(0);
    fQC_Event_Mult_Dist_DT      ->  Draw("HIST SAME");
    fQC_Event_Mult_Dist_MC      ->  Draw("HIST SAME");
    lLegend                     ->  Draw("SAME");
    cDrawMultDist               ->  SaveAs(Form("./%s/%s/%s/%s.pdf",fOutput_Directory.Data(),fOutput_Directory_Mult.Data(),"EVT","MultDist"));
    //
    //  --- Batch Mode
    gROOT->SetBatch(kFALSE);
}
