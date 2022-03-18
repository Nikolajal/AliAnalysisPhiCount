// File for 1-Dimensional Analysis:
// !TODO: 1. Move the results folder in the signal extraction folder 
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void
Production_SignalExtraction
 ( TString fOption = "", TString kFolder = "", Bool_t fSilent = true )  {
    //
    //-----------------------------//
    //  Setting general analysis   //
    //-----------------------------//
    //
    //  Verbose option
    if ( fSilent )  {
        gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    //
    //  Option chosing
    if ( !fChooseOption(fOption) ) return;
    //
    //  Generating the binning array--------------------------------------------------------------------------
    fSetAllBins();
    //
    //  Utility variables-------------------------------------------------------------------------------------
    Int_t fTotalCount, fProgrCount;
    //
    if ( kDoYield ) {
        // --- Retrieving PreProcessed Histograms
        TFile*      insFile_Data_YL         =   new TFile   ( Form(kAnalysis_InvMassHist,   (TString("Yield")   +kFolder).Data()) );
        TFile*      insFile_Resl_YL         =   new TFile   ( Form(kMassResolution_Anal,    (TString("Yield")   +kFolder).Data()) );
        //
        auto    h1D_ResolutionReference     =   uLoadHistograms<0,TH1F> ( insFile_Resl_YL, "hRes_RMS_3_1D" );
        auto    h2Db_ResolutionReference    =   uLoadHistograms<0,TH1F> ( insFile_Resl_YL, "hRes_RMS_3_2Db" );
        auto    fHEventCount                =   uLoadHistograms<0,TH1F> ( insFile_Data_YL, "fQC_Event_Enum_FLL" );
        auto    fHEvCountMlt                =   uLoadHistograms<0,TH1F> ( insFile_Data_YL, "fQC_Event_Enum_V0M" );
        auto    h1D_Nrec_PT                 =   uLoadHistograms<1,TH1F> ( insFile_Data_YL, "h1D_Nrec_PT_%i" );
        auto    h1D_Nrec_2Db_PT             =   uLoadHistograms<1,TH1F> ( insFile_Data_YL, "h1D_Nrec_2Db_PT_%i" );
        auto    h2D_Nrec_PT                 =   uLoadHistograms<2,TH2F> ( insFile_Data_YL, "h2D_Nrec_PT_%i_%i" );
        //
        // --- Building output and check plots directory
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,(TString("Yield")+kFolder+TString("/Systematics/")).Data())));
        TFile*  outFileRsltTot  =   new TFile   (( TString(Form(kAnalysis_SigExtr_Dir,(TString("Yield")+kFolder+TString("/Systematics/")).Data())) + TString("/FullResults.root") ).Data(),"recreate");
        //
        // --- Set the print progress utilities
        auto    kCurrent_Timer  =   "Analysis: YIELD -- Signal Extraction Systematics Evaluation 1D";
        fProgrCount = -1;
        fTotalCount =   kSyst_SEX_1D_Options.size();
        fStartTimer(kCurrent_Timer);
        //
        //  --- 1D Systematics Evaluation
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            //  --- Progress print
            fProgrCount++;
            fPrintLoopTimer(kCurrent_Timer,fProgrCount,fTotalCount,1);
            //
            gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,(TString("Yield")+kFolder+TString("/Systematics/")).Data()))+kCurrent_Syst);
            gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtrPlotDirectSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data()))+TString("/1D/"));
            TFile*  outFile_Check   =   new TFile   (Form(kASigExtr_FitChkPltSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data()),"recreate");
            auto    kPlotDirectory  =   TString(Form(kASigExtrPlotDirectSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data()))+TString("/1D/");
            //
            // --- Fit the model 1D
            auto    fFitResults_1DYield    =   FitModel    ( h1D_Nrec_PT, h1D_ResolutionReference, kPlotDirectory, kCurrent_Syst, kCurrent_Syst);
            //
            // --- Save to file
            TFile*  outFile_Result  =   new TFile   (Form(kASigExtr_FitChkRstSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),(kCurrent_Syst).Data(),"1D"),"recreate");
            //
            for ( auto hSave : fFitResults_1DYield )    hSave   ->  Write();
            outFileRsltTot  ->  cd();
            for ( auto hSave : fFitResults_1DYield )    hSave   ->  Write();
            //
            outFile_Result  ->  Close();
            outFile_Check   ->  Close();
        }
        //
        fStopTimer(kCurrent_Timer);
        //
        return;
        kCurrent_Timer  =   "Analysis: YIELD -- Signal Extraction Systematics Evaluation 2D";
        fProgrCount = -1;
        fTotalCount =   kSyst_SEX_1D_Options.size();
        fStartTimer(kCurrent_Timer);
        //
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            //  --- Progress print
            fProgrCount++;
            fPrintLoopTimer(kCurrent_Timer,fProgrCount,fTotalCount,1);
            //
            gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,(TString("Yield")+kFolder+TString("/Systematics/")).Data()))+kCurrent_Syst);
            gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtrPlotDirectSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data()))+TString("/2D/"));
            TFile*  outFile_Check   =   new TFile   (Form(kASigExtr_FitChkPltSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data()),"recreate");
            auto    kPlotDirectory  =   TString(Form(kASigExtrPlotDirectSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data()))+TString("/2D/");
            //
            // --- Fit the model 2D
            std::vector<TH1F*>  f1DCheck;
            auto    fFitResults_2DYield = FitModel( h1D_Nrec_2Db_PT, h2Db_ResolutionReference, h2D_Nrec_PT, f1DCheck, kPlotDirectory, kCurrent_Syst, kCurrent_Syst);
            //
            // --- Save to file
            TFile*  outFile_Result  =   new TFile   (Form(kASigExtr_FitChkRstSY,(TString("Yield")+kFolder+TString("/Systematics/")).Data(),kCurrent_Syst.Data(),"2D"),"recreate");
            //
            for ( auto hSave : f1DCheck )               hSave   ->  Write();
            for ( auto hSave : fFitResults_2DYield )    hSave   ->  Write();
            outFileRsltTot  ->  cd();
            for ( auto hSave : f1DCheck )               hSave   ->  Write();
            for ( auto hSave : fFitResults_2DYield )    hSave   ->  Write();
            //
            outFile_Result  ->  Close();
            outFile_Check   ->  Close();
        }
        //
        fStopTimer(kCurrent_Timer);
        //
        outFileRsltTot      ->  Close();
    }
}
