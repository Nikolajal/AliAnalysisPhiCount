//  Signal Extraction Systematic calculation
#include "../../inc/AliAnalysisPhiPair.h"

void
SY_SigExtraction
 ( TString fOption = "all", Bool_t fSilent = true, TString kFolder = "" )  {
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
        TFile*      insFile_Data_YL         =   new TFile   ( Form(kAnalysis_InvMassHist,   (TString("Yield")       +kFolder).Data()) );
        TFile*      insFile_Resl_YL         =   new TFile   ( Form(kMassResolution_Anal,    (TString("Yield")       +kFolder).Data()) );
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
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Yield/Systematics/")));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Yield/Systematics/"))+TString("1D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Yield/Systematics/"))+TString("2D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield/Systematics/"))+TString("1D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield/Systematics/"))+TString("2D/"));
        //
        //  --- Save FitResults
        std::vector<std::vector<TH1F*>> fFitResults_1DYield_Array;
        // --- Set the print progress utilities
        fTotalCount =   kSyst_SEX_1D_Options.size();
        fProgrCount =   0;
        fStartTimer("Signal Extraction Systematics Production 1D");
        //
        for ( auto kSystFit = 0; kSystFit < fTotalCount; kSystFit++ ) {
            auto    kCurrent_Option     =   kSyst_SEX_1D_Options.at(kSystFit);
            //
            // --- Building output and check plots directory
            gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield/Systematics/"))+TString("1D/") + kCurrent_Option);
            TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitChkPltSY,"Yield/Systematics/","1D",kCurrent_Option.Data()),"recreate");
            //
            // --- Fit the model 1D
            fFitResults_1DYield_Array   .push_back ( FitModel    ( h1D_Nrec_PT, h1D_ResolutionReference, Form( kASigExtr_Plot_Direct, "Yield/Systematics/" ) + TString( "1D/" ) + kCurrent_Option, kCurrent_Option, kCurrent_Option ) );
            //
            // --- Progressive Count
            fProgrCount++;
            fPrintLoopTimer("Signal Extraction Systematics Production 1D",fProgrCount,fTotalCount,1);
            //
            // --- Save to file
            TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitChkRstSY,"Yield/Systematics/","1D",kCurrent_Option.Data()),"recreate");
            //
            fHEventCount->Write();
            for ( auto hSave : fFitResults_1DYield_Array.at(kSystFit) )    hSave   ->  Write();
            //
            outFile_Check   ->  Close();
            outFile_Result  ->  Close();
        }
        fStopTimer("Signal Extraction Systematics Production 1D");
        //
        //  --- Save FitResults
        std::vector<std::vector<TH1F*>> f1DCheck_Array;
        std::vector<std::vector<TH2F*>> fFitResults_2DYield_Array;
        //
        // --- Set the print progress utilities
        fTotalCount =   kSyst_SEX_2D_Options.size();
        fProgrCount =   0;
        fStartTimer("Signal Extraction Systematics Production 2D");
        //
        for ( auto kSystFit = 0; kSystFit < fTotalCount; kSystFit++ ) {
            auto    kCurrent_Option     =   kSyst_SEX_2D_Options.at(kSystFit);
            //
            // --- Building output and check plots directory
            gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield/Systematics/"))+TString("2D/") + kCurrent_Option );
            TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitChkPltSY,"Yield/Systematics/","2D",kCurrent_Option.Data()),"recreate");
            //
            // --- Fit the model 2D
            std::vector<TH1F*> kUtility_Check1D;
            fFitResults_2DYield_Array   .push_back  ( FitModel   ( h1D_Nrec_2Db_PT, h2Db_ResolutionReference, h2D_Nrec_PT, kUtility_Check1D, Form( kASigExtr_Plot_Direct, "Yield/Systematics/" ) + TString( "2D/" ) + kCurrent_Option, kCurrent_Option, kCurrent_Option ) );
            f1DCheck_Array              .push_back  ( kUtility_Check1D );
            //
            // --- Progressive Count
            fProgrCount++;
            fPrintLoopTimer("Signal Extraction Systematics Production 2D",fProgrCount,fTotalCount,1);
            //
            // --- Save to file
            TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitChkRstSY,"Yield/Systematics/","2D",kCurrent_Option.Data()),"recreate");
            //
            fHEventCount->Write();
            for ( auto hSave : f1DCheck_Array.at(kSystFit) )            hSave   ->  Write();
            for ( auto hSave : fFitResults_2DYield_Array.at(kSystFit) ) hSave   ->  Write();
            //
            outFile_Check   ->  Close();
            outFile_Result  ->  Close();
        }
        fStopTimer("Signal Extraction Systematics Production 2D");
        //
    }
    if ( kDoMultiplicity ) {
        // --- Retrieving PreProcessed Histograms
        TFile*      insFile_Data_ML         =   new TFile   ( Form(kAnalysis_InvMassHist,   (TString("Multiplicity")    +kFolder).Data()) );
        TFile*      insFile_Resl_ML         =   new TFile   ( Form(kMassResolution_Anal,    (TString("Multiplicity")    +kFolder).Data()) );
        //
        auto    h1D_ResolutionReference     =   uLoadHistograms<0,TH1F> ( insFile_Resl_ML, "hRes_RMS_3_1D" );
        auto    h2Db_ResolutionReference    =   uLoadHistograms<0,TH1F> ( insFile_Resl_ML, "hRes_RMS_3_2Db" );
        auto    fHEventCount                =   uLoadHistograms<0,TH1F> ( insFile_Data_ML, "fQC_Event_Enum_FLL" );
        auto    fHEvCountMlt                =   uLoadHistograms<0,TH1F> ( insFile_Data_ML, "fQC_Event_Enum_V0M" );
        auto    h1D_Nrec_PT                 =   uLoadHistograms<1,TH1F> ( insFile_Data_ML, "h1D_Nrec_PT_%i" );
        auto    h1D_Nrec_2Db_PT             =   uLoadHistograms<1,TH1F> ( insFile_Data_ML, "h1D_Nrec_2Db_PT_%i" );
        auto    h2D_Nrec_PT                 =   uLoadHistograms<2,TH2F> ( insFile_Data_ML, "h2D_Nrec_PT_%i_%i" );
        //
        // --- Building output and check plots directory
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Multiplicity/Systematics/")));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Multiplicity/Systematics/"))+TString("1D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Multiplicity/Systematics/"))+TString("2D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Multiplicity/Systematics/"))+TString("1D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Multiplicity/Systematics/"))+TString("2D/"));
        //
        //  --- Save FitResults
        std::vector<std::vector<TH1F*>> fFitResults_1DYield_Array;
        // --- Set the print progress utilities
        fTotalCount =   kSyst_SEX_1D_Options.size()*nBinMult;
        fProgrCount =   0;
        fStartTimer("Signal Extraction Systematics Production 1D");
        //
        /*
         auto    iMult   =   -1;
         std::vector<std::vector<TH1F*>> fFitResults_1DYield_Array;
         for ( auto kTarget : h1D_Nrec_MT_PT )   {
             iMult++;
             fFitResults_1DYield_Array.push_back( FitModel    ( kTarget, h1D_ResolutionReference, Form( kASigExtr_Plot_Direct, "Multiplicity" ) + TString( "1D/" ), Form("MT_%i",iMult) ) );
             //
             // --- Progressive Count
             fProgrCount +=  nBinPT1D;
             fPrintLoopTimer("Multiplicity Analysis Signal Extraction",fProgrCount,fTotalCount,1);
         }
         */
        auto    iMult   =   -1;
        for ( Int_t iMlt = 0; iMlt < nBinMult; iMlt++ ) {
            iMult++;
            for ( auto kSystFit = 0; kSystFit < kSyst_SEX_1D_Options.size(); kSystFit++ ) {
                auto    kCurrent_Option     =   kSyst_SEX_1D_Options.at(kSystFit);
                //
                // --- Building output and check plots directory
                gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield/Systematics/"))+TString("1D/") + kCurrent_Option);
                TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitChkPltSY,"Yield/Systematics/","1D",kCurrent_Option.Data()),"recreate");
                //
                // --- Fit the model 1D
                fFitResults_1DYield_Array   .push_back ( FitModel    ( h1D_Nrec_PT, h1D_ResolutionReference, Form( kASigExtr_Plot_Direct, "Yield/Systematics/" ) + TString( "1D/" ) + kCurrent_Option, kCurrent_Option + TString(Form("MT_%i",iMult)), kCurrent_Option ) );
                //
                // --- Progressive Count
                fProgrCount++;
                fPrintLoopTimer("Signal Extraction Systematics Production 1D",fProgrCount,fTotalCount,1);
                //
                // --- Save to file
                TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitChkRstSY,"Yield/Systematics/","1D",kCurrent_Option.Data()),"recreate");
                //
                fHEventCount->Write();
                for ( auto hSave : fFitResults_1DYield_Array.at(kSystFit) )    hSave   ->  Write();
                //
                outFile_Check   ->  Close();
                outFile_Result  ->  Close();
            }
        }
        fStopTimer("Signal Extraction Systematics Production 1D");
        //
        /*
        //  --- Save FitResults
        std::vector<std::vector<TH1F*>> f1DCheck_Array;
        std::vector<std::vector<TH2F*>> fFitResults_2DYield_Array;
        //
        // --- Set the print progress utilities
        fTotalCount =   kSyst_SEX_2D_Options.size()*nBinMult;
        fProgrCount =   0;
        fStartTimer("Signal Extraction Systematics Production 2D");
        //
        for ( auto kSystFit = 0; kSystFit < kSyst_SEX_2D_Options.size(); kSystFit++ ) {
            auto    kCurrent_Option     =   kSyst_SEX_2D_Options.at(kSystFit);
            //
            // --- Building output and check plots directory
            gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield/Systematics/"))+TString("2D/") + kCurrent_Option );
            TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitChkPltSY,"Yield/Systematics/","2D",kCurrent_Option.Data()),"recreate");
            //
            // --- Fit the model 2D
            std::vector<TH1F*> kUtility_Check1D;
            fFitResults_2DYield_Array   .push_back  ( FitModel   ( h1D_Nrec_2Db_PT, h2Db_ResolutionReference, h2D_Nrec_PT, kUtility_Check1D, Form( kASigExtr_Plot_Direct, "Yield/Systematics/" ) + TString( "2D/" ) + kCurrent_Option, kCurrent_Option, kCurrent_Option ) );
            f1DCheck_Array              .push_back  ( kUtility_Check1D );
            //
            // --- Progressive Count
            fProgrCount++;
            fPrintLoopTimer("Signal Extraction Systematics Production 2D",fProgrCount,fTotalCount,1);
            //
            // --- Save to file
            TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitChkRstSY,"Yield/Systematics/","2D",kCurrent_Option.Data()),"recreate");
            //
            fHEventCount->Write();
            for ( auto hSave : f1DCheck_Array.at(kSystFit) )            hSave   ->  Write();
            for ( auto hSave : fFitResults_2DYield_Array.at(kSystFit) ) hSave   ->  Write();
            //
            outFile_Check   ->  Close();
            outFile_Result  ->  Close();
        }
        fStopTimer("Signal Extraction Systematics Production 2D");
        //
         */
    }
}
