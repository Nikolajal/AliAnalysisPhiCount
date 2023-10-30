#include "../../inc/AliAnalysisPhiPair.h"
//TODO: Add signal strength, SNR, Chi^2

void AN_SigExtraction   ( TString fOption = "yield", TString kFolder = "", Bool_t fSilent = true )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    // --- INFO on Set-up variables
    fChooseOption(fOption);
    //
    // --- Silencing warnings for smoother running
    if ( fSilent )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
        gErrorIgnoreLevel   =   kWarning;
    }
    //
    // --- Setting the input datastructure
    fSetAllBins();
    //
    // --- Utility variables
    Int_t fTotalCount, fProgrCount;
    //
    // --- YIELD ANALYSIS
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
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("1D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("2D/"));
        TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitCheckPlt,(TString("Yield")+kFolder).Data()),"recreate");
        //
        // --- Set the print progress utilities
        fTotalCount = nBinPT1D + nBinPT2D * ( nBinPT2D + 1 );
        fProgrCount = 0;
        fStartTimer("Yield Analysis Signal Extraction");
        //
        // --- Fit the model 1D
        auto    fFitResults_1DYield    =   FitModel    ( h1D_Nrec_PT, h1D_ResolutionReference, Form( kASigExtr_Plot_Direct, (TString("Yield")+kFolder).Data() ) + TString( "1D/" ) );
        //
        // --- Progressive Count
        fProgrCount +=  nBinPT1D;
        fPrintLoopTimer("Yield Analysis Signal Extraction",fProgrCount,fTotalCount,1);
        //
        // --- Fit the model 2D
        std::vector<TH1F*>  f1DCheck;
        auto    fFitResults_2DYield = FitModel( h1D_Nrec_2Db_PT, h2Db_ResolutionReference, h2D_Nrec_PT, f1DCheck, Form( kASigExtr_Plot_Direct, (TString("Yield")+kFolder).Data() ) + TString( "2D/" ) );
        //
        // --- Progressive Count & Stop timer
        fProgrCount +=  nBinPT2D * ( nBinPT2D + 1 );
        fPrintLoopTimer("Yield Analysis Signal Extraction",fProgrCount,fTotalCount,1);
        fStopTimer("Yield Analysis Signal Extraction");
        //
        // --- Save to file
        TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Yield")+kFolder).Data()),"recreate");
        //
        fHEventCount->Write();
        for ( auto hSave : fFitResults_1DYield )    hSave   ->  Write();
        for ( auto hSave : f1DCheck )               hSave   ->  Write();
        for ( auto hSave : fFitResults_2DYield )    hSave   ->  Write();
        //
        outFile_Check   ->  Close();
        outFile_Result  ->  Close();
    }
    // --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity ) {
        // --- Retrieving PreProcessed Histograms
        TFile*      insFile_Data_ML         =   new TFile   ( Form(kAnalysis_InvMassHist,   (TString("Multiplicity")    +kFolder).Data()) );
        TFile*      insFile_Resl_ML         =   new TFile   ( Form(kMassResolution_Anal,    (TString("Multiplicity")    +kFolder).Data()) );
        //
        auto    h1D_ResolutionReference     =   uLoadHistograms<0,TH1F> ( insFile_Resl_ML, "hRes_RMS_3_1D" );
        auto    h2Db_ResolutionReference    =   uLoadHistograms<0,TH1F> ( insFile_Resl_ML, "hRes_RMS_3_2Db" );
        auto    fHEventCount                =   uLoadHistograms<0,TH1F> ( insFile_Data_ML, "fQC_Event_Enum_FLL" );
        auto    fHEvCountMlt                =   uLoadHistograms<0,TH1F> ( insFile_Data_ML, "fQC_Event_Enum_V0M" );
        auto    h1D_Nrec_MT_PT              =   uLoadHistograms<2,TH1F> ( insFile_Data_ML, "h1D_Nrec_MT_%i_PT_%i" );
        auto    h1D_Nrec_2Db_MT_PT          =   uLoadHistograms<2,TH1F> ( insFile_Data_ML, "h1D_Nrec_2Db_MT_%i_PT_%i" );
        auto    h2D_Nrec_MT_PT              =   uLoadHistograms<3,TH2F> ( insFile_Data_ML, "h2D_Nrec_MT_%i_PT_%i_%i" );
        //
        // --- Building output and check plots directory
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",(Form(kAnalysis_SigExtr_Dir,(TString("Multiplicity")+kFolder).Data()))));
        Int_t   iMult = -1;
        for ( auto kTarget : h1D_Nrec_MT_PT )   {
            iMult++;
            gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kASigExtr_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/1D/",iMult))).Data()));
            gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kASigExtr_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/2D/",iMult))).Data()));
        }
        TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitCheckPlt,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        // --- Set the print progress utilities
        fTotalCount = ( nBinPT1D + nBinPT2D * ( nBinPT2D + 1 ) ) * ( nBinMult + 1 );
        fProgrCount = 0;
        fStartTimer("Multiplicity Analysis Signal Extraction");
        //
        // --- Fit the model 1D
        iMult   =   -1;
        std::vector<std::vector<TH1F*>> fFitResults_1DYield_Array;
        for ( auto kTarget : h1D_Nrec_MT_PT )   {
            iMult++;
            fFitResults_1DYield_Array.push_back( FitModel    ( kTarget, h1D_ResolutionReference, (TString(Form(kASigExtr_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/1D/",iMult))).Data(), Form("MT_%i",iMult) ) );
            //
            // --- Progressive Count
            fProgrCount +=  nBinPT1D;
            fPrintLoopTimer("Multiplicity Analysis Signal Extraction",fProgrCount,fTotalCount,1);
        }
        //
        // --- Fit the model 2D
        iMult   =   -1;
        std::vector<std::vector<TH1F*>> f1DCheck_Array;
        std::vector<std::vector<TH2F*>> fFitResults_2DYield_Array;
        for ( auto kTarget : h2D_Nrec_MT_PT )   {
            iMult++;
            std::vector<TH1F*>  f1DCheck;
            fFitResults_2DYield_Array.push_back(    FitModel    ( h1D_Nrec_2Db_MT_PT[iMult], h2Db_ResolutionReference, kTarget, f1DCheck, (TString(Form(kASigExtr_Plot_Direct,(TString("Multiplicity")+kFolder).Data()))+TString(Form("/MLT_%i/2D/",iMult))).Data(), Form("MT_%i",iMult) ) );
            f1DCheck_Array.push_back(  *(new std::vector<TH1F*> (f1DCheck))  );
            //
            // --- Progressive Count
            fProgrCount +=  nBinPT2D * ( nBinPT2D + 1 );
            fPrintLoopTimer("Multiplicity Analysis Signal Extraction",fProgrCount,fTotalCount,1);
        }
        //
        fStopTimer("Multiplicity Analysis Signal Extraction");
        //
        // --- Save to file
        TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        fHEventCount->Write();
        fHEvCountMlt->Write();
        for ( auto hVecSave : fFitResults_1DYield_Array )   for ( auto hSave : hVecSave )   hSave   ->  Write();
        for ( auto hVecSave : f1DCheck_Array )              for ( auto hSave : hVecSave )   hSave   ->  Write();
        for ( auto hVecSave : fFitResults_2DYield_Array )   for ( auto hSave : hVecSave )   hSave   ->  Write();
        //
        outFile_Check   ->  Close();
        outFile_Result  ->  Close();
    }
    // --- CORRELATION ANALYSIS
    if ( kDoCorrelation ) {
        // --- Retrieving PreProcessed Histograms
        TFile*      insFile_Data_CR         =   new TFile   ( Form(kAnalysis_InvMassHist,   (TString("Correlation") +kFolder).Data()) );
        TFile*      insFile_Resl_CR         =   new TFile   ( Form(kMassResolution_Anal,    (TString("Correlation") +kFolder).Data()) );
        //
        auto    h1D_ResolutionReference     =   uLoadHistograms<0,TH1F> ( insFile_Resl_CR, "hRes_RMS_3_1D" );
        auto    h2Db_ResolutionReference    =   uLoadHistograms<0,TH1F> ( insFile_Resl_CR, "hRes_RMS_3_2Db" );
        auto    fHEventCount                =   uLoadHistograms<0,TH1F> ( insFile_Data_CR, "fQC_Event_Enum_FLL" );
        auto    fHEvCountMlt                =   uLoadHistograms<0,TH1F> ( insFile_Data_CR, "fQC_Event_Enum_V0M" );
        auto    h1D_Nrec_2Db_PT             =   uLoadHistograms<1,TH1F> ( insFile_Data_CR, "h1D_Nrec_2Db_PT_%i" );
        auto    h2D_Nrec_CR_PT              =   uLoadHistograms<3,TH2F> ( insFile_Data_CR, "h2D_Nrec_CR_%i_PT_%i_%i" );
        //
        // --- Building output and check plots directory
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,(TString("Correlation")+kFolder).Data())));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString("1D/"));
        gROOT                               ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,(TString("Correlation")+kFolder).Data()))+TString("2D/"));
        TFile*  outFile_Check               =   new TFile   (Form(kASigExtr_FitCheckPlt,(TString("Correlation")+kFolder).Data()),"recreate");
        //
        // --- Set the print progress utilities
        fTotalCount = ( nBinPT1D + nBinPT2D * ( nBinPT2D + 1 ) ) * ( nBinCrPh );
        fProgrCount = 0;
        fStartTimer("Correlation Analysis Signal Extraction");
        //
        // --- Fit the model 2D
        auto    iCrPh   =   -1;
        std::vector<std::vector<TH1F*>> f1DCheck_Array;
        std::vector<std::vector<TH2F*>> fFitResults_2DYield_Array;
        for ( auto kTarget : h2D_Nrec_CR_PT )   {
            iCrPh++;
            std::vector<TH1F*>  f1DCheck;
            fFitResults_2DYield_Array.push_back(    FitModel    ( h1D_Nrec_2Db_PT, h2Db_ResolutionReference, kTarget, f1DCheck, Form( kASigExtr_Plot_Direct,(TString("Correlation")+kFolder).Data()) + TString( "2D/" ), Form("CR_%i",iCrPh) ) );
            f1DCheck_Array.push_back(  *(new std::vector<TH1F*> (f1DCheck))  );
            //
            // --- Progressive Count
            fProgrCount +=  nBinPT2D * ( nBinPT2D + 1 );
            fPrintLoopTimer("Correlation Analysis Signal Extraction",fProgrCount,fTotalCount,1);
        }
        //
        fStopTimer("Correlation Analysis Signal Extraction");
        //
        // --- Save to file
        TFile*      outFile_Result  =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Correlation")+kFolder).Data()),"recreate");
        //
        fHEventCount->Write();
        fHEvCountMlt->Write();
        for ( auto hVecSave : f1DCheck_Array )              for ( auto hSave : hVecSave )   hSave   ->  Write();
        for ( auto hVecSave : fFitResults_2DYield_Array )   for ( auto hSave : hVecSave )   hSave   ->  Write();
        //
        outFile_Check   ->  Close();
        outFile_Result  ->  Close();
    }
}
