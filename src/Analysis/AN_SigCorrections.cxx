#include "../../inc/AliAnalysisPhiPair.h"

void AN_SigCorrections   ( TString fOption = "all", TString kFolder = "", Bool_t fSilent = true )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- INFO on Set-up variables
    fChooseOption(fOption);
    //
    //  --- Silencing warnings for smoother running
    if ( fSilent )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
        gErrorIgnoreLevel   =   kWarning;
    }
    //
    //  --- Setting the input datastructure
    fSetAllBins();
    //
    //  --- Utility variables
    Int_t fTotalCount, fProgrCount;
    //
    // --- YIELD ANALYSIS
    if ( kDoYield ) {
        //  --- Load Files
        TFile*  insFile_Data_YL     =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("Yield")+kFolder).Data()));
        TFile*  insFile_Effc_YL     =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()));
        //
        //  --- Build output structure
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data())));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("/1D"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("/2D"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()))+TString("/Full"));
        auto        kPlotFolder1D   =   TString( Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()) ) + TString( "/1D/" );
        auto        kPlotFolder2D   =   TString( Form(kASigExtp_Plot_Direct,(TString("Yield")+kFolder).Data()) ) + TString( "/2D/" );
        //
        //  --- Load Histograms
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "fQC_Event_Enum_FLL" );
        auto        h1D_Nraw        =   uLoadHistograms<0,TH1F> ( insFile_Data_YL,  "h1D_Nraw_" );
        auto        h1D_Nrec        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec" );
        auto        h1D_Ngen        =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen" );
        auto        h2D_Nraw        =   uLoadHistograms<0,TH2F> ( insFile_Data_YL,  "anSS2D_" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_YL,  "h1D_Ngen_2Db" );
        //
        TFile* kTemp = new TFile("/Users/nikolajal/alice/AliAnalysisPhiCount/result/__HISTORY/Yield_7TeV/PreProcessing/IM_MonteCarloTruth.root");
        auto        hGEN_INELVTX_1D =   uLoadHistograms<0,TH1F> ( kTemp,  "hGEN_INELVTX_1D" );
        auto        hGEN_1D         =   uLoadHistograms<0,TH1F> ( kTemp,  "hGEN_1D" );
        auto        hSigLoss        =   (TH1F*)(hGEN_1D->Clone());
        hSigLoss->Divide(hGEN_1D,hGEN_INELVTX_1D);
        //
        //  --- Output File for Fit Check
        TFile*  outFile_Chck_YL     =   new TFile(Form(kASigExtp_FitCheckPlt,"Correlation"),"recreate");
        //
        //  --- Minimum Bias Normalisation
        auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger));
        auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex));
        auto        kN_MB           =   (fHEventCount->GetBinContent(kEventCount::kVertex10));
        Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        if ( kEnergy == 5 ) {   f1DCorrection   *=  kTriggerEff15n17pq; f2DCorrection   *=  kTriggerEff15n17pq; }
        if ( kEnergy == 7 ) {   f1DCorrection   *=  kTriggerEff10bcdef; f2DCorrection   *=  kTriggerEff10bcdef; }
        //
        fSetAllCustomFunctions();
        //
                    h1D_Nraw        ->  Scale(1.,"width");
        //
        TH1F*       h1D_Nraw_stat   =   uEfficiencyCorrection1D ( h1D_Nraw, h1D_Nrec, h1D_Ngen, f1DCorrection );
        TH1F*       h1D_Nraw_syst   =   fSetSystErrors ( h1D_Nraw_stat );
        SetAxis(h1D_Nraw_stat,"PT1D");
        SetAxis(h1D_Nraw_syst,"PT1D");
        h1D_Nraw_stat   ->  SetName("h1D_Nraw_stat");
        h1D_Nraw_syst   ->  SetName("h1D_Nraw_syst");
        //
        auto        k1D_Results     =   uMeasureFullYield(h1D_Nraw_stat,h1D_Nraw_syst,kStandardSystematicFitFunctions,100,kPlotFolder1D,"1D");
        //
                    h2D_Nraw        ->  Scale(1.,"width");
        auto        h2D_Nraw_stat   =   uEfficiencyCorrection2D ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection );
        TH2F*       h2D_Nraw_syst   =   fSetSystErrors ( h2D_Nraw_stat );
        //
        auto        k2D_Results     =   uMeasureFullYield2D( h2D_Nraw_stat, h2D_Nraw_stat, kStandardSystematicFitFunctions, 100, kPlotFolder2D, "2D_%i" );
        //
        //  --- Output
        TFile *outFile_Rslt_YL      =   new TFile   (Form(kASigExtp_FitCheckRst,(TString("Yield")+kFolder).Data()),"recreate");
        //
        h1D_Nraw_stat               ->  Write();
        h1D_Nraw_syst               ->  Write();
        //
        outFile_Rslt_YL->Close();
        outFile_Chck_YL->Close();
        insFile_Effc_YL->Close();
        insFile_Data_YL->Close();
    }
    // --- CORRELATION ANALYSIS
    if ( kDoCorrelation && false ) {
        //  --- Load Files
        TFile*  insFile_Data_CR     =   new TFile   (Form(kASigExtr_FitCheckRst,"Correlation"));
        TFile*  insFile_Effc_CR     =   new TFile   (Form(kAnalysis_MCTruthHist,"Correlation"));
        //
        //  --- Build output structure
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,"Correlation")));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Correlation")));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Correlation"))+TString("/1D"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Correlation"))+TString("/2D"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Correlation"))+TString("/Full"));
        //
        //  --- Load Histograms
        auto        fHEventCount    =   uLoadHistograms<0,TH1F> ( insFile_Data_CR, "fQC_Event_Enum_FLL" );
        auto        h1D_Nrec_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_CR, "h1D_Nrec_2Db" );
        auto        h1D_Ngen_2Db    =   uLoadHistograms<0,TH1F> ( insFile_Effc_CR, "h1D_Ngen_2Db" );
        auto        h2D_Nraw_CR     =   uLoadHistograms<1,TH2F> ( insFile_Data_CR, "anSS2D_CR_%i" );
        //
        //  --- Minimum Bias Normalisation
        auto        kN_Trg          =   (fHEventCount->GetBinContent(kEventCount::kTrigger) );
        auto        kN_Vtx          =   (fHEventCount->GetBinContent(kEventCount::kVertex)  );
        auto        kN_MB           =   (fHEventCount->GetBinContent(kEventCount::kVertex10));
        Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        //
        //  --- Output File for Fit Check
        TFile*  outFile_Chck_CR     =   new TFile(Form(kASigExtp_FitCheckPlt,"Correlation"),"recreate");
        //
        //  --- Efficiency correction
        std::vector<std::vector<TH1F*>> h2D_Nres_Cd1_Stat;
        std::vector<std::vector<TH1F*>> h2D_Nres_Cd1_Syst;
        //
        for ( auto h2D_Nraw : h2D_Nraw_CR ) {
            // TODO: CHECK FIX  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            h2D_Nraw->Scale(1.,"width");
            // TODO: CHECK FIX  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            h2D_Nres_Cd1_Stat.push_back( fEfficiencycorrection ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection ) );
            // TODO: CHECK FIX  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            for ( Int_t iBin = 0; iBin < h2D_Nraw->GetNbinsX(); ++iBin ) for ( Int_t jBin = 0; jBin < h2D_Nraw->GetNbinsY(); ++jBin ) h2D_Nraw->SetBinError(iBin,jBin,0);
            // TODO: CHECK FIX  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            h2D_Nres_Cd1_Syst.push_back( fEfficiencycorrection ( h2D_Nraw, h1D_Nrec_2Db, h1D_Ngen_2Db, f2DCorrection ) );
        }
        //
        //  --- Final Result Histograms
        //
        hName   =   Form( "h1D_PhiCorr_Stat" );
        hTitle  =   Form( "h1D_PhiCorr_Stat" );
        TH1F*   h1D_PhiCorr_Stat    =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
        SetAxis( h1D_PhiCorr_Stat, "CR" );
        //
        hName   =   Form( "h1D_PhiCorr_Syst" );
        hTitle  =   Form( "h1D_PhiCorr_Syst" );
        TH1F*   h1D_PhiCorr_Syst    =   new TH1F ( hName, hTitle, nBinCrPh, fArrCrPh );
        SetAxis( h1D_PhiCorr_Syst, "CR" );
        //
        // --- Set the print progress utilities
        fTotalCount = nBinCrPh;
        fProgrCount = 0;
        fStartTimer("Correlation Analysis Signal Correction");
        //
        /*
        for ( auto iTer = 0; iTer < h2D_Nres_Cd1_Stat.size(); ++iTer ) {
            hName   =   Form( "h2D_Nres_Cd1_Stat_0_CR_%i", iTer );
            hTitle  =   Form( "h2D_Nres_Cd1_Stat_0_CR_%i", iTer );
            auto    h2D_Nres_Cd1_Stat_0 =   new TH1F( hName, hTitle, nBinPT2D, fArrPT2D );
            //
            hName   =   Form( "h2D_Nres_Cd1_Syst_0_CR_%i", iTer );
            hTitle  =   Form( "h2D_Nres_Cd1_Syst_0_CR_%i", iTer );
            auto    h2D_Nres_Cd1_Syst_0 =   new TH1F( hName, hTitle, nBinPT2D, fArrPT2D );
            //
            for ( auto jTer = 0; jTer < h2D_Nres_Cd1_Stat.at(iTer).size(); ++jTer ) {
                auto    kCurrent_FullYield  =   fMeasureFullYield(h2D_Nres_Cd1_Stat.at(iTer).at(jTer),h2D_Nres_Cd1_Syst.at(iTer).at(jTer),Form("1D_2D_CR_%i_%i",iTer,jTer+1),Form(kASigExtp_Plot_Direct,"Correlation"));
                h2D_Nres_Cd1_Stat_0 ->  SetBinContent   ( jTer+1, kCurrent_FullYield[10] );
                h2D_Nres_Cd1_Stat_0 ->  SetBinError     ( jTer+1, kCurrent_FullYield[11] );
                h2D_Nres_Cd1_Syst_0 ->  SetBinContent   ( jTer+1, kCurrent_FullYield[10] );
                h2D_Nres_Cd1_Syst_0 ->  SetBinError     ( jTer+1, kCurrent_FullYield[12] );
            }
            //
            auto fIntegral_Stat =   0.;
            auto fIntegral_Syst =   0.;
            auto fIntegral_Val_ =   uHistoIntegralAndError(h2D_Nres_Cd1_Stat.at(iTer),fIntegral_Stat);
                 fIntegral_Val_ =   uHistoIntegralAndError(h2D_Nres_Cd1_Syst.at(iTer),fIntegral_Syst);
            //
            push_to_front( h2D_Nres_Cd1_Stat.at(iTer), h2D_Nres_Cd1_Stat_0 );
            push_to_front( h2D_Nres_Cd1_Syst.at(iTer), h2D_Nres_Cd1_Syst_0 );
            //
            auto    kCurrent_FullYield  =   fMeasureFullYield( h2D_Nres_Cd1_Stat.at(iTer).at(0),h2D_Nres_Cd1_Syst.at(iTer).at(0),Form("1D_2D_CR_%i_%i",iTer,0),Form(kASigExtp_Plot_Direct,"Correlation"));
            //
            auto    fExtrapol_Val_ =   kCurrent_FullYield[10];
            auto    fExtrapol_Stat =   kCurrent_FullYield[11];
            auto    fExtrapol_Syst =   kCurrent_FullYield[12];
            //
            auto    fSecondIntErr   =   1.;
            auto fTotalYield    =   fExtrapol_Val_ +    fIntegral_Val_ +    h2D_Nres_Cd1_Stat.at(iTer).at(0)->IntegralAndError(-1,1000,fSecondIntErr,"width");
            auto fTotalEStat    =   SquareSum( { fIntegral_Stat, fExtrapol_Stat, fSecondIntErr } );
                                                                            h2D_Nres_Cd1_Syst.at(iTer).at(0)->IntegralAndError(-1,1000,fSecondIntErr,"width");
            auto fTotalESyst    =   SquareSum( { fIntegral_Syst, fExtrapol_Syst, fSecondIntErr } );
            h1D_PhiCorr_Stat    ->  SetBinContent   ( iTer+1, fTotalYield );
            h1D_PhiCorr_Stat    ->  SetBinError     ( iTer+1, fTotalEStat );
            h1D_PhiCorr_Syst    ->  SetBinContent   ( iTer+1, fTotalYield );
            h1D_PhiCorr_Syst    ->  SetBinError     ( iTer+1, fTotalESyst );
            //
            // --- Progressive Count
            fProgrCount++;
            fPrintLoopTimer("Correlation Analysis Signal Correction",fProgrCount,fTotalCount,1);
        }
         */
        //
        fStopTimer("Correlation Analysis Signal Correction");
        //
        //  --- Output
        TFile *outFile_Rslt_CR      =   new TFile   (Form(kASigExtp_FitCheckRst,"Correlation"),"recreate");
        //
        h1D_PhiCorr_Stat->Scale(1.,"width");
        h1D_PhiCorr_Syst->Scale(1.,"width");
        //
        h1D_PhiCorr_Stat->Write();
        h1D_PhiCorr_Syst->Write();
        //
        outFile_Rslt_CR->Close();
        outFile_Chck_CR->Close();
        insFile_Effc_CR->Close();
        insFile_Data_CR->Close();
    }
}
