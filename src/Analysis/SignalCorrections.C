// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void
testfunc
( TH1F* h ) {
    h->GetXaxis()->SetNdivisions(6);
    h->GetXaxis()->SetBinLabel(6,"[0-100]");
    h->GetXaxis()->SetBinLabel(5,"[0-5]");
    h->GetXaxis()->SetBinLabel(4,"[5-15]");
    h->GetXaxis()->SetBinLabel(3,"[15-30]");
    h->GetXaxis()->SetBinLabel(2,"[30-50]");
    h->GetXaxis()->SetBinLabel(1,"[50-100]");
    h->SetMarkerStyle(markers[1]);
    h->SetMarkerColor(colors[3]);
    h->SetMaximum(1.2*h->GetMaximum());
}

void
SignalCorrections
 ( TString fOption = "", bool fSilent = true, TString kFolder = "" )    {
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //-// OPTIONS
    
    // Silencing warnings for smoother
    if ( fSilent )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    fChooseOption(fOption);
    gErrorIgnoreLevel   =   kFatal;
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (Form(kASigExtr_FitCheckRst,"Yield"));
    TFile*  insFile_EF_Yield            =   new TFile   (Form(kAnalysis_MCTruthHist,"Yield"));
    TFile*  insFile_DT_Multi            =   new TFile   (Form(kASigExtr_FitCheckRst,"Multiplicity"));
    TFile*  insFile_EF_Multi            =   new TFile   (Form(kAnalysis_MCTruthHist,"Multiplicity"));
    
    
    // Recovering the histograms-------------------------------------------------------------------------------
    
    // >-> GENERAL ANALYSIS //
    //
    TH1D       *hEvntEff;
    TH1D       *hEvntMlt;
    //
    hName       =   "fQC_Event_Enum_FLL";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "fQC_Event_Enum_V0M";
    hEvntMlt    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    //
    // >-> YIELD ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hRAW_1D;
    TH1F       *hREC_1D;
    TH1F       *hGEN_1D;
    TH1F       *hGEN_INELVTX_1D;
    TH1F       *hGEN_INELFLL_1D;
    TH1F       *hREC_1D_in_2D_bin;
    TH1F       *hGEN_1D_in_2D_bin;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hRAW_1D";
    hRAW_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "hREC_1D";
    hREC_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_1D";
    hGEN_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_INELFLL_1D";
    hGEN_INELFLL_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_INELVTX_1D";
    hGEN_INELVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hREC_1D_in_2D_bin";
    hREC_1D_in_2D_bin     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_1D_in_2D_bin";
    hGEN_1D_in_2D_bin     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hRAW_2D;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hRAW_2D";
    hRAW_2D     =   (TH2F*)(insFile_DT_Yield->Get(hName));
    //
    // >-> MULTIPLICITY ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //---------------------//
    //  Setting up output  //
    //---------------------//
    //
    // Generating the binning array--------------------------------------------------------------------------
    //
    fSetAllBins();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept;
    //
    // Creating the histograms-------------------------------------------------------------------------------
    //
    //---------------------//
    // Preprocessing input //
    //---------------------//
    //
    //                 N_raw            f_norm X f_vtx X f_SL
    // N_res = --------------------- X -----------------------
    //          EXA X DpT X Dy X BR             N_MB
    //
    // Scaling in pT [Done in PreProcessing]
    //
    // Scaling for efficiencies
    //
    auto        kN_Trg          =   (hEvntEff->GetBinContent(kEventCount::kTrigger));
    auto        kN_Vtx          =   (hEvntEff->GetBinContent(kEventCount::kVertex));
    auto        kN_MB           =   (hEvntEff->GetBinContent(kEventCount::kVertex10));
    Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(1./kSignalMiss1D) *(kN_Vtx/kN_Trg);
    Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(1./kSignalMiss2D) *(kN_Vtx/kN_Trg);
    //
    // >-> YIELD ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F   *hRES_1D_Stat    =   fEfficiencycorrection   ( fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection), hGEN_1D, hGEN_INELVTX_1D );
    TH1F   *hRES_1D_Syst    =   fSetSystErrors          ( hRES_1D_Stat );
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    std::vector<TH1F*>  hRES_2D_Cond1_Stat  =   fEfficiencycorrection   ( hRAW_2D,hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
    std::vector<TH1F*>  hRES_2D_Cond1_Syst  =   fSetSystErrors          ( hRES_2D_Cond1_Stat );
    //
    hName   =   Form("hRES_2D_Cond2_Stat");
    hTitle  =   Form("hRES_2D_Cond2_Stat");
    TH1F       *hRES_2D_Cond2_Stat              =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("hRES_2D_Cond2_Syst");
    hTitle  =   Form("hRES_2D_Cond2_Syst");
    TH1F       *hRES_2D_Cond2_Syst              =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("hRES_2D_Cond3_Stat");
    hTitle  =   Form("hRES_2D_Cond3_Stat");
    TH1F       *hRES_2D_Cond3_Stat              =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("hRES_2D_Cond3_Syst");
    hTitle  =   Form("hRES_2D_Cond3_Syst");
    TH1F       *hRES_2D_Cond3_Syst              =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("gConditional_Mean_PT");
    hTitle  =   Form("gConditional_Mean_PT");
    TGraphMultiErrors  *gConditional_Mean_PT    =   new TGraphMultiErrors(hName,hTitle,nBinPT2D,2);
    //
    //
    hName   =   Form("hMultipleResults_Stat");
    hTitle  =   Form("hMultipleResults_Stat");
    TH1F*   hMultipleResults_Stat               =   new TH1F(hName,hTitle,6,0,6);
    //
    hName   =   Form("hMultipleResults_Syst");
    hTitle  =   Form("hMultipleResults_Syst");
    TH1F*   hMultipleResults_Syst               =   new TH1F(hName,hTitle,6,0,6);
    //
    hName   =   Form("hMeanPT_2D_Stat");
    hTitle  =   Form("hMeanPT_2D_Stat");
    //TH1F*   hMeanPT_2D_Stat                     =   new TH1F(hName,hTitle,fArrPT2D_Comp,nBinPT2D+1);
    //
    hName   =   Form("hMeanPT_2D_Syst");
    hTitle  =   Form("hMeanPT_2D_Syst");
    //TH1F*   hMeanPT_2D_Syst                     =   new TH1F(hName,hTitle,fArrPT2D_Comp,nBinPT2D+1);
    //
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    if ( kDoYield ) {
        //
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,"Yield")));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield")));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield"))+TString("/1D"));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield"))+TString("/2D"));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield"))+TString("/Full"));
        //
        fStartTimer("Fit_for_extrapolation");
        //
        // Output File for Fit Check
        TFile*  outCheckFitYld  =   new TFile(Form(kASigExtp_FitCheckPlt,"Yield"),"recreate");
        //
        // Total Fit number and progressive
        //
        Int_t   fTotalCount = 2+nBinPT2D;
        Int_t   fProgrCount = 0;
        //
        auto    fResults = fMeasureFullYield(hRES_1D_Stat,hRES_1D_Syst,"1D",Form(kASigExtp_Plot_Direct,"Yield"));
        //  Progressive Count
        fProgrCount++;
        //  Print loop Timer
        fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
        //
        hMultipleResults_Stat   ->  SetBinContent   (1, fResults[0]);
        hMultipleResults_Stat   ->  SetBinError     (1, fResults[1]);
        hMultipleResults_Syst   ->  SetBinContent   (1, fResults[0]);
        hMultipleResults_Syst   ->  SetBinError     (1, fResults[3]);
        //
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            fResults = fMeasureFullYield(hRES_2D_Cond1_Stat.at(iFit),hRES_2D_Cond1_Syst.at(iFit),Form("1D_2D_%i",iFit),Form(kASigExtp_Plot_Direct,"Yield"));
            //
            hRES_2D_Cond2_Stat  ->  SetBinContent   ( iFit+1, fResults[10] );
            hRES_2D_Cond2_Stat  ->  SetBinError     ( iFit+1, fResults[11] );
            hRES_2D_Cond2_Syst  ->  SetBinContent   ( iFit+1, fResults[10] );
            hRES_2D_Cond2_Syst  ->  SetBinError     ( iFit+1, fResults[12] );
            //
            //hMeanPT_2D_Stat     ->  SetBinContent   ( iFit+1, fResults[5]);
            //hMeanPT_2D_Stat     ->  SetBinError     ( iFit+1, fResults[6]);
            //hMeanPT_2D_Syst     ->  SetBinContent   ( iFit+1, fResults[5]);
            //hMeanPT_2D_Syst     ->  SetBinError     ( iFit+1, fResults[8]);
        }
        //
        auto fResult2       = fMeasureFullYield(hRES_2D_Cond2_Stat,hRES_2D_Cond2_Syst,"2D_-1",Form(kASigExtp_Plot_Direct,"Yield"));
        //
        auto fExtrap_1_Stat =   0.;
        auto fExtrap_1_Syst =   0.;
        auto fExtrap_1_Val_ =   hRES_2D_Cond2_Stat->IntegralAndError(-1.,10000.,fExtrap_1_Stat,"width");
             fExtrap_1_Val_ =   hRES_2D_Cond2_Syst->IntegralAndError(-1.,10000.,fExtrap_1_Syst,"width");
        //
        auto fExtrap_2_Val_ =   fResult2[10];
        auto fExtrap_2_Stat =   fResult2[11];
        auto fExtrap_2_Syst =   fResult2[12];
        //
        auto fIntegral_Stat =   0.;
        auto fIntegral_Syst =   0.;
        auto fIntegral_Val_ =   uHistoIntegralAndError(hRES_2D_Cond1_Stat,fIntegral_Stat);
             fIntegral_Val_ =   uHistoIntegralAndError(hRES_2D_Cond1_Syst,fIntegral_Syst);
        //
        hMultipleResults_Stat   ->  SetBinContent   (2, fIntegral_Val_ + fExtrap_1_Val_ + fExtrap_2_Val_                 );
        hMultipleResults_Stat   ->  SetBinError     (2, SquareSum({ fIntegral_Stat, fExtrap_1_Stat, fExtrap_2_Stat })    );
        hMultipleResults_Syst   ->  SetBinContent   (2, fIntegral_Val_ + fExtrap_1_Val_ + fExtrap_2_Val_                 );
        hMultipleResults_Syst   ->  SetBinError     (2, SquareSum({ fIntegral_Syst, fExtrap_1_Syst, fExtrap_2_Syst })    );
        //
        hMultipleResults_Stat   ->  SetBinContent   (3, (hMultipleResults_Stat->GetBinContent(2))/(hMultipleResults_Stat->GetBinContent(1)) );
        hMultipleResults_Stat   ->  SetBinError     (3, hMultipleResults_Stat->GetBinContent(3)*SquareSum( { hMultipleResults_Stat->GetBinError(2)/hMultipleResults_Stat->GetBinContent(2) , hMultipleResults_Stat->GetBinError(1)/hMultipleResults_Stat->GetBinContent(1) } ) );
        hMultipleResults_Syst   ->  SetBinContent   (3, (hMultipleResults_Stat->GetBinContent(2))/(hMultipleResults_Stat->GetBinContent(1)) );
        hMultipleResults_Syst   ->  SetBinError     (3, 0 );
        //
        hMultipleResults_Stat   ->  SetBinContent   (4, (hMultipleResults_Stat->GetBinContent(2))/((hMultipleResults_Stat->GetBinContent(1))*(hMultipleResults_Stat->GetBinContent(1))) );
        hMultipleResults_Stat   ->  SetBinError     (4, hMultipleResults_Stat->GetBinContent(4)*SquareSum( { hMultipleResults_Stat->GetBinError(2)/hMultipleResults_Stat->GetBinContent(2) , 4*hMultipleResults_Stat->GetBinError(1)/hMultipleResults_Stat->GetBinContent(1) } ) );
        hMultipleResults_Syst   ->  SetBinContent   (4, (hMultipleResults_Stat->GetBinContent(2))/((hMultipleResults_Stat->GetBinContent(1))*(hMultipleResults_Stat->GetBinContent(1))) );
        hMultipleResults_Syst   ->  SetBinError     (4, 0 );
        //
        hMultipleResults_Stat   ->  SetBinContent   (5, fSigmaPhiValue(hMultipleResults_Stat->GetBinContent(1),hMultipleResults_Stat->GetBinContent(2)) );
        hMultipleResults_Stat   ->  SetBinError     (5, fSigmaPhiError(hMultipleResults_Stat->GetBinContent(1),hMultipleResults_Stat->GetBinContent(2),hMultipleResults_Stat->GetBinError(1),hMultipleResults_Stat->GetBinError(2)) );
        hMultipleResults_Syst   ->  SetBinContent   (5, fSigmaPhiValue(hMultipleResults_Stat->GetBinContent(1),hMultipleResults_Stat->GetBinContent(2)) );
        hMultipleResults_Syst   ->  SetBinError     (5, 0 );
        //
        
        hMultipleResults_Stat   ->  SetBinContent   (6, fGammaPhiValue(hMultipleResults_Stat->GetBinContent(1),hMultipleResults_Stat->GetBinContent(2)) );
        hMultipleResults_Stat   ->  SetBinError     (6, fGammaPhiError(hMultipleResults_Stat->GetBinContent(1),hMultipleResults_Stat->GetBinContent(2),hMultipleResults_Stat->GetBinError(1),hMultipleResults_Stat->GetBinError(2)) );
        hMultipleResults_Syst   ->  SetBinContent   (6, fGammaPhiValue(hMultipleResults_Stat->GetBinContent(1),hMultipleResults_Stat->GetBinContent(2)) );
        hMultipleResults_Syst   ->  SetBinError     (6, 0 );
        //
        fStopTimer("Fit_for_extrapolation");
        //
        outCheckFitYld  ->  Close();
    }
    //
    if ( kDoMultiplicity )  {
        //
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,"Multiplicity")));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Multiplicity")));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Multiplicity"))+TString("/1D"));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Multiplicity"))+TString("/2D"));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Multiplicity"))+TString("/Full"));
        //
        //  Event Count for Normalisation
        TH1D*   hUtilEventMultiplicity  =   (TH1D*)(insFile_DT_Multi->Get("fQC_Event_Enum_V0M"));
        //
        //  REC & GEN for Efficiencies
        TH1F*   hREC_1D_for_MT          =   (TH1F*)(insFile_EF_Multi->Get("hREC_1D"));
        TH1F*   hGEN_1D_for_MT          =   (TH1F*)(insFile_EF_Multi->Get("hGEN_1D"));
        TH1F*   hREC_2D_for_MT          =   (TH1F*)(insFile_EF_Multi->Get("hREC_1D_in_2D_bin"));
        TH1F*   hGEN_2D_for_MT          =   (TH1F*)(insFile_EF_Multi->Get("hGEN_1D_in_2D_bin"));
        //
        //
        fStartTimer("Fit_for_extrapolation");
        //
        // Output File for Fit Check
        TFile*  outCheckFitYld  =   new TFile(Form(kASigExtp_FitCheckPlt,"Mutliplicity"),"recreate");
        //
        // Total Fit number and progressive
        //
        Int_t   fTotalCount = ( nBinMult + 1 ) * ( nBinPT2D + 2 );
        Int_t   fProgrCount = 0;
        //
        std::vector<TH1F*> f1DSpectra;
        std::vector<Double_t> fYield1D;
        std::vector<Double_t> fYield2D;
        for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
            //  Normalisation Factor for Multiplicity bin
            auto    kNormalisation1D    =   (1.)/(fEvaluateINELgt0(iMlt-1,hUtilEventMultiplicity) * kBR );
            auto    kNormalisation2D    =   (1.)/(fEvaluateINELgt0(iMlt-1,hUtilEventMultiplicity) * kBR * kBR );
            //
            //  Histogram Retrieving
            hName                       =   Form("hRAW_1D_in_Mlt_%i",iMlt);
            TH1F*   hRAW_1D_in_MT       =   (TH1F*)(insFile_DT_Multi->Get(hName));
            //
            //  Efficiency and Normalisation Correction
            f1DSpectra.push_back(fEfficiencycorrection( hRAW_1D_in_MT, hREC_1D_for_MT, hGEN_1D_for_MT, kNormalisation1D ));
            // *************************************************************
            // TODO: Implement Syst Errors
            // f1DSpectra_Stat (..)
            // f1DSpectra_Syst (..)
            // *************************************************************
            //
            //  Yield evaluation
            auto    fResults    =   fMeasureFullYield(f1DSpectra.at(iMlt),f1DSpectra.at(iMlt),Form("1D_in_Mlt_%i",iMlt),Form(kASigExtp_Plot_Direct,"Multiplicity"));
            fYield1D.push_back(fResults[0]);   // !TEMP
            // *************************************************************
            // TODO: Save results to a histo
            // *************************************************************
            //
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            //  Histogram Retrieving
            hName                           =   Form("hRAW_2D_in_Mlt_%i",iMlt);
            TH2F*   hRAW_2D_in_MT           =   (TH2F*)(insFile_DT_Multi->Get(hName));
            TH1F*   hRES_2D_Cond2_in_MT     =   new TH1F ("hRES_2D_Cond1_in_MT","hRES_2D_Cond1_in_MT",nBinPT2D,fArrPT2D);
            //
            //  Efficiency and Normalisation Correction
            std::vector<TH1F*>  hRES_2D_Cond1_in_MT  =   fEfficiencycorrection   ( hRAW_2D_in_MT, hREC_2D_for_MT, hGEN_2D_for_MT, kNormalisation2D );
            //
            //  Yield Evaluation
            for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
                //  Progressive Count
                fProgrCount++;
                //  Print loop Timer
                fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
                //
                fResults = fMeasureFullYield(hRES_2D_Cond1_in_MT.at(iFit),hRES_2D_Cond1_in_MT.at(iFit),Form("1D_2D_%i",iFit+1),Form(kASigExtp_Plot_Direct,"Multiplicity"));
                //
                hRES_2D_Cond2_in_MT  ->  SetBinContent   ( iFit+1, fResults[10] );
                hRES_2D_Cond2_in_MT  ->  SetBinError     ( iFit+1, fResults[10]*0.04 );
            }
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            fResults = fMeasureFullYield(hRES_2D_Cond2_in_MT,hRES_2D_Cond2_in_MT,Form("1D_2D_%i",0),Form(kASigExtp_Plot_Direct,"Multiplicity"));
            //
            auto fExtrap_1_Stat =   0.;
            auto fExtrap_1_Syst =   0.;
            auto fExtrap_1_Val_ =   hRES_2D_Cond2_in_MT->IntegralAndError(-1.,10000.,fExtrap_1_Stat,"width");
                 fExtrap_1_Val_ =   hRES_2D_Cond2_in_MT->IntegralAndError(-1.,10000.,fExtrap_1_Syst,"width");
            //
            auto fExtrap_2_Val_ =   fResults[10];
            auto fExtrap_2_Stat =   fResults[11];
            auto fExtrap_2_Syst =   fResults[12];
            //
            auto fIntegral_Stat =   0.;
            auto fIntegral_Syst =   0.;
            auto fIntegral_Val_ =   uHistoIntegralAndError(hRES_2D_Cond1_in_MT,fIntegral_Stat);
                 fIntegral_Val_ =   uHistoIntegralAndError(hRES_2D_Cond1_in_MT,fIntegral_Syst);
            //
            fYield2D.push_back( fIntegral_Val_ + fExtrap_1_Val_ + fExtrap_2_Val_ );
            
            
            
            
            
            // TODO: Clean this FIX
            // *************************************************************
            //
            /*
            TH1F**      hRAW_1D_in_MT   =   new TH1F*   [nBinMult+1];
            //
            //  Defining cumulative histogram over measurable pT
            //
            for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
                hName                   =   Form("hRAW_1D_in_Mlt_%i",iMlt);
                hRAW_1D_in_MT[iMlt]     =   (TH1F*)(insFile_DT_Multi->Get(hName));
            }
            TFile*  kReference  =   new TFile("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Multiplicity/phi_pp5_mul_24Feb2021.root");
            std::vector<TH1F*> fRef_Array;
            std::vector<Double_t> fTest;
            std::vector<Double_t> fTest2D;
            if ( iMlt >  0 )    fRef_Array.push_back((TH1F*)(kReference->Get(Form("h%i_%i",(int)fArrMult[iMlt-1],(int)fArrMult[iMlt]))));
            if ( iMlt == 0 )    fRef_Array.push_back((TH1F*)(kReference->Get(Form("h%i_%i",0,100))));
            // *************************************************************
            
            auto    fResults    =   fMeasureFullYield(fSaveArray.at(iMlt),fSaveArray.at(iMlt),Form("1D_in_Mlt_%i",iMlt),Form(kASigExtp_Plot_Direct,"Multiplicity"));
            fTest.push_back(fResults[0]);
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            //
            */
        }
        
        // TODO: This is a fix to see first results
        // *************************************************************
        
        gROOT->SetBatch();
        TCanvas*    c1  =   new TCanvas();
        c1->Divide(3,2);
        TH1F* hShow1D   =   new TH1F("hShow1D","hShow1D",6,0,6);
        TH1F* hShow2D   =   new TH1F("hShow2D","hShow2D",6,0,6);
        TH1F* hShowR1   =   new TH1F("hShowR1","hShowR1",6,0,6);
        TH1F* hShowR2   =   new TH1F("hShowR2","hShowR2",6,0,6);
        TH1F* hShowP1   =   new TH1F("hShowP1","hShowP1",6,0,6);
        TH1F* hShowP2   =   new TH1F("hShowP2","hShowP2",6,0,6);
        int    iTer    =   1;
        for ( auto kYield : fYield1D ) {
            auto k1D    =   kYield;
            auto k2D    =   fYield2D.at(iTer-1);
            //
            hShow1D->SetBinContent  ( 7-iTer, k1D );
            hShow1D->SetBinError    ( 7-iTer, 0 );
            hShow2D->SetBinContent  ( 7-iTer, k2D );
            hShow2D->SetBinError    ( 7-iTer, 0 );
            hShowR1->SetBinContent  ( 7-iTer, k2D/k1D );
            hShowR1->SetBinError    ( 7-iTer, 0 );
            hShowR2->SetBinContent  ( 7-iTer, k2D/(k1D*k1D) );
            hShowR2->SetBinError    ( 7-iTer, 0 );
            hShowP1->SetBinContent  ( 7-iTer, fSigmaPhiValue(k1D,k2D) );
            hShowP1->SetBinError    ( 7-iTer, 0 );
            hShowP2->SetBinContent  ( 7-iTer, fGammaPhiValue(k1D,k2D) );
            hShowP2->SetBinError    ( 7-iTer, 0 );
            //
            iTer++;
        }
        testfunc(hShow1D);
        testfunc(hShow2D);
        testfunc(hShowR1);
        testfunc(hShowR2);
        testfunc(hShowP1);
        testfunc(hShowP2);
        //
        c1->cd(1);
        hShow1D->Draw("EP");
        uLatex->DrawLatexNDC(0.18,0.80,"#frac{dN_{#phi}}{dy}");
        c1->cd(4);
        hShow2D->Draw("EP");
        uLatex->DrawLatexNDC(0.18,0.80,"#frac{dN_{#phi#phi}}{dy}");
        c1->cd(2);
        hShowR1->Draw("EP");
        uLatex->DrawLatexNDC(0.18,0.80,"#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT}");
        c1->cd(5);
        hShowR2->Draw("EP");
        uLatex->DrawLatexNDC(0.18,0.80,"#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT^{2}}");
        c1->cd(3);
        hShowP1->Draw("EP");
        uLatex->DrawLatexNDC(0.18,0.80,"#sigma^{2}_{#phi}");
        c1->cd(6);
        hShowP2->Draw("EP");
        uLatex->DrawLatexNDC(0.18,0.80,"#gamma_{#phi}");
        c1->SaveAs(Form("%s/TEST.pdf",Form(kASigExtp_Plot_Direct,"Multiplicity")));
        delete c1;
        gROOT->SetBatch(kFALSE);
        //
        // *************************************************************
        //
        //--------------------------//
        //  Printing output objects //
        //--------------------------//
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,(TString("Multiplicity")+kFolder).Data())));
        TFile *outFil2  =   new TFile   (Form(kASigExtp_FitCheckRst,"Multiplicity"),"recreate");
        for ( auto kYield : f1DSpectra ) {
            kYield->Write();
        }
        hShow1D->Write();
        hShow2D->Write();
        hShowR1->Write();
        hShowR2->Write();
        hShowP1->Write();
        hShowP2->Write();
        //
        fStopTimer("Fit_for_extrapolation");
    }
    //
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Yield Analysis
    //
    gROOT->SetBatch();
    if ( kDoYield )  {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,(TString("Yield")+kFolder).Data())));
        TFile *outFil2  =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield"),"recreate");
        hEvntEff                ->Write();
        //
        hName   =   Form("hRES_1D_Stat");
        hTitle  =   Form("1D Spectrum Stat Err");
        hRES_1D_Stat            ->SetNameTitle(hName,hTitle);
        hRES_1D_Stat            ->Write();
        //
        hName   =   Form("hRES_1D_Syst");
        hTitle  =   Form("1D Spectrum Syst Err");
        hRES_1D_Syst            ->SetNameTitle(hName,hTitle);
        hRES_1D_Syst            ->Write();
        //
        gConditional_Mean_PT    ->Write();
        hRES_2D_Cond2_Stat      ->Write();
        hRES_2D_Cond2_Syst      ->Write();
        hMultipleResults_Stat   ->Write();
        hMultipleResults_Syst   ->Write();
        //
        gROOT->SetBatch(kTRUE);
        //
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            hName   =   Form("hRES_2D_Cond1_Stat_%i",iFit);
            hTitle  =   Form("Conditional Spectrum in p_{T} [%.2f-%.2f] Stat Err",fArrPT2D[iFit],fArrPT2D[iFit+1]);
            hRES_2D_Cond1_Stat.at(iFit)->SetNameTitle(hName,hTitle);
            hRES_2D_Cond1_Stat.at(iFit)->Write();
            hName   =   Form("hRES_2D_Cond1_Syst_%i",iFit);
            hTitle  =   Form("Conditional Spectrum in p_{T} [%.2f-%.2f] Syst Err",fArrPT2D[iFit],fArrPT2D[iFit+1]);
            hRES_2D_Cond1_Syst.at(iFit)->SetNameTitle(hName,hTitle);
            hRES_2D_Cond1_Syst.at(iFit)->Write();
        }
        auto cDrawResult = uPlotSpectrum(hRES_1D_Stat,hRES_1D_Syst,"1D");
        cDrawResult ->  SetLogx();
        cDrawResult ->  SetLeftMargin(0.16);
        cDrawResult ->  SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),"/1D/Yield_1D.pdf"));
        delete      cDrawResult;
        //
        TCanvas    *cDrawFullResults    =   new TCanvas("","",900,1200);
        cDrawFullResults->Divide(3,4);
        cDrawResult =   uPlotSpectrum(hRES_2D_Cond2_Stat,hRES_2D_Cond2_Syst,"12D");
        uLatex->DrawLatexNDC(0.2,0.53,Form("p_{T,#phi_{2}} #in [%.1f;%.1f]",0.,fArrPT2D[0]));
        cDrawFullResults    ->cd    (1);
        cDrawResult ->  SetLogx();
        cDrawResult ->  SetLeftMargin(0.16);
        cDrawResult ->  DrawClonePad();
        cDrawResult ->  SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),Form("/2D/Yield_2D_%i.pdf",-1)));
        delete          cDrawResult;
        for ( Int_t iHisto = 0; iHisto < hRES_2D_Cond1_Stat.size(); iHisto++ )  {
            cDrawResult =   uPlotSpectrum(hRES_2D_Cond1_Stat.at(iHisto),hRES_2D_Cond1_Syst.at(iHisto),"12D");
            uLatex->DrawLatexNDC(0.2,0.53,Form("p_{T,#phi_{2}} #in [%.1f;%.1f]",fArrPT2D[iHisto],fArrPT2D[iHisto+1]));
            cDrawFullResults    ->cd    (iHisto+2);
            cDrawResult ->  SetLogx();
            cDrawResult ->  SetLeftMargin(0.16);
            cDrawResult ->  DrawClonePad();
            cDrawResult ->  SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),Form("/2D/Yield_2D_%i.pdf",iHisto)));
            delete          cDrawResult;
        }
        cDrawFullResults->SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),Form("/2D/Yield_2D.pdf")));
        delete      cDrawFullResults;
        //
        cDrawResult = new TCanvas();
        uPlotSpectrum(hRES_1D_Stat,hRES_1D_Syst,"1D");
        
        cDrawResult ->  SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),"/Full/Production.pdf"));
        delete      cDrawResult;
        //
        gROOT->SetBatch(kFALSE);
        //
        outFil2->Close();
    }
    gROOT->SetBatch(kFALSE);
    //
    // >-> Close input File
    //
    insFile_DT_Yield    ->Close();
    insFile_EF_Yield    ->Close();
    //
    gErrorIgnoreLevel   =   kInfo;
}
