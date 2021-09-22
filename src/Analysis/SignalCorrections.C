// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void SignalCorrections ( TString fOption = "", bool fSilent = true, TString kFolder = "" )
{
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
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (Form(kASigExtr_FitCheckRst,"Yield"));
    TFile*  insFile_EF_Yield            =   new TFile   (Form(kAnalysis_MCTruthHist,"Yield"));
    //
    //  Generating Folder structure
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtp_Dir,"Yield")));
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield")));
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield"))+TString("/1D"));
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield"))+TString("/2D"));
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtp_Plot_Direct,"Yield"))+TString("/Full"));
    
    
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
    TGraphMultiErrors  *gConditional_Mean_PT    =   new TGraphMultiErrors(hName,hTitle,nBinPT2D+2,2);
    //
    hName   =   Form("gYieldResult");
    hTitle  =   Form("gYieldResult");
    TGraphMultiErrors  *gYieldResult            =   new TGraphMultiErrors(hName,hTitle,2,2);
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
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
    gYieldResult            ->  SetPoint        ( 0,    1,      fResults[0]  );
    gYieldResult            ->  SetPointEX      ( 0,    0.1,    0.1          );
    gYieldResult            ->  SetPointEY      ( 0,    0,      fResults[1],    fResults[2] );
    gYieldResult            ->  SetPointEY      ( 0,    1,      fResults[3],    fResults[4] );
    gConditional_Mean_PT    ->  SetPoint        ( 0,    -2,     fResults[5] );
    gConditional_Mean_PT    ->  SetPointEX      ( 0,    0.5,    0.5         );
    gConditional_Mean_PT    ->  SetPointEY      ( 0,    0,      fResults[6],    fResults[7] );
    gConditional_Mean_PT    ->  SetPointEY      ( 0,    1,      fResults[8],    fResults[9] );
    //
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        //  Progressive Count
        fProgrCount++;
        //  Print loop Timer
        fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
        //
        fResults = fMeasureFullYield(hRES_2D_Cond1_Stat.at(iFit),hRES_2D_Cond1_Syst.at(iFit),Form("1D_2D_%i",iFit),Form(kASigExtp_Plot_Direct,"Yield"));
        //
        auto    fCenter =   hRES_2D_Cond2_Stat->    GetBinCenter(iFit+1);
        auto    fXError =   hRES_2D_Cond2_Stat->    GetBinLowEdge(iFit+2) - fCenter;
        gConditional_Mean_PT    ->  SetPoint        ( iFit+1,   fCenter,    fResults[5] );
        gConditional_Mean_PT    ->  SetPointEX      ( iFit+1,   fXError,    fXError     );
        gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   0,          fResults[6],    fResults[7] );
        gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   1,          fResults[8],    fResults[9] );
        //
        hRES_2D_Cond2_Stat  ->  SetBinContent   ( iFit+1, fResults[10] );
        hRES_2D_Cond2_Stat  ->  SetBinError     ( iFit+1, fResults[11] );
        hRES_2D_Cond2_Syst  ->  SetBinContent   ( iFit+1, fResults[10] );
        hRES_2D_Cond2_Syst  ->  SetBinError     ( iFit+1, fResults[12] );
        //
        hRES_2D_Cond3_Stat  ->  SetBinContent   ( iFit+1, fResults[13] );
        hRES_2D_Cond3_Stat  ->  SetBinError     ( iFit+1, fResults[14] );
        hRES_2D_Cond3_Syst  ->  SetBinContent   ( iFit+1, fResults[13] );
        hRES_2D_Cond3_Syst  ->  SetBinError     ( iFit+1, fResults[15] );
    }
    //
    auto fResult2       = fMeasureFullYield(hRES_2D_Cond2_Stat,hRES_2D_Cond2_Syst,"2D");
    //
    auto fExtrap_1_Stat =   0.;
    auto fExtrap_1_Syst =   0.;
    auto fExtrap_1_Val_ =   hRES_2D_Cond2_Stat->IntegralAndError(-1.,10000.,fExtrap_1_Stat,"width");
         fExtrap_1_Val_ =   hRES_2D_Cond2_Syst->IntegralAndError(-1.,10000.,fExtrap_1_Syst,"width");
    //
    cout << fExtrap_1_Val_ << endl;
    //
    auto fExtrap_2_Val_ =   fResult2[10];
    auto fExtrap_2_Stat =   fResult2[11];
    auto fExtrap_2_Syst =   fResult2[12];
    //
    cout << fExtrap_2_Val_ << endl;
    //
    auto fIntegral_Stat =   0.;
    auto fIntegral_Syst =   0.;
    auto fIntegral_Val_ =   uHistoIntegralAndError(hRES_2D_Cond1_Stat,fIntegral_Stat);
         fIntegral_Val_ =   uHistoIntegralAndError(hRES_2D_Cond1_Syst,fIntegral_Syst);
    //
    cout << fIntegral_Val_ << endl;
    cout << fIntegral_Val_ + fExtrap_1_Val_ + fExtrap_2_Val_/2. << endl;
    //
    gYieldResult            ->  SetPoint        ( 1,    2,      fIntegral_Val_ + fExtrap_1_Val_ + fExtrap_2_Val_/2. );
    gYieldResult            ->  SetPointEX      ( 1,    0.1,    0.1         );
    gYieldResult            ->  SetPointEY      ( 1,    0,      SquareSum({ fIntegral_Stat, fExtrap_1_Stat, fExtrap_2_Stat/2. }),   SquareSum({ fIntegral_Stat, fExtrap_1_Stat, fExtrap_2_Stat/2. }));
    gYieldResult            ->  SetPointEY      ( 1,    1,      SquareSum({ fIntegral_Syst, fExtrap_1_Syst, fExtrap_2_Syst/2. }),   SquareSum({ fIntegral_Syst, fExtrap_1_Syst, fExtrap_2_Syst/2. }));
    //
    fStopTimer("Fit_for_extrapolation");
    //
    outCheckFitYld  ->  Close();
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
        //
        gYieldResult            ->Write();
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
        for ( Int_t iHisto = 0; iHisto < hRES_2D_Cond1_Stat.size(); iHisto++ )  {
            cDrawResult =   uPlotSpectrum(hRES_2D_Cond1_Stat.at(iHisto),hRES_2D_Cond1_Syst.at(iHisto),"12D");
            cDrawFullResults    ->cd    (iHisto+1);
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
        gPad->SetLogy();
        uSetHisto(gYieldResult);
        TLegend    *lLegend =   new TLegend();
        lLegend->AddEntry(gYieldResult,"Point","P");
        gYieldResult->Draw("APS ; 2 ; 5");
        lLegend->Draw("same");
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
}
