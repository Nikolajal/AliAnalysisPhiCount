// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void SignalCorrections ( bool fSilent = true, TString fOption = "" )
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
    
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> GENERAL ANALYSIS //
    //
    TH1D       *hEvntEff;
    TH1D       *hEvntMlt;
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "fQC_Event_Enum_Mult";
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
    TH1F       *hTRU_RECVTX_1D;
    TH1F       *hTRU_ALLVTX_1D;
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
    hName       =   "hTRU_INELFLL_1D";
    hTRU_RECVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hTRU_INELVTX_1D";
    hTRU_ALLVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
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
    /*
    // >-> MULTIPLICITY ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hRAW_1D_in_Mlt;
    TH1F      **hREC_1D_in_Mlt;
    TH1F      **hGEN_1D_in_Mlt;
    TH1F      **hTRU_RECVTX_1D_in_Mlt;
    TH1F      **hTRU_ALLVTX_1D_in_Mlt;
    TH1F      **hREC_1D_in_2D_bin_in_Mlt;
    TH1F      **hGEN_1D_in_2D_bin_in_Mlt;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hRAW_1D_in_Mlt              =   new TH1F*[nBinMult+1];
    hREC_1D_in_Mlt              =   new TH1F*[nBinMult+1];
    hGEN_1D_in_Mlt              =   new TH1F*[nBinMult+1];
    hTRU_RECVTX_1D_in_Mlt       =   new TH1F*[nBinMult+1];
    hTRU_ALLVTX_1D_in_Mlt       =   new TH1F*[nBinMult+1];
    hREC_1D_in_2D_bin_in_Mlt    =   new TH1F*[nBinMult+1];
    hGEN_1D_in_2D_bin_in_Mlt    =   new TH1F*[nBinMult+1];
    for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
        hName       =   Form("hRAW_1D_in_Mlt_%i",iMlt);
        hRAW_1D_in_Mlt[iMlt]     =   (TH1F*)(insFile_DT_Mltty->Get(hName));
        hRAW_1D_in_Mlt[iMlt]->GetEntries();
        //
        hName       =   Form("hREC_1D_in_Mlt_%i",iMlt);
        hREC_1D_in_Mlt[iMlt]     =   (TH1F*)(insFile_EF_Mltty->Get(hName));
        //
        hName       =   Form("hGEN_1D_in_Mlt_%i",iMlt);
        hGEN_1D_in_Mlt[iMlt]     =   (TH1F*)(insFile_EF_Mltty->Get(hName));
        //
        hName       =   Form("hREC_1D_in_Mlt_in_2Dbin_%i",iMlt);
        hREC_1D_in_2D_bin_in_Mlt[iMlt]  =   (TH1F*)(insFile_EF_Mltty->Get(hName));
        //
        hName       =   Form("hGEN_1D_in_Mlt_in_2Dbin_%i",iMlt);
        hGEN_1D_in_2D_bin_in_Mlt[iMlt]  =   (TH1F*)(insFile_EF_Mltty->Get(hName));
    }
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hRAW_2D_in_Mlt;
    //
    //  Defining cumulative histogram over measurable pT
    //
    /*
    hRAW_2D_in_Mlt                     =   new TH2F*[nBinMult];
    for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
        hName               =   Form("hRAW_2D_in_Mlt_%i",iMlt);
        hRAW_2D_in_Mlt[iMlt]=   (TH2F*)(insFile_DT_Mltty->Get(hName));
    }
     */
    /*
    //
    // >-> RAPIDITY ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hRAW_1D_in_Rap;
    TH1F      **hREC_1D_in_Rap;
    TH1F      **hGEN_1D_in_Rap;
    TH1F      **hTRU_RECVTX_1D_in_Rap;
    TH1F      **hTRU_ALLVTX_1D_in_Rap;
    TH1F      **hREC_1D_in_2D_bin_in_Rap;
    TH1F      **hGEN_1D_in_2D_bin_in_Rap;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hRAW_1D_in_Rap              =   new TH1F*[nBinRap_];
    hREC_1D_in_Rap              =   new TH1F*[nBinRap_];
    hGEN_1D_in_Rap              =   new TH1F*[nBinRap_];
    hTRU_RECVTX_1D_in_Rap       =   new TH1F*[nBinRap_];
    hTRU_ALLVTX_1D_in_Rap       =   new TH1F*[nBinRap_];
    hREC_1D_in_2D_bin_in_Rap    =   new TH1F*[nBinRap_];
    hGEN_1D_in_2D_bin_in_Rap    =   new TH1F*[nBinRap_];
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ ) {
        hName       =   Form("hRAW_1D_in_Rap_%i",iRap);
        hRAW_1D_in_Rap[iRap]     =   (TH1F*)(insFile_DT_Rapty->Get(hName));
        //
        hName       =   Form("hREC_1D_in_RP_%i",iRap);
        hREC_1D_in_Rap[iRap]     =   (TH1F*)(insFile_EF_Rapty->Get(hName));
        //
        hName       =   Form("hGEN_1D_in_RP_%i",iRap);
        hGEN_1D_in_Rap[iRap]     =   (TH1F*)(insFile_EF_Rapty->Get(hName));
        //
        hName       =   Form("hREC_1D_in_RP_in_2Dbin_%i",iRap);
        hREC_1D_in_2D_bin_in_Rap[iRap]  =   (TH1F*)(insFile_EF_Rapty->Get(hName));
        //
        hName       =   Form("hGEN_1D_in_RP_in_2Dbin_%i",iRap);
        hGEN_1D_in_2D_bin_in_Rap[iRap]  =   (TH1F*)(insFile_EF_Rapty->Get(hName));
    }
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hRAW_2D_in_Rap;
    //
    //  Defining cumulative histogram over measurable pT
    //
    /*
    hRAW_2D_in_Rap                     =   new TH2F*[nBinRap_];
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ ) {
        hName               =   Form("hRAW_2D_in_Rap_%i",iRap);
        hRAW_2D_in_Rap[iRap]=   (TH2F*)(insFile_DT_Rapty->Get(hName));
    }
     */
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
    auto        kN_Trg          =   (hEvntEff->GetBinContent(1));//kEventCount::kTrigger));
    auto        kN_Vtx          =   (hEvntEff->GetBinContent(7));//kEventCount::kVertex));
    auto        kN_MB           =   (hEvntEff->GetBinContent(8));//kEventCount::kVertex10));
    Double_t    f1DCorrection   =   (1./kBR)        *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    //
    // >-> YIELD ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F   *hRES_1D_Stat    =  /* fEfficiencycorrection   (*/ fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection)/*,hTRU_ALLVTX_1D,hTRU_RECVTX_1D )*/;
    TH1F   *hRES_1D_Syst    =   fSetSystErrors          ( hRES_1D_Stat );
    //
    auto cDrawResult = uPlotSpectrum(hRES_1D_Stat,hRES_1D_Syst);
    cDrawResult ->  SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),"/Yield_1D.pdf"));
    delete      cDrawResult;
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    std::vector<TH1F*>  hRES_2D_Cond1_Stat  =   fEfficiencycorrection   ( hRAW_2D,hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
    std::vector<TH1F*>  hRES_2D_Cond1_Syst  =   fSetSystErrors          ( hRES_2D_Cond1_Stat );
    //
    TCanvas    *cDrawFullResults    =   new TCanvas("","",900,1200);
    cDrawFullResults->Divide(3,4);
    for ( Int_t iHisto = 0; iHisto < hRES_2D_Cond1_Stat.size(); iHisto++ )  {
        cDrawResult = uPlotSpectrum(hRES_2D_Cond1_Stat.at(iHisto),hRES_2D_Cond1_Syst.at(iHisto));
        cDrawFullResults->cd(iHisto+1);
        cDrawResult->DrawClonePad();
        cDrawResult ->  SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),Form("/Yield_2D_%i.pdf",iHisto)));
        delete      cDrawResult;
    }
    cDrawFullResults->SaveAs(Form("%s%s",Form(kASigExtp_Plot_Direct,"Yield"),Form("/Yield_2D.pdf")));
    delete      cDrawFullResults;
    //
    hName   =   Form("hRES_2D_Cond2_Stat");
    hTitle  =   Form("hRES_2D_Cond2_Stat");
    TH1F       *hRES_2D_Cond2_Stat  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("hRES_2D_Cond2_Syst");
    hTitle  =   Form("hRES_2D_Cond2_Syst");
    TH1F       *hRES_2D_Cond2_Syst  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("gConditional_Mean_PT");
    hTitle  =   Form("gConditional_Mean_PT");
    TGraphMultiErrors  *gConditional_Mean_PT    =   new TGraphMultiErrors(hName,hTitle,nBinPT2D+2,2);
    //
    // >-> RAPIDITY ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    /*
    TH1F  **hRES_1D_Stat_in_Rap    =   new TH1F*[nBinRap_];//fEfficiencycorrection   ( fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection),hTRU_RECVTX_1D,hTRU_ALLVTX_1D );
    TH1F  **hRES_1D_Syst_in_Rap    =   new TH1F*[nBinRap_];//fSetSystErrors          ( hRES_1D_Stat );
    //
    TH1F  **hRES_1D_Stat_in_Mlt    =   new TH1F*[nBinMult+1];//fEfficiencycorrection   ( fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection),hTRU_RECVTX_1D,hTRU_ALLVTX_1D );
    TH1F  **hRES_1D_Syst_in_Mlt    =   new TH1F*[nBinMult+1];//fSetSystErrors          ( hRES_1D_Stat );
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    std::vector<TH1F*>* hRES_2D_Cond1_Stat_in_Rap  =   new std::vector<TH1F*> [nBinRap_];//fEfficiencycorrection   ( hRAW_2D,hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
    std::vector<TH1F*>* hRES_2D_Cond1_Syst_in_Rap  =   new std::vector<TH1F*> [nBinRap_];//fSetSystErrors          ( hRES_2D_Cond1_Stat );
    TH1F      **hRES_2D_Cond2_Stat_in_Rap           =   new TH1F*[nBinRap_];
    TH1F      **hRES_2D_Cond2_Syst_in_Rap           =   new TH1F*[nBinRap_];
    //
    std::vector<TH1F*>* hRES_2D_Cond1_Stat_in_Mlt  =   new std::vector<TH1F*> [nBinMult+1];//fEfficiencycorrection   ( hRAW_2D,hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
    std::vector<TH1F*>* hRES_2D_Cond1_Syst_in_Mlt  =   new std::vector<TH1F*> [nBinMult+1];//fSetSystErrors          ( hRES_2D_Cond1_Stat );
    TH1F      **hRES_2D_Cond2_Stat_in_Mlt           =   new TH1F*[nBinMult+1];
    TH1F      **hRES_2D_Cond2_Syst_in_Mlt           =   new TH1F*[nBinMult+1];

    //
    /*
    hName   =   Form("");
    hTitle  =   Form("");
    TH1F       *hRES_2D_Cond2_Stat  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("");
    hTitle  =   Form("");
    TH1F       *hRES_2D_Cond2_Syst  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName   =   Form("");
    hTitle  =   Form("");
    TGraphMultiErrors  *gConditional_Mean_PT    =   new TGraphMultiErrors(hName,hTitle,nBinPT2D+2,2);
    //
     */
    /*
    for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
        auto    kN_INELgt0  =   (kN_MB);//fEvaluateINELgt0(iMlt-1,hEvntMlt);
        cout << "INEL0 : " << kN_INELgt0 << endl;
        if ( iMlt != 0 )   { kN_INELgt0 *=   (fArrMult[iMlt]-fArrMult[iMlt-1])/(100.*kMultTrgEff[iMlt-1]); }
        cout << "INEL0 : " << kN_INELgt0 << endl;
        cout << "INELS : " << fEvaluateINELgt0(iMlt-1,hEvntMlt) << endl;
        f1DCorrection   =   (1./kBR)        *(1./kN_INELgt0);
        f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_INELgt0);
        cout << " - - " << f1DCorrection << " - - " << endl;
        //
        hRES_1D_Stat_in_Mlt[iMlt] =   fEfficiencycorrection ( hRAW_1D_in_Mlt[iMlt],hREC_1D,hGEN_1D,f1DCorrection);
        hRES_1D_Syst_in_Mlt[iMlt] =   fSetSystErrors        ( hRES_1D_Stat_in_Mlt[iMlt] );
        /*
        hRES_2D_Cond1_Stat_in_Rap[iMlt] =   fEfficiencycorrection   ( hRAW_2D_in_Mlt[iMlt],hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
        hRES_2D_Cond1_Syst_in_Rap[iMlt] =   fSetSystErrors          ( hRES_2D_Cond1_Stat_in_Mlt[iMlt] );
        */
    /*
        hName   =   Form("");
        hTitle  =   Form("");
        hRES_2D_Cond2_Stat_in_Mlt[iMlt]  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
        hRES_2D_Cond2_Syst_in_Mlt[iMlt]  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    }
    TH1F      **hPubMult =   new TH1F*[nMltTrgECls+2];
    //
    for ( Int_t iMult = 0; iMult <  nMltTrgECls; iMult++ )   {
        hName           =   Form("h%i_%i",(int)kMltTrgECls[iMult],(int)kMltTrgECls[iMult+1]);
        hPubMult[iMult] =   (TH1F*)(insPublishedRslt->Get(hName));
    }
    //
    /*
    TH1F      **hPubAdpt  =  new TH1F*[nBinMult+1];
    //
    hName           =   Form("h0_100");
    hPubAdpt[0] =   new TH1F(*((TH1F*)(insPublishedRslt->Get(hName))));
    hPubAdpt[1] =   new TH1F(*hPubMult[0]);
    for ( Int_t iMult = 1; iMult < nBinMult; iMult++  )   {
        auto    h1 = *hPubMult[2*iMult-1];
        auto    h2 = *hPubMult[2*iMult];
        h1.Scale(kMltTrgECls[2*iMult+1]-kMltTrgECls[2*iMult]);
        h2.Scale(kMltTrgECls[2*iMult]-kMltTrgECls[2*iMult-1]);
        auto h3 = h1+h2;
        h3.Scale(1./(kMltTrgECls[2*iMult+1]-kMltTrgECls[2*iMult-1]));
        hPubAdpt[iMult+1] =   new TH1F(h3);
    }
    */
    //
    // Check singular Integral and extrapolate
    // Make a full Yield in multiplicity
    //
    //float fSResult [] = {0.0341,0.1052,0.078950000,0.057233333,0.038350000,0.019580000};
    /*
    for ( Int_t iMult = 0; iMult < nBinMult-1; iMult++ )   {
        
        hPubMult[iMult]->SetMarkerColor(kRed);
        hPubMult[iMult]->SetMarkerStyle(22);
        
        auto    fCheckTru   =   new TH1F(*hPubMult[iMult]);
        auto    fCheckTr2   =   new TH1F(*hPubMult[iMult]);
        fSetFunction(fLevyTsallis);
        fCheckTru->Fit(fLevyTsallis);
        auto    fURaFailr   =   new TH1F(*hRES_1D_Stat_in_Mlt[iMult+1]);
        auto    fURaFail2   =   new TH1F(*hRES_1D_Stat_in_Mlt[iMult+1]);
        fURaFail2->Divide(fCheckTru);
        fCheckTr2->Divide(fCheckTru);
        
        fURaFailr->SetMarkerColor(kBlue);
        fURaFailr->SetMarkerStyle(23);
        fURaFailr->SetMarkerColor(kBlue);
        fURaFailr->SetMarkerStyle(23);
        
        TCanvas *fCanvasCheck =   new TCanvas();
        fCanvasCheck->Divide(2,1);
        fCanvasCheck->cd(1);
        gPad->SetLogy();
        fCheckTru->Draw();
        fURaFailr->Draw("same");
        gPad->BuildLegend();
        fCanvasCheck->cd(2);
        fURaFail2->Draw();
        fCheckTr2->Draw("same");
        fCanvasCheck->SaveAs(Form("pppp_%i.pdf",iMult));
        delete fCanvasCheck;
    }
    /*
     hName   =   Form("");
     hTitle  =   Form("");
     TH1F       *hRES_2D_Cond2_Stat  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
     //
     hName   =   Form("");
     hTitle  =   Form("");
     TH1F       *hRES_2D_Cond2_Syst  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
     //
     hName   =   Form("");
     hTitle  =   Form("");
     TGraphMultiErrors  *gConditional_Mean_PT    =   new TGraphMultiErrors(hName,hTitle,nBinPT2D+2,2);
     //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ ) {
        //hRES_1D_Stat_in_Rap[iRap] =   fEfficiencycorrection ( hRAW_1D_in_Rap[iRap],hREC_1D,hGEN_1D,f1DCorrection);
        //hRES_1D_Syst_in_Rap[iRap] =   fSetSystErrors        ( hRES_1D_Stat_in_Rap[iRap] );
        
        hRES_2D_Cond1_Stat_in_Rap[iRap] =   fEfficiencycorrection   ( hRAW_2D_in_Rap[iRap],hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
        hRES_2D_Cond1_Syst_in_Rap[iRap] =   fSetSystErrors          ( hRES_2D_Cond1_Stat_in_Rap[iRap] );
        
        hName   =   Form("");
        hTitle  =   Form("");
        hRES_2D_Cond2_Stat_in_Rap[iRap]  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
        hRES_2D_Cond2_Syst_in_Rap[iRap]  =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    }
     */
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    /*
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
    auto    fResults = fMeasureFullYield(hRES_1D_Stat,hRES_1D_Syst,"1D");
    //  Progressive Count
    fProgrCount++;
    //  Print loop Timer
    fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
    //
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
        fResults = fMeasureFullYield(hRES_2D_Cond1_Stat.at(iFit),hRES_2D_Cond1_Syst.at(iFit),Form("1D_2D_%i",iFit));
        //
        auto    fCenter =   hRES_2D_Cond2_Stat->    GetBinCenter(iFit+1);
        auto    fXError =   hRES_2D_Cond2_Stat->    GetBinLowEdge(iFit+2) - fCenter;
        gConditional_Mean_PT    ->  SetPoint        ( iFit+1,   fCenter,    fResults[5] );
        gConditional_Mean_PT    ->  SetPointEX      ( iFit+1,   fXError,    fXError     );
        gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   0,          fResults[6],    fResults[7] );
        gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   1,          fResults[8],    fResults[9] );
        //
        hRES_2D_Cond2_Stat  ->  SetBinContent   ( iFit+1, fResults[0] );
        hRES_2D_Cond2_Stat  ->  SetBinError     ( iFit+1, fResults[1] );
        hRES_2D_Cond2_Syst  ->  SetBinContent   ( iFit+1, fResults[0] );
        hRES_2D_Cond2_Syst  ->  SetBinError     ( iFit+1, fResults[3] );
    }
    //
    fResults = fMeasureFullYield(hRES_2D_Cond2_Stat,hRES_2D_Cond2_Syst,"2D");
    gConditional_Mean_PT    ->  SetPoint        ( nBinPT2D+1,   -1,     fResults[5] );
    gConditional_Mean_PT    ->  SetPointEX      ( nBinPT2D+1,   0.5,    0.5         );
    gConditional_Mean_PT    ->  SetPointEY      ( nBinPT2D+1,   0,      fResults[6],    fResults[7] );
    gConditional_Mean_PT    ->  SetPointEY      ( nBinPT2D+1,   1,      fResults[8],    fResults[9] );
    //
    fStopTimer("Fit_for_extrapolation");
    //
    
    //  -------------------------------------------------------- //
    /*
    //
    fStartTimer("Fit_for_extrapolation");
    //
    // Output File for Fit Check
    TFile*  outCheckFitMlt  =   new TFile(fMltSigCh2k,"recreate");
    //
    // Total Fit number and progressive
    //
    fTotalCount = nBinMult * ( 2 + nBinPT2D );
    fProgrCount = 0;
    //
    for ( Int_t iMlt = 0; iMlt < nBinMult; iMlt++ ) {
        //auto fResults = fMeasureFullYield(hRES_1D_Stat_in_Mlt[iMlt],hRES_1D_Syst_in_Mlt[iMlt],Form("1D_in_Mlt_%i",iMlt));
        //  Progressive Count
        fProgrCount++;
        //  Print loop Timer
        fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
        //
        /*
        gConditional_Mean_PT    ->  SetPoint        ( 0,    -2,     fResults[5] );
        gConditional_Mean_PT    ->  SetPointEX      ( 0,    0.5,    0.5         );
        gConditional_Mean_PT    ->  SetPointEY      ( 0,    0,      fResults[6],    fResults[7] );
        gConditional_Mean_PT    ->  SetPointEY      ( 0,    1,      fResults[8],    fResults[9] );
         */
        /*
        fYield_Stat ->  SetPoint        (0,1,fResults[0]);
        fYield_Syst ->  SetPoint        (0,1,fResults[0]);
        fYield_Stat ->  SetPointError   (0,0,0,fResults[1],fResults[2]);
        fYield_Syst ->  SetPointError   (0,0,0,fResults[3],fResults[4]);
         */
        /*
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            fResults = fMeasureFullYield(hRES_2D_Cond1_Stat_in_Rap[iRap].at(iFit),hRES_2D_Cond1_Syst_in_Rap[iRap].at(iFit),Form("1D_2D_%i_in_Rap_%i",iFit,iRap));
            //
            auto    fCenter =   hRES_2D_Cond2_Stat->    GetBinCenter(iFit+1);
            auto    fXError =   hRES_2D_Cond2_Stat->    GetBinLowEdge(iFit+2) - fCenter;
            gConditional_Mean_PT    ->  SetPoint        ( iFit+1,   fCenter,    fResults[5] );
            gConditional_Mean_PT    ->  SetPointEX      ( iFit+1,   fXError,    fXError     );
            gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   0,          fResults[6],    fResults[7] );
            gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   1,          fResults[8],    fResults[9] );
            //
            hRES_2D_Cond2_Stat_in_Rap[iRap]  ->  SetBinContent   ( iFit+1, fResults[0] );
            hRES_2D_Cond2_Stat_in_Rap[iRap]  ->  SetBinError     ( iFit+1, fResults[1] );
            hRES_2D_Cond2_Syst_in_Rap[iRap]  ->  SetBinContent   ( iFit+1, fResults[0] );
            hRES_2D_Cond2_Syst_in_Rap[iRap]  ->  SetBinError     ( iFit+1, fResults[3] );
        }
        //
        fResults = fMeasureFullYield(hRES_2D_Cond2_Stat_in_Rap[iRap],hRES_2D_Cond2_Syst_in_Rap[iRap],Form("2D_in_Rap_%i",iRap));
        
        hTest_Stat->SetBinContent(  iRap+1, fResults[0]);
        hTest_Stat->SetBinError(    iRap+1, fResults[1]);
        hTest_Syst->SetBinContent(  iRap+1, fResults[0]);
        hTest_Syst->SetBinError(    iRap+1, fResults[3]);
        
        gConditional_Mean_PT    ->  SetPoint        ( nBinPT2D+1,   -1,     fResults[5] );
        gConditional_Mean_PT    ->  SetPointEX      ( nBinPT2D+1,   0.5,    0.5         );
        gConditional_Mean_PT    ->  SetPointEY      ( nBinPT2D+1,   0,      fResults[6],    fResults[7] );
        gConditional_Mean_PT    ->  SetPointEY      ( nBinPT2D+1,   1,      fResults[8],    fResults[9] );
         *//*
    }
    //
    fStopTimer("Fit_for_extrapolation");
    //
    fStartTimer("Fit_for_extrapolation");
    //
    // Output File for Fit Check
    TFile*  outCheckFitRap  =   new TFile(fRapSigCh2k,"recreate");
    //
    // Total Fit number and progressive
    //
    fTotalCount = nBinRap_ * ( 2 + nBinPT2D );
    fProgrCount = 0;
    //
    TH1F   *hTest_Stat  =   new TH1F("testStat","",nBinRap_,fArrRap_);
    TH1F   *hTest_Syst  =   new TH1F("testSyst","",nBinRap_,fArrRap_);
    
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ ) {
        //fResults = fMeasureFullYield(hRES_1D_Stat_in_Rap[iRap],hRES_1D_Syst_in_Rap[iRap],Form("1D_in_Rap_%i",iRap));
        //  Progressive Count
        //fProgrCount++;
        //  Print loop Timer
        //fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
        //
        /*
        gConditional_Mean_PT    ->  SetPoint        ( 0,    -2,     fResults[5] );
        gConditional_Mean_PT    ->  SetPointEX      ( 0,    0.5,    0.5         );
        gConditional_Mean_PT    ->  SetPointEY      ( 0,    0,      fResults[6],    fResults[7] );
        gConditional_Mean_PT    ->  SetPointEY      ( 0,    1,      fResults[8],    fResults[9] );
         */
        /*
        fYield_Stat ->  SetPoint        (0,1,fResults[0]);
        fYield_Syst ->  SetPoint        (0,1,fResults[0]);
        fYield_Stat ->  SetPointError   (0,0,0,fResults[1],fResults[2]);
        fYield_Syst ->  SetPointError   (0,0,0,fResults[3],fResults[4]);
         */
        /*
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            fResults = fMeasureFullYield(hRES_2D_Cond1_Stat_in_Rap[iRap].at(iFit),hRES_2D_Cond1_Syst_in_Rap[iRap].at(iFit),Form("1D_2D_%i_in_Rap_%i",iFit,iRap));
            //
            auto    fCenter =   hRES_2D_Cond2_Stat->    GetBinCenter(iFit+1);
            auto    fXError =   hRES_2D_Cond2_Stat->    GetBinLowEdge(iFit+2) - fCenter;
            gConditional_Mean_PT    ->  SetPoint        ( iFit+1,   fCenter,    fResults[5] );
            gConditional_Mean_PT    ->  SetPointEX      ( iFit+1,   fXError,    fXError     );
            gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   0,          fResults[6],    fResults[7] );
            gConditional_Mean_PT    ->  SetPointEY      ( iFit+1,   1,          fResults[8],    fResults[9] );
            //
            hRES_2D_Cond2_Stat_in_Rap[iRap]  ->  SetBinContent   ( iFit+1, fResults[0] );
            hRES_2D_Cond2_Stat_in_Rap[iRap]  ->  SetBinError     ( iFit+1, fResults[1] );
            hRES_2D_Cond2_Syst_in_Rap[iRap]  ->  SetBinContent   ( iFit+1, fResults[0] );
            hRES_2D_Cond2_Syst_in_Rap[iRap]  ->  SetBinError     ( iFit+1, fResults[3] );
        }
        //
        fResults = fMeasureFullYield(hRES_2D_Cond2_Stat_in_Rap[iRap],hRES_2D_Cond2_Syst_in_Rap[iRap],Form("2D_in_Rap_%i",iRap));
        
        hTest_Stat->SetBinContent(  iRap+1, fResults[0]);
        hTest_Stat->SetBinError(    iRap+1, fResults[1]);
        hTest_Syst->SetBinContent(  iRap+1, fResults[0]);
        hTest_Syst->SetBinError(    iRap+1, fResults[3]);
        
        gConditional_Mean_PT    ->  SetPoint        ( nBinPT2D+1,   -1,     fResults[5] );
        gConditional_Mean_PT    ->  SetPointEX      ( nBinPT2D+1,   0.5,    0.5         );
        gConditional_Mean_PT    ->  SetPointEY      ( nBinPT2D+1,   0,      fResults[6],    fResults[7] );
        gConditional_Mean_PT    ->  SetPointEY      ( nBinPT2D+1,   1,      fResults[8],    fResults[9] );
         */
    /*}*/
    //
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Trigger Analysis
    //
    //TFile *outFil1  =   new TFile   (fTrgSigCorr,"recreate");
    //
    //outFil1->Close();
    //
    // >> Yield Analysis
    //
    TFile *outFil2  =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield"),"recreate");
    //
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
    //
    outFil2->Close();
    //
    // >-> Close input File
    //
    //outCheckFitYld      ->Close();
    //outCheckFitRap      ->Close();
    //outCheckFitMlt      ->Close();
    insFile_DT_Yield    ->Close();
    insFile_EF_Yield    ->Close();
}


/*
TCanvas * fDraw = new TCanvas("","",800,600);
fDraw->Divide(4,3);
fDraw->cd(1);
gPad->SetLogy();
gPad->SetLogx();
fGraphMultiErrors(hRES_1D_Stat,hRES_1D_Syst)->Draw("APS ; PE5 ; PE2 ; PE");
fDraw->cd(2);
gPad->SetLogy();
gPad->SetLogx();
fGraphMultiErrors(hRES_2D_Cond2_Stat,hRES_2D_Cond2_Syst)->Draw("APS ; PE5 ; PE2 ; PE");
for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
    fDraw->cd(iFit+3);
    gPad->SetLogy();
    gPad->SetLogx();
    fGraphMultiErrors(hRES_2D_Cond1_Stat.at(iFit),hRES_2D_Cond1_Syst.at(iFit))->Draw("APS ; PE5 ; PE2 ; PE");
}
 */

/*
 // File for 1-Dimensional Analysis:
 // !TODO: All Set!
 #include "../../inc/AliAnalysisPhiPair.h"
 #include "RooMsgService.h"

 void SignalCorrections ( bool fSilent = true, TString fOption = "" )
 {
     //---------------------//
     //  Setting up input   //
     //---------------------//
     
     //-// OPTIONS
     
     // Silencing warnings for smoother
     if ( fSilent )
     {
         RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
         RooMsgService::instance().setSilentMode(fSilent);
     }
     fChooseOption(fOption);
     
     // Retrieving PreProcessed data histograms
     TFile*  insFile_DT_Yield            =   new TFile   (fYldSigExtr);
     TFile*  insFile_DT_Mult             =   new TFile   (fMltSigExtr);
     TFile*  insFile_EF_Yield            =   new TFile   (fYldPrePrMC);
     TFile*  insFile_EF_Mult             =   new TFile   (fMltPrePrMC);
     TFile*  insPublishedRslt            =   new TFile   ("./result/phi_pp5_mul_24Feb2021.root");
     
     // Recovering the histograms-------------------------------------------------------------------------------

     // >-> GENERAL ANALYSIS //
     //
     TH1D       *hEvntEff;
     TH1D       *hEvntMlt;
     //
     hName       =   "fQC_Event_Enumerate";
     hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
     //
     hName       =   "fQC_Event_Enum_Mult";
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
     TH1F       *hTRU_RECVTX_1D;
     TH1F       *hTRU_ALLVTX_1D;
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
     hName       =   "hTRU_RECVTX_1D";
     hTRU_RECVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
     //
     hName       =   "hTRU_ALLVTX_1D";
     hTRU_ALLVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
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
     hRAW_2D->GetEntries();
     //
     /*
     // >-> MULTIPLICITY ANALYSIS //

     // >->-->-> 1-Dimension analysis //
     //
     //  Declaring all histograms
     //
     TH1F      **hRAW_1D_in_MT               = new TH1F     *[nBinMult+1];
     TH1F       *hREC_1D_in_MT_0;
     TH1F       *hGEN_1D_in_MT_0;
     TH1F       *hREC_1D_in_MT_in_2D_bin_0;
     TH1F       *hGEN_1D_in_MT_in_2D_bin_0;
     //
     //  Defining cumulative histogram over measurable pT
     //
     hName       =   "hREC_1D_in_MT_0";
     hREC_1D_in_MT_0     =   (TH1F*)(insFile_EF_Mult->Get(hName));
     //
     hName       =   "hGEN_1D_in_MT_0";
     hGEN_1D_in_MT_0     =   (TH1F*)(insFile_EF_Mult->Get(hName));
     //
     hName       =   "hREC_1D_in_MT_in_2D_bin_0";
     hREC_1D_in_MT_in_2D_bin_0     =   (TH1F*)(insFile_EF_Mult->Get(hName));
     //
     hName       =   "hGEN_1D_in_MT_in_2D_bin_0";
     hGEN_1D_in_MT_in_2D_bin_0     =   (TH1F*)(insFile_EF_Mult->Get(hName));
     
     //  Defining MT-Differential histograms over measurable pT
     //
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )   {
         hName = Form("hRAW_1D_in_MT_%i",iMult);
         hRAW_1D_in_MT[iMult]   =   (TH1F*)(insFile_DT_Mult->Get(hName));
         hRAW_1D_in_MT[iMult]->Draw();
     }

     // >->-->-> 2-Dimension analysis //
     //
     //  Declaring all histograms
     //
     TH2F      **hRAW_2D_in_MT               = new TH2F     *[nBinMult+1];
     //
     //  Defining MT-Differential histograms over measurable pT
     //
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )   {
         hName   =   Form("hRAW_2D_in_MT_%i",iMult);
         hRAW_2D_in_MT[iMult]   =   (TH2F*)(insFile_DT_Mult->Get(hName));
     }
     //
     //---------------------//
     //  Setting up output  //
     //---------------------//
     //
     // Generating the binning array--------------------------------------------------------------------------
     //
     fSetBinPT1D();
     fSetBinIM1D();
     fSetBinPT2D();
     fSetBinIM2D();
     fSetBinRap_();
     fSetBinMult();
     fSetBinNTup();
     Int_t       U_AccCand[1024];
     Int_t       U_nAccept;
     //
     // Creating the histograms-------------------------------------------------------------------------------
     //
     // >-> YIELD ANALYSIS //
     //
     TGraphAsymmErrors   *fYield_Stat        =   new TGraphAsymmErrors();
     hName       =   Form("fYield_Stat");
     hTitle      =   Form("fYield_Stat");
     fYield_Stat->SetNameTitle(hName,hTitle);
     //
     TGraphAsymmErrors   *fYield_Syst        =   new TGraphAsymmErrors();
     hName       =   Form("fYield_Syst");
     hTitle      =   Form("fYield_Syst");
     fYield_Syst->SetNameTitle(hName,hTitle);
     //
     // >->-->-> 1-Dimension analysis //
     //
     //  Declaring all histograms
     //
     TGraphAsymmErrors   *gRES_1D_Stat       =   new TGraphAsymmErrors();
     hName       =   Form("gRES_1D_Stat");
     hTitle      =   Form("gRES_1D_Stat");
     gRES_1D_Stat->SetNameTitle(hName,hTitle);
     //
     TGraphAsymmErrors   *gRES_1D_Syst       =   new TGraphAsymmErrors();
     hName       =   Form("gRES_1D_Syst");
     hTitle      =   Form("gRES_1D_Syst");
     gRES_1D_Syst->SetNameTitle(hName,hTitle);
     //
     TGraphAsymmErrors   *gMPT_1D_Stat       =   new TGraphAsymmErrors();
     hName       =   Form("gMPT_1D_Stat");
     hTitle      =   Form("gMPT_1D_Stat");
     gMPT_1D_Stat->SetNameTitle(hName,hTitle);
     //
     TGraphAsymmErrors   *gMPT_1D_Syst       =   new TGraphAsymmErrors();
     hName       =   Form("gMPT_1D_Syst");
     hTitle      =   Form("gMPT_1D_Syst");
     gMPT_1D_Syst->SetNameTitle(hName,hTitle);
     //
     // >->-->-> 2-Dimension analysis //
     //
     //  Declaring all histograms
     //
     std::vector<TGraphAsymmErrors*> gRES_2D_Stat;
     std::vector<TGraphAsymmErrors*> gRES_2D_Syst;
     //
     TGraphAsymmErrors  *gYield_Profile_Stat =   new TGraphAsymmErrors();
     TGraphAsymmErrors  *gYield_Profile_Syst =   new TGraphAsymmErrors();
     //
     hName       =   Form("gYield_Profile_Stat");
     hTitle      =   Form("<Y_{#phi#phi}> in p_{T} [%.2f#;%.2f], Stat. Errors",fArrPT2D[0],fArrPT2D[nBinPT2D]);
     gYield_Profile_Stat->SetNameTitle(hName,hTitle);
     //
     hName       =   Form("gYield_Profile_Syst");
     hTitle      =   Form("<Y_{#phi#phi}> in p_{T} [%.2f#;%.2f], Syst. Errors",fArrPT2D[0],fArrPT2D[nBinPT2D]);
     gYield_Profile_Syst->SetNameTitle(hName,hTitle);
     //
     TGraphAsymmErrors  *gMPT_2D_Stat =   new TGraphAsymmErrors();
     TGraphAsymmErrors  *gMPT_2D_Syst =   new TGraphAsymmErrors();
     //
     hName       =   Form("gMPT_2D_Stat");
     hTitle      =   Form("<Y_{#phi#phi}> in p_{T} [%.2f#;%.2f], Stat. Errors",fArrPT2D[0],fArrPT2D[nBinPT2D]);
     gMPT_2D_Stat->SetNameTitle(hName,hTitle);
     //
     hName       =   Form("gMPT_2D_Syst");
     hTitle      =   Form("<Y_{#phi#phi}> in p_{T} [%.2f#;%.2f], Syst. Errors",fArrPT2D[0],fArrPT2D[nBinPT2D]);
     gMPT_2D_Syst->SetNameTitle(hName,hTitle);
     //
     /*
     // >-> MULTIPLICITY ANALYSIS //
     //
     // >->-->-> 1-Dimension analysis //
     //
     TGraphAsymmErrors  **fYield_Stat_in_MT  =   new TGraphAsymmErrors  *[nBinMult+1];
     TGraphAsymmErrors  **fYield_Syst_in_MT  =   new TGraphAsymmErrors  *[nBinMult+1];
     //
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )  {
         fYield_Stat_in_MT[iMult]            =   new TGraphAsymmErrors();
         hName                               =   Form("fYield_Stat_%i",iMult);
         hTitle                              =   Form("fYield_Stat_%i",iMult);
         fYield_Stat_in_MT[iMult]            ->  SetNameTitle(hName,hTitle);
         //
         fYield_Syst_in_MT[iMult]            =   new TGraphAsymmErrors();
         hName                               =   Form("fYield_Syst_%i",iMult);
         hTitle                              =   Form("fYield_Syst_%i",iMult);
         fYield_Syst_in_MT[iMult]            ->  SetNameTitle(hName,hTitle);
     }
     //
     //  Declaring all histograms
     //
     TGraphAsymmErrors  **gRES_1D_Stat_in_MT =   new TGraphAsymmErrors  *[nBinMult+1];
     TGraphAsymmErrors  **gRES_1D_Syst_in_MT =   new TGraphAsymmErrors  *[nBinMult+1];
     //
     //  Defining MT-Differential histograms over measurable pT
     //
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
     {
         gRES_1D_Stat_in_MT[iMult]           =   new TGraphAsymmErrors();
         hName                               =   Form("gRES_1D_Stat_in_MT_%i",iMult);
         hTitle                              =   Form("gRES_1D_Stat_in_MT_%i",iMult);
         gRES_1D_Stat_in_MT[iMult]           ->  SetNameTitle(hName,hTitle);
         //
         gRES_1D_Syst_in_MT[iMult]           =   new TGraphAsymmErrors();
         hName                               =   Form("gRES_1D_Syst_in_MT_%i",iMult);
         hTitle                              =   Form("gRES_1D_Syst_in_MT_%i",iMult);
         gRES_1D_Syst_in_MT[iMult]           ->  SetNameTitle(hName,hTitle);
     }
     // >->-->-> 2-Dimension analysis //
     //
     //  Declaring all histograms
     //
     TGraphAsymmErrors  **fYield_Profile_Stat_in_MT  =   new TGraphAsymmErrors  *[nBinMult+1];
     TGraphAsymmErrors  **fYield_Profile_Syst_in_MT  =   new TGraphAsymmErrors  *[nBinMult+1];
     //
     //  Defining MT-Differential histograms over measurable pT
     //
     fYield_Profile_Stat_in_MT[0]    =   new TGraphAsymmErrors();
     hName       =   Form("fYield_Profile_Stat_in_MT%i",0);
     hTitle      =   Form("<Y_{#phi#phi}> in Mult. [%.2f#;%.2f] in p_{T} [%.2f#;%.2f], Stat. Errors",fArrMult[0],fArrMult[nBinMult],fArrPT2D[0],fArrPT2D[nBinPT2D]);
     fYield_Profile_Stat_in_MT[0]->SetNameTitle(hName,hTitle);
     //
     fYield_Profile_Syst_in_MT[0]    =   new TGraphAsymmErrors();
     hName       =   Form("fYield_Profile_Syst_in_MT%i",0);
     hTitle      =   Form("<Y_{#phi#phi}> in Mult. [%.2f#;%.2f] in p_{T} [%.2f#;%.2f], Stat. Errors",fArrMult[0],fArrMult[nBinMult],fArrPT2D[0],fArrPT2D[nBinPT2D]);
     fYield_Profile_Syst_in_MT[0]->SetNameTitle(hName,hTitle);
     
     for ( Int_t iMult = 1; iMult <= nBinMult; iMult++ )
     {
         fYield_Profile_Stat_in_MT[iMult]    =   new TGraphAsymmErrors();
         hName       =   Form("fYield_Profile_Stat_in_MT%i",iMult);
         hTitle      =   Form("<Y_{#phi#phi}> in Mult. [%.2f#;%.2f] in p_{T} [%.2f#;%.2f], Stat. Errors",fArrMult[iMult-1],fArrMult[iMult],fArrPT2D[0],fArrPT2D[nBinPT2D]);
         fYield_Profile_Stat_in_MT[iMult]->SetNameTitle(hName,hTitle);
         //
         fYield_Profile_Syst_in_MT[iMult]    =   new TGraphAsymmErrors();
         hName       =   Form("fYield_Profile_Syst_in_MT%i",iMult);
         hTitle      =   Form("<Y_{#phi#phi}> in Mult. [%.2f#;%.2f] in p_{T} [%.2f#;%.2f], Stat. Errors",fArrMult[iMult-1],fArrMult[iMult],fArrPT2D[0],fArrPT2D[nBinPT2D]);
         fYield_Profile_Syst_in_MT[iMult]->SetNameTitle(hName,hTitle);
     }
     //
     std::vector<TGraphAsymmErrors*>*gRES_2D_Stat_in_MT  =   new std::vector<TGraphAsymmErrors*> [nBinMult+1];
     std::vector<TGraphAsymmErrors*>*gRES_2D_Syst_in_MT  =   new std::vector<TGraphAsymmErrors*> [nBinMult+1];
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
     Double_t    f1DCorrection   =   (1./kBR)        *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
     Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
     //
     gRES_1D_Stat                =   fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection);
     gRES_1D_Stat                =   fEfficiencycorrection(gRES_1D_Stat,hTRU_RECVTX_1D,hTRU_ALLVTX_1D);
     hName       =   Form("gRES_1D_Stat");
     hTitle      =   Form("<Y_{#phi}> Inclusive, Stat. Errors");
     gRES_1D_Stat->SetNameTitle(hName,hTitle);
     //
     gRES_1D_Syst                =   fSetSystErrors(gRES_1D_Stat);
     hName       =   Form("gRES_1D_Syst");
     hTitle      =   Form("<Y_{#phi}> Inclusive, Syst. Errors");
     gRES_1D_Syst->SetNameTitle(hName,hTitle);
     //
     gRES_2D_Stat                =   fEfficiencycorrection(hRAW_2D,hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection);
     for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )  {
         hName       =   Form("gRES_2D_Stat_%i",iHisto);
         hTitle      =   Form("<Y_{#phi#phi}> in p_{T} [%.2f#;%.2f], Stat. Errors",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
         gRES_2D_Stat.at(iHisto)->SetNameTitle(hName,hTitle);
     }
     //
     gRES_2D_Syst                =   fSetSystErrors(gRES_2D_Stat);
     for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )  {
         hName       =   Form("gRES_2D_Syst_%i",iHisto);
         hTitle      =   Form("<Y_{#phi#phi}> in p_{T} [%.2f#;%.2f], Syst. Errors",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
         gRES_2D_Syst.at(iHisto)->SetNameTitle(hName,hTitle);
         
     }
     //
     /*
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
     {
         auto    kN_INELgt0  =   fEvaluateINELgt0(iMult-1,hEvntMlt);
         f1DCorrection   =   (1./kBR)        *(1./kN_INELgt0);
         f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_INELgt0);
         //
         if ( iMult == 0 )   {
             gRES_1D_Stat_in_MT[0]   =   fEfficiencycorrection(hRAW_1D_in_MT[0],hREC_1D_in_MT_0,hGEN_1D_in_MT_0,f1DCorrection);
             hName       =   Form("gRES_1D_Stat_%i",0);
             hTitle      =   Form("<Y_{#phi}> in Mult. [%.2f#;%.2f], Stat. Errors",fArrMult[0],fArrMult[nBinMult]);
             gRES_1D_Stat_in_MT[0]   ->SetNameTitle(hName,hTitle);
             //
             gRES_1D_Syst_in_MT[0]   =   fSetSystErrors(gRES_1D_Stat_in_MT[0]);
             hName       =   Form("gRES_1D_Syst_%i",0);
             hTitle      =   Form("<Y_{#phi}> in Mult. [%.2f#;%.2f], Syst. Errors",fArrMult[0],fArrMult[nBinMult]);
             gRES_1D_Syst_in_MT[0] ->SetNameTitle(hName,hTitle);
         }   else    {
             gRES_1D_Stat_in_MT[iMult]=   fEfficiencycorrection(hRAW_1D_in_MT[iMult],hREC_1D_in_MT_0,hGEN_1D_in_MT_0,f1DCorrection);
             hName       =   Form("gRES_1D_Stat_%i",iMult);
             hTitle      =   Form("<Y_{#phi}> in Mult. [%.2f#;%.2f], Stat. Errors",fArrMult[iMult-1],fArrMult[iMult]);
             gRES_1D_Stat_in_MT[iMult] ->SetNameTitle(hName,hTitle);
             //
             gRES_1D_Syst_in_MT[iMult]=   fSetSystErrors(gRES_1D_Stat_in_MT[iMult]);
             hName       =   Form("gRES_1D_Syst_%i",iMult);
             hTitle      =   Form("<Y_{#phi}> in Mult. [%.2f#;%.2f], Syst. Errors",fArrMult[iMult-1],fArrMult[iMult]);
             gRES_1D_Syst_in_MT[iMult] ->SetNameTitle(hName,hTitle);
         }
         //
         gRES_2D_Stat_in_MT[iMult]   =   fEfficiencycorrection(hRAW_2D_in_MT[iMult],hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection);
         for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )  {
             hName       =   Form("gRES_2D_Stat_%i_%i",iHisto,iMult);
             hTitle      =   Form("<Y_{#phi#phi}> in Mult. [%.2f#;%.2f] in p_{T} [%.2f#;%.2f], Stat. Errors",fArrMult[iMult],fArrMult[iMult+1],fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
             gRES_2D_Stat_in_MT[iMult].at(iHisto)->SetNameTitle(hName,hTitle);
         }
         //
         gRES_2D_Syst_in_MT[iMult]   =   fSetSystErrors(gRES_2D_Stat_in_MT[iMult]);
         for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )  {
             hName       =   Form("gRES_2D_Syst_%i_%i",iHisto,iMult);
             hTitle      =   Form("<Y_{#phi#phi}> in Mult. [%.2f#;%.2f] in p_{T} [%.2f#;%.2f], Syst. Errors",fArrMult[iMult],fArrMult[iMult+1],fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
             gRES_2D_Syst_in_MT[iMult].at(iHisto)->SetNameTitle(hName,hTitle);
         }
         
     }
     //
     //-------------------------//
     //  Filling output objects //
     //-------------------------//
     //
     fStartTimer("Fit_for_extrapolation");
     //
     // Output File for Fit Check
     TFile*  outCheckFitYld  =   new TFile(fYldSigCh2k,"recreate");
     //
     // Total Fit number and progressive
     //
     Int_t   fTotalCount = nBinPT1D+(1+nBinPT2D)*nBinPT2D;
     Int_t   fProgrCount = 0;
     //
     //  Progressive Count
     //
     fProgrCount++;
     //
     //  Print loop Timer
     //
     fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
     //
     auto fResults = fMeasureFullYield(gRES_1D_Stat,gRES_1D_Syst,"1D");
     fYield_Stat ->  SetPoint        (0,1,fResults[0]);
     fYield_Syst ->  SetPoint        (0,1,fResults[0]);
     fYield_Stat ->  SetPointError   (0,0,0,fResults[1],fResults[2]);
     fYield_Syst ->  SetPointError   (0,0,0,fResults[3],fResults[4]);
     gMPT_1D_Stat->  SetPoint        (0,1,fResults[5]);
     gMPT_1D_Syst->  SetPoint        (0,1,fResults[5]);
     gMPT_1D_Stat->  SetPointError   (0,.5,.5,fResults[6],fResults[7]);
     gMPT_1D_Syst->  SetPointError   (0,.5,.5,fResults[8],fResults[9]);
     //
     for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
         //  Progressive Count
         fProgrCount++;
         //  Print loop Timer
         fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
         //
         fResults = fMeasureFullYield(gRES_2D_Stat.at(iFit),gRES_2D_Syst.at(iFit),Form("1D_2D_%i",iFit));
         //
         auto binwidth   =   (fArrPT2D[iFit+1]-fArrPT2D[iFit])*0.5;
         auto bincenter  =   (fArrPT2D[iFit+1]+fArrPT2D[iFit])*0.5;
         gYield_Profile_Stat ->  SetPoint        (iFit,bincenter,fResults[0]);
         gYield_Profile_Syst ->  SetPoint        (iFit,bincenter,fResults[0]);
         gYield_Profile_Stat ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[2]);
         gYield_Profile_Syst ->  SetPointError   (iFit,binwidth,binwidth,fResults[3],fResults[4]);
         gMPT_2D_Stat        ->  SetPoint        (iFit+1,bincenter,fResults[5]);
         gMPT_2D_Syst        ->  SetPoint        (iFit+1,bincenter,fResults[5]);
         gMPT_2D_Stat        ->  SetPointError   (iFit+1,binwidth,binwidth,fResults[6],fResults[7]);
         gMPT_2D_Syst        ->  SetPointError   (iFit+1,binwidth,binwidth,fResults[8],fResults[9]);
     }
     //
     fResults = fMeasureFullYield(gYield_Profile_Stat,gYield_Profile_Syst,"2D");
     fYield_Stat ->  SetPoint        (1,2,fResults[0]);
     fYield_Syst ->  SetPoint        (1,2,fResults[0]);
     fYield_Stat ->  SetPointError   (1,0,0,fResults[1],fResults[2]);
     fYield_Syst ->  SetPointError   (1,0,0,fResults[3],fResults[4]);
     gMPT_2D_Stat->  SetPoint        (0,0,fResults[5]);
     gMPT_2D_Syst->  SetPoint        (0,0,fResults[5]);
     gMPT_2D_Stat->  SetPointError   (0,.5,.5,fResults[6],fResults[7]);
     gMPT_2D_Syst->  SetPointError   (0,.5,.5,fResults[8],fResults[9]);
     //
     outCheckFitYld->Close();
     //
     /*
     // Output File for Fit Check
     TFile*  outCheckFitMlt  =   new TFile(fMltSigCh2k,"recreate");
     //
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ ) {
         if ( !kDoMultiplicity ) break;
         //  Progressive Count
         fProgrCount++;
         //  Print loop Timer
         fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
         //
         fResults = fMeasureFullYield(gRES_1D_Stat_in_MT[iMult],gRES_1D_Syst_in_MT[iMult],Form("1D_in_MT_%i",iMult));
         //
         fYield_Stat_in_MT[iMult] ->  SetPoint        (0,1,fResults[0]);
         fYield_Syst_in_MT[iMult] ->  SetPoint        (0,1,fResults[0]);
         fYield_Stat_in_MT[iMult] ->  SetPointError   (0,0,0,fResults[1],fResults[2]);
         fYield_Syst_in_MT[iMult] ->  SetPointError   (0,0,0,fResults[3],fResults[4]);
         //*/
         /*
         for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
             //  Progressive Count
             fProgrCount++;
             //  Print loop Timer
             fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
             //
             fResults = fMeasureFullYield(gRES_2D_Stat_in_MT[iMult].at(iFit),gRES_2D_Syst_in_MT[iMult].at(iFit),Form("2D_MT_%i_PT_%i",iMult,iFit));
             //
             auto binwidth   =   (fArrPT2D[iFit+1]-fArrPT2D[iFit])*0.5;
             auto bincenter  =   (fArrPT2D[iFit+1]+fArrPT2D[iFit])*0.5;
             fYield_Profile_Stat_in_MT[iMult] ->  SetPoint        (iFit,bincenter,fResults[0]);
             fYield_Profile_Syst_in_MT[iMult] ->  SetPoint        (iFit,bincenter,fResults[0]);
             fYield_Profile_Stat_in_MT[iMult] ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[1]);
             fYield_Profile_Syst_in_MT[iMult] ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[1]);
         }
         fResults = fMeasureFullYield(fYield_Profile_Stat_in_MT[iMult],fYield_Profile_Syst_in_MT[iMult],Form("2D_%i",iMult));
         fYield_Stat_in_MT[iMult]    ->  SetPoint        (1,2,fResults[0]);
         fYield_Syst_in_MT[iMult]    ->  SetPoint        (1,2,fResults[0]);
         fYield_Stat_in_MT[iMult]    ->  SetPointError   (1,0,0,fResults[1],fResults[1]);
         fYield_Syst_in_MT[iMult]    ->  SetPointError   (1,0,0,fResults[1],fResults[1]);
        /
     //}
     //
     fStopTimer("Fit_for_extrapolation");
     //
     //outCheckFitMlt->Close();
     //
     //--------------------------//
     //
     /*
     TH1F      **hPubMult =   new TH1F*[nMltTrgECls];
     //
     for ( Int_t iMult = 0; iMult <  nMltTrgECls; iMult++ )   {
         hName           =   Form("h%i_%i",(int)kMltTrgECls[iMult],(int)kMltTrgECls[iMult+1]);
         hPubMult[iMult] =   (TH1F*)(insPublishedRslt->Get(hName));
         hPubMult[iMult]->GetEntries();
     }
     //
     TH1F      **hPubAdpt  =  new TH1F*[nBinMult+1];
     //
     hName           =   Form("h0_100");
     hPubAdpt[0] =   new TH1F(*((TH1F*)(insPublishedRslt->Get(hName))));
     hPubAdpt[1] =   new TH1F(*hPubMult[0]);
     for ( Int_t iMult = 1; iMult < nBinMult; iMult++  )   {
         auto    h1 = *hPubMult[2*iMult-1];
         auto    h2 = *hPubMult[2*iMult];
         h1.Scale(kMltTrgECls[2*iMult+1]-kMltTrgECls[2*iMult]);
         h2.Scale(kMltTrgECls[2*iMult]-kMltTrgECls[2*iMult-1]);
         auto h3 = h1+h2;
         h3.Scale(1./(kMltTrgECls[2*iMult+1]-kMltTrgECls[2*iMult-1]));
         hPubAdpt[iMult+1] =   new TH1F(h3);
     }
     //
     // Check singular Integral and extrapolate
     // Make a full Yield in multiplicity
     //
     float fSResult [] = {0.0341,0.1052,0.078950000,0.057233333,0.038350000,0.019580000};
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )   {
         
         auto    fCheckTru   =   new TH1F(*hPubAdpt[iMult]);
         auto    hHisto      =   fMakeMeTH1F(gRES_1D_Stat_in_MT[iMult]);
         auto    hHist2      =   fMakeMeTH1F(gRES_1D_Stat_in_MT[iMult]);
         hHisto->Divide(hHist2,fCheckTru);
         
         TCanvas *fCanvasCheck =   new TCanvas();
         fCanvasCheck->SetLogy();
         fCheckTru->Fit(fLevyFit1D,"IMREQ0S","");
         fCheckTru->Draw("same");
         fLevyFit1D->Draw("same");
         hHisto->Draw("same");
         fCanvasCheck->SaveAs(Form("Check_%i.pdf",iMult));
         delete fCanvasCheck;
         
         TCanvas *fCanvasCheck2  =   new TCanvas();
         //fCanvasCheck2->SetLogy();
         hHisto->Fit("pol0");
         hHisto->Draw();
         fCanvasCheck2->SaveAs(Form("Check2_%i.pdf",iMult));
         delete fCanvasCheck2;
     }
     auto test = fPublishSpectrum(gMPT_2D_Stat,gMPT_2D_Syst);
     test->SaveAs("test.pdf");
     delete test;
     
     auto cDraw = new TCanvas();
     gStyle->SetOptStat(0);
     cDraw->SetLogy();
     TMultiGraph* mdraw = new TMultiGraph();
     mdraw->Add( gRES_1D_Stat, "EP*");
     mdraw->Add(gYield_Profile_Stat , "EP*");
     mdraw->Draw("ALP");
     cDraw->SaveAs("test2.pdf");
     delete cDraw;
     
     auto cDraw2 = new TCanvas();
     gStyle->SetOptStat(0);
     auto rrr = fMakeMeTH1F(gYield_Profile_Stat);
     gRES_1D_Stat->Fit(fLevyTsallis);
     auto fff = new TH1F(*rrr);
     fff->Divide(fLevyTsallis);
     fff->Draw();
     cDraw2->SaveAs("test3.pdf");
     delete cDraw2;
     
     cout << "------------------------------------------------" << endl;
     cout << "I                                              I" << endl;
     cout << "I       SUMMARY -- Inclusive Analysis          I" << endl;
     cout << "I                                              I" << endl;
     cout << "------------------------------------------------" << endl;
     cout << Form("-                  syst.+   %.5f -> %.1f",fYield_Syst->GetErrorYhigh(0),100*(fYield_Syst->GetErrorYhigh(0))/fYield_Syst->GetPointY(0))<<  endl;
     cout << Form("-                  stat.+   %.5f -> %.1f",fYield_Stat->GetErrorYhigh(0),100*(fYield_Stat->GetErrorYhigh(0))/fYield_Stat->GetPointY(0))<<  endl;
     cout << "- Y_{phi}          =   " << fYield_Stat->GetPointY(0) << endl;
     cout << Form("-                  stat.-   %.5f -> %.1f",fYield_Stat->GetErrorYlow(0),100*(fYield_Stat->GetErrorYlow(0))/fYield_Stat->GetPointY(0))<<  endl;
     cout << Form("-                  syst.-   %.5f -> %.1f",fYield_Syst->GetErrorYlow(0),100*(fYield_Syst->GetErrorYlow(0))/fYield_Syst->GetPointY(0))<<  endl;
     cout << "-                                               " << endl;
     cout << Form("-                  syst.+   %.5f -> %.1f",gMPT_1D_Syst->GetErrorYhigh(0),100*(gMPT_1D_Syst->GetErrorYhigh(0))/gMPT_1D_Syst->GetPointY(0))<<  endl;
     cout << Form("-                  stat.+   %.5f -> %.1f",gMPT_1D_Stat->GetErrorYhigh(0),100*(gMPT_1D_Stat->GetErrorYhigh(0))/gMPT_1D_Stat->GetPointY(0))<<  endl;
     cout << "- mean pT          =   " << gMPT_1D_Stat->GetPointY(0) << endl;
     cout << Form("-                  stat.-   %.5f -> %.1f",gMPT_1D_Stat->GetErrorYlow(0),100*(gMPT_1D_Stat->GetErrorYlow(0))/gMPT_1D_Stat->GetPointY(0))<<  endl;
     cout << Form("-                  syst.-   %.5f -> %.1f",gMPT_1D_Syst->GetErrorYlow(0),100*(gMPT_1D_Syst->GetErrorYlow(0))/gMPT_1D_Syst->GetPointY(0))<<  endl;
     cout << "-                                               " << endl;
     cout << Form("-                  syst.+   %.5f -> %.1f",fYield_Syst->GetErrorYhigh(1),100*(fYield_Syst->GetErrorYhigh(1))/fYield_Syst->GetPointY(1))<<  endl;
     cout << Form("-                  stat.+   %.5f -> %.1f",fYield_Stat->GetErrorYhigh(1),100*(fYield_Stat->GetErrorYhigh(1))/fYield_Stat->GetPointY(1))<<  endl;
     cout << "- Y_{phi,phi}      =   " << fYield_Stat->GetPointY(1) << endl;
     cout << Form("-                  stat.-   %.5f -> %.1f",fYield_Stat->GetErrorYlow(1),100*(fYield_Stat->GetErrorYlow(1))/fYield_Stat->GetPointY(1))<<  endl;
     cout << Form("-                  syst.-   %.5f -> %.1f",fYield_Syst->GetErrorYlow(1),100*(fYield_Syst->GetErrorYlow(1))/fYield_Syst->GetPointY(1))<<  endl;
     cout << "-                                               " << endl;
     auto fGamma         =   fGammaPhiValue(fYield_Stat->GetPointY(0),fYield_Stat->GetPointY(1));
     auto fGammaStatHig  =   fGammaPhiError(fYield_Stat->GetPointY(0),fYield_Stat->GetPointY(1),fYield_Stat->GetErrorYhigh(0),fYield_Stat->GetErrorYhigh(1));
     auto fGammaStatLow  =   fGammaPhiError(fYield_Stat->GetPointY(0),fYield_Stat->GetPointY(1),fYield_Stat->GetErrorYlow(0),fYield_Stat->GetErrorYlow(1));
     auto fGammaSystHig  =   fGammaPhiError(fYield_Stat->GetPointY(0),fYield_Stat->GetPointY(1),fYield_Syst->GetErrorYhigh(0),fYield_Syst->GetErrorYhigh(1));
     auto fGammaSystLow  =   fGammaPhiError(fYield_Stat->GetPointY(0),fYield_Stat->GetPointY(1),fYield_Syst->GetErrorYlow(0),fYield_Syst->GetErrorYlow(1));
     cout << Form("-                  syst.+   %.5f -> %.1f",fGammaSystHig,100*(fGammaSystHig/(fGamma))) <<  endl;
     cout << Form("-                  stat.+   %.5f -> %.1f",fGammaStatHig,100*(fGammaStatHig/(fGamma))) <<  endl;
     cout << "- gamma_{phi}      =   " << fGamma << endl;
     cout << Form("-                  stat.-   %.5f -> %.1f",fGammaStatLow,100*(fGammaStatLow/(fGamma))) <<  endl;
     cout << Form("-                  syst.-   %.5f -> %.1f",fGammaSystLow,100*(fGammaSystLow/(fGamma))) <<  endl;
     cout << "-                                               " << endl;
     cout << "------------------------------------------------" << endl;
     cout << "I       COMPARISON -- 7TeV pp                  I" << endl;
     cout << "------------------------------------------------" << endl;
     cout << "-                                               " << endl;
     cout << Form("- Y_{phi}           %.5f / %.5f -> %.5f",fYield_Stat->GetPointY(0),0.0318,fYield_Stat->GetPointY(0)/0.0318) << endl;
     cout << Form("- mean pT           %.5f / %.5f -> %.5f",gMPT_1D_Stat->GetPointY(0),1.132,gMPT_1D_Stat->GetPointY(0)/1.132) << endl;
     cout << "-                                               " << endl;
     cout << "------------------------------------------------" << endl;
     cout << "I       COMPARISON -- 5TeV pp                  I" << endl;
     cout << "------------------------------------------------" << endl;
     cout << "-                                               " << endl;
     cout << Form("- Y_{phi}           %.5f / %.5f -> %.5f",fYield_Stat->GetPointY(0),0.0301,fYield_Stat->GetPointY(0)/0.0301) << endl;
     cout << Form("- mean pT           %.5f / %.5f -> %.5f",fResults[5],1.132,fResults[5]/1.132) << endl;
     cout << "-                                               " << endl;
     cout << "------------------------------------------------" << endl;
     cout << "I                                              I" << endl;
     cout << "I       SUMMARY -- Multiplicity Analysis       I" << endl;
     cout << "I                                              I" << endl;
     cout << "------------------------------------------------" << endl;
     /*
     for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )   {
         cout << "Mult" << iMult << " phi:    " <<    fYield_Stat_in_MT[iMult]->GetPointY(0) << "   fResult: " << fSResult[iMult] << " Ratio: " << fYield_Stat_in_MT[iMult]->GetPointY(0)/fSResult[iMult] << endl;
         cout << "Mult" << iMult << " phiphi: " <<    fYield_Stat_in_MT[iMult]->GetPointY(1) << endl;
     }*/
     
     /*
     TFile *fCheck   =   new TFile("./result/HEPData-ins1762364-v1-Table_4.root");
     hName       =   "Table 4/Graph1D_y1";
     TGraphAsymmErrors  *gCheck   =   (TGraphAsymmErrors*)(fCheck->Get(hName));
     hName       =   "Table 4/Hist1D_y1";
     TH1F  *hCheck   =   (TH1F*)(fCheck->Get(hName));
     
     
     TFile * TheCheck = new TFile("eee.root","recreate");
     TCanvas *fCanvasCheck3  =   new TCanvas();
     gStyle->SetOptStat(0);
     fCanvasCheck3->SetLogy();
     
     auto gRES_1D = fSumGraphErrors(gRES_1D_Stat,gRES_1D_Syst);
     
     // X-Axis formatting
     gRES_1D->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
     gRES_1D->GetXaxis()->SetTitleOffset(1.15);
     
     // Y-Axis formatting
     gRES_1D->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
     gRES_1D->GetYaxis()->SetTitleOffset(1.15);
     
     gRES_1D->SetMarkerStyle(26);
     gRES_1D->SetMarkerColor(2);
     gRES_1D->SetLineColor(2);
     gRES_1D->SetFillColor(2);
     gRES_1D->SetMarkerSize(1);
     
     gCheck->SetMarkerStyle(27);
     gCheck->SetMarkerColor(4);
     gCheck->SetLineColor(4);
     gCheck->SetFillColor(4);
     gCheck->SetMarkerSize(1);
     
     TMultiGraph * m1 = new TMultiGraph();
     m1->Add(gRES_1D);
     m1->Add(gCheck);
     m1->Draw("APE");
     
     TLegend *fLegend = new TLegend(0.6,0.7,0.9,0.9);
     fLegend->AddEntry(gRES_1D,"My Results","PE");
     fLegend->AddEntry(gCheck,"Pub. Results","PE");
     
     fLegend->Draw("same");
     fCanvasCheck3->SaveAs(Form("Check2_%i.pdf",-1));
     fCanvasCheck3->Write();
     delete fCanvasCheck3;
     
     //--------------------------
     
     TCanvas *fCanvasCheck4  =   new TCanvas();
     gStyle->SetOptStat(0);
     fCanvasCheck4->SetLogx();
     
     TH1F* hEFF = new TH1F(*hREC_1D);
     TH1F* hRES = new TH1F(*hREC_1D);
     TH1F* hFRC = new TH1F(*hREC_1D);
     
     hEFF->Divide(hREC_1D,hGEN_1D,1.,1.,"b");
     hRES->Divide(hRAW_1D,hEFF,f1DCorrection);
     
     auto hDraw = fCheckPublishedData(hRES,hCheck,gCheck);
     auto gErrors = fSumGraphErrors(gRES_1D,gCheck);
     
     for (Int_t i = 0; i < gCheck->GetN(); i++)  {
         auto fError = gErrors->GetErrorYhigh(i);
         auto fValY = gErrors->GetPointY(i);
         hDraw->SetBinError(i+1,fError/fValY);
     }
     hDraw->SetMaximum(1.5);
     hDraw->SetMinimum(0.5);
     hDraw->SetMarkerStyle(27);
     hDraw->Draw();
     fCheckPublishedData(hRES,hCheck,gCheck)->Write();
     
     fCanvasCheck4->SaveAs(Form("Check2_%i.pdf",-2));
     delete fCanvasCheck4;
     
     TheCheck->Close();
      
     //
     //--------------------------//
     //
     //--------------------------//
     //  Printing output objects //
     //--------------------------//
     //
     // >> Trigger Analysis
     //
     //TFile *outFil1  =   new TFile   (fTrgSigCorr,"recreate");
     //
     //outFil1->Close();
     //
     // >> Yield Analysis
     //
     TFile *outFil2  =   new TFile   (fYldSigCorr,"recreate");
     //
     hEvntEff                ->Write();
     //
     fYield_Stat             ->Write();
     fYield_Syst             ->Write();
     gMPT_1D_Stat            ->Write();
     gMPT_1D_Syst            ->Write();
     gMPT_2D_Stat            ->Write();
     gMPT_2D_Syst            ->Write();
     gRES_1D_Stat->Write();
     gRES_1D_Syst->Write();
     for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
         gRES_2D_Stat.at(iFit)->Write();
         gRES_2D_Syst.at(iFit)->Write();
     }
     gYield_Profile_Stat     ->Write();
     gYield_Profile_Syst     ->Write();
     //
     outFil2->Close();
     //
     // >> Multiplicity Analysis
     //
     /*
     if  ( kDoMultiplicity ) {
         TFile *outFil3  =   new TFile   (fMltSigCorr,"recreate");
         //
         hEvntEff->Write();
         for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
         {
             if ( !kDoMultiplicity ) break;
             fYield_Stat_in_MT[iMult]    ->Write();
             fYield_Syst_in_MT[iMult]    ->Write();
             gRES_1D_Stat_in_MT[iMult]   ->Write();
             gRES_1D_Syst_in_MT[iMult]   ->Write();
             /*
             for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
                 gRES_2D_Stat_in_MT[iMult].at(iFit)->Write();
                 gRES_2D_Syst_in_MT[iMult].at(iFit)->Write();
             }*//*
         }
         //
         outFil3->Close();
     }
     //
     // >-> Close input File
     //
     insFile_DT_Yield    ->Close();
     insFile_DT_Mult     ->Close();
     insFile_EF_Yield    ->Close();
     insFile_EF_Mult     ->Close();
 }


 */
