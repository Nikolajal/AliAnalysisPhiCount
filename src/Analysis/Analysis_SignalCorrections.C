// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void Analysis_SignalCorrections ( bool fSilent = true, TString fOption = "" )
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
    TH1F       *hEFF_1D;
    TH1F       *hGEN_1D;
    TH1F       *hREC_1D_in_2D_bin;
    TH1F       *hGEN_1D_in_2D_bin;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hRAW_1D";
    hRAW_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "hEFF_1D";
    hEFF_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hREC_1D";
    hREC_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_1D";
    hGEN_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
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
    TH2F       *hEFF_2D;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hRAW_2D";
    hRAW_2D     =   (TH2F*)(insFile_DT_Yield->Get(hName));
    hRAW_2D->GetEntries();
    //

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
    Double_t    f1DCorrection   =   (1./kBR)        *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg)*(1./0.99);
    Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    //
    gRES_1D_Stat                =   fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection);
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
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        if ( !kDoMultiplicity ) break;
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
        /*
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
         */
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
        //
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
         */
    }
    //
    fStopTimer("Fit_for_extrapolation");
    //
    outCheckFitMlt->Close();
    //
    //--------------------------//
    //
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
    cout << Form("- mean pT           %.5f / %.5f -> %.5f",fResults[5],1.132,fResults[5]/1.132) << endl;
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
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )   {
        cout << "Mult" << iMult << " phi:    " <<    fYield_Stat_in_MT[iMult]->GetPointY(0) << "   fResult: " << fSResult[iMult] << " Ratio: " << fYield_Stat_in_MT[iMult]->GetPointY(0)/fSResult[iMult] << endl;
        cout << "Mult" << iMult << " phiphi: " <<    fYield_Stat_in_MT[iMult]->GetPointY(1) << endl;
    }
    
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
     */
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
            }*/
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

