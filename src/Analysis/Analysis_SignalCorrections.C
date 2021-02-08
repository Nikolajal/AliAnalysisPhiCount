// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void Analysis_SignalCorrections ( bool fSilent = true )
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
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (fYldSigExtr);
    TFile*  insFile_DT_Mult             =   new TFile   (fMltSigExtr);
    TFile*  insFile_EF_Yield            =   new TFile   (fYldPrePrMC);
    TFile*  insFile_EF_Mult             =   new TFile   (fMltPrePrMC);
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> GENERAL ANALYSIS //
    //
    TH1D       *hEvntEff;
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    //
    // >-> YIELD ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hRAW_1D;
    TH1F       *hEFF_1D;
    TH1F       *hREC_1D;
    TH1F       *hGEN_1D;
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
    //
    hName       =   "hEFF_2D";
    hEFF_2D     =   (TH2F*)(insFile_EF_Yield->Get(hName));
    //

    // >-> MULTIPLICITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hRAW_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F      **hEFF_1D_in_MT               = new TH1F     *[nBinMult];
    //
    //  Defining MT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        //
        hName = Form("hRAW_1D_in_MT_%i",iHisto);
        hRAW_1D_in_MT[iHisto]   =   (TH1F*)(insFile_DT_Mult->Get(hName));
        //
        hName = Form("hEFF_1D_in_MT_%i",iHisto);
        hEFF_1D_in_MT[iHisto]   =   (TH1F*)(insFile_EF_Mult->Get(hName));
        //
    }

    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hRAW_2D_in_MT               = new TH2F     *[nBinMult];
    TH2F      **hEFF_2D_in_MT               = new TH2F     *[nBinMult];
    //
    //  Defining MT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        //
        hName   =   Form("hRAW_2D_in_MT_%i",iHisto);
        hRAW_2D_in_MT[iHisto]   =   (TH2F*)(insFile_DT_Mult->Get(hName));
        //
        hName   =   Form("hEFF_2D_in_MT_%i",iHisto);
        hEFF_2D_in_MT[iHisto]   =   (TH2F*)(insFile_EF_Mult->Get(hName));
        //
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
    TGraphAsymmErrors   *gRES_1D_Syst       =   new TGraphAsymmErrors();
    TH1F       *hRES_1D;
    //
    hName       =   Form("hRES_1D");
    hTitle      =   Form("hRES_1D");
    hRES_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRES_1D,"PT 1D");
    //
    hName       =   Form("gRES_1D_Stat");
    hTitle      =   Form("gRES_1D_Stat");
    gRES_1D_Stat->SetNameTitle(hName,hTitle);
    //
    hName       =   Form("gRES_1D_Syst");
    hTitle      =   Form("gRES_1D_Syst");
    gRES_1D_Syst->SetNameTitle(hName,hTitle);
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hRES_2D;
    //
    hName       =   Form("hRES_2D");
    hTitle      =   Form("hRES_2D");
    hRES_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hRES_2D,"PT 2D");
    //
    TGraphAsymmErrors   *fYield_Profile_Stat    =   new TGraphAsymmErrors();
    hName       =   Form("fYield_Stat");
    hTitle      =   Form("fYield_Stat");
    fYield_Stat->SetNameTitle(hName,hTitle);
    //
    TGraphAsymmErrors   *fYield_Profile_Syst    =   new TGraphAsymmErrors();
    hName       =   Form("fYield_Syst");
    hTitle      =   Form("fYield_Syst");
    fYield_Syst->SetNameTitle(hName,hTitle);
    //
    
    // >-> MULTIPLICITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    TGraphAsymmErrors  **fYield_Stat_in_MT  =   new TGraphAsymmErrors  *[nBinMult];
    TGraphAsymmErrors  **fYield_Syst_in_MT  =   new TGraphAsymmErrors  *[nBinMult];
    //
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )   {
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
    TGraphAsymmErrors  **gRES_1D_Stat_in_MT =   new TGraphAsymmErrors  *[nBinMult];
    TGraphAsymmErrors  **gRES_1D_Syst_in_MT =   new TGraphAsymmErrors  *[nBinMult];
    TH1F      **hRES_1D_in_MT               = new TH1F     *[nBinMult];
    //
    //  Defining MT-Differential histograms over measurable pT
    //
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )
    {
        hName = Form("hRES_1D_in_MT_%i",iMult);
        hTitle = Form("hRES_1D_in_MT_%i",iMult);
        hRES_1D_in_MT[iMult]   =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hRES_1D_in_MT[iMult],"PT 1D");
        //
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
    TH2F      **hRES_2D_in_MT               = new TH2F     *[nBinMult];
    TGraphAsymmErrors  **fYield_Profile_Stat_in_MT  =   new TGraphAsymmErrors  *[nBinMult];
    TGraphAsymmErrors  **fYield_Profile_Syst_in_MT  =   new TGraphAsymmErrors  *[nBinMult];
    //
    //  Defining MT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName       =   Form("hRES_2D_in_MT_%i",iHisto);
        hTitle      =   Form("hRES_2D_in_MT_%i",iHisto);
        hRES_2D_in_MT[iHisto]   =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hRES_2D_in_MT[iHisto],"PT 2D");
        //
        fYield_Profile_Stat_in_MT[iHisto]    =   new TGraphAsymmErrors();
        hName       =   Form("fYield_Stat_in_MT_%i",iHisto);
        hTitle      =   Form("fYield_Stat_in_MT_%i",iHisto);
        fYield_Profile_Stat_in_MT[iHisto]->SetNameTitle(hName,hTitle);
        //
        fYield_Profile_Syst_in_MT[iHisto]    =   new TGraphAsymmErrors();
        hName       =   Form("fYield_Syst_in_MT_%i",iHisto);
        hTitle      =   Form("fYield_Syst_in_MT_%i",iHisto);
        fYield_Profile_Syst_in_MT[iHisto]->SetNameTitle(hName,hTitle);
    }
    //
    //
    //
    //---------------------//
    // Preprocessing input //
    //---------------------//
    //
    //         N_raw    1     1    Eff    1     1
    // N_res = ----- X --- X --- X --- X --- X ---
    //          eff    DpT   Dy    Evn   BR    Vrt
    //
    // Error propagation
    hRES_1D->Sumw2();
    hRES_2D->Sumw2();
    //
    // Scaling in pT [Done in PreProcessing]
    //
    // Scaling for efficiencies
    //
    gRES_1D_Stat    =   fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,1./((hEvntEff->GetBinContent(9))*kRapidityIntvl*kBranchingRtio*kVertexEfficnc));
    //
    hRES_1D->Divide(hRAW_1D,hEFF_1D,1.,1.,"");
    hRES_2D->Divide(hRAW_2D,hEFF_2D,1.,1.,"");
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRES_1D_in_MT[iHisto]->Divide(hRAW_1D_in_MT[iHisto],hEFF_1D_in_MT[iHisto],1.,1.,"");
        hRES_2D_in_MT[iHisto]->Divide(hRAW_2D_in_MT[iHisto],hEFF_2D_in_MT[iHisto],1.,1.,"");
    }
    //
    // Scaling in Rapidity Interval
    hRES_1D->Scale(1./kRapidityIntvl);
    hRES_2D->Scale(1./kRapidityIntvl);
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRES_1D_in_MT[iHisto]->Scale(1./kRapidityIntvl);
        hRES_2D_in_MT[iHisto]->Scale(1./kRapidityIntvl);
    }
    //
    // Scaling for events
    hRES_1D->Scale(kTriggerEfficnc/(hEvntEff->GetBinContent(1)));
    hRES_2D->Scale(kTriggerEfficnc/(hEvntEff->GetBinContent(1)));
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRES_1D_in_MT[iHisto]->Scale(kTriggerEfficnc/(hEvntEff->GetBinContent(1)));
        hRES_2D_in_MT[iHisto]->Scale(kTriggerEfficnc/(hEvntEff->GetBinContent(1)));
    }
    //
    // Scaling for Branching Ratio
    hRES_1D->Scale(1./kBranchingRtio);
    hRES_2D->Scale(1./(kBranchingRtio*kBranchingRtio));
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRES_1D_in_MT[iHisto]->Scale(1./kBranchingRtio);
        hRES_2D_in_MT[iHisto]->Scale(1./(kBranchingRtio*kBranchingRtio));
    }
    //
    // Scaling for Vertex Selection Efficiency
    //
    hRES_1D->Scale(1./kVertexEfficnc);
    hRES_2D->Scale(1./kVertexEfficnc);
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRES_1D_in_MT[iHisto]->Scale(1./kVertexEfficnc);
        hRES_2D_in_MT[iHisto]->Scale(1./kVertexEfficnc);
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
    TGraphAsymmErrors   *fCheck = new TGraphAsymmErrors(hRES_1D);
    auto fResults = fMeasureFullYield(fCheck,fCheck,"1D");
    fYield_Stat ->  SetPoint        (0,1,fResults[0]);
    fYield_Syst ->  SetPoint        (0,1,fResults[0]);
    fYield_Stat ->  SetPointError   (0,0,0,fResults[1],fResults[1]);
    fYield_Syst ->  SetPointError   (0,0,0,fResults[1],fResults[1]);
    //
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        //  Progressive Count
        fProgrCount++;
        //  Print loop Timer
        fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
        //
        fResults = fMeasureFullYield(hRES_2D->ProjectionX(Form("dd_%i",iFit),iFit+1,iFit+1),Form("2D_PT_%i",iFit));
        //
        auto binwidth   =   (fArrPT2D[iFit+1]-fArrPT2D[iFit])*0.5;
        auto bincenter  =   (fArrPT2D[iFit+1]+fArrPT2D[iFit])*0.5;
        fYield_Profile_Stat ->  SetPoint        (iFit,bincenter,fResults[0]);
        fYield_Profile_Syst ->  SetPoint        (iFit,bincenter,fResults[0]);
        fYield_Profile_Stat ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[1]);
        fYield_Profile_Syst ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[1]);
    }
    //
    /*
    fResults = fMeasureFullYield(fYield_Profile_Stat,fYield_Profile_Syst,"2D");
    fYield_Stat ->  SetPoint        (1,2,fResults[0]);
    fYield_Syst ->  SetPoint        (1,2,fResults[0]);
    fYield_Stat ->  SetPointError   (1,0,0,fResults[1],fResults[1]);
    fYield_Syst ->  SetPointError   (1,0,0,fResults[1],fResults[1]);
    //
    outCheckFitYld->Close();
    //
    // Output File for Fit Check
    TFile*  outCheckFitMlt  =   new TFile(fMltSigCh2k,"recreate");
    //
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ ) {
        //  Progressive Count
        fProgrCount++;
        //  Print loop Timer
        fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
        //
        fResults = fMeasureFullYield(hRES_1D_in_MT[iMult],Form("1D_in_MT_%i",iMult));
        //
        fYield_Stat_in_MT[iMult] ->  SetPoint        (0,1,fResults[0]);
        fYield_Syst_in_MT[iMult] ->  SetPoint        (0,1,fResults[0]);
        fYield_Stat_in_MT[iMult] ->  SetPointError   (0,0,0,fResults[1],fResults[1]);
        fYield_Syst_in_MT[iMult] ->  SetPointError   (0,0,0,fResults[1],fResults[1]);
        //
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
            //  Progressive Count
            fProgrCount++;
            //  Print loop Timer
            fPrintLoopTimer("Fit_for_extrapolation",fProgrCount,fTotalCount,1);
            //
            fResults = fMeasureFullYield(hRES_2D_in_MT[iMult]->ProjectionX(Form("dd_%i",iFit),iFit+1,iFit+1),Form("2D_MT_%i_PT_%i",iMult,iFit));
            //
            auto binwidth   =   (fArrPT2D[iFit+1]-fArrPT2D[iFit])*0.5;
            auto bincenter  =   (fArrPT2D[iFit+1]+fArrPT2D[iFit])*0.5;
            fYield_Profile_Stat_in_MT[iMult] ->  SetPoint        (iFit,bincenter,fResults[0]);
            fYield_Profile_Syst_in_MT[iMult] ->  SetPoint        (iFit,bincenter,fResults[0]);
            fYield_Profile_Stat_in_MT[iMult] ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[1]);
            fYield_Profile_Syst_in_MT[iMult] ->  SetPointError   (iFit,binwidth,binwidth,fResults[1],fResults[1]);
        }
        fResults = fMeasureFullYield(fYield_Profile_Stat_in_MT[iMult],fYield_Profile_Syst_in_MT[iMult],"2D");
        fYield_Stat_in_MT[iMult]    ->  SetPoint        (1,2,fResults[0]);
        fYield_Syst_in_MT[iMult]    ->  SetPoint        (1,2,fResults[0]);
        fYield_Stat_in_MT[iMult]    ->  SetPointError   (1,0,0,fResults[1],fResults[1]);
        fYield_Syst_in_MT[iMult]    ->  SetPointError   (1,0,0,fResults[1],fResults[1]);
    }
    //
    fStopTimer("Fit_for_extrapolation");
    //
    outCheckFitMlt->Close();
     */
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
    hRES_1D                 ->Write();
    hRES_2D                 ->Write();
    gRES_1D_Stat            ->Write();
    gRES_1D_Syst            ->Write();
    fYield_Stat             ->Write();
    fYield_Syst             ->Write();
    fYield_Profile_Stat     ->Write();
    fYield_Profile_Syst     ->Write();
    //
    outFil2->Close();
    //
    // >> Multiplicity Analysis
    //
    TFile *outFil3  =   new TFile   (fMltSigCorr,"recreate");
    //
    hEvntEff->Write();
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )
    {
        hRES_1D_in_MT[iMult]        ->Write();
        hRES_2D_in_MT[iMult]        ->Write();
        fYield_Stat_in_MT[iMult]    ->Write();
        fYield_Syst_in_MT[iMult]    ->Write();
    }
    //
    outFil3->Close();
    //
    // >-> Close input File
    //
    insFile_DT_Yield    ->Close();
    insFile_DT_Mult     ->Close();
    insFile_EF_Yield    ->Close();
    insFile_EF_Mult     ->Close();
}
