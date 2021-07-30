// File for 1-Dimensional Analysis:
// !TODO: 1. Move the results folder in the signal extraction folder 
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void Production_SignalExtraction ( bool fSilent = true, TString fOption = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //>-> OPTIONS
    
    // Silencing warnings for smoother running
    if ( fSilent )
    {
        gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    fChooseOption(fOption);
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (Form(kAnalysis_InvMassHist,"yield"));
    TFile*  insFile_DT_Slp              =   new TFile   (Form(kMassResolution_Anal,"yield"));
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hREC_1D;
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    TH1F       *hSlop_Reference;
    TH1F       *hSlop_Referen2D;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hREC_1D";
    hREC_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
    hName                       =   "hSigmaCnt_1D";
    hSlop_Reference             =   (TH1F*)(insFile_DT_Slp->Get(hName));
    hName                       =   "hSigmaCnt_1D_in_2D_bin";
    hSlop_Referen2D             =   (TH1F*)(insFile_DT_Slp->Get(hName));
    //
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        hName                   =   Form("hREC_1D_in_PT_%i",iHisto);
        hREC_1D_in_PT[iHisto]   =   (TH1F*)(insFile_DT_Yield->Get(hName));
    }

    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hREC_2D;
    TH1F      **hREC_1D_in_PT_2D_bin        = new TH1F     *[nBinPT2D];
    TH2F     ***hREC_2D_in_PT               = new TH2F    **[nBinPT2D];
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hREC_2D";
    hREC_2D     =   (TH2F*)(insFile_DT_Yield->Get(hName));
    //
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        hName = Form("hREC_1D_in_PT_2D_bin_%i",iHisto);
        hREC_1D_in_PT_2D_bin[iHisto]    =   (TH1F*)(insFile_DT_Yield->Get(hName));
        hREC_2D_in_PT[iHisto]           =   new TH2F       *[nBinPT2D];
        
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i",iHisto,jHisto);
            hREC_2D_in_PT[iHisto][jHisto]    = (TH2F*)(insFile_DT_Yield->Get(hName));
        }
    }
    //
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Creating the histograms-------------------------------------------------------------------------------

    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    fSetBinRap_();
    fSetBinMult();
    fSetBinNTup();
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hRAW_1D;
    TH1F       *hRAW_1D_in_2D_bin;
    TH1F       *hRAW_Mss_1D;
    TH1F       *hRAW_Mss_1D_in_2D_bin;
    TH1F       *hRAW_Wdt_1D;
    TH1F       *hRAW_Wdt_1D_in_2D_bin;
    TH1F       *hRAW_Slp_1D;
    TH1F       *hRAW_Slp_1D_in_2D_bin;
    //
    //  Defining yield histogram over measurable pT
    //
    hName       =   Form("hRAW_1D");
    hTitle      =   Form("hRAW_1D");
    hRAW_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_1D,"PT 1D");
    //
    hName       =   Form("hRAW_1D_in_2D_bin");
    hTitle      =   Form("hRAW_1D_in_2D_bin");
    hRAW_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hRAW_Mss_1D");
    hTitle      =   Form("hRAW_Mss_1D");
    hRAW_Mss_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Mss_1D,"PT 1D");
    //
    hName       =   Form("hRAW_Mss_1D_in_2D_bin");
    hTitle      =   Form("hRAW_Mss_1D_in_2D_bin");
    hRAW_Mss_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hRAW_Mss_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hRAW_Wdt_1D");
    hTitle      =   Form("hRAW_Wdt_1D");
    hRAW_Wdt_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Wdt_1D,"PT 1D");
    //
    hName       =   Form("hRAW_Wdt_1D_in_2D_bin");
    hTitle      =   Form("hRAW_Wdt_1D_in_2D_bin");
    hRAW_Wdt_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hRAW_Wdt_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hRAW_Slp_1D");
    hTitle      =   Form("hRAW_Slp_1D");
    hRAW_Slp_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Slp_1D,"PT 1D");
    //
    hName       =   Form("hRAW_Slp_1D_in_2D_bin");
    hTitle      =   Form("hRAW_Slp_1D_in_2D_bin");
    hRAW_Slp_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hRAW_Slp_1D_in_2D_bin,"PT 1D");
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hRAW_2D;
    //
    //  Defining yield histogram over measurable pT
    //
    hName       =   Form("hRAW_2D");
    hTitle      =   Form("hRAW_2D");
    hRAW_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hRAW_2D,"PT 2D");
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    //  Making utility variables
    Int_t fTotalCount, fProgrCount;
    //
    // Total Fit number abnd progressive
    fTotalCount = nOptions;
    fProgrCount = 0.;
    //
    // Starting cycle
    //
    fStartTimer("Systematics Yield Extraction Evaluation 1D");
    //
    for ( Int_t iSys = 0; iSys < nOptions; iSys++ ) {
        //
        gROOT                       ->  ProcessLine (Form(".! mkdir -p ./result/yield/ExtractionSystematics/ExtractionCheck/%s/1D/",sOptions.at(iSys).Data()));
        TFile      *fCheckFit       =   new TFile   (Form("./result/yield/ExtractionSystematics/ExtractionCheck/%s/1D/CheckFitResults_%s.root",sOptions.at(iSys).Data(),sOptions.at(iSys).Data()),"recreate");
        //
        //>>    Fit the Model
        auto    fFitResults_1DYield = FitModel        (hREC_1D_in_PT,hSlop_Reference,Form("./result/yield/ExtractionSystematics/ExtractionCheck/%s/1D/",sOptions.at(iSys).Data()),sOptions.at(iSys).Data(),sOptions.at(iSys).Data());
        //
        //>>    Progressive Count
        fProgrCount++;
        //
        TFile      *fResultFit       =   new TFile   (Form("./result/yield/ExtractionSystematics/ExtractionCheck/%s/1D/FitResults_%s.root",sOptions.at(iSys).Data(),sOptions.at(iSys).Data()),"recreate");
        for ( auto hSave : fFitResults_1DYield )    {
            if ( strncmp(hSave->GetName(),"hRAW_1D",7) == 0 )   hSave->Scale(1.,"width");
            if ( strncmp(hSave->GetName(),"h1D_anBB",7) == 0 )  hSave->Scale(1.,"width");
            hSave   ->  Write();
        }
        //
        fResultFit->Close();
        fCheckFit->Close();
        //
        //>>    Print Progress
        fPrintLoopTimer("Systematics Yield Extraction Evaluation 1D",fProgrCount,fTotalCount,1);
    }
    //
    fStopTimer("Systematics Yield Extraction Evaluation 1D");
    //
    // Total Fit number abnd progressive
    fTotalCount = nOption2;
    fProgrCount = 0.;
    //
    fStartTimer("Systematics Yield Extraction Evaluation 2D");
    //
    for ( Int_t iSys = 0; iSys < nOption2; iSys++ ) {
        //
        gROOT                       ->  ProcessLine (Form(".! mkdir -p ./result/yield/ExtractionSystematics/ExtractionCheck/%s/2D/",sOption2.at(iSys).Data()));
        TFile      *fCheckFit       =   new TFile   (Form("./result/yield/ExtractionSystematics/ExtractionCheck/%s/2D/CheckFitResults_%s.root",sOption2.at(iSys).Data(),sOption2.at(iSys).Data()),"recreate");
        //
        //>>    Fit
        std::vector<TH1F*>  f1DCheck;
        auto    fFitResults_2DYield = FitModel(hREC_1D_in_PT_2D_bin,hSlop_Referen2D,hREC_2D_in_PT,f1DCheck,Form("./result/yield/ExtractionSystematics/ExtractionCheck/%s/2D/",sOption2.at(iSys).Data()),sOption2.at(iSys).Data(),sOption2.at(iSys).Data());
        //
        //Progressive Count
        fProgrCount++;
        //
        //
        TFile      *fResultFit       =   new TFile   (Form("./result/yield/ExtractionSystematics/ExtractionCheck/%s/2D/FitResults_%s.root",sOption2.at(iSys).Data(),sOption2.at(iSys).Data()),"recreate");
        for ( auto hSave : f1DCheck )    {
            if ( strncmp(hSave->GetName(),"hRAW_1D_in_2D_bin",17) == 0 )    hSave->Scale(1.,"width");
            if ( strncmp(hSave->GetName(),"2Dbin_anBB",17) == 0 )           hSave->Scale(1.,"width");
            hSave   ->  Write();
        }
        for ( auto hSave : fFitResults_2DYield )    {
            if ( strncmp(hSave->GetName(),"hRAW_2D",7) == 0 )   hSave->Scale(1.,"width");
            if ( strncmp(hSave->GetName(),"anSB2D",6) == 0 )    hSave->Scale(1.,"width");
            if ( strncmp(hSave->GetName(),"anBS2D",6) == 0 )    hSave->Scale(1.,"width");
            if ( strncmp(hSave->GetName(),"anBB2D",6) == 0 )    hSave->Scale(1.,"width");
            hSave   ->  Write();
        }
        //
        fResultFit->Close();
        fCheckFit->Close();
        //
        //>>    Print Progress
        fPrintLoopTimer("Systematics Yield Extraction Evaluation 2D",fProgrCount,fTotalCount,1);
    }
    fStopTimer("Systematics Yield Extraction Evaluation 2D");
    insFile_DT_Yield->Close();
}
