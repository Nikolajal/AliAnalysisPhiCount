#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void SignalExtraction ( TString fOption = "", bool fSilent = true )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //>-> OPTIONS
    
    bool kDosss = true;
    
    // Silencing warnings for smoother running
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    fChooseOption(fOption);
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (Form(kAnalysis_InvMassHist,"Yield"));
    TFile*  insFile_DT_Mult             =   new TFile   (Form(kAnalysis_InvMassHist,"Multiplicity"));
    TFile*  insFile_DT_Rap              =   new TFile   (Form(kAnalysis_InvMassHist,"Rapidity"));
    TFile*  insFile_DT_Slp              =   new TFile   (Form(kMassResolution_Anal,"Yield"));
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1D       *hEvntEff;
    TH1D       *hEvntMlt;
    TH1F       *hREC_1D;
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    TH1F       *hSlop_Reference;
    TH1F       *hSlop_Referen2D;
    //
    //  Utility
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "fQC_Event_Enum_Mult";
    hEvntMlt    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hREC_1D";
    hREC_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
    //
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        hName                   =   Form("hREC_1D_in_PT_%i",iHisto);
        hREC_1D_in_PT[iHisto]   =   (TH1F*)(insFile_DT_Yield->Get(hName));
    }
    //
    hName                       =   "hSigmaCnt_1D";
    hSlop_Reference             =   (TH1F*)(insFile_DT_Slp->Get(hName));
    hName                       =   "hSigmaCnt_1D_in_2D_bin";
    hSlop_Referen2D             =   (TH1F*)(insFile_DT_Slp->Get(hName));

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
    
    // >-> MULTIPLICITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_MT               = new TH1F     *[nBinMult+1];
    TH1F     ***hREC_1D_in_MT_in_PT         = new TH1F    **[nBinMult+1];
    //
    //  Defining MT-Differential histograms over measurable pT
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iMult);
        hREC_1D_in_MT[iMult]   =   (TH1F*)(insFile_DT_Mult->Get(hName));
        
        hREC_1D_in_MT_in_PT[iMult] = new TH1F     *[nBinPT1D];

        for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
        {
            hName = Form("hREC_1D_in_MT_PT_%i_%i",iMult,jHisto);
            hREC_1D_in_MT_in_PT[iMult][jHisto]   = (TH1F*)(insFile_DT_Mult->Get(hName));
        }
    }
    
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hREC_2D_in_MT               = new TH2F     *[nBinMult+1];
    TH1F     ***hREC_1D_in_MT_in_PT_2D_bin  = new TH1F    **[nBinMult+1];
    TH2F    ****hREC_2D_in_MT_in_PT         = new TH2F   ***[nBinMult+1];
    //
    //  Defining MT-Differential histograms over measurable pT
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hName   =   Form("hREC_2D_in_MT_%i",iMult);
        hREC_2D_in_MT[iMult]   =   (TH2F*)(insFile_DT_Mult->Get(hName));

        hREC_1D_in_MT_in_PT_2D_bin[iMult]  = new
        TH1F     *[nBinPT2D];
        hREC_2D_in_MT_in_PT[iMult]         = new TH2F    **[nBinPT2D];
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_in_MT_%i",iMult,jHisto);
            hREC_1D_in_MT_in_PT_2D_bin[iMult][jHisto]  =   (TH1F*)(insFile_DT_Mult->Get(hName));
            
            hREC_2D_in_MT_in_PT[iMult][jHisto]         = new TH2F     *[nBinPT2D];
            for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_in_MT_%i",jHisto,kHisto,iMult);
                hREC_2D_in_MT_in_PT[iMult][jHisto][kHisto]    =    (TH2F*)(insFile_DT_Mult->Get(hName));
            }
        }
    }
    
    // >-> RAPIDITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_RP               = new TH1F     *[nBinRap_];
    TH1F     ***hREC_1D_in_RP_in_PT         = new TH1F    **[nBinRap_];
    //
    //  Defining RP-Differential histograms over measurable pT
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName = Form("hREC_1D_in_RP_%i",iRap);
        hREC_1D_in_RP[iRap]   =   (TH1F*)(insFile_DT_Rap->Get(hName));
        
        hREC_1D_in_RP_in_PT[iRap] = new TH1F     *[nBinPT1D];

        for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
        {
            hName = Form("hREC_1D_in_RP_PT_%i_%i",iRap,jHisto);
            hREC_1D_in_RP_in_PT[iRap][jHisto]   = (TH1F*)(insFile_DT_Rap->Get(hName));
        }
    }

    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hREC_2D_in_RP               = new TH2F     *[nBinRap_];
    TH1F     ***hREC_1D_in_RP_in_PT_2D_bin  = new TH1F    **[nBinRap_];
    TH2F    ****hREC_2D_in_RP_in_PT         = new TH2F   ***[nBinRap_];
    //
    //  Defining RP-Differential histograms over measurable pT
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName   =   Form("hREC_2D_in_RP_%i",iRap);
        hREC_2D_in_RP[iRap]   =   (TH2F*)(insFile_DT_Rap->Get(hName));

        hREC_1D_in_RP_in_PT_2D_bin[iRap]  = new
        TH1F     *[nBinPT2D];
        hREC_2D_in_RP_in_PT[iRap]         = new TH2F    **[nBinPT2D];
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_in_RP_%i",jHisto,iRap);
            hREC_1D_in_RP_in_PT_2D_bin[iRap][jHisto]  =   (TH1F*)(insFile_DT_Rap->Get(hName));
            hREC_2D_in_RP_in_PT[iRap][jHisto]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_in_RP_%i",jHisto,kHisto,iRap);
                hREC_2D_in_RP_in_PT[iRap][jHisto][kHisto]    =    (TH2F*)(insFile_DT_Rap->Get(hName));
            }
        }
    }
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    fSetBinRap_();
    fSetBinMult();
    fSetBinNTup();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept;
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    //>>    Making utility variables
    Int_t fTotalCount, fProgrCount;
    //
    if ( kDoYield )   {
        //
        //>>    Timer
        fStartTimer("Yield Analysis Signal Extrapolation");
        //
        //>>    Total Fit number abnd progressive
        fTotalCount = nBinPT1D + nBinPT2D * ( nBinPT2D + 1 );
        fProgrCount = 0;
        //
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,"Yield")));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield"))+TString("1D/"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kASigExtr_Plot_Direct,"Yield"))+TString("2D/"));
        TFile      *outCheckFitYld  =   new TFile   (Form(kASigExtr_FitCheckPlt,"Yield"),"recreate");
        //
        //>>    Fit the Model
        auto    fFitResults_1DYield = FitModel        (hREC_1D_in_PT,hSlop_Reference,Form(kASigExtr_Plot_Direct,"Yield")+TString("1D/"));
        //
        //>>    Progressive Count
        fProgrCount +=  nBinPT1D;
        //
        //>>    Print Progress
        fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
        //
        //>>    Fit
        std::vector<TH1F*>  f1DCheck;
        auto    fFitResults_2DYield = FitModel(hREC_1D_in_PT_2D_bin,hSlop_Referen2D,hREC_2D_in_PT,f1DCheck,Form(kASigExtr_Plot_Direct,"Yield")+TString("2D/"));
        //
        //Progressive Count
        fProgrCount +=  nBinPT2D;
        //
        fStopTimer("Yield Analysis Signal Extrapolation");
        //
        // >> YIELD ANALYSIS //
        TFile *outFil2  =   new TFile   (Form(kASigExtr_FitCheckRst,"Yield"),"recreate");
        //
        hEvntEff->Write();
        hEvntMlt->Write();
        for ( auto hSave : fFitResults_1DYield )    {
            if ( strncmp(hSave->GetName(),"hRAW_1D",7) == 0 )   hSave->Scale(1.,"width");
            if ( strncmp(hSave->GetName(),"h1D_anBB",7) == 0 )  hSave->Scale(1.,"width");
            hSave   ->  Write();
            if ( strncmp(hSave->GetName(),"h1D_bWidt",7) == 0 ) uPlotReferenceValue(hSave,kPhiMesonWidth,kPhiMesonWdErr)->SaveAs("./result/wid.pdf");
            if ( strncmp(hSave->GetName(),"h1D_bSlop",7) == 0 ) uPlotReferenceValue(hSave,0.001,0.0001)->SaveAs("./result/res.pdf");
            if ( strncmp(hSave->GetName(),"h1D_bMass",7) == 0 ) uPlotReferenceValue(hSave,kPhiMesonMass_,kPhiMesonMsErr)->SaveAs("./result/mas.pdf");
        }
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
        outCheckFitYld->Close();
        //
        outFil2->Close();
    }
    //
    /*
    if ( kDoRapidity )   {
        //
        //>>    Timer
        fStartTimer("Rapidity Analysis Signal Extrapolation");
        //
        //>>    Total Fit number abnd progressive
        fTotalCount = nBinRap_ * ( nBinPT1D + nBinPT2D * ( nBinPT2D + 1 ) );
        fProgrCount = 0;
        //
        TFile      *outCheckFit     =   new TFile(fRapSigChek,"recreate");
        gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/rapidity/ExtractionCheck/1D/"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/rapidity/ExtractionCheck/2D/"));
        //
        std::vector<std::vector<TH1F*>> fFitResults_1DYield;
        std::vector<std::vector<TH2F*>> fFitResults_2DYield;
        std::vector<std::vector<TH1F*>> f1DCheck;
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ ) {
            gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/rapidity/ExtractionCheck/1D/%i",iRap));
            //
            //>>    Fit the Model
            //fFitResults_1DYield.push_back ( FitModel        (hREC_1D_in_RP_in_PT[iRap],hSlop_Reference,Form("./result/rapidity/ExtractionCheck/1D/%i",iRap),Form("RAP_%i",iRap)) );
            //
            //>>    Progressive Count
            fProgrCount +=  nBinPT1D;
            //
            //>>    Print Progress
            fPrintLoopTimer("Rapidity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/rapidity/ExtractionCheck/2D/%i",iRap));
            //
            //>>    Fit
            std::vector<TH1F*>  f1DCheck_Temp;
            fFitResults_2DYield.push_back( FitModel(hREC_1D_in_PT_2D_bin,hSlop_Referen2D,hREC_2D_in_RP_in_PT[iRap],f1DCheck_Temp,Form("./result/rapidity/ExtractionCheck/2D/%i",iRap),Form("RAP_%i",iRap)) );
            f1DCheck.push_back(f1DCheck_Temp);
            //
            //Progressive Count
            fProgrCount +=  nBinPT2D;
            fPrintLoopTimer("Rapidity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
        }
        //
        fStopTimer("Rapidity Analysis Signal Extrapolation");
        //
        // >> YIELD ANALYSIS //
        TFile *outFil2  =   new TFile   (fRapSigExtr,"recreate");
        //
        hEvntEff->Write();
        hEvntMlt->Write();

        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ ) {
            
            for ( auto hSave : fFitResults_1DYield.at(iRap) )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"h1D_anBB",7) == 0 )  hSave->Scale(1.,"width");
                hSave   ->  SetName(Form("%s_in_Rap_%i",hSave->GetName(),iRap));
                hSave   ->  Write();
            }
            for ( auto hSave : f1DCheck.at(iRap) )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D_in_2D_bin",17) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"2Dbin_anBB",17) == 0 )           hSave->Scale(1.,"width");
                hSave   ->  SetName(Form("%s_in_Rap_%i",hSave->GetName(),iRap));
                hSave   ->  Write();
            }
            for ( auto hSave : fFitResults_2DYield.at(iRap) )    {
                if ( strncmp(hSave->GetName(),"hRAW_2D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anSB2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBS2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBB2D",6) == 0 )    hSave->Scale(1.,"width");
                hSave   ->  SetName(Form("%s_in_Rap_%i",hSave->GetName(),iRap));
                hSave   ->  Write();
            }
        }
        //
        outCheckFit->Close();
        //
        outFil2->Close();
    }
    //
    if ( kDoMultiplicity )   {
        //
        //>>    Timer
        fStartTimer("Multiplicity Analysis Signal Extrapolation");
        //
        //>>    Total Fit number abnd progressive
        fTotalCount = nBinMult * ( nBinPT1D + nBinPT2D * ( nBinPT2D + 1 ) );
        fProgrCount = 0;
        //
        TFile      *outCheckFit     =   new TFile(fMltSigChek,"recreate");
        gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/multiplicity/ExtractionCheck/1D/"));
        gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/multiplicity/ExtractionCheck/2D/"));
        //
        std::vector<std::vector<TH1F*>> fFitResults_1DYield;
        std::vector<std::vector<TH2F*>> fFitResults_2DYield;
        std::vector<std::vector<TH1F*>> f1DCheck;
        for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
            gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/multiplicity/ExtractionCheck/1D/%i",iMlt));
            //
            //>>    Fit the Model
            fFitResults_1DYield.push_back ( FitModel        (hREC_1D_in_MT_in_PT[iMlt],hSlop_Reference,Form("./result/multiplicity/ExtractionCheck/1D/%i",iMlt),Form("MLT_%i",iMlt)) );
            //
            //>>    Progressive Count
            fProgrCount +=  nBinPT1D;
            //
            //>>    Print Progress
            fPrintLoopTimer("Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            gROOT                       ->  ProcessLine(Form(".! mkdir -p ./result/multiplicity/ExtractionCheck/2D/%i",iMlt));
            //
            //>>    Fit
            std::vector<TH1F*>  f1DCheck_Temp;
            //fFitResults_2DYield.push_back( FitModel(hREC_1D_in_PT_2D_bin,hSlop_Referen2D,hREC_2D_in_MT_in_PT[iMlt],f1DCheck_Temp,Form("./result/multiplicity/ExtractionCheck/2D/%i",iMlt),Form("MLT_%i",iMlt)) );
            //f1DCheck.push_back(f1DCheck_Temp);
            //
            //Progressive Count
            fProgrCount +=  nBinPT2D;
            fPrintLoopTimer("Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
        }
        //
        fStopTimer("Multiplicity Analysis Signal Extrapolation");
        //
        // >> YIELD ANALYSIS //
        TFile *outFil2  =   new TFile   (fMltSigExtr,"recreate");
        //
        hEvntEff->Write();
        hEvntMlt->Write();

        for ( Int_t iMlt = 0; iMlt <= nBinMult; iMlt++ ) {
            for ( auto hSave : fFitResults_1DYield.at(iMlt) )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"h1D_anBB",7) == 0 )  hSave->Scale(1.,"width");
                hSave   ->  SetName(Form("%s_in_Mlt_%i",hSave->GetName(),iMlt));
                hSave   ->  Write();
            }
            for ( auto hSave : f1DCheck.at(iMlt) )    {
                if ( strncmp(hSave->GetName(),"hRAW_1D_in_2D_bin",17) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"2Dbin_anBB",17) == 0 )           hSave->Scale(1.,"width");
                hSave   ->  SetName(Form("%s_in_Mlt_%i",hSave->GetName(),iMlt));
                hSave   ->  Write();
            }
            for ( auto hSave : fFitResults_2DYield.at(iMlt) )    {
                if ( strncmp(hSave->GetName(),"hRAW_2D",7) == 0 )   hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anSB2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBS2D",6) == 0 )    hSave->Scale(1.,"width");
                if ( strncmp(hSave->GetName(),"anBB2D",6) == 0 )    hSave->Scale(1.,"width");
                hSave   ->  SetName(Form("%s_in_Mlt_%i",hSave->GetName(),iMlt));
                hSave   ->  Write();
            }
        }
        //
        outCheckFit->Close();
        //
        outFil2->Close();
    }
     */
    //
    // >-> Close input File
    //
    insFile_DT_Yield    ->Close();
    insFile_DT_Rap      ->Close();
    insFile_DT_Mult     ->Close();
}
