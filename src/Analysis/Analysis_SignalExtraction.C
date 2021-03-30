#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void Analysis_SignalExtraction ( bool fSilent = true, TString fOption = "" )
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
    TFile*  insFile_DT_Yield            =   new TFile   (fYldPreProc);
    TFile*  insFile_DT_Mult             =   new TFile   (fMltPreProc);
    TFile*  insFile_DT_Rap              =   new TFile   (fRapPreProc);
    
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
    //
    //  Utility
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "fQC_Event_Enum_Mult";
    hEvntMlt    =   (TH1D*)(insFile_DT_Mult->Get(hName));
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
    TH1F      **hREC_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F     ***hREC_1D_in_MT_in_PT         = new TH1F    **[nBinMult];
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
    TH2F      **hREC_2D_in_MT               = new TH2F     *[nBinMult];
    TH1F     ***hREC_1D_in_MT_in_PT_2D_bin  = new TH1F    **[nBinMult];
    TH2F    ****hREC_2D_in_MT_in_PT         = new TH2F   ***[nBinMult];
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
            hName = Form("hREC_1D_in_PT_2D_bin_%i_%i",iMult,jHisto);
            hREC_1D_in_MT_in_PT_2D_bin[iMult][jHisto]  =   (TH1F*)(insFile_DT_Mult->Get(hName));
            
            hREC_2D_in_MT_in_PT[iMult][jHisto]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_%i",iMult,jHisto,kHisto);
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
            hName = Form("hREC_1D_in_PT_2D_bin_%i_%i",iRap,jHisto);
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
    hRAW_Mss_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Mss_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hRAW_Wdt_1D");
    hTitle      =   Form("hRAW_Wdt_1D");
    hRAW_Wdt_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Wdt_1D,"PT 1D");
    //
    hName       =   Form("hRAW_Wdt_1D_in_2D_bin");
    hTitle      =   Form("hRAW_Wdt_1D_in_2D_bin");
    hRAW_Wdt_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Wdt_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hRAW_Slp_1D");
    hTitle      =   Form("hRAW_Slp_1D");
    hRAW_Slp_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_Slp_1D,"PT 1D");
    //
    hName       =   Form("hRAW_Slp_1D_in_2D_bin");
    hTitle      =   Form("hRAW_Slp_1D_in_2D_bin");
    hRAW_Slp_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
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
    
    // >-> MULTIPLICITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hRAW_1D_in_MT               = new TH1F     *[nBinMult];
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hName       =   Form("hRAW_1D_in_MT_%i",iMult);
        hTitle      =   Form("hRAW_1D_in_MT_%i",iMult);
        hRAW_1D_in_MT[iMult]     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hRAW_1D_in_MT[iMult],"PT 1D");
    }
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hRAW_2D_in_MT               = new TH2F     *[nBinMult];
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hName       =   Form("hRAW_2D_in_MT_%i",iMult);
        hTitle      =   Form("hRAW_2D_in_MT_%i",iMult);
        hRAW_2D_in_MT[iMult]     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hRAW_2D_in_MT[iMult],"PT 2D");
    }
    
    // >-> RAPIDITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hRAW_1D_in_RP               = new TH1F     *[nBinRap_];
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName       =   Form("hRAW_1D_in_RP_%i",iRap);
        hTitle      =   Form("hRAW_1D_in_RP_%i",iRap);
        hRAW_1D_in_RP[iRap]     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hRAW_1D_in_RP[iRap],"PT 1D");
    }
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hRAW_2D_in_RP               = new TH2F     *[nBinRap_];
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName       =   Form("hRAW_2D_in_RP_%i",iRap);
        hTitle      =   Form("hRAW_2D_in_RP_%i",iRap);
        hRAW_2D_in_RP[iRap]     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hRAW_2D_in_RP[iRap],"PT 2D");
    }
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    // Storing Results for shape
    //
    RooFitResult  **fShapeStore =   new RooFitResult   *[nBinPT2D];
    //
    //  Making utility variables
    Int_t fTotalCount, fProgrCount;
    //
    if ( kDoYield || kDoRapidity )   {
        fStartTimer("Yield Analysis Signal Extrapolation");
        
        // Output File for Fit Check
        TFile*  outCheckFitYld  =   new TFile(fYldSigChek,"recreate");
         
        // Total Fit number abnd progressive
        fTotalCount = nBinPT1D+(1+nBinPT2D)*nBinPT2D;
        fProgrCount = 0;
        //
        // Starting cycle
        cout << Form("[INFO] Starting yield analysis in 1D") << endl;
        for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
        {
            if ( kDoRapidity && !kDoYield ) { cout << Form("[INFO] Only recovering 2D Shapes for Rapidity analysis, skipping this one... ") << endl; break; }
            //Progressive Count
            fProgrCount++;
            //
            // Print Progress
            fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            //
            // Fit
            auto fResults       =   FitModel(hREC_1D_in_PT[iFit],"1D",true,iFit,1);
            //
            if ( !fResults ) continue;
            //
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
            hRAW_1D->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_1D->SetBinError        (iFit+1,N_Raw->getError());
            //
            // Building N_Raw mass histogram
            N_Raw           = static_cast<RooRealVar*>(fResults->floatParsFinal().at(kMass));
            hRAW_Mss_1D->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_Mss_1D->SetBinError        (iFit+1,N_Raw->getError());
            //
            // Building N_Raw width histogram
            N_Raw           = static_cast<RooRealVar*>(fResults->floatParsFinal().at(kWidt));
            hRAW_Wdt_1D->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_Wdt_1D->SetBinError        (iFit+1,N_Raw->getError());
            //
            // Building N_Raw width histogram
            N_Raw           = static_cast<RooRealVar*>(fResults->floatParsFinal().at(kSlop));
            hRAW_Slp_1D->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_Slp_1D->SetBinError        (iFit+1,N_Raw->getError());
        }
        cout << Form("[INFO] Starting yield analysis in 1D in 2D bins") << endl;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            //Progressive Count
            fProgrCount++;
            
            // Print Progress
            fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            
            // Fit
            fShapeStore[iFit]   =   FitModel(hREC_1D_in_PT_2D_bin[iFit],"2D_bin",true,iFit,2);
            
            if ( !fShapeStore[iFit]  ) continue;
            //
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fShapeStore[iFit] ->floatParsFinal().at(Signal));
            hRAW_1D_in_2D_bin->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_1D_in_2D_bin->SetBinError        (iFit+1,N_Raw->getError());
            //
            // Building N_Raw histogram
            N_Raw      = static_cast<RooRealVar*>(fShapeStore[iFit] ->floatParsFinal().at(kMass));
            hRAW_Mss_1D_in_2D_bin->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_Mss_1D_in_2D_bin->SetBinError        (iFit+1,N_Raw->getError());
            //
            // Building N_Raw histogram
            N_Raw      = static_cast<RooRealVar*>(fShapeStore[iFit] ->floatParsFinal().at(kWidt));
            hRAW_Wdt_1D_in_2D_bin->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_Wdt_1D_in_2D_bin->SetBinError        (iFit+1,N_Raw->getError());
            //
            // Building N_Raw histogram
            N_Raw      = static_cast<RooRealVar*>(fShapeStore[iFit] ->floatParsFinal().at(kSlop));
            hRAW_Slp_1D_in_2D_bin->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_Slp_1D_in_2D_bin->SetBinError        (iFit+1,N_Raw->getError());
        }
        //
        cout << Form("[INFO] Starting yield analysis in 2D") << endl;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            if ( kDoRapidity && !kDoYield ) { cout << Form("[INFO] Only recovering 2D Shapes for Rapidity analysis, skipping this one... ") << endl; break; }
            for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ )
            {
                //Progressive Count
                fProgrCount++;
                
                // Print Progress
                fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                
                if ( !fShapeStore[iFit] ) continue;
                if ( !fShapeStore[jFit] ) continue;
                
                // Fit
                auto fResults       =   FitModel(hREC_2D_in_PT[iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],"2D",true,iFit,jFit);
                
                if ( !fResults ) continue;
                
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(SignlSignl));
                hRAW_2D->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                hRAW_2D->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
            }
        }
        //
        if ( !kDoRapidity ) { cout<< "IM HERE " << endl; delete[] fShapeStore; }
        fStopTimer("Yield Analysis Signal Extrapolation");
        
        outCheckFitYld->Close();
    }
    //
    if ( kDoRapidity )   {
        fStartTimer("Yield in Rapidity Analysis Signal Extrapolation");
        
        // Total Fit number and progressive
        fTotalCount = (nBinPT1D+(1+nBinPT2D)*nBinPT2D)*(nBinRap_);
        fProgrCount = 0;
        
        // Output File for Fit Check
        TFile*  outCheckFitRap  =   new TFile(fRapSigChek,"recreate");
        
        // Starting cycle
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
        {
            //
            cout << Form("[INFO] Starting Rapidity yield analysis in bin [%.2f;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]) << endl;
            //
            cout << Form("[INFO] Starting Rapidity yield analysis in 1D") << endl;
            for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
            {
                //Progressive Count
                fProgrCount++;
                
                // Print Progress
                fPrintLoopTimer("Yield in Rapidity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                
                // Fit
                auto fResults = FitModel(hREC_1D_in_RP_in_PT[iRap][iFit],Form("1D_RP_%i",iRap),true,iFit,1);
                
                if ( !fResults ) continue;
                
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
                hRAW_1D_in_RP[iRap]->SetBinContent     (iFit+1,N_Raw->getVal());
                hRAW_1D_in_RP[iRap]->SetBinError       (iFit+1,N_Raw->getError());
            }
            //
            cout << Form("[INFO] Starting Rapidity yield analysis in 2D") << endl;
            for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
            {
                for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ )
                {
                    //Progressive Count
                    fProgrCount++;
                    
                    // Print Progress
                    fPrintLoopTimer("Yield in Rapidity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                    
                    if ( !fShapeStore[iFit] ) break;
                    if ( !fShapeStore[jFit] ) continue;
                    
                    // Fit
                    auto fResults       =   FitModel(hREC_2D_in_RP_in_PT[iRap][iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],Form("2D_RP_%i",iRap),true,iFit,jFit,"");
                    
                    if ( !fResults ) continue;
                    
                    // Building N_Raw histogram
                    auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(SignlSignl));
                    hRAW_2D_in_RP[iRap]->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                    hRAW_2D_in_RP[iRap]->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
                }
            }
        }
        outCheckFitRap->Close();
        //
        delete[] fShapeStore;
        //
        fStopTimer("Yield in Rapidity Analysis Signal Extrapolation");
    }
    //
    if ( kDoMultiplicity )   {
        fStartTimer("Yield in Multiplicity Analysis Signal Extrapolation");
        
        // Total Fit number and progressive
        fTotalCount = (nBinPT1D+(1+nBinPT2D)*nBinPT2D)*(nBinMult+1);
        fProgrCount = 0;
        
        // Output File for Fit Check
        TFile*  outCheckFitMlt  =   new TFile(fMltSigChek,"recreate");
        
        // Starting cycle
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
        {
            // Storing Results for shape
            RooFitResult  **fShapeStore =   new RooFitResult   *[nBinPT2D];
            //
            if ( iMult != 0 )   cout << Form("[INFO] Starting Multiplicity yield analysis in bin [%.2f;%.2f]",fArrMult[iMult-1],fArrMult[iMult]) << endl;
            else                cout << Form("[INFO] Starting Multiplicity yield analysis in bin [%.2f;%.2f]",fArrMult[0],fArrMult[nBinMult]) << endl;
            //
            cout << Form("[INFO] Starting Multiplicity yield analysis in 1D") << endl;
            for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
            {
                //Progressive Count
                fProgrCount++;
                
                // Print Progress
                fPrintLoopTimer("Yield in Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                
                // Fit
                auto fResults = FitModel(hREC_1D_in_MT_in_PT[iMult][iFit],Form("1D_MT_%i",iMult),true,iFit,1);
                
                if ( !fResults ) continue;
                
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
                hRAW_1D_in_MT[iMult]->SetBinContent     (iFit+1,N_Raw->getVal());
                hRAW_1D_in_MT[iMult]->SetBinError       (iFit+1,N_Raw->getError());
            }
            //
            cout << Form("[INFO] Starting Multiplicity yield analysis in 1D in 2D bins") << endl;
            for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
            {
                //Progressive Count
                fProgrCount++;
                
                // Print Progress
                fPrintLoopTimer("Yield in Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                
                // Fit
                fShapeStore[iFit]   =   FitModel(hREC_1D_in_MT_in_PT_2D_bin[iMult][iFit],Form("2D_bin_MT_%i",iMult),true,iFit,2);
            }
            //
            cout << Form("[INFO] Starting Multiplicity yield analysis in 2D") << endl;
            for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
            {
                for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ )
                {
                    //Progressive Count
                    fProgrCount++;
                    
                    // Print Progress
                    fPrintLoopTimer("Yield in Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                    
                    if ( !fShapeStore[iFit] ) break;
                    if ( !fShapeStore[jFit] ) continue;
                    
                    // Fit
                    auto fResults       =   FitModel(hREC_2D_in_MT_in_PT[iMult][iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],Form("2D_MT_%i",iMult),true,iFit,jFit);
                    
                    if ( !fResults ) continue;
                    
                    // Building N_Raw histogram
                    auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(SignlSignl));
                    hRAW_2D_in_MT[iMult]->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                    hRAW_2D_in_MT[iMult]->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
                }
            }
            delete[] fShapeStore;
        }
        
        outCheckFitMlt->Close();
        
        fStopTimer("Yield in Multiplicity Analysis Signal Extrapolation");
    }
    //
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    
    // >> YIELD ANALYSIS //
    //
    hRAW_1D->Scale(1.,"width");
    hRAW_2D->Scale(1.,"width");
    //
    // >> MULTIPLICITY ANALYSIS //
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hRAW_1D_in_MT[iMult]->Scale(1.,"width");
        hRAW_2D_in_MT[iMult]->Scale(1.,"width");
    }
    //
    // >> RAPIDITY ANALYSIS //
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hRAW_1D_in_RP[iRap]->Scale(1.,"width");
        hRAW_2D_in_RP[iRap]->Scale(1.,"width");
    }
    //
    //
    // >> TRIGGER ANALYSIS //
    //
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Trigger Analysis
    //
    //TFile *outFil1  =   new TFile   (fTrgSigExtr,"recreate");
    //
    //outFil1->Close();
    //
    // >> Yield Analysis
    //
    if ( kDoYield ) {
        TFile *outFil2  =   new TFile   (fYldSigExtr,"recreate");
        //
        hEvntEff->Write();
        hEvntMlt->Write();
        hRAW_1D->Write();
        hRAW_1D_in_2D_bin->Write();
        hRAW_2D->Write();
        //
        outFil2->Close();
    }
    //
    // >> Multiplicity Analysis
    //
    if ( kDoMultiplicity ) {
        TFile *outFil3  =   new TFile   (fMltSigExtr,"recreate");
        //
        hEvntEff->Write();
        hEvntMlt->Write();
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
        {
            hRAW_1D_in_MT[iMult]->Write();
            hRAW_2D_in_MT[iMult]->Write();
        }
        //
        outFil3->Close();
    }
    //
    // >> Rapidity Analysis
    //
    if ( kDoRapidity ) {
        TFile *outFil4  =   new TFile   (fRapSigExtr,"recreate");
        //
        hEvntEff->Write();
        hEvntMlt->Write();
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
        {
            hRAW_1D_in_RP[iRap]->Write();
            hRAW_2D_in_RP[iRap]->Write();
        }
        //
        outFil4->Close();
    }
    //
    if ( kDosss )   {
        TFile *outFil5  =   new TFile   ("./result/test.root","recreate");
        //
        hRAW_Mss_1D             ->Write();
        hRAW_Mss_1D_in_2D_bin   ->Write();
        hRAW_Wdt_1D             ->Write();
        hRAW_Wdt_1D_in_2D_bin   ->Write();
        hRAW_Slp_1D             ->Write();
        hRAW_Slp_1D_in_2D_bin   ->Write();
        //
        outFil5->Close();
    }
    // >-> Close input File
    //
    insFile_DT_Yield    ->Close();
    insFile_DT_Mult     ->Close();
}
