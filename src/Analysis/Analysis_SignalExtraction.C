#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void Analysis_SignalExtraction ( bool fSilent = true )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //>-> OPTIONS
    
    // Silencing warnings for smoother running
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (fYldPreProc);
    TFile*  insFile_DT_Mult             =   new TFile   (fMltPreProc);
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> YIELD ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1D       *hEvntEff;
    TH1F       *hREC_1D;
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    //
    //  Utility
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
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
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iHisto);
        hREC_1D_in_MT[iHisto]   =   (TH1F*)(insFile_DT_Mult->Get(hName));
        
        hREC_1D_in_MT_in_PT[iHisto] = new TH1F     *[nBinPT1D];

        for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
        {
            hName = Form("hREC_1D_in_MT_PT_%i_%i",iHisto,jHisto);
            hREC_1D_in_MT_in_PT[iHisto][jHisto]   = (TH1F*)(insFile_DT_Mult->Get(hName));
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
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName   =   Form("hREC_2D_in_MT_%i",iHisto);
        hREC_2D_in_MT[iHisto]   =   (TH2F*)(insFile_DT_Mult->Get(hName));

        hREC_1D_in_MT_in_PT_2D_bin[iHisto]  = new TH1F     *[nBinPT2D];
        hREC_2D_in_MT_in_PT[iHisto]         = new TH2F    **[nBinPT2D];
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_%i",iHisto,jHisto);
            hREC_1D_in_MT_in_PT_2D_bin[iHisto][jHisto]  =   (TH1F*)(insFile_DT_Mult->Get(hName));
            
            hREC_2D_in_MT_in_PT[iHisto][jHisto]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_%i",iHisto,jHisto,kHisto);
                hREC_2D_in_MT_in_PT[iHisto][jHisto][kHisto]    =    (TH2F*)(insFile_DT_Mult->Get(hName));
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
    //
    //  Defining yield histogram over measurable pT
    //
    hName       =   Form("hRAW_1D");
    hTitle      =   Form("hRAW_1D");
    hRAW_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hRAW_1D,"PT 1D");
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
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName       =   Form("hRAW_1D_in_MT_%i",iHisto);
        hTitle      =   Form("hRAW_1D_in_MT_%i",iHisto);
        hRAW_1D_in_MT[iHisto]     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hRAW_1D_in_MT[iHisto],"PT 1D");
    }
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hRAW_2D_in_MT               = new TH2F     *[nBinMult];
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName       =   Form("hRAW_2D_in_MT_%i",iHisto);
        hTitle      =   Form("hRAW_2D_in_MT_%i",iHisto);
        hRAW_2D_in_MT[iHisto]     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hRAW_2D_in_MT[iHisto],"PT 2D");
    }
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Yield Analysis Signal Extrapolation");
    
    // Output File for Fit Check
    TFile*  outCheckFitYld  =   new TFile(fYldSigChek,"recreate");
     
    // Storing Results for shape
    RooFitResult  **fShapeStore =   new RooFitResult   *[nBinPT2D];
     
    // Total Fit number abnd progressive
    Int_t   fTotalCount = nBinPT1D+(1+nBinPT2D)*nBinPT2D;
    Int_t   fProgrCount = 0;
    
    // Starting cycle
    cout << Form("[INFO] Starting yield analysis in 1D") << endl;
    for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
    {
        //Progressive Count
        fProgrCount++;
        //
        // Print Progress
        fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
        //
        // Fit
        auto fResults       =   FitModel(hREC_1D_in_PT[iFit],"CH3",true,iFit,1);
        
        if ( !fResults ) continue;
        
        // Building N_Raw histogram
        auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
        hRAW_1D->SetBinContent      (iFit+1,N_Raw->getVal());
        hRAW_1D->SetBinError        (iFit+1,N_Raw->getError());
    }
    
    cout << Form("[INFO] Starting yield analysis in 1D in 2D bins") << endl;
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        //Progressive Count
        fProgrCount++;
        
        // Print Progress
        fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
        
        // Fit
        fShapeStore[iFit]   =   FitModel(hREC_1D_in_PT_2D_bin[iFit],"CH3",true,iFit,2);
    }
     
    cout << Form("[INFO] Starting yield analysis in 2D") << endl;
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            //Progressive Count
            fProgrCount++;
            
            // Print Progress
            fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            
            if ( !fShapeStore[iFit] ) break;
            if ( !fShapeStore[jFit] ) continue;
            
            // Fit
            auto fResults       =   FitModel(hREC_2D_in_PT[iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],"CH3",true,iFit,jFit);
            
            if ( !fResults ) continue;
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
            hRAW_2D->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            hRAW_2D->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
     
    delete[] fShapeStore;
    fStopTimer("Yield Analysis Signal Extrapolation");
    
    outCheckFitYld->Close();
    
    fStartTimer("Yield in Multiplicity Analysis Signal Extrapolation");
    
    // Total Fit number and progressive
    fTotalCount = (nBinPT1D+(1+nBinPT2D)*nBinPT2D)*nBinMult;
    fProgrCount = 0;
    
    // Output File for Fit Check
    TFile*  outCheckFitMlt  =   new TFile(fMltSigChek,"recreate");
    
    // Starting cycle
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )
    {
        // Storing Results for shape
        RooFitResult  **fShapeStore =   new RooFitResult   *[nBinPT2D];
        
        cout << Form("[INFO] Starting Multiplicity yield analysis in bin [%1.2f-%1.2f]",fArrMult[iMult],fArrMult[iMult+1]) << endl;
        
        cout << Form("[INFO] Starting Multiplicity yield analysis in 1D") << endl;
        for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
        {
            //Progressive Count
            fProgrCount++;
            
            // Print Progress
            fPrintLoopTimer("Yield in Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            
            // Fit
            auto fResults = FitModel(hREC_1D_in_MT_in_PT[iMult][iFit],"CH3",true,iFit,1);
            
            if ( !fResults ) continue;
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
            hRAW_1D_in_MT[iMult]->SetBinContent     (iFit+1,N_Raw->getVal());
            hRAW_1D_in_MT[iMult]->SetBinError       (iFit+1,N_Raw->getError());
        }
        
        cout << Form("[INFO] Starting Multiplicity yield analysis in 1D in 2D bins") << endl;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            //Progressive Count
            fProgrCount++;
            
            // Print Progress
            fPrintLoopTimer("Yield in Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            
            // Fit
            fShapeStore[iFit]   =   FitModel(hREC_1D_in_MT_in_PT_2D_bin[iMult][iFit],"",true,iFit,2);
        }
        
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
                auto fResults       =   FitModel(hREC_2D_in_MT_in_PT[iMult][iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],"",true,iFit,jFit);
                
                if ( !fResults ) continue;
                
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
                hRAW_2D_in_MT[iMult]->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                hRAW_2D_in_MT[iMult]->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
            }
        }
        delete[] fShapeStore;
    }
    
    outCheckFitMlt->Close();
    
    fStopTimer("Yield in Multiplicity Analysis Signal Extrapolation");
    
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
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRAW_1D_in_MT[iHisto]->Scale(1.,"width");
        hRAW_2D_in_MT[iHisto]->Scale(1.,"width");
    }
    //
    
    // >> TRIGGER ANALYSIS //
    
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
    TFile *outFil2  =   new TFile   (fYldSigExtr,"recreate");
    //
    hEvntEff->Write();
    hRAW_1D->Write();
    hRAW_2D->Write();
    //
    outFil2->Close();
    //
    // >> Multiplicity Analysis
    //
    TFile *outFil3  =   new TFile   (fMltSigExtr,"recreate");
    //
    hEvntEff->Write();
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hRAW_1D_in_MT[iHisto]->Write();
        hRAW_2D_in_MT[iHisto]->Write();
    }
    //
    outFil3->Close();
    //
    // >-> Close input File
    //
    insFile_DT_Yield    ->Close();
    insFile_DT_Mult     ->Close();
}

    
    /*
    //---------------------//
    //      Analysis       //-------------------------------------------------------------------------------
    //---------------------//
    
    // Output File for Fit Check
    TFile*  outFile_FT  =   new TFile(fFitResHist,"recreate");
    
    // Fit Results and PlotOn object
    RooFitResult *** Results = new RooFitResult **  [nBinPT2D];
    RooFitResult **  utility = new RooFitResult *   [nBinPT2D];
    
    //------ 1D Histogram of N_raw ------//
    
    for (Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
    {
        // Not considering pT < 0.4 GeV
        if ( fArrPT1D[iFit+1] <= 0.41 ) continue;
        
        // Fit
        auto fResults = FitModel(hIM_1D_Rec_PT_S[iFit],"",bSave,iFit,1);        // FitModel with Save enabled will write every fit performed on a Canvas
        
        // Building N_Raw histogram
        auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
        h1D_Raw->SetBinContent      (iFit+1,N_Raw->getVal());
        h1D_Raw->SetBinError        (iFit+1,N_Raw->getError());
    }
    
    //------ 2D Histogram of N_raw ------//
    
    // Preparing 1D fit for shape in pT bins
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        Results[iFit]   = new RooFitResult * [nBinPT2D];
        utility[iFit]   = FitModel(hdM_dpT_Tot_Rec[iFit],"CH3",bSave,iFit,2);
    }
    
    // 2D Fits
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
        for (int jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            // Not considering pT < 0.4 GeV
            if ( fArrPT2D[jFit+1] <= 0.41 ) continue;
            if ( iFit == 11 && jFit == 2  ) continue;
            if ( iFit == 2  && jFit == 11 ) continue;
            
            // Fit
            Results[iFit][jFit] = FitModel(hdM_dpT_Tot_Rec2D[iFit][jFit],utility[iFit],utility[jFit],"CH3",bSave,iFit,jFit);
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(Results[iFit][jFit]->floatParsFinal().at(SignlSignl));
            h2D_Raw->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            h2D_Raw->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
    */
