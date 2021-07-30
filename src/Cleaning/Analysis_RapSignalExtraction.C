#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void Analysis_RapSignalExtraction ( bool fSilent = true, TString fOption = "" )
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
    fChooseOption(fOption);
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (fRapPreProc);
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> RAPIDITY ANALYSIS //
    //
    // Rapidity Analysis
    //
    hName   =   "hBKG_BKG_Rap";
    TH1D    *hBKG_BKG_Rap           =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName   =   "hBKG_SIG_Rap";
    TH1D    *hBKG_SIG_Rap           =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName   =   "hBKG_BKG_BKG_SIG_Rap";
    TH1D    *hBKG_BKG_BKG_SIG_Rap   =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
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

    // >-> RAPIDITY ANALYSIS //
    //
    // Rapidity Analysis
    //
    hName   =   "hSIG_Rap";
    hTitle  =   "Rapidity Phi hSIG_Rap";
    TH1D    *hSIG_Rap   =   new TH1D (hName,hTitle,nBinRap_,fArrRap_);
    SetAxis(hSIG_Rap,"PT 2D");
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    // Storing Results for shape
    //
    TFile*  outCheckFitYld  =   new TFile("example.root","recreate");
    //
    FitModelRap(hBKG_SIG_Rap,hBKG_BKG_Rap,hBKG_BKG_BKG_SIG_Rap);
    
    outCheckFitYld->Close();
    /*
    //  Making utility variables
    Int_t fTotalCount, fProgrCount;
    //
    if ( kDoYield )   {
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
            //Progressive Count
            fProgrCount++;
            //
            // Print Progress
            fPrintLoopTimer("Yield Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
            //
            // Fit
            auto fResults       =   FitModel(hREC_1D_in_PT[iFit],"1D",true,iFit,1);
            
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
            fShapeStore[iFit]   =   FitModel(hREC_1D_in_PT_2D_bin[iFit],"2D_bin",true,iFit,2);
            
            if ( !fShapeStore[iFit]  ) continue;
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fShapeStore[iFit] ->floatParsFinal().at(Signal));
            hRAW_1D_in_2D_bin->SetBinContent      (iFit+1,N_Raw->getVal());
            hRAW_1D_in_2D_bin->SetBinError        (iFit+1,N_Raw->getError());
        }
        //
        cout << Form("[INFO] Starting yield analysis in 2D") << endl;
        for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ )
        {
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
        delete[] fShapeStore;
        fStopTimer("Yield Analysis Signal Extrapolation");
        
        outCheckFitYld->Close();
    }
    //
    if ( kDoMultiplicity )   {
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
            //
            cout << Form("[INFO] Starting Multiplicity yield analysis in bin [%1.2f-%1.2f]",fArrMult[iMult],fArrMult[iMult+1]) << endl;
            //
            cout << Form("[INFO] Starting Multiplicity yield analysis in 1D") << endl;
            for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
            {
                //Progressive Count
                fProgrCount++;
                
                // Print Progress
                fPrintLoopTimer("Yield in Multiplicity Analysis Signal Extrapolation",fProgrCount,fTotalCount,1);
                
                // Fit
                auto fResults = FitModel(hREC_1D_in_MT_in_PT[iMult][iFit],Form("1D_%i",iMult),true,iFit,1);
                
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
                fShapeStore[iFit]   =   FitModel(hREC_1D_in_MT_in_PT_2D_bin[iMult][iFit],Form("2D_bin_%i",iMult),true,iFit,2);
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
                    auto fResults       =   FitModel(hREC_2D_in_MT_in_PT[iMult][iFit][jFit],fShapeStore[iFit],fShapeStore[jFit],Form("2D_%i",iMult),true,iFit,jFit);
                    
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
     */
    //
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    //
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Trigger Analysis
    //
    TFile *outFil1  =   new TFile   ("example.root","recreate");
    //
    outFil1->Close();
    //
    //
    // >-> Close input File
    //
    insFile_DT_Yield    ->Close();
}
