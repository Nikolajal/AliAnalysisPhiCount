#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void PreProcessing_data ( string fFileName = "" )
{
    if ( fFileName == "" )
    {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing_data(\"Root_file_name.root\")" << endl;
        return;
    }
    
    //Retrieving Event data
    TFile *insFileDT        =   new TFile   (fFileName.c_str());
    
    //Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileDT->Get(fPhiCandidate_Tree);
    TTree   *TKaonCandidate =   (TTree*)insFileDT->Get(fKaonCandidate_Tree);
    
    // Define tree data structures
    Struct_PhiCandidate     evPhiCandidate;
    Struct_KaonCandidate    evKaonCandidate;

    if ( !TPhiCandidate && !TKaonCandidate )
    {
        cout << "Input Data Tree not found!" << endl;
        return;
    }
    if ( !TPhiCandidate )
    {
        cout << "[INFO] No PhiCandidate Tree, switching to Kaon Analysis" << endl;
        TKaonCandidate-> SetBranchAddress   ("fMultiplicity",    &evKaonCandidate.Multiplicity);
        TKaonCandidate-> SetBranchAddress   ("nKaon",           &evKaonCandidate.nKaon);
        TKaonCandidate-> SetBranchAddress   ("Px",              &evKaonCandidate.Px);
        TKaonCandidate-> SetBranchAddress   ("Py",              &evKaonCandidate.Py);
        TKaonCandidate-> SetBranchAddress   ("Pz",              &evKaonCandidate.Pz);
        TKaonCandidate-> SetBranchAddress   ("Charge",          &evKaonCandidate.Charge);
        TKaonCandidate-> SetBranchAddress   ("TOFSigma",        &evKaonCandidate.SigmaTOF);
        TKaonCandidate-> SetBranchAddress   ("TPCSigma",        &evKaonCandidate.SigmaTPC);
    }
    else if ( !TKaonCandidate )
    {
        cout << "[INFO] No KaonCandidate Tree, switching to Phi Analysis" << endl;
        TPhiCandidate-> SetBranchAddress    ("fMultiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
    }
    else
    {
        TPhiCandidate-> SetBranchAddress    ("fMultiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
        
        TKaonCandidate-> SetBranchAddress   ("fMultiplicity",    &evKaonCandidate.Multiplicity);
        TKaonCandidate-> SetBranchAddress   ("nKaon",           &evKaonCandidate.nKaon);
        TKaonCandidate-> SetBranchAddress   ("Px",              &evKaonCandidate.Px);
        TKaonCandidate-> SetBranchAddress   ("Py",              &evKaonCandidate.Py);
        TKaonCandidate-> SetBranchAddress   ("Pz",              &evKaonCandidate.Pz);
        TKaonCandidate-> SetBranchAddress   ("Charge",          &evKaonCandidate.Charge);
        TKaonCandidate-> SetBranchAddress   ("TOFSigma",        &evKaonCandidate.SigmaTOF);
        TKaonCandidate-> SetBranchAddress   ("TPCSigma",        &evKaonCandidate.SigmaTPC);
    
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
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept;
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >>->> YIELD ANALYSIS //
    // 1D
    TH1F       *hREC_1D;
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    TH1F       *hREC_1D_in_Rap, *hREC_1D_in_Rap_All, *hREF_1D_in_Rap;
    
    // 2D
    TH2F       *hREC_2D;
    TH1F      **hREC_1D_in_PT_2D_bin        = new TH1F     *[nBinPT2D];
    TH2F      **hREC_1D_in_Rap_2D_Bin       = new TH2F     *[nBinRap_];
    TH2F      **hREC_2D_in_Rap              = new TH2F     *[nBinRap_];
    TH2F     ***hREC_2D_in_PT               = new TH2F    **[nBinPT2D];

    hName = "hREC_1D_in_Rap";
    hTitle= "Rapidity difference for #phi meson candidates";
    hREC_1D_in_Rap  =   new TH1F (hName,hTitle,100,-1.,1.);
    
    hName = "hREC_1D_in_Rap_All";
    hTitle= "Rapidity difference for #phi meson candidates";
    hREC_1D_in_Rap_All  =   new TH1F (hName,hTitle,100,-1.,1.);
    
    hName = "hREF_1D_in_Rap";
    hTitle= "Rapidity distribution for #phi meson candidates";
    hREF_1D_in_Rap  =   new TH1F (hName,hTitle,100,-.5,.5);
    
    hName = "hREC_1D";
    hTitle= "--";
    hREC_1D  =   new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
    
    hName = "hREC_2D";
    hTitle= "--";
    hREC_2D  =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
    
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        hName = Form("hREC_1D_in_PT_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT1D[iHisto],fArrPT1D[iHisto+1]);
        hREC_1D_in_PT[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT[iHisto],"IM 1D");
    }
    
    for ( Int_t iHisto = 0; iHisto < nBinRap_; iHisto++ )
    {
        hName = Form("hREC_2D_in_Rap_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in |#Delta y| range [%.1f-%.1f]",fArrRap_[iHisto],fArrRap_[iHisto+1]);
        hREC_2D_in_Rap[iHisto]   = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
        SetAxis(hREC_2D_in_Rap[iHisto],"IM 2D");
    }
    
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        hName = Form("hREC_1D_in_PT_2D_bin_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        hREC_1D_in_PT_2D_bin[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT_2D_bin[iHisto],"IM 1D");
        
        hREC_2D_in_PT[iHisto]    = new TH2F *    [nBinPT2D];
        
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hREC_2D_in_PT[iHisto][jHisto]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            SetAxis(hREC_2D_in_PT[iHisto][jHisto],"IM 2D");
        }
    }

    // >>->> CORRELATION ANALYSIS //

    // >>->> TRIGGER ANALYSIS //
    hName   =   "";
    hTitle  =   "";
    TH1F   *hTriggerEvt =   new TH1F (hName,hTitle,5,-0.5,4.5);
    
    // Creating the Number of events Result Histogram------------------------------------------------------------------
    hName                       = "Entry_DT";
    hTitle                      = "Events in DT";
    TH1F *          hUtlEntry   = new TH1F (hName,hTitle,1,0.5,1.5);
    hUtlEntry                   ->GetXaxis()->SetTitle("");
    hUtlEntry                   ->GetYaxis()->SetTitle("Events");
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Analysis");
    
    // Evaluating entries and saving them for later
    Int_t nEvents = TPhiCandidate->GetEntries();
    hUtlEntry     ->SetBinContent(1,nEvents);
    Int_t nOverflow = 0;

    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents*0+3; iEvent++ )
    {
        gROOT->SetBatch(true);
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        gROOT->SetBatch(false);
        
        if ( iEvent%1000000 == 0 && iEvent != 0) fPrintLoopTimer("Analysis",iEvent,nEvents);
        
        // Skipping overflow
        if ( (int)(evPhiCandidate.nPhi) >= 153 )
        {
            cout << "[INFO] Skipping overflow event" << endl;
            nOverflow++;
            continue;
        }
        // Skipping non candidate events
        if ( (int)(evPhiCandidate.nPhi) == 0 )
        {
            hTriggerEvt->Fill(0);
            continue;
        } 

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        Bool_t  fCheckFill1 =   false;
        Bool_t  fCheckFill2 =   false;
        Bool_t  fCheckFill3 =   false;
        Bool_t  fCheckFill4 =   false;

        for ( Int_t iPhi = 0; iPhi < (int)(evPhiCandidate.nPhi); iPhi++ )
        {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fAcceptCandidate(LPhi_candidate1.Rapidity(),LPhi_candidate1.Pt(),evPhiCandidate.Multiplicity) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )
        {
            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[U_AccCand[iPhi]],evPhiCandidate.Py[U_AccCand[iPhi]],evPhiCandidate.Pz[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[iPhi]]);

            // >>->> 1-Dimensional Analysis Fill
            // Trigger
            if ( !fCheckFill1 )
            {
                hTriggerEvt->Fill(1);
                fCheckFill1 = true;
            }

            // Full Yield
            hREC_1D                                                     ->  Fill(evPhiCandidate.InvMass[iPhi]);

            // pT-Differential Yield
            hREC_1D_in_PT[fGetBinPT1D(LPhi_candidate1.Pt())]            ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]]);
            hREC_1D_in_PT_2D_bin[fGetBinPT2D(LPhi_candidate1.Pt())]     ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]]);

            // Correlation
            hREF_1D_in_Rap                                              ->  Fill(LPhi_candidate1.Rapidity());

            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )
            {
                // Building Second Candidate
                LPhi_candidate2.SetXYZM(evPhiCandidate.Px[U_AccCand[jPhi]],evPhiCandidate.Py[U_AccCand[jPhi]],evPhiCandidate.Pz[U_AccCand[jPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]]);

                // Selecting valid candidates
                if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi) ) continue;
                
                // >>->> 2-Dimensional Analysis Fill
                // Trigger
                if ( !fCheckFill2 )
                {
                    hTriggerEvt->Fill(2);
                    fCheckFill2 = true;
                }

                // Full Yield
                hREC_2D_in_PT[fGetBinPT2D(LPhi_candidate1.Pt())][fGetBinPT2D(LPhi_candidate2.Pt())] ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]],0.5);
                hREC_2D                                                                             ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]],0.5);
                
                hREC_1D_in_Rap->    Fill(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity(),0.5);
                hREC_1D_in_Rap->    Fill( fabs(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity()) ,0.5);
                hREC_1D_in_Rap_All->Fill( (LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity()) ,0.5);
                 /*
                if ( fGetBinRap_(fabs(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity())) == -1 )
                {
                    cout << endl;
                    cout << fGetBinRap_(fabs(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity())) << endl;
                    cout << "Rap1:" << LPhi_candidate1.Rapidity() << endl;
                    cout << "Rap2:" << LPhi_candidate2.Rapidity() << endl;
                    cout << "DRap:" << fabs(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity()) << endl;
                    cout << endl;
                    continue;
                }
                hREC_2D_in_Rap[fGetBinRap_(fabs(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity()))]->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]],0.5);
                */
                for ( Int_t kPhi = 0; kPhi < U_nAccept; kPhi++ )
                {
                    // Selecting valid candidates
                    if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi, kPhi) ) continue;
                    
                    // >>->> 3-Dimensional Analysis Fill
                    // Trigger
                    if ( !fCheckFill3 )
                    {
                        hTriggerEvt->Fill(3);
                        fCheckFill3 = true;
                    }

                    for ( Int_t lPhi = 0; lPhi < U_nAccept; lPhi++ )
                    {
                        // Selecting valid candidates
                        if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi, kPhi, lPhi) ) continue;

                        // >>->> 4-Dimensional Analysis Fill
                        // Trigger
                        if ( !fCheckFill4 )
                        {
                            hTriggerEvt->Fill(4);
                            fCheckFill4 = true;
                        }
                    }
                }
            }
        }
    }

    fStopTimer("Analysis");
    cout << "[INFO] The overflow events were: " << nOverflow << endl;
    cout << "[INFO] The overflow events were the " << 100*(nOverflow*1.)/(1.*nEvents) << "% of total events" << endl;

    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    
    TFile *outFile  =   new TFile   (fInvMasHist,"recreate");
    
    hUtlEntry->Write();
    hREF_1D_in_Rap->Write();
    hREC_1D_in_Rap->Write();
    hREC_1D_in_Rap_All->Write();
    hREC_1D->Write();
    hREC_2D->Write();
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hREC_1D_in_PT[iHisto]   ->Write();
    }
    for (int iHisto = 0; iHisto < nBinRap_; iHisto++)
    {
        hREC_2D_in_Rap[iHisto]   ->Write();
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hREC_1D_in_PT_2D_bin[iHisto]   ->Write();
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hREC_2D_in_PT[iHisto][jHisto]->Write();
        }
    }
    
    outFile->Close();
    insFileDT->Close();
}
