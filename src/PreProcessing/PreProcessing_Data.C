#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void PreProcessing_Data ( string fFileName = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    if ( fFileName == "" )
    {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing_Data(\"Root_file_name.root\")" << endl;
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
        TKaonCandidate-> SetBranchAddress   ("Multiplicity",   &evKaonCandidate.Multiplicity);
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
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",   &evPhiCandidate.Multiplicity);
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
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",   &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
        
        TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
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
    TH1F       *hREC_1D;
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   Form("hREC_1D");
    hTitle      =   Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fMinPT1D,fMaxPT1D);
    hREC_1D     =   new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
    //
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        hName = Form("hREC_1D_in_PT_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT1D[iHisto],fArrPT1D[iHisto+1]);
        hREC_1D_in_PT[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT[iHisto],"IM 1D");
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
    hName       =   Form("hREC_2D");
    hTitle      =   Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D);
    hREC_2D     =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
    //
    //  Defining pT-Differential histograms over measurable pT
    //
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

    // >-> MULTIPLICITY ANALYSIS //

    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F     ***hREC_1D_in_MT_in_PT         = new TH1F    **[nBinMult];
    //
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fMinPT1D,fMaxPT1D,fArrMult[iHisto],fArrMult[iHisto+1]);
        hREC_1D_in_MT[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_MT[iHisto],"IM 1D");

        hREC_1D_in_MT_in_PT[iHisto] = new TH1F     *[nBinPT1D];

        for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
        {
            hName = Form("hREC_1D_in_MT_PT_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fArrPT1D[jHisto],fArrPT1D[jHisto+1],fArrMult[iHisto],fArrMult[iHisto+1]);
            hREC_1D_in_MT_in_PT[iHisto][jHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_MT_in_PT[iHisto][jHisto],"IM 1D");
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
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName   =   Form("hREC_2D_in_MT_%i",iHisto);
        hTitle  =   Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,fArrMult[iHisto],fArrMult[iHisto+1]);
        hREC_2D_in_MT[iHisto]   =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);

        hREC_1D_in_MT_in_PT_2D_bin[iHisto]  = new TH1F     *[nBinPT2D];
        hREC_2D_in_MT_in_PT[iHisto]         = new TH2F    **[nBinPT2D];
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrMult[iHisto],fArrMult[iHisto+1]);
            hREC_1D_in_MT_in_PT_2D_bin[iHisto][jHisto]  = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_MT_in_PT_2D_bin[iHisto][jHisto],"IM 1D");
            
            hREC_2D_in_MT_in_PT[iHisto][jHisto]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_%i",iHisto,jHisto,kHisto);
                hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1],fArrMult[iHisto],fArrMult[iHisto+1]);
                hREC_2D_in_MT_in_PT[iHisto][jHisto][kHisto]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
                SetAxis(hREC_2D_in_MT_in_PT[iHisto][jHisto][kHisto],"IM 2D");
            }
        }
    }

    // >-> TRIGGER ANALYSIS //
    //
    //  Overall Triggering in Events
    //
    hName   =   "hTrigger";
    hTitle  =   "Triggered events for NTuples candidates";
    TH1D   *hTriggerEvt =   new TH1D (hName,hTitle,5,-0.5,4.5);
    hTriggerEvt->GetYaxis()->SetTitle("% of events accepted");
    //
    //  Triggering in Events in pT
    //
    hName   =   "hTrigger1D";
    hTitle  =   "Triggered events for Single Phis in PT";
    TH1D   *hTriggerEvt1D =   new TH1D (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hTriggerEvt1D,"PT 1D");
    hTriggerEvt1D->GetYaxis()->SetTitle("% of events accepted");
    //
    //  Triggering in Events in pT
    //
    hName   =   "hTrigger2D";
    hTitle  =   "Triggered events for Double Phis in PT";
    TH2D   *hTriggerEvt2D =   new TH2D (hName,hTitle,nBinPT1D,fArrPT1D,nBinPT1D,fArrPT1D);
    SetAxis(hTriggerEvt2D,"PT 2D");
    hTriggerEvt2D->GetZaxis()->SetTitle("% of events accepted");

    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Analysis");
    
    // Evaluating entries and saving them for later
    Int_t nEvents = TPhiCandidate->GetEntries();

    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Analysis",iEvent,nEvents,1000000);

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

        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; iPhi++ )  {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fAcceptCandidate(LPhi_candidate1.Rapidity(),evPhiCandidate.InvMass[iPhi],LPhi_candidate1.Pt(),evPhiCandidate.Multiplicity) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        
        if ( U_nAccept == 0 )   hTriggerEvt->Fill(0);
        
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;
        
            // Selecting valid candidates
            if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi) ) continue;

            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[U_AccCand[iPhi]],evPhiCandidate.Py[U_AccCand[iPhi]],evPhiCandidate.Pz[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[iPhi]]);

            // >-> 1-Dimensional Analysis Fill   //
            //
            // >->-->-> Utilities
            //
            Int_t   indexPT1D = fGetBinPT1D(LPhi_candidate1.Pt());
            Int_t   indexPT2D = fGetBinPT2D(LPhi_candidate1.Pt());
            Int_t   indexMult = fGetBinMult(evPhiCandidate.Multiplicity);
            Float_t iInvMass_ = evPhiCandidate.InvMass[U_AccCand[iPhi]];
            //
            // >->-->-> Trigger
            //
            if ( !fCheckFill1 ) {
                hTriggerEvt->Fill(1);
                fCheckFill1 = true;
            }
            hTriggerEvt1D                                       ->Fill(LPhi_candidate1.Pt());
            // 
            // >->-->-> Yield
            //
            hREC_1D                                             ->  Fill(iInvMass_);
            hREC_1D_in_PT[indexPT1D]                            ->  Fill(iInvMass_);
            hREC_1D_in_PT_2D_bin[indexPT2D]                     ->  Fill(iInvMass_);
            //
            // >->-->-> Multiplicity
            //
            hREC_1D_in_MT[indexMult]                            ->  Fill(iInvMass_);
            hREC_1D_in_MT_in_PT[indexMult][indexPT1D]           ->  Fill(iInvMass_);
            hREC_1D_in_MT_in_PT_2D_bin[indexMult][indexPT2D]    ->  Fill(iInvMass_);
            
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )    {
                // Must have at least 2 candidates
                if ( U_nAccept < 2 ) break;

                // Selecting valid candidates
                if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi) ) continue;
                
                // Building Second Candidate
                LPhi_candidate2.SetXYZM(evPhiCandidate.Px[U_AccCand[jPhi]],evPhiCandidate.Py[U_AccCand[jPhi]],evPhiCandidate.Pz[U_AccCand[jPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]]);

                // >-> 2-Dimensional Analysis Fill
                //
                // >->-->-> Utilities
                //
                Int_t   jndexPT1D = fGetBinPT1D(LPhi_candidate2.Pt());
                Int_t   jndexPT2D = fGetBinPT2D(LPhi_candidate2.Pt());
                Float_t jInvMass_ = evPhiCandidate.InvMass[U_AccCand[jPhi]];
                //
                // >->-->-> Trigger
                //
                if ( !fCheckFill2 ) {
                    hTriggerEvt->Fill(2);
                    fCheckFill2 = true;
                }
                hTriggerEvt2D                                           ->Fill(LPhi_candidate1.Pt(),LPhi_candidate2.Pt(),0.5);
                // 
                // >->-->-> Yield
                //
                hREC_2D                                                 ->  Fill(iInvMass_,jInvMass_,0.5);
                hREC_2D_in_PT[indexPT2D][jndexPT2D]                     ->  Fill(iInvMass_,jInvMass_,0.5);
                //
                // >->-->-> Multiplicity
                //
                hREC_2D_in_MT[indexMult]                                ->  Fill(iInvMass_,jInvMass_,0.5);
                hREC_2D_in_MT_in_PT[indexMult][indexPT2D][jndexPT2D]    ->  Fill(iInvMass_,jInvMass_,0.5);
                
                for ( Int_t kPhi = 0; kPhi < U_nAccept; kPhi++ )    {
                    // Must have at least 3 candidates
                    if ( U_nAccept < 3 ) break;

                    // Selecting valid candidates
                    if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi, kPhi) ) continue;
                    
                    // >->->-> 3-Dimensional Analysis Fill
                    // Trigger
                    if ( !fCheckFill3 ) {
                        hTriggerEvt->Fill(3);
                        fCheckFill3 = true;
                    }

                    for ( Int_t lPhi = 0; lPhi < U_nAccept; lPhi++ )    {
                        // Must have at least 4 candidates
                        if ( U_nAccept < 4 ) break;

                        // Selecting valid candidates
                        if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi, kPhi, lPhi) ) continue;

                        // >->->-> 4-Dimensional Analysis Fill
                        // Trigger
                        if ( !fCheckFill4 ) {
                            hTriggerEvt->Fill(4);
                            fCheckFill4 = true;
                        }
                    }
                }
            }
        }
    }

    fStopTimer("Analysis");

    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    
    // >-> Trigger Analysis
    //
    // >->-> Bin Normalisation
    hTriggerEvt1D->Scale(1.,"width");
    hTriggerEvt2D->Scale(1.,"width");
    //
    // >->-> Event Normalisation
    //
    hTriggerEvt->Scale(100./nEvents);
    hTriggerEvt1D->Scale(100./nEvents);
    hTriggerEvt2D->Scale(100./nEvents);

    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >-> Trigger Analysis
    //
    TFile *outFil1  =   new TFile   (fTrgPreProc,"recreate");
    //
    hTriggerEvt->Write();
    hTriggerEvt1D->Write();
    hTriggerEvt2D->Write();
    //
    outFil1->Close();
    //
    // >-> Yield Analysis
    //
    TFile *outFil2  =   new TFile   (fYldPreProc,"recreate");
    //
    hREC_1D->Write();
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hREC_1D_in_PT[iHisto]   ->Write();
    }
    //
    hREC_2D->Write();
    //
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hREC_1D_in_PT_2D_bin[iHisto]   ->Write();
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hREC_2D_in_PT[iHisto][jHisto]->Write();
        }
    }
    //
    outFil2->Close();
    //
    // >-> Multiplicity Analysis
    //
    TFile *outFil3  =   new TFile   (fMltPreProc,"recreate");
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hREC_1D_in_MT[iHisto]->Write();
        hREC_2D_in_MT[iHisto]->Write();

        for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
        {
            hREC_1D_in_MT_in_PT[iHisto][jHisto]->Write();
        }
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hREC_1D_in_MT_in_PT_2D_bin[iHisto][jHisto]->Write();
            
            for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
            {
                hREC_2D_in_MT_in_PT[iHisto][jHisto][kHisto]->Write();
            }
        }
    }
    //
    outFil3->Close();
    //
    // >-> Close input File
    //
    insFileDT->Close();
    //
}
