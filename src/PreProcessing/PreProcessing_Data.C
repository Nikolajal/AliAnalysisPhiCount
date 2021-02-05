#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void PreProcessing_Data ( string fFileName = "", Int_t nEventsCut = -1, string fOption = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    //
    // >-> Initialisation warnings
    //
    if ( fFileName == "" )  {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing_Data.C(\"Root_file_name.root\")" << endl;
        return;
    }
    if ( nEventsCut != -1 ) cout << "[WARNING] Choosing to limit the datasample to " << nEventsCut << " events" <<endl;
    
    // Retrieving Event data
    TFile *insFileDT        =   new TFile   (fFileName.c_str());
    
    // Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileDT       ->Get(fPhiCandidate_Tree);
    TTree   *TKaonCandidate =   (TTree*)insFileDT       ->Get(fKaonCandidate_Tree);
    
    // Retrieving Event Count Histogram
    TList  *fQCOutputList   =   (TList*)insFileDT       ->Get("fQCOutputList");
    TH1D   *fHEventCount    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enumerate");
    
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
        TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
        TKaonCandidate-> SetBranchAddress   ("nPhi",            &evKaonCandidate.nKaon);
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
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
        TPhiCandidate-> SetBranchAddress    ("Nature",          &evPhiCandidate.Nature);
    }
    else
    {
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
        TPhiCandidate-> SetBranchAddress    ("Nature",          &evPhiCandidate.Nature);
        
        TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
        TKaonCandidate-> SetBranchAddress   ("nPhi",           &evKaonCandidate.nKaon);
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
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
    {
        hName = Form("hREC_1D_in_PT_%i",iPT1D);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]);
        hREC_1D_in_PT[iPT1D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT[iPT1D],"IM 1D");
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
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )
    {
        hName = Form("hREC_1D_in_PT_2D_bin_%i",iPT2D);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]);
        hREC_1D_in_PT_2D_bin[iPT2D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT_2D_bin[iPT2D],"IM 1D");
        
        hREC_2D_in_PT[iPT2D]    = new TH2F *    [nBinPT2D];
        
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i",iPT2D,jPT2D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrPT2D[jPT2D],fArrPT2D[jPT2D+1]);
            hREC_2D_in_PT[iPT2D][jPT2D]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            SetAxis(hREC_2D_in_PT[iPT2D][jPT2D],"IM 2D");
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
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iMult);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fMinPT1D,fMaxPT1D,fArrMult[iMult],fArrMult[iMult+1]);
        hREC_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_MT[iMult],"IM 1D");

        hREC_1D_in_MT_in_PT[iMult] = new TH1F     *[nBinPT1D];

        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
        {
            hName = Form("hREC_1D_in_MT_PT_%i_%i",iMult,iPT1D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1],fArrMult[iMult],fArrMult[iMult+1]);
            hREC_1D_in_MT_in_PT[iMult][iPT1D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_MT_in_PT[iMult][iPT1D],"IM 1D");
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
    for ( Int_t iMult = 0; iMult < nBinMult; iMult++ )
    {
        hName   =   Form("hREC_2D_in_MT_%i",iMult);
        hTitle  =   Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,fArrMult[iMult],fArrMult[iMult+1]);
        hREC_2D_in_MT[iMult]   =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);

        hREC_1D_in_MT_in_PT_2D_bin[iMult]  = new TH1F     *[nBinPT2D];
        hREC_2D_in_MT_in_PT[iMult]         = new TH2F    **[nBinPT2D];
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_%i",iMult,iPT2D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrMult[iMult],fArrMult[iMult+1]);
            hREC_1D_in_MT_in_PT_2D_bin[iMult][iPT2D]  = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_MT_in_PT_2D_bin[iMult][iPT2D],"IM 1D");
            
            hREC_2D_in_MT_in_PT[iMult][iPT2D]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_%i",iMult,iPT2D,jPT2D);
                hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c, Multiplicity [%.1f-%.1f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrPT2D[jPT2D],fArrPT2D[jPT2D+1],fArrMult[iMult],fArrMult[iMult+1]);
                hREC_2D_in_MT_in_PT[iMult][iPT2D][jPT2D]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
                SetAxis(hREC_2D_in_MT_in_PT[iMult][iPT2D][jPT2D],"IM 2D");
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
    
    // Evaluating entries and saving them for later
    Int_t nEvents = (!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : nEventsCut);

    if ( nEvents > 0 )  fStartTimer("Phi Analysis");
    
    TH1F * check = new TH1F ("TEST","TEST",nBinPT1D,fArrPT1D);
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Phi Analysis",iEvent,nEvents,1000000);

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        Bool_t  fCheckFill1 =   false;
        Bool_t  fCheckFill2 =   false;
        Bool_t  fCheckFill3 =   false;
        Bool_t  fCheckFill4 =   false;
        
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; iPhi++ )  {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            evPhiCandidate.pT[iPhi]     =   LPhi_candidate1.Pt();
            evPhiCandidate.Rap[iPhi]    =   LPhi_candidate1.Rapidity();
            if ( !fAcceptCandidate(evPhiCandidate.Rap[iPhi],evPhiCandidate.InvMass[iPhi],evPhiCandidate.pT[iPhi],evPhiCandidate.Multiplicity) ) continue;
            if ( evPhiCandidate.Nature[iPhi] == 1 ) check->Fill(LPhi_candidate1.Pt());
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        
        if ( U_nAccept == 0 )   hTriggerEvt->Fill(0);
        
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;
        
            // Selecting valid candidates
            if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi) ) continue;
            
            // >-> 1-Dimensional Analysis Fill   //
            //
            // >->-->-> Utilities
            //
            Int_t   iPT1D       = fGetBinPT1D(evPhiCandidate.pT[U_AccCand[iPhi]]);
            Int_t   iPT2D       = fGetBinPT2D(evPhiCandidate.pT[U_AccCand[iPhi]]);
            Int_t   iMult       = fGetBinMult(evPhiCandidate.Multiplicity);
            Float_t iInvMass_   = evPhiCandidate.InvMass[U_AccCand[iPhi]];
            //
            // >->-->-> Trigger
            //
            if ( !fCheckFill1 ) {
                hTriggerEvt->Fill(1);
                fCheckFill1 = true;
            }
            hTriggerEvt1D                                   ->  Fill(LPhi_candidate1.Pt());
            // 
            // >->-->-> Yield
            //
            hREC_1D                                         ->  Fill(iInvMass_);
            hREC_1D_in_PT[iPT1D]                            ->  Fill(iInvMass_);
            hREC_1D_in_PT_2D_bin[iPT2D]                     ->  Fill(iInvMass_);
            //
            // >->-->-> Multiplicity
            //
            hREC_1D_in_MT[iMult]                            ->  Fill(iInvMass_);
            hREC_1D_in_MT_in_PT[iMult][iPT1D]               ->  Fill(iInvMass_);
            hREC_1D_in_MT_in_PT_2D_bin[iMult][iPT2D]        ->  Fill(iInvMass_);
            
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )    {
                // Must have at least 2 candidates
                if ( U_nAccept < 2 ) break;

                // Selecting valid candidates
                if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi) ) continue;
                
                // >-> 2-Dimensional Analysis Fill
                //
                // >->-->-> Utilities
                //
                Int_t   jPT1D       = fGetBinPT1D(evPhiCandidate.pT[U_AccCand[jPhi]]);
                Int_t   jPT2D       = fGetBinPT2D(evPhiCandidate.pT[U_AccCand[jPhi]]);
                Float_t jInvMass_   = evPhiCandidate.InvMass[U_AccCand[jPhi]];
                //
                // >->-->-> Trigger
                //
                if ( !fCheckFill2 ) {
                    hTriggerEvt->Fill(2);
                    fCheckFill2 = true;
                }
                hTriggerEvt2D                                           ->  Fill(LPhi_candidate1.Pt(),LPhi_candidate2.Pt(),0.5);
                // 
                // >->-->-> Yield
                //
                hREC_2D                                                 ->  Fill(iInvMass_,jInvMass_,0.5);
                hREC_2D_in_PT[iPT2D][jPT2D]                             ->  Fill(iInvMass_,jInvMass_,0.5);
                //
                // >->-->-> Multiplicity
                //
                hREC_2D_in_MT[iMult]                                    ->  Fill(iInvMass_,jInvMass_,0.5);
                hREC_2D_in_MT_in_PT[iMult][iPT2D][jPT2D]                ->  Fill(iInvMass_,jInvMass_,0.5);
                
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
    
    if ( nEvents > 0 )  fStopTimer("Phi Analysis");
    
    // Evaluating entries and saving them for later
    nEvents = 0;//(!TKaonCandidate) ? 0 : ( nEventsCut == -1.? TKaonCandidate->GetEntries() : nEventsCut);
    /*
    if ( nEvents > 0 )  fStartTimer("Kaon Analysis");

    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TKaonCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Kaon Analysis",iEvent,nEvents,1000000);
    }
    
    if ( nEvents > 0 )  fStopTimer("Kaon Analysis");
    */
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    //
    fStartTimer("post-processing output objects");
    //
    // >-> Event count recovery
    //
    nEvents =   fHEventCount->GetBinContent(9);
    //for ( int i = 0; i <= fHEventCount->GetBinContent(10); i++ ) { hTriggerEvt ->  Fill(0); }
    //
    // >-> Trigger Analysis
    //
    // >->-> Bin Normalisation
    hTriggerEvt1D   ->Scale(1.,"width");
    hTriggerEvt2D   ->Scale(1.,"width");
    //
    // >->-> Event Normalisation
    //
    hTriggerEvt     ->Scale(100./fHEventCount->GetBinContent(1));
    hTriggerEvt1D   ->Scale(100./nEvents);
    hTriggerEvt2D   ->Scale(100./nEvents);
    //
    fStopTimer("post-processing output objects");
    //
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    fStartTimer("writing output objects to file(s)");
    //
    // >-> Trigger Analysis
    //
    TFile *outFil1  =   new TFile   (fTrgPreProc,"recreate");
    //
    fHEventCount    ->Write();
    hTriggerEvt     ->Write();
    hTriggerEvt1D   ->Write();
    hTriggerEvt2D   ->Write();
    //
    outFil1->Close();
    //
    // >-> Yield Analysis
    //
    TFile *outFil2  =   new TFile   (fYldPreProc,"recreate");
    //
    check->Write();
    fHEventCount    ->Write();
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
    fHEventCount->Write();
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
    fStopTimer("writing output objects to file(s)");
    // >-> Close input File
    //
    insFileDT->Close();
    //
}
