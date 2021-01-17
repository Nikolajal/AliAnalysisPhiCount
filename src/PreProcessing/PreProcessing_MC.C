#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: All Set!

void PreProcessing_MC ( string fFileName = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    // >-> OPTIONS
    
    if ( fFileName == "" )
    {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing_MC(\"Root_file_name.root\")" << endl;
        return;
    }
    
    //Retrieving Event data
    TFile *insFileMC        =   new TFile   (fFileName.c_str());
    
    //Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileMC->Get(fPhiCandidateEff_Tree);
    TTree   *TKaonCandidate =   (TTree*)insFileMC->Get(fKaonCandidateEff_Tree);
    TTree   *TPhi_Multref   =   (TTree*)insFileMC->Get(fPhiCandidate_Tree);
    TTree   *TKaon_Multref  =   (TTree*)insFileMC->Get(fKaonCandidate_Tree);
    
    // Define tree data structures
    Struct_PhiEfficiency    evPhiEfficiency;
    Struct_KaonEfficiency   evKaonEfficiency;

    if ( !TPhiCandidate && !TKaonCandidate )
    {
        cout << "Input Data Tree not found!" << endl;
        return;
    }
    if ( !TPhiCandidate )
    {
        cout << "[INFO] No PhiCandidate Tree, switching to Kaon Analysis" << endl;
        TKaonCandidate->SetBranchAddress    ("nKaon",           &evKaonEfficiency.nKaon);
        TKaonCandidate->SetBranchAddress    ("Px",              &evKaonEfficiency.Px);
        TKaonCandidate->SetBranchAddress    ("Py",              &evKaonEfficiency.Py);
        TKaonCandidate->SetBranchAddress    ("Pz",              &evKaonEfficiency.Pz);
        TKaonCandidate->SetBranchAddress    ("Selection",       &evKaonEfficiency.Selection);
        TKaon_Multref ->SetBranchAddress    ("Multiplicity",    &evKaonEfficiency.Multiplicity);
    }
    else if ( !TKaonCandidate )
    {
        cout << "[INFO] No KaonCandidate Tree, switching to Phi Analysis" << endl;
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiEfficiency.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiEfficiency.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiEfficiency.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiEfficiency.Pz);
        TPhiCandidate-> SetBranchAddress    ("Selection",       &evPhiEfficiency.Selection);
        TPhi_Multref -> SetBranchAddress    ("Multiplicity",    &evPhiEfficiency.Multiplicity);
    }
    else
    {
        TKaonCandidate->SetBranchAddress    ("nKaon",           &evKaonEfficiency.nKaon);
        TKaonCandidate->SetBranchAddress    ("Px",              &evKaonEfficiency.Px);
        TKaonCandidate->SetBranchAddress    ("Py",              &evKaonEfficiency.Py);
        TKaonCandidate->SetBranchAddress    ("Pz",              &evKaonEfficiency.Pz);
        TKaonCandidate->SetBranchAddress    ("Selection",       &evKaonEfficiency.Selection);
        TKaon_Multref ->SetBranchAddress    ("Multiplicity",    &evKaonEfficiency.Multiplicity);
        
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiEfficiency.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiEfficiency.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiEfficiency.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiEfficiency.Pz);
        TPhiCandidate-> SetBranchAddress    ("Selection",       &evPhiEfficiency.Selection);
        TPhi_Multref -> SetBranchAddress    ("Multiplicity",    &evPhiEfficiency.Multiplicity);
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

    // >> YIELD ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hREC_1D;
    TH1F       *hGEN_1D;
    TH1F       *hTRU_1D;
    TH1F       *hREC_1D_in_2Dbin;
    TH1F       *hGEN_1D_in_2Dbin;
    TH1F       *hTRU_1D_in_2Dbin;
    TH1F       *hEFF_1D;
    TH1F       *hEFF_1D_in_2Dbin;
    //
    //  Defining Efficiency and check utilities
    //
    hName       =   Form("hREC_1D");
    hTitle      =   Form("hREC_1D");
    hREC_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hREC_1D,"PT 1D");
    //
    hName       =   Form("hGEN_1D");
    hTitle      =   Form("hGEN_1D");
    hGEN_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hGEN_1D,"PT 1D");
    //
    hName       =   Form("hTRU_1D");
    hTitle      =   Form("hTRU_1D");
    hTRU_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hTRU_1D,"PT 1D");
    //
    hName       =   Form("hEFF_1D");
    hTitle      =   Form("hEFF_1D");
    hEFF_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hEFF_1D,"PT 1D");
    //
    hName       =   Form("hREC_1D_in_2Dbin");
    hTitle      =   Form("hREC_1D_in_2Dbin");
    hREC_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hREC_1D_in_2Dbin,"PT 1D");
    //
    hName       =   Form("hGEN_1D_in_2Dbin");
    hTitle      =   Form("hGEN_1D_in_2Dbin");
    hGEN_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_1D_in_2Dbin,"PT 1D");
    //
    hName       =   Form("hTRU_1D_in_2Dbin");
    hTitle      =   Form("hTRU_1D_in_2Dbin");
    hTRU_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_1D_in_2Dbin,"PT 1D");
    //
    hName       =   Form("hEFF_1D_in_2Dbin");
    hTitle      =   Form("hEFF_1D_in_2Dbin");
    hEFF_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_1D_in_2Dbin,"PT 1D");
    //
    // >>-->> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hREC_2D;
    TH2F       *hGEN_2D;
    TH2F       *hTRU_2D;
    TH2F       *hEFF_2D;
    TH2F       *hEFF_2D_fr_1D;
    //
    //  Defining Efficiency and check utilities
    //
    hName       =   Form("hREC_2D");
    hTitle      =   Form("hREC_2D");
    hREC_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hREC_2D,"PT 2D");
    //
    hName       =   Form("hGEN_2D");
    hTitle      =   Form("hGEN_2D");
    hGEN_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_2D,"PT 1D");
    //
    hName       =   Form("hTRU_2D");
    hTitle      =   Form("hTRU_2D");
    hTRU_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_2D,"PT 2D");
    //
    hName       =   Form("hEFF_2D");
    hTitle      =   Form("hEFF_2D");
    hEFF_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_2D,"PT 2D");
    //
    hName       =   Form("hEFF_2D_fr_1D");
    hTitle      =   Form("hEFF_2D_fr_1D");
    hEFF_2D_fr_1D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_2D_fr_1D,"PT 2D");
    //

    // >> MULTIPLICITY ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F      **hGEN_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F      **hTRU_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F      **hEFF_1D_in_MT               = new TH1F     *[nBinMult];
    TH1F      **hREC_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult];
    TH1F      **hGEN_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult];
    TH1F      **hTRU_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult];
    TH1F      **hEFF_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult];
    //
    //  Defining MT-Differential histograms
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iHisto);
        hTitle= Form("hREC_1D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hREC_1D_in_MT[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hREC_1D_in_MT[iHisto],"PT 1D");
        
        hName = Form("hGEN_1D_in_MT_%i",iHisto);
        hTitle= Form("hGEN_1D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hGEN_1D_in_MT[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hGEN_1D_in_MT[iHisto],"PT 1D");
        
        hName = Form("hTRU_1D_in_MT_%i",iHisto);
        hTitle= Form("hTRU_1D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hTRU_1D_in_MT[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hTRU_1D_in_MT[iHisto],"PT 1D");
        
        hName = Form("hEFF_1D_in_MT_%i",iHisto);
        hTitle= Form("hEFF_1D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hEFF_1D_in_MT[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hEFF_1D_in_MT[iHisto],"PT 1D");
        
        hName = Form("hREC_1D_in_MT_in_2Dbin_%i",iHisto);
        hTitle= Form("hREC_1D_in_MT_in_2Dbin Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hREC_1D_in_MT_in_2Dbin[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hREC_1D_in_MT_in_2Dbin[iHisto],"PT 1D");
            
        hName = Form("hGEN_1D_in_MT_in_2Dbin_%i",iHisto);
        hTitle= Form("hGEN_1D_in_MT_in_2Dbin Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hGEN_1D_in_MT_in_2Dbin[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hGEN_1D_in_MT_in_2Dbin[iHisto],"PT 1D");
        
        hName = Form("hTRU_1D_in_MT_in_2Dbin_%i",iHisto);
        hTitle= Form("hTRU_1D_in_MT_in_2Dbin Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hTRU_1D_in_MT_in_2Dbin[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hTRU_1D_in_MT_in_2Dbin[iHisto],"PT 1D");
        
        hName = Form("hEFF_1D_in_MT_in_2Dbin_%i",iHisto);
        hTitle= Form("hEFF_1D_in_MT_in_2Dbin Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hEFF_1D_in_MT_in_2Dbin[iHisto]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hEFF_1D_in_MT_in_2Dbin[iHisto],"PT 1D");
    }
    //
    // >>-->> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hREC_2D_in_MT               = new TH2F     *[nBinMult];
    TH2F      **hGEN_2D_in_MT               = new TH2F     *[nBinMult];
    TH2F      **hTRU_2D_in_MT               = new TH2F     *[nBinMult];
    TH2F      **hEFF_2D_in_MT               = new TH2F     *[nBinMult];
    TH2F      **hEFF_2D_in_MT_fr_1D         = new TH2F     *[nBinMult];
    //
    //  Defining MT-Differential histograms
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hName = Form("hREC_2D_in_MT_%i",iHisto);
        hTitle= Form("hREC_2D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hREC_2D_in_MT[iHisto]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hREC_2D_in_MT[iHisto],"PT 2D");
        
        hName = Form("hGEN_2D_in_MT_%i",iHisto);
        hTitle= Form("hGEN_2D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hGEN_2D_in_MT[iHisto]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_2D_in_MT[iHisto],"PT 2D");
        
        hName = Form("hTRU_2D_in_MT_%i",iHisto);
        hTitle= Form("hTRU_2D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hTRU_2D_in_MT[iHisto]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_2D_in_MT[iHisto],"PT 2D");
        
        hName = Form("hEFF_2D_in_MT_%i",iHisto);
        hTitle= Form("hEFF_2D_in_MT Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hEFF_2D_in_MT[iHisto]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_2D_in_MT[iHisto],"PT 2D");
        
        hName = Form("hEFF_2D_in_MT_fr_1D_%i",iHisto);
        hTitle= Form("hEFF_2D_in_MT_fr_1D Multiplicity [%.1f-%.1f]",fArrMult[iHisto],fArrMult[iHisto+1]);
        hEFF_2D_in_MT_fr_1D[iHisto]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_2D_in_MT_fr_1D[iHisto],"PT 2D");
    }
    //

    // >> TRIGGER ANALYSIS //

    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Analysis");
    
    // Evaluating entries
    Int_t nEvents = 1000000; //TPhiCandidate->GetEntries();
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        TPhi_Multref ->GetEntry(iEvent);
        
        evKaonEfficiency.Multiplicity    *= 1./4.;
        evPhiEfficiency.Multiplicity     *= 1./4.;
        
        fPrintLoopTimer("Analysis",iEvent,nEvents,1000000);
        
        cout << evPhiEfficiency.Selection[iPhi] << endl;

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;

        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPhi; iPhi++ )
        {
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],evPhiEfficiency.InvMass[iPhi]);
            if ( !fAcceptCandidate(LPhi_candidate1.Rapidity(),evPhiEfficiency.InvMass[iPhi],LPhi_candidate1.Pt(),evPhiEfficiency.Multiplicity) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            cout << evPhiEfficiency.Selection[iPhi] << endl;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )
        {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;

            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[U_AccCand[iPhi]],evPhiEfficiency.Py[U_AccCand[iPhi]],evPhiEfficiency.Pz[U_AccCand[iPhi]],evPhiEfficiency.InvMass[U_AccCand[iPhi]]);

            // >> 1-Dimensional Analysis Fill   //
            //
            // >>-->> Utilities
            //
            Int_t   indexMult = fGetBinMult(evPhiEfficiency.Multiplicity);
            Int_t   indexSele = evPhiEfficiency.Selection[U_AccCand[iPhi]];
            Float_t indexTMom = LPhi_candidate1.Pt();
            //
            // >>-->> Yield
            //
            if ( indexSele >= 0 ) hTRU_1D                           ->  Fill(indexTMom);
            if ( indexSele >= 1 ) hGEN_1D                           ->  Fill(indexTMom);
            if ( indexSele >= 2 ) hREC_1D                           ->  Fill(indexTMom);
            if ( indexSele >= 0 ) hTRU_1D_in_2Dbin                  ->  Fill(indexTMom);
            if ( indexSele >= 1 ) hGEN_1D_in_2Dbin                  ->  Fill(indexTMom);
            if ( indexSele >= 2 ) hREC_1D_in_2Dbin                  ->  Fill(indexTMom);
            //
            // >>-->> Multiplicity
            //
            if ( indexSele >= 0 ) hTRU_1D_in_MT[indexMult]          ->  Fill(indexTMom);
            if ( indexSele >= 1 ) hGEN_1D_in_MT[indexMult]          ->  Fill(indexTMom);
            if ( indexSele >= 2 ) hREC_1D_in_MT[indexMult]          ->  Fill(indexTMom);
            if ( indexSele >= 0 ) hTRU_1D_in_MT_in_2Dbin[indexMult] ->  Fill(indexTMom);
            if ( indexSele >= 1 ) hGEN_1D_in_MT_in_2Dbin[indexMult] ->  Fill(indexTMom);
            if ( indexSele >= 2 ) hREC_1D_in_MT_in_2Dbin[indexMult] ->  Fill(indexTMom);
            
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )
            {
                // Must have at least 2 candidates
                if ( U_nAccept < 2 ) break;

                // Building Second Candidate
                LPhi_candidate2.SetXYZM(evPhiEfficiency.Px[U_AccCand[jPhi]],evPhiEfficiency.Py[U_AccCand[jPhi]],evPhiEfficiency.Pz[U_AccCand[jPhi]],evPhiEfficiency.InvMass[U_AccCand[jPhi]]);

                // >> 2-Dimensional Analysis Fill   //
                //
                // >>-->> Utilities
                //
                Int_t   jndexSele = evPhiEfficiency.Selection[U_AccCand[jPhi]];
                Float_t jndexTMom = LPhi_candidate2.Pt();
                //
                // >>-->> Yield
                //
                if ( indexSele >= 0 && jndexSele >= 0 ) hTRU_2D                     ->  Fill(indexTMom,jndexTMom,0.5);
                if ( indexSele >= 1 && jndexSele >= 1 ) hGEN_2D                     ->  Fill(indexTMom,jndexTMom,0.5);
                if ( indexSele >= 2 && jndexSele >= 2 ) hREC_2D                     ->  Fill(indexTMom,jndexTMom,0.5);
                //
                // >>-->> Multiplicity
                //
                if ( indexSele >= 0 && jndexSele >= 0 ) hTRU_2D_in_MT[indexMult]    ->  Fill(indexTMom,jndexTMom,0.5);
                if ( indexSele >= 1 && jndexSele >= 1 ) hGEN_2D_in_MT[indexMult]    ->  Fill(indexTMom,jndexTMom,0.5);
                if ( indexSele >= 2 && jndexSele >= 2 ) hREC_2D_in_MT[indexMult]    ->  Fill(indexTMom,jndexTMom,0.5);
                
                for ( Int_t kPhi = 0; kPhi < U_nAccept; kPhi++ )
                {
                    // Must have at least 3 candidates
                    if ( U_nAccept < 3 ) break;
                    
                    // >>->> 3-Dimensional Analysis Fill

                    for ( Int_t lPhi = 0; lPhi < U_nAccept; lPhi++ )
                    {
                        // Must have at least 4 candidates
                        if ( U_nAccept < 4 ) break;

                        // >>->> 4-Dimensional Analysis Fill
                        
                    }
                }
            }
        }
    }
    
    fStopTimer("Analysis");
    
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    
    // >> YIELD ANALYSIS //
    //
    hEFF_1D                             ->Divide(hREC_1D,           hGEN_1D,            1.,1.,"b");
    hEFF_1D_in_2Dbin                    ->Divide(hREC_1D_in_2Dbin,  hGEN_1D_in_2Dbin,   1.,1.,"b");
    hEFF_2D                             ->Divide(hREC_2D,           hGEN_2D,            1.,1.,"b");
    hREC_1D->Scale(1.,"width");
    hGEN_1D->Scale(1.,"width");
    hTRU_1D->Scale(1.,"width");
    hREC_1D_in_2Dbin->Scale(1.,"width");
    hGEN_1D_in_2Dbin->Scale(1.,"width");
    hTRU_1D_in_2Dbin->Scale(1.,"width");
    hREC_2D->Scale(1.,"width");
    hGEN_2D->Scale(1.,"width");
    hTRU_2D->Scale(1.,"width");
    //
    
    // >> MULTIPLICITY ANALYSIS //
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hEFF_1D_in_MT[iHisto]               ->Divide(hREC_1D_in_MT[iHisto],         hGEN_1D_in_MT[iHisto],          1.,1.,"b");
        hEFF_1D_in_MT_in_2Dbin[iHisto]      ->Divide(hREC_1D_in_MT_in_2Dbin[iHisto],hGEN_1D_in_MT_in_2Dbin[iHisto], 1.,1.,"b");
        hEFF_2D_in_MT[iHisto]               ->Divide(hREC_2D_in_MT[iHisto],         hGEN_2D_in_MT[iHisto],          1.,1.,"b");
        hREC_1D_in_MT[iHisto]->Scale(1.,"width");
        hGEN_1D_in_MT[iHisto]->Scale(1.,"width");
        hTRU_1D_in_MT[iHisto]->Scale(1.,"width");
        hREC_1D_in_MT_in_2Dbin[iHisto]->Scale(1.,"width");
        hGEN_1D_in_MT_in_2Dbin[iHisto]->Scale(1.,"width");
        hTRU_1D_in_MT_in_2Dbin[iHisto]->Scale(1.,"width");
        hREC_2D_in_MT[iHisto]->Scale(1.,"width");
        hGEN_2D_in_MT[iHisto]->Scale(1.,"width");
        hTRU_2D_in_MT[iHisto]->Scale(1.,"width");
    }
    //
    
    // >> TRIGGER ANALYSIS //
    
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Trigger Analysis
    //
    TFile *outFil1  =   new TFile   (fTrgPrePrMC,"recreate");
    //
    outFil1->Close();
    //
    // >> Yield Analysis
    //
    TFile *outFil2  =   new TFile   (fYldPrePrMC,"recreate");
    //
    hGEN_1D->Write();
    hTRU_1D->Write();
    hEFF_1D->Write();
    hREC_1D_in_2Dbin->Write();
    hGEN_1D_in_2Dbin->Write();
    hTRU_1D_in_2Dbin->Write();
    hEFF_1D_in_2Dbin->Write();
    hREC_2D->Write();
    hGEN_2D->Write();
    hTRU_2D->Write();
    hEFF_2D->Write();
    //
    outFil2->Close();
    //
    // >> Multiplicity Analysis
    //
    TFile *outFil3  =   new TFile   (fMltPrePrMC,"recreate");
    //
    for ( Int_t iHisto = 0; iHisto < nBinMult; iHisto++ )
    {
        hREC_1D_in_MT[iHisto]->Write();
        hGEN_1D_in_MT[iHisto]->Write();
        hTRU_1D_in_MT[iHisto]->Write();
        hEFF_1D_in_MT[iHisto]->Write();
        hREC_1D_in_MT_in_2Dbin[iHisto]->Write();
        hGEN_1D_in_MT_in_2Dbin[iHisto]->Write();
        hTRU_1D_in_MT_in_2Dbin[iHisto]->Write();
        hEFF_1D_in_MT_in_2Dbin[iHisto]->Write();
        hREC_2D_in_MT[iHisto]->Write();
        hGEN_2D_in_MT[iHisto]->Write();
        hTRU_2D_in_MT[iHisto]->Write();
        hEFF_2D_in_MT[iHisto]->Write();
    }
    //
    outFil3->Close();
    //
    // >-> Close input File
    //
    insFileMC->Close();
    //
}
