#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: All Set!

void Production_MeasureMass ( string fFileName = "/Volumes/[HD][Nikolajal]_Toshiba/Dataset/2010/Sim/2021_05_12/LHC14j4_STD.root", Int_t nEventsCut = -1., TString fOption = "yield" )   {
    //---------------------//
    //  Setting up input   //
    //---------------------//
    //
    // >-> Initialisation warnings
    //
    if ( fFileName == "" )  {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage PreProcessing_MC.C(\"Root_file_name.root\")" << endl;
        return;
    }
    if ( nEventsCut != -1 ) cout << "[WARNING] Choosing to limit the datasample to " << nEventsCut << " events" <<endl;
    fChooseOption(fOption);
    
    // Retrieving Event data
    TFile *insFileDT        =   new TFile   (fFileName.c_str());
    
    // Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileDT       ->Get(Form("%s%s",fPhiCandidate_Tree,""));
    TTree   *TKaonCandidate =   (TTree*)insFileDT       ->Get(Form("%s%s",fKaonCandidate_Tree,""));
    
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
        TKaonCandidate-> SetBranchAddress   ("EventMask",       &evKaonCandidate.EventMask);
        TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
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
        TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiCandidate.EventMask);
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("TrueInvMass",     &evPhiCandidate.TrueInvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
    }
    else
    {
        TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiCandidate.EventMask);
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
        TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
        TPhiCandidate-> SetBranchAddress    ("TrueInvMass",     &evPhiCandidate.TrueInvMass);
        TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
        TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
        
        TKaonCandidate-> SetBranchAddress   ("EventMask",       &evKaonCandidate.EventMask);
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
    fSetAllBins();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept,  U_nAccep2;
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >> YIELD ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //>>    Declaring all histograms
    TH1F  **hMassResolution_1D;
    TH1F  **hMassResolution_1D_in_2D_bin;
    TH2F ***hMassResolution_2D;
    TH1F  **hMassDistribution_1D;
    TH1F  **hMassDistribution_1D_in_2D_bin;
    TH2F ***hMassDistribution_2D;
    TH1F  **hMasTDistribution_1D;
    TH1F  **hMasTDistribution_1D_in_2D_bin;
    TH2F ***hMasTDistribution_2D;
    //
    TH1F   *hSlop_1D;
    TH1F   *hSlop_1D_in_2D_bin;
    //
    //>>    Defining all histograms
    hMassResolution_1D              =   new TH1F   *[nBinPT1D];
    hMassDistribution_1D            =   new TH1F   *[nBinPT1D];
    hMasTDistribution_1D            =   new TH1F   *[nBinPT1D];
    for ( Int_t iTer = 0; iTer < nBinPT1D; iTer++ )  {
        hName   =   Form("hMassResolution_1D_%i",iTer);
        hTitle  =   Form("hMassResolution_1D_%i",iTer);
        hMassResolution_1D[iTer]            =   new TH1F(hName,hTitle,nBinIMRs,fArrIMRs);
        hName   =   Form("hMassDistribution_1D_%i",iTer);
        hTitle  =   Form("hMassDistribution_1D_%i",iTer);
        hMassDistribution_1D[iTer]            =   new TH1F(hName,hTitle,nBinIM1D,fArrIM1D);
        hName   =   Form("hMasTDistribution_1D_%i",iTer);
        hTitle  =   Form("hMasTDistribution_1D_%i",iTer);
        hMasTDistribution_1D[iTer]            =   new TH1F(hName,hTitle,nBinIM1D,fArrIM1D);
    }
    hMassResolution_1D_in_2D_bin    =   new TH1F   *[nBinPT2D];
    hMassResolution_2D              =   new TH2F  **[nBinPT2D];
    hMassDistribution_1D_in_2D_bin  =   new TH1F   *[nBinPT2D];
    hMassDistribution_2D            =   new TH2F  **[nBinPT2D];
    hMasTDistribution_1D_in_2D_bin  =   new TH1F   *[nBinPT2D];
    hMasTDistribution_2D            =   new TH2F  **[nBinPT2D];
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )  {
        hName   =   Form("hMassResolution_1D_in_2D_bin_%i",iTer);
        hTitle  =   Form("hMassResolution_1D_in_2D_bin_%i",iTer);
        hMassResolution_1D_in_2D_bin[iTer]  =   new TH1F(hName,hTitle,nBinIMRs,fArrIMRs);
        hName   =   Form("hMassDistribution_1D_in_2D_bin_%i",iTer);
        hTitle  =   Form("hMassDistribution_1D_in_2D_bin_%i",iTer);
        hMassDistribution_1D_in_2D_bin[iTer]  =   new TH1F(hName,hTitle,nBinIM2D,fArrIM2D);
        hName   =   Form("hMasTDistribution_1D_in_2D_bin_%i",iTer);
        hTitle  =   Form("hMasTDistribution_1D_in_2D_bin_%i",iTer);
        hMasTDistribution_1D_in_2D_bin[iTer]  =   new TH1F(hName,hTitle,nBinIM2D,fArrIM2D);
        //
        hMassResolution_2D[iTer]            =   new TH2F   *[nBinPT2D];
        hMassDistribution_2D[iTer]          =   new TH2F   *[nBinPT2D];
        hMasTDistribution_2D[iTer]          =   new TH2F   *[nBinPT2D];
        for ( Int_t jTer = 0; jTer < nBinPT2D; jTer++ )  {
            hName   =   Form("hMassResolution_2D_%i_%i",iTer,jTer);
            hTitle  =   Form("hMassResolution_2D_%i_%i",iTer,jTer);
            hMassResolution_2D[iTer][jTer]  =   new TH2F(hName,hTitle,nBinIMR2,fArrIMR2,nBinIMR2,fArrIMR2);
            hName   =   Form("hMassDistribution_2D_%i_%i",iTer,jTer);
            hTitle  =   Form("hMassDistribution_2D_%i_%i",iTer,jTer);
            hMassDistribution_2D[iTer][jTer]  =   new TH2F(hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            hName   =   Form("hMasTDistribution_2D_%i_%i",iTer,jTer);
            hTitle  =   Form("hMasTDistribution_2D_%i_%i",iTer,jTer);
            hMasTDistribution_2D[iTer][jTer]  =   new TH2F(hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
        }
    }
    //
    hName   =   Form("hSlop_1D");
    hTitle  =   Form("hSlop_1D");
    hSlop_1D    =   new TH1F(hName,hTitle,nBinPT1D,fArrPT1D);
    hName   =   Form("hSlop_1D_in_2D_bin");
    hTitle  =   Form("hSlop_1D_in_2D_bin");
    hSlop_1D_in_2D_bin    =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    //>>    Evaluating entries and saving them for later
    Int_t nEvents = (!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : nEventsCut);
    //
    //>>    If there are actually entries start timer
    if ( nEvents > 0 )  fStartTimer("Resolution Histogram Production");
    //
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Resolution Histogram Production",iEvent,nEvents,kPrintIntervalPP);
        
        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        // Discarding Pile-up events
        if ( fCheckMask(evPhiCandidate.EventMask,1) ) continue;
        
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; iPhi++ )  {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fAcceptCandidate(evPhiCandidate.InvMass[iPhi],LPhi_candidate1.Pt()) ) continue;
            if ( !kDoRapidity && !fCutRapidity(LPhi_candidate1.Rapidity()) ) continue;
            if (  kOnlyTrue && !evPhiCandidate.Nature[iPhi] ) continue;
            U_AccCand[U_nAccept] = iPhi;
            evPhiCandidate.pT[iPhi]     =   LPhi_candidate1.Pt();
            evPhiCandidate.Rap[iPhi]    =   LPhi_candidate1.Rapidity();
            evPhiCandidate.iRap[iPhi]   =   fGetBinRap_(evPhiCandidate.Rap[iPhi]);
            evPhiCandidate.iPT1D[iPhi]  =   fGetBinPT1D(evPhiCandidate.pT[iPhi]);
            evPhiCandidate.iPT2D[iPhi]  =   fGetBinPT2D(evPhiCandidate.pT[iPhi]);
            evPhiCandidate.kHasRap[iPhi]=   evPhiCandidate.iRap[iPhi] != -1;
            U_nAccept++;
        }
        //
        // >>-->> Utilities
        //
        // >>-->>-->> Multiplicity
        //
        Int_t   iMult               =   fGetBinMult(evPhiCandidate.Multiplicity);
        evPhiCandidate.kHasMult     =   evPhiCandidate.iMult != -1;
        //
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;

            // Selecting valid candidates
            if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi) ) continue;
            //
            // >> 1-Dimensional Analysis Fill
            //
            // >>-->> Utilities
            //
            // >>-->>-->> True Phis
            //
            Int_t   iPT1D               =   evPhiCandidate.iPT1D[U_AccCand[iPhi]];
            Int_t   iPT2D               =   evPhiCandidate.iPT2D[U_AccCand[iPhi]];
            Float_t iInvarMass          =   evPhiCandidate.InvMass[U_AccCand[iPhi]];
            Float_t iTrueIMass          =   evPhiCandidate.TrueInvMass[U_AccCand[iPhi]];
            Float_t iRapidity           =   evPhiCandidate.Rap[U_AccCand[iPhi]];
            Int_t   iRap                =   evPhiCandidate.iRap[U_AccCand[iPhi]];
            Bool_t  fHasRapidity        =   evPhiCandidate.kHasRap[U_AccCand[iPhi]];
            //
            // >->-->-> Yield
            //
            if ( kDoYield && fHasRapidity && iTrueIMass !=0 ) {
                hMassResolution_1D[iPT1D]           ->Fill(-iTrueIMass+iInvarMass);
                hMassResolution_1D_in_2D_bin[iPT2D] ->Fill(-iTrueIMass+iInvarMass);
                hMassDistribution_1D[iPT1D]           ->Fill(iInvarMass);
                hMassDistribution_1D_in_2D_bin[iPT2D] ->Fill(iInvarMass);
                hMasTDistribution_1D[iPT1D]           ->Fill(iTrueIMass);
                hMasTDistribution_1D_in_2D_bin[iPT2D] ->Fill(iTrueIMass);
            }
            //
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )    {
                // Must have at least 2 candidates
                if ( U_nAccept < 2 ) break;

                // Selecting valid candidates
                if ( !fAcceptCandidate( evPhiCandidate, U_AccCand, iPhi, jPhi) ) continue;
                
                // >-> 2-Dimensional Analysis Fill
                //
                // >->-->-> Utilities
                //
                // >>-->>-->> True Phis
                //
                Int_t   jPT1D               =   evPhiCandidate.iPT1D[U_AccCand[jPhi]];
                Int_t   jPT2D               =   evPhiCandidate.iPT2D[U_AccCand[jPhi]];
                Float_t jInvarMass          =   evPhiCandidate.InvMass[U_AccCand[jPhi]];
                Float_t jTrueIMass          =   evPhiCandidate.TrueInvMass[U_AccCand[jPhi]];
                Float_t jRapidity           =   evPhiCandidate.Rap[U_AccCand[jPhi]];
                Int_t   jRap                =   evPhiCandidate.iRap[U_AccCand[jPhi]];
                        fHasRapidity        =   evPhiCandidate.kHasRap[U_AccCand[iPhi]] && evPhiCandidate.kHasRap[U_AccCand[jPhi]];
                Int_t   ijRap               =   fGetBinRap_(evPhiCandidate.Rap[U_AccCand[iPhi]]-evPhiCandidate.Rap[U_AccCand[jPhi]]);
                Bool_t  fHasDiffRap         =   ijRap != -1;
                //
                if  ( fHasRapidity && kDoYield && iTrueIMass !=0  && jTrueIMass !=0 )    {
                //
                // >->-->-> Yield
                //
                    hMassResolution_2D[iPT2D][jPT2D]    ->Fill(-iTrueIMass+iInvarMass,-jTrueIMass+jInvarMass);
                    hMassDistribution_2D[iPT2D][jPT2D]  ->Fill(iInvarMass,jInvarMass);
                    hMasTDistribution_2D[iPT2D][jPT2D]  ->Fill(iTrueIMass,jTrueIMass);
                //
                }
                //
            }
        }
    }
    //>>    If there are actually entries start timer
    if ( nEvents > 0 )  fStopTimer("Resolution Histogram Production");
    //
    //--------------------------//
    //  Saving output objects   //
    //--------------------------//
    //
    // >> Yield Analysis
    //
    gROOT->ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Dir_.Data(),"yield")));
    TFile  *fOutTest    =   new TFile(Form(kMassResolution_Prod,"yield"),"recreate");
    for ( Int_t iTer = 0; iTer < nBinPT1D; iTer++ )  {
        hMassResolution_1D[iTer]            ->  Write();
        hMassDistribution_1D[iTer]          ->  Write();
        hMasTDistribution_1D[iTer]          ->  Write();
    }
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )  {
        hMassDistribution_1D_in_2D_bin[iTer]  ->  Write();
        hMasTDistribution_1D_in_2D_bin[iTer]  ->  Write();
        hMassResolution_1D_in_2D_bin[iTer]            ->  Write();
        for ( Int_t jTer = 0; jTer < nBinPT2D; jTer++ )  {
            hMassDistribution_2D[iTer][jTer]  ->  Write();
            hMasTDistribution_2D[iTer][jTer]  ->  Write();
            hMassResolution_2D[iTer][jTer]           ->  Write();
        }
    }
    fOutTest->Close();
    //
    // >-> Close input File
    //
    insFileDT->Close();
    //
}
