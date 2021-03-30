#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: All Set!

void PreProcessing_MC ( string fFileName = "", Int_t nEventsCut = -1., TString fOption = "" )
{
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
    
    //Retrieving Event data
    TFile *insFileMC        =   new TFile   (fFileName.c_str());
    
    //Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileMC->Get(fPhiCandidateEff_Tree);
    TTree   *TKaonCandidate =   nullptr;//(TTree*)insFileMC->Get(fKaonCandidateEff_Tree);
    
    // Retrieving Event Count Histogram
    TList  *fQCOutputList   =   (TList*)insFileMC       ->Get("fQCOutputList");
    TH1D   *fHEventCount    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enumerate");
    
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
        TKaonCandidate->SetBranchAddress    ("EventMask",       &evKaonEfficiency.EventMask);
        TKaonCandidate->SetBranchAddress    ("TrueEventMask",   &evKaonEfficiency.TrueEventMask);
        TKaonCandidate->SetBranchAddress    ("Multiplicity",    &evKaonEfficiency.Multiplicity);
        TKaonCandidate->SetBranchAddress    ("nKaon",           &evKaonEfficiency.nKaon);
        TKaonCandidate->SetBranchAddress    ("Px",              &evKaonEfficiency.Px);
        TKaonCandidate->SetBranchAddress    ("Py",              &evKaonEfficiency.Py);
        TKaonCandidate->SetBranchAddress    ("Pz",              &evKaonEfficiency.Pz);
        TKaonCandidate->SetBranchAddress    ("Selection",       &evKaonEfficiency.Selection);
    }
    else if ( !TKaonCandidate )
    {
        cout << "[INFO] No KaonCandidate Tree, switching to Phi Analysis" << endl;
        TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiEfficiency.EventMask);
        TPhiCandidate-> SetBranchAddress    ("TrueEventMask",   &evPhiEfficiency.TrueEventMask);
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiEfficiency.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiEfficiency.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiEfficiency.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiEfficiency.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiEfficiency.Pz);
        TPhiCandidate-> SetBranchAddress    ("Selection",       &evPhiEfficiency.Selection);
    }
    else
    {
        TKaonCandidate->SetBranchAddress    ("EventMask",       &evKaonEfficiency.EventMask);
        TKaonCandidate->SetBranchAddress    ("TrueEventMask",   &evKaonEfficiency.TrueEventMask);
        TKaonCandidate->SetBranchAddress    ("Multiplicity",    &evKaonEfficiency.Multiplicity);
        TKaonCandidate->SetBranchAddress    ("nKaon",           &evKaonEfficiency.nKaon);
        TKaonCandidate->SetBranchAddress    ("Px",              &evKaonEfficiency.Px);
        TKaonCandidate->SetBranchAddress    ("Py",              &evKaonEfficiency.Py);
        TKaonCandidate->SetBranchAddress    ("Pz",              &evKaonEfficiency.Pz);
        TKaonCandidate->SetBranchAddress    ("Selection",       &evKaonEfficiency.Selection);
        
        TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiEfficiency.EventMask);
        TPhiCandidate-> SetBranchAddress    ("TrueEventMask",   &evPhiEfficiency.TrueEventMask);
        TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiEfficiency.Multiplicity);
        TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiEfficiency.nPhi);
        TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiEfficiency.Px);
        TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiEfficiency.Py);
        TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiEfficiency.Pz);
        TPhiCandidate-> SetBranchAddress    ("Selection",       &evPhiEfficiency.Selection);
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
    Int_t       U_nAccept,  U_nAccep2;
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >> YIELD ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hREC_1D;
    TH1F       *hGEN_1D;
    TH1F       *hTRU_1D;
    TH1F       *hTRU_ALL_1D;
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
    hName       =   Form("hTRU_ALL_1D");
    hTitle      =   Form("hTRU_ALL_1D");
    hTRU_ALL_1D =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hTRU_ALL_1D,"PT 1D");
    //
    hName       =   Form("hEFF_1D");
    hTitle      =   Form("hEFF_1D");
    hEFF_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hEFF_1D,"PT 1D");
    //
    hName       =   Form("hREC_1D_in_2D_bin");
    hTitle      =   Form("hREC_1D_in_2D_bin");
    hREC_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hREC_1D_in_2Dbin,"PT 1D");
    //
    hName       =   Form("hGEN_1D_in_2D_bin");
    hTitle      =   Form("hGEN_1D_in_2D_bin");
    hGEN_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_1D_in_2Dbin,"PT 1D");
    //
    hName       =   Form("hTRU_1D_in_2D_bin");
    hTitle      =   Form("hTRU_1D_in_2D_bin");
    hTRU_1D_in_2Dbin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_1D_in_2Dbin,"PT 1D");
    //
    hName       =   Form("hEFF_1D_in_2D_bin");
    hTitle      =   Form("hEFF_1D_in_2D_bin");
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
    TH2F       *hTRU_ALL_2D;
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
    hName       =   Form("hTRU_ALL_2D");
    hTitle      =   Form("hTRU_ALL_2D");
    hTRU_ALL_2D =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_ALL_2D,"PT 2D");
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
    TH1F      **hREC_1D_in_MT               = new TH1F     *[nBinMult+1];
    TH1F      **hGEN_1D_in_MT               = new TH1F     *[nBinMult+1];
    TH1F      **hTRU_1D_in_MT               = new TH1F     *[nBinMult+1];
    TH1F      **hEFF_1D_in_MT               = new TH1F     *[nBinMult+1];
    TH1F      **hREC_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult+1];
    TH1F      **hGEN_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult+1];
    TH1F      **hTRU_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult+1];
    TH1F      **hEFF_1D_in_MT_in_2Dbin      = new TH1F     *[nBinMult+1];
    //
    //  Defining MT-Differential histograms
    //
    hName = Form("hREC_1D_in_MT_%i",0);
    hTitle= Form("hREC_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hREC_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hREC_1D_in_MT[0],"PT 1D");
    
    hName = Form("hGEN_1D_in_MT_%i",0);
    hTitle= Form("hGEN_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hGEN_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hGEN_1D_in_MT[0],"PT 1D");
    
    hName = Form("hTRU_1D_in_MT_%i",0);
    hTitle= Form("hTRU_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hTRU_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hTRU_1D_in_MT[0],"PT 1D");
    
    hName = Form("hEFF_1D_in_MT_%i",0);
    hTitle= Form("hEFF_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hEFF_1D_in_MT[0],"PT 1D");
    
    hName = Form("hREC_1D_in_MT_in_2Dbin_%i",0);
    hTitle= Form("hREC_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hREC_1D_in_MT_in_2Dbin[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hREC_1D_in_MT_in_2Dbin[0],"PT 1D");
        
    hName = Form("hGEN_1D_in_MT_in_2Dbin_%i",0);
    hTitle= Form("hGEN_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hGEN_1D_in_MT_in_2Dbin[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_1D_in_MT_in_2Dbin[0],"PT 1D");
    
    hName = Form("hTRU_1D_in_MT_in_2Dbin_%i",0);
    hTitle= Form("hTRU_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hTRU_1D_in_MT_in_2Dbin[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_1D_in_MT_in_2Dbin[0],"PT 1D");
    
    hName = Form("hEFF_1D_in_MT_in_2Dbin_%i",0);
    hTitle= Form("hEFF_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_1D_in_MT_in_2Dbin[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_1D_in_MT_in_2Dbin[0],"PT 1D");
    
    for ( Int_t iMult = 1; iMult <= nBinMult; iMult++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iMult);
        hTitle= Form("hREC_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hREC_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hREC_1D_in_MT[iMult],"PT 1D");
        
        hName = Form("hGEN_1D_in_MT_%i",iMult);
        hTitle= Form("hGEN_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hGEN_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hGEN_1D_in_MT[iMult],"PT 1D");
        
        hName = Form("hTRU_1D_in_MT_%i",iMult);
        hTitle= Form("hTRU_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hTRU_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hTRU_1D_in_MT[iMult],"PT 1D");
        
        hName = Form("hEFF_1D_in_MT_%i",iMult);
        hTitle= Form("hEFF_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hEFF_1D_in_MT[iMult],"PT 1D");
        
        hName = Form("hREC_1D_in_MT_in_2Dbin_%i",iMult);
        hTitle= Form("hREC_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hREC_1D_in_MT_in_2Dbin[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hREC_1D_in_MT_in_2Dbin[iMult],"PT 1D");
            
        hName = Form("hGEN_1D_in_MT_in_2Dbin_%i",iMult);
        hTitle= Form("hGEN_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hGEN_1D_in_MT_in_2Dbin[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_1D_in_MT_in_2Dbin[iMult],"PT 1D");
        
        hName = Form("hTRU_1D_in_MT_in_2Dbin_%i",iMult);
        hTitle= Form("hTRU_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hTRU_1D_in_MT_in_2Dbin[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_1D_in_MT_in_2Dbin[iMult],"PT 1D");
        
        hName = Form("hEFF_1D_in_MT_in_2Dbin_%i",iMult);
        hTitle= Form("hEFF_1D_in_MT_in_2Dbin Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_1D_in_MT_in_2Dbin[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_1D_in_MT_in_2Dbin[iMult],"PT 1D");
    }
    //
    // >>-->> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hREC_2D_in_MT               = new TH2F     *[nBinMult+1];
    TH2F      **hGEN_2D_in_MT               = new TH2F     *[nBinMult+1];
    TH2F      **hTRU_2D_in_MT               = new TH2F     *[nBinMult+1];
    TH2F      **hEFF_2D_in_MT               = new TH2F     *[nBinMult+1];
    TH2F      **hEFF_2D_in_MT_fr_1D         = new TH2F     *[nBinMult+1];
    //
    //  Defining MT-Differential histograms
    //
    hName = Form("hREC_2D_in_MT_%i",0);
    hTitle= Form("hREC_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hREC_2D_in_MT[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hREC_2D_in_MT[0],"PT 2D");
    
    hName = Form("hGEN_2D_in_MT_%i",0);
    hTitle= Form("hGEN_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hGEN_2D_in_MT[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_2D_in_MT[0],"PT 2D");
    
    hName = Form("hTRU_2D_in_MT_%i",0);
    hTitle= Form("hTRU_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hTRU_2D_in_MT[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_2D_in_MT[0],"PT 2D");
    
    hName = Form("hEFF_2D_in_MT_%i",0);
    hTitle= Form("hEFF_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_2D_in_MT[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_2D_in_MT[0],"PT 2D");
    
    hName = Form("hEFF_2D_in_MT_fr_1D_%i",0);
    hTitle= Form("hEFF_2D_in_MT_fr_1D Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_2D_in_MT_fr_1D[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_2D_in_MT_fr_1D[0],"PT 2D");
    
    for ( Int_t iMult = 1; iMult <= nBinMult; iMult++ )
    {
        hName = Form("hREC_2D_in_MT_%i",iMult);
        hTitle= Form("hREC_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hREC_2D_in_MT[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hREC_2D_in_MT[iMult],"PT 2D");
        
        hName = Form("hGEN_2D_in_MT_%i",iMult);
        hTitle= Form("hGEN_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hGEN_2D_in_MT[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_2D_in_MT[iMult],"PT 2D");
        
        hName = Form("hTRU_2D_in_MT_%i",iMult);
        hTitle= Form("hTRU_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hTRU_2D_in_MT[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_2D_in_MT[iMult],"PT 2D");
        
        hName = Form("hEFF_2D_in_MT_%i",iMult);
        hTitle= Form("hEFF_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_2D_in_MT[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_2D_in_MT[iMult],"PT 2D");
        
        hName = Form("hEFF_2D_in_MT_fr_1D_%i",iMult);
        hTitle= Form("hEFF_2D_in_MT_fr_1D Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_2D_in_MT_fr_1D[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_2D_in_MT_fr_1D[iMult],"PT 2D");
    }
    //
    
    // >> RAPIDITY ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_RP               = new TH1F     *[nBinRap_];
    TH1F      **hGEN_1D_in_RP               = new TH1F     *[nBinRap_];
    TH1F      **hTRU_1D_in_RP               = new TH1F     *[nBinRap_];
    TH1F      **hEFF_1D_in_RP               = new TH1F     *[nBinRap_];
    TH1F      **hREC_1D_in_RP_in_2Dbin      = new TH1F     *[nBinRap_];
    TH1F      **hGEN_1D_in_RP_in_2Dbin      = new TH1F     *[nBinRap_];
    TH1F      **hTRU_1D_in_RP_in_2Dbin      = new TH1F     *[nBinRap_];
    TH1F      **hEFF_1D_in_RP_in_2Dbin      = new TH1F     *[nBinRap_];
    //
    //  Defining RP-Differential histograms
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName = Form("hREC_1D_in_RP_%i",iRap);
        hTitle= Form("hREC_1D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_1D_in_RP[iRap]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hREC_1D_in_RP[iRap],"PT 1D");
        
        hName = Form("hGEN_1D_in_RP_%i",iRap);
        hTitle= Form("hGEN_1D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hGEN_1D_in_RP[iRap]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hGEN_1D_in_RP[iRap],"PT 1D");
        
        hName = Form("hTRU_1D_in_RP_%i",iRap);
        hTitle= Form("hTRU_1D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hTRU_1D_in_RP[iRap]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hTRU_1D_in_RP[iRap],"PT 1D");
        
        hName = Form("hEFF_1D_in_RP_%i",iRap);
        hTitle= Form("hEFF_1D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hEFF_1D_in_RP[iRap]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hEFF_1D_in_RP[iRap],"PT 1D");
        
        hName = Form("hREC_1D_in_RP_in_2Dbin_%i",iRap);
        hTitle= Form("hREC_1D_in_RP_in_2Dbin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_1D_in_RP_in_2Dbin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hREC_1D_in_RP_in_2Dbin[iRap],"PT 1D");
            
        hName = Form("hGEN_1D_in_RP_in_2Dbin_%i",iRap);
        hTitle= Form("hGEN_1D_in_RP_in_2Dbin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hGEN_1D_in_RP_in_2Dbin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_1D_in_RP_in_2Dbin[iRap],"PT 1D");
        
        hName = Form("hTRU_1D_in_RP_in_2Dbin_%i",iRap);
        hTitle= Form("hTRU_1D_in_RP_in_2Dbin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hTRU_1D_in_RP_in_2Dbin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_1D_in_RP_in_2Dbin[iRap],"PT 1D");
        
        hName = Form("hEFF_1D_in_RP_in_2Dbin_%i",iRap);
        hTitle= Form("hEFF_1D_in_RP_in_2Dbin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hEFF_1D_in_RP_in_2Dbin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_1D_in_RP_in_2Dbin[iRap],"PT 1D");
    }
    //
    // >>-->> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F      **hREC_2D_in_RP               = new TH2F     *[nBinRap_];
    TH2F      **hGEN_2D_in_RP               = new TH2F     *[nBinRap_];
    TH2F      **hTRU_2D_in_RP               = new TH2F     *[nBinRap_];
    TH2F      **hEFF_2D_in_RP               = new TH2F     *[nBinRap_];
    TH2F      **hEFF_2D_in_RP_fr_1D         = new TH2F     *[nBinRap_];
    //
    //  Defining RP-Differential histograms
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName = Form("hREC_2D_in_RP_%i",iRap);
        hTitle= Form("hREC_2D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_2D_in_RP[iRap]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hREC_2D_in_RP[iRap],"PT 2D");
        
        hName = Form("hGEN_2D_in_RP_%i",iRap);
        hTitle= Form("hGEN_2D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hGEN_2D_in_RP[iRap]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_2D_in_RP[iRap],"PT 2D");
        
        hName = Form("hTRU_2D_in_RP_%i",iRap);
        hTitle= Form("hTRU_2D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hTRU_2D_in_RP[iRap]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_2D_in_RP[iRap],"PT 2D");
        
        hName = Form("hEFF_2D_in_RP_%i",iRap);
        hTitle= Form("hEFF_2D_in_RP Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hEFF_2D_in_RP[iRap]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_2D_in_RP[iRap],"PT 2D");
        
        hName = Form("hEFF_2D_in_RP_fr_1D_%i",iRap);
        hTitle= Form("hEFF_2D_in_RP_fr_1D Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hEFF_2D_in_RP_fr_1D[iRap]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_2D_in_RP_fr_1D[iRap],"PT 2D");
    }
    
    // >> TRIGGER ANALYSIS //

    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Analysis");
    
    // Evaluating entries
    Int_t nEvents = (!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : nEventsCut);
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Analysis",iEvent,nEvents,kPrintIntervalPP);

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPhi; iPhi++ ) {
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],kPhiMesonMass_);
            if ( !fAcceptCandidate(kPhiMesonMass_,LPhi_candidate1.Pt()) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;

            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[U_AccCand[iPhi]],evPhiEfficiency.Py[U_AccCand[iPhi]],evPhiEfficiency.Pz[U_AccCand[iPhi]],kPhiMesonMass_);

            // >> 1-Dimensional Analysis Fill
            //
            // >>-->> Utilities
            //
            // >>-->>-->> Event
            //
            Bool_t  fNoVtxRec           =   (fCheckMask(evPhiEfficiency.TrueEventMask,0) || fCheckMask(evPhiEfficiency.TrueEventMask,1));
            Bool_t  fIsMBevent          =   (evPhiEfficiency.TrueEventMask == 0);
            //
            // >>-->>-->> Multiplicity
            //
            Int_t   iMult               =   fGetBinMult(evPhiEfficiency.Multiplicity);
            Bool_t  fHasMultiplicity    =   iMult != -1;
            //
            // >>-->>-->> True Phis
            //
            Int_t   iSelection          =   (int)evPhiEfficiency.Selection[U_AccCand[iPhi]];
            Bool_t  iIsGen              =   (iSelection >= 1);
            Bool_t  iIsRec              =   (iSelection >= 2);
            Float_t iTransMom           =   LPhi_candidate1.Pt();
            Float_t iRapidity           =   LPhi_candidate1.Rapidity();
            Int_t   iRap                =   fGetBinRap_(LPhi_candidate1.Rapidity());
            Bool_t  fHasRapidity        =   iRap != -1;
            //
            if ( (fNoVtxRec || fIsMBevent) && fHasRapidity )  {
                hTRU_ALL_1D                                 ->  Fill(iTransMom);
            }
            if ( fIsMBevent && fHasRapidity )   {
            //
            // >>-->> Rapidity
            //
                hTRU_1D_in_RP[iRap]                         ->  Fill(iTransMom);
                hTRU_1D_in_RP_in_2Dbin[iRap]                ->  Fill(iTransMom);
                if ( iIsGen )   {
                    hGEN_1D_in_RP[iRap]                     ->  Fill(iTransMom);
                    hGEN_1D_in_RP_in_2Dbin[iRap]            ->  Fill(iTransMom);
                } if ( iIsRec )  {
                    hREC_1D_in_RP[iRap]                     ->  Fill(iTransMom);
                    hREC_1D_in_RP_in_2Dbin[iRap]            ->  Fill(iTransMom);
                }
            //
            // >>-->> Yield
            //
            
                hTRU_1D                                     ->  Fill(iTransMom);
                hTRU_1D_in_2Dbin                            ->  Fill(iTransMom);
                if ( iIsGen )   {
                    hGEN_1D                                 ->  Fill(iTransMom);
                    hGEN_1D_in_2Dbin                        ->  Fill(iTransMom);
                } if ( iIsRec )  {
                    hREC_1D                                 ->  Fill(iTransMom);
                    hREC_1D_in_2Dbin                        ->  Fill(iTransMom);
                }
            //
            // >>-->> Multiplicity
            //
                if ( fHasMultiplicity ) {
                    hTRU_1D_in_MT[iMult+1]                  ->  Fill(iTransMom);
                    hTRU_1D_in_MT[0]                        ->  Fill(iTransMom);
                    hTRU_1D_in_MT_in_2Dbin[iMult+1]         ->  Fill(iTransMom);
                    hTRU_1D_in_MT_in_2Dbin[0]               ->  Fill(iTransMom);
                    if ( iIsGen )   {
                        hGEN_1D_in_MT[iMult+1]              ->  Fill(iTransMom);
                        hGEN_1D_in_MT[0]                    ->  Fill(iTransMom);
                        hGEN_1D_in_MT_in_2Dbin[iMult+1]     ->  Fill(iTransMom);
                        hGEN_1D_in_MT_in_2Dbin[0]           ->  Fill(iTransMom);
                    } if ( iIsRec )  {
                        hREC_1D_in_MT[iMult+1]              ->  Fill(iTransMom);
                        hREC_1D_in_MT[0]                    ->  Fill(iTransMom);
                        hREC_1D_in_MT_in_2Dbin[iMult+1]     ->  Fill(iTransMom);
                        hREC_1D_in_MT_in_2Dbin[0]           ->  Fill(iTransMom);
                    }
                }
            }
            //
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )    {
                // Must have at least 2 candidates
                if ( U_nAccept < 2 ) break;
                
                // Protection against auto-correlation
                if ( iPhi == jPhi ) continue;

                // Building Second Candidate
                LPhi_candidate2.SetXYZM(evPhiEfficiency.Px[U_AccCand[jPhi]],evPhiEfficiency.Py[U_AccCand[jPhi]],evPhiEfficiency.Pz[U_AccCand[jPhi]],kPhiMesonMass_);

                // >> 2-Dimensional Analysis Fill
                //
                // >>-->> Utilities
                //
                // >>-->>-->> True Phis
                //
                Int_t   jSelection          =   (int)evPhiEfficiency.Selection[U_AccCand[jPhi]];
                Bool_t  jIsGen              =   (jSelection >= 1);
                Bool_t  jIsRec              =   (jSelection >= 2);
                Float_t jTransMom           =   LPhi_candidate2.Pt();
                Float_t jRapidity           =   LPhi_candidate2.Rapidity();
                Int_t   jRap                =   fGetBinRap_(LPhi_candidate2.Rapidity());
                        fHasRapidity        =   jRap != -1 && iRap != -1;
                Int_t   ijRap               =   fGetBinRap_(LPhi_candidate2.Rapidity()-LPhi_candidate1.Rapidity());
                //
                if ( (fNoVtxRec || fIsMBevent) && fHasRapidity )  {
                    hTRU_ALL_2D                                     ->  Fill(iTransMom,jTransMom,0.5);
                }
                if ( fIsMBevent )   {
                //
                // >>-->> Rapidity
                //
                    if ( ijRap != -1 )   {
                        hTRU_2D_in_RP[ijRap]                            ->  Fill(iTransMom,jTransMom,0.5);
                        if ( iIsGen && jIsGen )   {
                            hGEN_2D_in_RP[ijRap]                        ->  Fill(iTransMom,jTransMom,0.5);
                        } if ( iIsRec && jIsRec)  {
                            hREC_2D_in_RP[ijRap]                        ->  Fill(iTransMom,jTransMom,0.5);
                        }
                    }
                //
                // >>-->> Yield
                //
                    if ( fHasRapidity ) {
                        hTRU_2D                                     ->  Fill(iTransMom,jTransMom,0.5);
                        if ( iIsGen && jIsGen )   {
                            hGEN_2D                                 ->  Fill(iTransMom,jTransMom,0.5);
                        } if ( iIsRec && jIsRec )  {
                            hREC_2D                                 ->  Fill(iTransMom,jTransMom,0.5);
                        }
                //
                // >>-->> Multiplicity
                //
                        if ( fHasMultiplicity ) {
                            hTRU_2D_in_MT[iMult+1]                  ->  Fill(iTransMom,jTransMom,0.5);
                            hTRU_2D_in_MT[0]                        ->  Fill(iTransMom,jTransMom,0.5);
                            if ( iIsGen && jIsGen )   {
                                hGEN_2D_in_MT[iMult+1]              ->  Fill(iTransMom,jTransMom,0.5);
                                hGEN_2D_in_MT[0]                    ->  Fill(iTransMom,jTransMom,0.5);
                            } if ( iIsRec && jIsRec )  {
                                hREC_2D_in_MT[iMult+1]              ->  Fill(iTransMom,jTransMom,0.5);
                                hREC_2D_in_MT[0]                    ->  Fill(iTransMom,jTransMom,0.5);
                            }
                        }
                    }
                }
                //
                for ( Int_t kPhi = 0; kPhi < U_nAccept; kPhi++ )    {
                    // Must have at least 3 candidates
                    if ( U_nAccept < 3 ) break;
                    //
                    // >>->> 3-Dimensional Analysis Fill
                    //
                    for ( Int_t lPhi = 0; lPhi < U_nAccept; lPhi++ )    {
                        // Must have at least 4 candidates
                        if ( U_nAccept < 4 ) break;
                        //
                        // >>->> 4-Dimensional Analysis Fill
                        //
                    }
                }
            }
        }
    }
    
    fStopTimer("Analysis");
    
    
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    //
    // >> YIELD ANALYSIS //
    //
    auto fNormEvent = fHEventCount->GetBinContent(1);
    hEFF_1D                             ->Divide(hREC_1D,           hGEN_1D,            1.,1.,"b");
    hEFF_1D_in_2Dbin                    ->Divide(hREC_1D_in_2Dbin,  hGEN_1D_in_2Dbin,   1.,1.,"b");
    hEFF_2D                             ->Divide(hREC_2D,           hGEN_2D,            1.,1.,"b");
    hREC_1D->Scale(1.,"width");
    hGEN_1D->Scale(1.,"width");
    hTRU_1D->Scale(1.,"width");
    hTRU_ALL_1D->Scale(1.,"width");
    hREC_1D_in_2Dbin->Scale(1.,"width");
    hGEN_1D_in_2Dbin->Scale(1.,"width");
    hTRU_1D_in_2Dbin->Scale(1.,"width");
    hREC_2D->Scale(1.,"width");
    hGEN_2D->Scale(1.,"width");
    hTRU_2D->Scale(1.,"width");
    hTRU_ALL_2D->Scale(1.,"width");
    hREC_1D->Scale(1./fNormEvent);
    hGEN_1D->Scale(1./fNormEvent);
    hTRU_1D->Scale(1./fNormEvent);
    hTRU_ALL_1D->Scale(1./fNormEvent);
    hREC_1D_in_2Dbin->Scale(1./fNormEvent);
    hGEN_1D_in_2Dbin->Scale(1./fNormEvent);
    hTRU_1D_in_2Dbin->Scale(1./fNormEvent);
    hREC_2D->Scale(1./fNormEvent);
    hGEN_2D->Scale(1./fNormEvent);
    hTRU_2D->Scale(1./fNormEvent);
    hTRU_ALL_2D->Scale(1./fNormEvent);
    //
    
    // >> MULTIPLICITY ANALYSIS //
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hEFF_1D_in_MT[iMult]               ->Divide(hREC_1D_in_MT[iMult],         hGEN_1D_in_MT[iMult],          1.,1.,"b");
        hEFF_1D_in_MT_in_2Dbin[iMult]      ->Divide(hREC_1D_in_MT_in_2Dbin[iMult],hGEN_1D_in_MT_in_2Dbin[iMult], 1.,1.,"b");
        hEFF_2D_in_MT[iMult]               ->Divide(hREC_2D_in_MT[iMult],         hGEN_2D_in_MT[iMult],          1.,1.,"b");
        hREC_1D_in_MT[iMult]->Scale(1.,"width");
        hGEN_1D_in_MT[iMult]->Scale(1.,"width");
        hTRU_1D_in_MT[iMult]->Scale(1.,"width");
        hREC_1D_in_MT_in_2Dbin[iMult]->Scale(1.,"width");
        hGEN_1D_in_MT_in_2Dbin[iMult]->Scale(1.,"width");
        hTRU_1D_in_MT_in_2Dbin[iMult]->Scale(1.,"width");
        hREC_2D_in_MT[iMult]->Scale(1.,"width");
        hGEN_2D_in_MT[iMult]->Scale(1.,"width");
        hTRU_2D_in_MT[iMult]->Scale(1.,"width");
        hREC_1D_in_MT[iMult]->Scale(1./fNormEvent);
        hGEN_1D_in_MT[iMult]->Scale(1./fNormEvent);
        hTRU_1D_in_MT[iMult]->Scale(1./fNormEvent);
        hREC_1D_in_MT_in_2Dbin[iMult]->Scale(1./fNormEvent);
        hGEN_1D_in_MT_in_2Dbin[iMult]->Scale(1./fNormEvent);
        hTRU_1D_in_MT_in_2Dbin[iMult]->Scale(1./fNormEvent);
        hREC_2D_in_MT[iMult]->Scale(1./fNormEvent);
        hGEN_2D_in_MT[iMult]->Scale(1./fNormEvent);
        hTRU_2D_in_MT[iMult]->Scale(1./fNormEvent);
    }
    //
    
    // >> RAPIDITY ANALYSIS //
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hEFF_1D_in_RP[iRap]               ->Divide(hREC_1D_in_RP[iRap],         hGEN_1D_in_RP[iRap],          1.,1.,"b");
        hEFF_1D_in_RP_in_2Dbin[iRap]      ->Divide(hREC_1D_in_RP_in_2Dbin[iRap],hGEN_1D_in_RP_in_2Dbin[iRap], 1.,1.,"b");
        hEFF_2D_in_RP[iRap]               ->Divide(hREC_2D_in_RP[iRap],         hGEN_2D_in_RP[iRap],          1.,1.,"b");
        hREC_1D_in_RP[iRap]->Scale(1.,"width");
        hGEN_1D_in_RP[iRap]->Scale(1.,"width");
        hTRU_1D_in_RP[iRap]->Scale(1.,"width");
        hREC_1D_in_RP_in_2Dbin[iRap]->Scale(1.,"width");
        hGEN_1D_in_RP_in_2Dbin[iRap]->Scale(1.,"width");
        hTRU_1D_in_RP_in_2Dbin[iRap]->Scale(1.,"width");
        hREC_2D_in_RP[iRap]->Scale(1.,"width");
        hGEN_2D_in_RP[iRap]->Scale(1.,"width");
        hTRU_2D_in_RP[iRap]->Scale(1.,"width");
        hREC_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hGEN_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hTRU_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hREC_1D_in_RP_in_2Dbin[iRap]->Scale(1./fNormEvent);
        hGEN_1D_in_RP_in_2Dbin[iRap]->Scale(1./fNormEvent);
        hTRU_1D_in_RP_in_2Dbin[iRap]->Scale(1./fNormEvent);
        hREC_2D_in_RP[iRap]->Scale(1./fNormEvent);
        hGEN_2D_in_RP[iRap]->Scale(1./fNormEvent);
        hTRU_2D_in_RP[iRap]->Scale(1./fNormEvent);
    }
    //
    
    // >> TRIGGER ANALYSIS //
    
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> Trigger Analysis
    //
    if ( kDoTrigger )   {
        TFile *outFil1  =   new TFile   (fTrgPrePrMC,"recreate");
        //
        outFil1->Close();
    }
    //
    // >> Yield Analysis
    //
    if ( kDoYield || kDoRapidity ) {
        TFile *outFil2  =   new TFile   (fYldPrePrMC,"recreate");
        //
        hREC_1D->Write();
        hGEN_1D->Write();
        hTRU_1D->Write();
        hTRU_ALL_1D->Write();
        hEFF_1D->Write();
        hREC_1D_in_2Dbin->Write();
        hGEN_1D_in_2Dbin->Write();
        hTRU_1D_in_2Dbin->Write();
        hEFF_1D_in_2Dbin->Write();
        hREC_2D->Write();
        hGEN_2D->Write();
        hTRU_2D->Write();
        hTRU_ALL_2D->Write();
        hEFF_2D->Write();
        //
        outFil2->Close();
    }
    //
    // >> Multiplicity Analysis
    //
    if ( kDoMultiplicity )  {
        TFile *outFil3  =   new TFile   (fMltPrePrMC,"recreate");
        //
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
        {
            hREC_1D_in_MT[iMult]->Write();
            hGEN_1D_in_MT[iMult]->Write();
            hTRU_1D_in_MT[iMult]->Write();
            hEFF_1D_in_MT[iMult]->Write();
            hREC_1D_in_MT_in_2Dbin[iMult]->Write();
            hGEN_1D_in_MT_in_2Dbin[iMult]->Write();
            hTRU_1D_in_MT_in_2Dbin[iMult]->Write();
            hEFF_1D_in_MT_in_2Dbin[iMult]->Write();
            hREC_2D_in_MT[iMult]->Write();
            hGEN_2D_in_MT[iMult]->Write();
            hTRU_2D_in_MT[iMult]->Write();
            hEFF_2D_in_MT[iMult]->Write();
        }
        //
        outFil3->Close();
    }
    //
    if ( kDoRapidity )  {
        TFile *outFil4  =   new TFile   (fRapPrePrMC,"recreate");
        //
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
        {
            hREC_1D_in_RP[iRap]->Write();
            hGEN_1D_in_RP[iRap]->Write();
            hTRU_1D_in_RP[iRap]->Write();
            hEFF_1D_in_RP[iRap]->Write();
            hREC_1D_in_RP_in_2Dbin[iRap]->Write();
            hGEN_1D_in_RP_in_2Dbin[iRap]->Write();
            hTRU_1D_in_RP_in_2Dbin[iRap]->Write();
            hEFF_1D_in_RP_in_2Dbin[iRap]->Write();
            hREC_2D_in_RP[iRap]->Write();
            hGEN_2D_in_RP[iRap]->Write();
            hTRU_2D_in_RP[iRap]->Write();
            hEFF_2D_in_RP[iRap]->Write();
        }
        //
        outFil4->Close();
    }
    //
    // >-> Close input File
    //
    insFileMC->Close();
    //
}
