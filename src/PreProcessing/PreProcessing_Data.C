#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void PreProcessing_Data ( string fFileName = "", TString fOption = "", Int_t nEventsCut = -1, TString kFolder = "" )
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
    fChooseOption(fOption);
    
    // Retrieving Event data
    TFile *insFileDT        =   new TFile   (fFileName.c_str());
    
    // Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileDT       ->Get(Form("%s%s",fPhiCandidate_Tree,""));
    TTree   *TKaonCandidate =   (TTree*)insFileDT       ->Get(Form("%s%s",fKaonCandidate_Tree,""));
    
    // Retrieving Event Count Histogram
    TList  *fQCOutputList   =   (TList*)insFileDT       ->Get("fQCOutputList");
    TH1D   *fHEventCount    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enum_FLL");
    TH1D   *fHEvCountMlt    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enum_V0M");
    
    // Define tree data structures
    Struct_PhiCandidate     evPhiCandidate;
    Struct_KaonCandidate    evKaonCandidate;
    
    // Setting the input Candidates in the Trees
    if ( !fSetCandidates(TPhiCandidate,evPhiCandidate,TKaonCandidate,evKaonCandidate) )     return;
    

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
    Int_t       U_AccCand[1024],    U_AccCan2[1024];
    Int_t       U_nAccept,          U_nAccep2;
    
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
    hTitle      =   Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c",fMinPT1D,fMaxPT1D);
    hREC_1D     =   new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
    //
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
    {
        hName = Form("hREC_1D_in_PT_%i",iPT1D);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]);
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
    hTitle      =   Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D);
    hREC_2D     =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
    //
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )
    {
        hName = Form("hREC_1D_in_PT_2D_bin_%i",iPT2D);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]);
        hREC_1D_in_PT_2D_bin[iPT2D]   = new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
        SetAxis(hREC_1D_in_PT_2D_bin[iPT2D],"IM 1D");
        
        hREC_2D_in_PT[iPT2D]    = new TH2F *    [nBinPT2D];
        
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i",iPT2D,jPT2D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrPT2D[jPT2D],fArrPT2D[jPT2D+1]);
            hREC_2D_in_PT[iPT2D][jPT2D]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            SetAxis(hREC_2D_in_PT[iPT2D][jPT2D],"IM 2D");
        }
    }

    // >-> MULTIPLICITY ANALYSIS //
    
    // >->-->-> General analysis //
    //
    TH1D       *fHEventCount_in_MT;
    TH1D       *fHEventCount_PhiCandidate_General_in_MT;
    TH1D       *fHEventCount_PhiCandidate_Accepted_in_MT;
    //
    fHEventCount_in_MT                          =   new TH1D("fHEventCount_in_MT",                          "fHEventCount_in_MT",                       nBinMult,   fArrMult);
    fHEventCount_PhiCandidate_General_in_MT     =   new TH1D("fHEventCount_PhiCandidate_General_in_MT",     "fHEventCount_PhiCandidate_General_in_MT",  nBinMult,   fArrMult);
    fHEventCount_PhiCandidate_Accepted_in_MT    =   new TH1D("fHEventCount_PhiCandidate_Accepted_in_MT",    "fHEventCount_PhiCandidate_Accepted_in_MT", nBinMult,   fArrMult);
    //
    
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_MT               = new TH1F     *[nBinMult+1];
    TH1F     ***hREC_1D_in_MT_in_PT         = new TH1F    **[nBinMult+1];
    //
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    hName = Form("hREC_1D_in_MT_%i",0);
    hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fMinPT1D,fMaxPT1D,fArrMult[0],fArrMult[nBinMult]);
    hREC_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
    SetAxis(hREC_1D_in_MT[0],"IM 1D");

    hREC_1D_in_MT_in_PT[0] = new TH1F     *[nBinPT1D];

    for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
    {
        hName = Form("hREC_1D_in_MT_PT_%i_%i",0,iPT1D);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1],fArrMult[0],fArrMult[nBinMult]);
        hREC_1D_in_MT_in_PT[0][iPT1D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_MT_in_PT[0][iPT1D],"IM 1D");
    }
    
    for ( Int_t iMult = 1; iMult <= nBinMult; iMult++ )
    {
        hName = Form("hREC_1D_in_MT_%i",iMult);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fMinPT1D,fMaxPT1D,fArrMult[iMult-1],fArrMult[iMult]);
        hREC_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_MT[iMult],"IM 1D");

        hREC_1D_in_MT_in_PT[iMult] = new TH1F     *[nBinPT1D];

        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
        {
            hName = Form("hREC_1D_in_MT_PT_%i_%i",iMult,iPT1D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1],fArrMult[iMult-1],fArrMult[iMult]);
            hREC_1D_in_MT_in_PT[iMult][iPT1D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_MT_in_PT[iMult][iPT1D],"IM 1D");
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
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    hName   =   Form("hREC_2D_in_MT_%i",0);
    hTitle  =   Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,fArrMult[0],fArrMult[nBinMult]);
    hREC_2D_in_MT[0]   =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);

    hREC_1D_in_MT_in_PT_2D_bin[0]  = new TH1F     *[nBinPT2D];
    hREC_2D_in_MT_in_PT[0]         = new TH2F    **[nBinPT2D];
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )
    {
        hName = Form("hREC_1D_in_PT_2D_bin_%i_in_MT_%i",iPT2D,0);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrMult[0],fArrMult[nBinMult]);
        hREC_1D_in_MT_in_PT_2D_bin[0][iPT2D]  = new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
        SetAxis(hREC_1D_in_MT_in_PT_2D_bin[0][iPT2D],"IM 1D");
        
        hREC_2D_in_MT_in_PT[0][iPT2D]         = new TH2F     *[nBinPT2D];
        
        for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i_in_MT_%i",iPT2D,jPT2D,0);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrPT2D[jPT2D],fArrPT2D[jPT2D+1],fArrMult[0],fArrMult[nBinMult]);
            hREC_2D_in_MT_in_PT[0][iPT2D][jPT2D]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            SetAxis(hREC_2D_in_MT_in_PT[0][iPT2D][jPT2D],"IM 2D");
        }
    }
    
    for ( Int_t iMult = 1; iMult <= nBinMult; iMult++ )
    {
        hName   =   Form("hREC_2D_in_MT_%i",iMult);
        hTitle  =   Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,fArrMult[iMult-1],fArrMult[iMult]);
        hREC_2D_in_MT[iMult]   =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);

        hREC_1D_in_MT_in_PT_2D_bin[iMult]  = new TH1F     *[nBinPT2D];
        hREC_2D_in_MT_in_PT[iMult]         = new TH2F    **[nBinPT2D];
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_in_MT_%i",iPT2D,iMult);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrMult[iMult-1],fArrMult[iMult]);
            hREC_1D_in_MT_in_PT_2D_bin[iMult][iPT2D]  = new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
            SetAxis(hREC_1D_in_MT_in_PT_2D_bin[iMult][iPT2D],"IM 1D");
            
            hREC_2D_in_MT_in_PT[iMult][iPT2D]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_in_MT_%i",iPT2D,jPT2D,iMult);
                hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c, Multiplicity [%.2f#;%.2f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrPT2D[jPT2D],fArrPT2D[jPT2D+1],fArrMult[iMult-1],fArrMult[iMult]);
                hREC_2D_in_MT_in_PT[iMult][iPT2D][jPT2D]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
                SetAxis(hREC_2D_in_MT_in_PT[iMult][iPT2D][jPT2D],"IM 2D");
            }
        }
    }
    //
    // >-> RAPIDITY ANALYSIS //
    
    // >->-->-> General analysis //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_RP               = new TH1F     *[nBinRap_];
    TH1F     ***hREC_1D_in_RP_in_PT         = new TH1F    **[nBinRap_];
    //
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName = Form("hREC_1D_in_RP_%i",iRap);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fMinPT1D,fMaxPT1D,fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_1D_in_RP[iRap]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_RP[iRap],"IM 1D");

        hREC_1D_in_RP_in_PT[iRap] = new TH1F     *[nBinPT1D];

        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
        {
            hName = Form("hREC_1D_in_RP_PT_%i_%i",iRap,iPT1D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1],fArrRap_[iRap],fArrRap_[iRap+1]);
            hREC_1D_in_RP_in_PT[iRap][iPT1D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_RP_in_PT[iRap][iPT1D],"IM 1D");
        }
    }
    
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_KP               = new TH1F     *[nBinRap_];
    TH1F     ***hREC_1D_in_KP_in_PT         = new TH1F    **[nBinRap_];
    //
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName = Form("hREC_1D_in_KP_%i",iRap);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fMinPT1D,fMaxPT1D,fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_1D_in_KP[iRap]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_KP[iRap],"IM 1D");

        hREC_1D_in_KP_in_PT[iRap] = new TH1F     *[nBinPT1D];

        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )
        {
            hName = Form("hREC_1D_in_KP_PT_%i_%i",iRap,iPT1D);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fArrPT1D[iPT1D],fArrPT1D[iPT1D+1],fArrRap_[iRap],fArrRap_[iRap+1]);
            hREC_1D_in_KP_in_PT[iRap][iPT1D]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
            SetAxis(hREC_1D_in_KP_in_PT[iRap][iPT1D],"IM 1D");
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
    //  Defining cumulative histogram over measurable pT
    //  Defining pT-Differential histograms over measurable pT
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hName   =   Form("hREC_2D_in_RP_%i",iRap);
        hTitle  =   Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_2D_in_RP[iRap]   =   new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);

        hREC_1D_in_RP_in_PT_2D_bin[iRap]  = new TH1F     *[nBinPT2D];
        hREC_2D_in_RP_in_PT[iRap]         = new TH2F    **[nBinPT2D];
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )
        {
            hName = Form("hREC_1D_in_PT_2D_bin_%i_in_RP_%i",iPT2D,iRap);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrRap_[iRap],fArrRap_[iRap+1]);
            hREC_1D_in_RP_in_PT_2D_bin[iRap][iPT2D]  = new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
            SetAxis(hREC_1D_in_RP_in_PT_2D_bin[iRap][iPT2D],"IM 1D");
            
            hREC_2D_in_RP_in_PT[iRap][iPT2D]         = new TH2F     *[nBinPT2D];
            
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ )
            {
                hName = Form("hREC_2D_in_PT_%i_%i_in_RP_%i",iPT2D,jPT2D,iRap);
                hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.2f#;%.2f] GeV/c and [%.2f#;%.2f] GeV/c, Rapidity [%.2f#;%.2f]",fArrPT2D[iPT2D],fArrPT2D[iPT2D+1],fArrPT2D[jPT2D],fArrPT2D[jPT2D+1],fArrRap_[iRap],fArrRap_[iRap+1]);
                hREC_2D_in_RP_in_PT[iRap][iPT2D][jPT2D]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
                SetAxis(hREC_2D_in_RP_in_PT[iRap][iPT2D][jPT2D],"IM 2D");
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
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    // Evaluating entries and saving them for later
    Int_t nEvents = (!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : nEventsCut);
    
    if ( nEvents > 0 )  fStartTimer("Phi Analysis");
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Phi Analysis",iEvent,nEvents,kPrintIntervalPP);
        
        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        Bool_t  fCheckFill1 =   false;
        Bool_t  fCheckFill2 =   false;
        Bool_t  fCheckFill3 =   false;
        Bool_t  fCheckFill4 =   false;
        
        // Discarding Pile-up events
        //if ( kDoYield           &&  fCheckMask(evPhiCandidate.EventMask,1) ) continue;
        //if ( kDoMultiplicity    &&  ( fCheckMask(evPhiCandidate.EventMask,1) || fCheckMask(evPhiCandidate.EventMask,2) )) continue;
        
        fHEventCount_PhiCandidate_General_in_MT->Fill(evPhiCandidate.Multiplicity);
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
            evPhiCandidate.kHasRap[iPhi]=   fabs(LPhi_candidate1.Rapidity()) <0.5;
            U_nAccept++;
        }
        //
        if ( U_nAccept == 0 )       hTriggerEvt->Fill(0);
        //
        //__________ORDERING by PT
        auto kContinue = kTRUE;
        if ( U_nAccept > 1 )    {
            while ( kContinue ) {
                auto fN = 0;
                for ( Int_t iPhi = 0; iPhi < U_nAccept-1; iPhi++ )    {
                    auto fPT1 = evPhiCandidate.pT[U_AccCand[iPhi]];
                    auto fPT2 = evPhiCandidate.pT[U_AccCand[iPhi+1]];
                    if ( fPT2 >= fPT1 )  {
                        fN++;
                        continue;
                    }
                    auto fN1            =   U_AccCand[iPhi];
                    auto fN2            =   U_AccCand[iPhi+1];
                    U_AccCand[iPhi]     =   fN2;
                    U_AccCand[iPhi+1]   =   fN1;
                }
                if ( fN ==  U_nAccept-1 ) kContinue = kFALSE;
            }
        }
        //
        // >>-->> Utilities
        //
        // >>-->>-->> Multiplicity
        //
        Int_t   iMult               =   fGetBinMult(evPhiCandidate.Multiplicity);
        evPhiCandidate.kHasMult     =   iMult != -1;
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
            Float_t iRapidity           =   evPhiCandidate.Rap[U_AccCand[iPhi]];
            Int_t   iRap                =   evPhiCandidate.iRap[U_AccCand[iPhi]];
            Bool_t  fHasRapidity        =   evPhiCandidate.kHasRap[U_AccCand[iPhi]];
            //
            // >->-->-> Trigger
            //
            if ( !fCheckFill1 ) {
                fHEventCount_PhiCandidate_Accepted_in_MT->Fill(evPhiCandidate.Multiplicity);
                hTriggerEvt->Fill(1);
                fCheckFill1 = true;
            }
            hTriggerEvt1D                                       ->  Fill(evPhiCandidate.pT[U_AccCand[iPhi]]);
            //
            if  ( fabs(iRapidity) < .5 )    {
            //
            // >->-->-> Rapidity
            //
                if ( kDoRapidity )  {
                    hREC_1D_in_RP[iRap]                             ->  Fill(iInvarMass);
                    hREC_1D_in_RP_in_PT[iRap][iPT1D]                ->  Fill(iInvarMass);
                    hREC_1D_in_RP_in_PT_2D_bin[iRap][iPT2D]         ->  Fill(iInvarMass);
                }
            //
            // >->-->-> Yield
            //
                if ( kDoYield ) {
                    hREC_1D                                         ->  Fill(iInvarMass);
                    hREC_1D_in_PT[iPT1D]                            ->  Fill(iInvarMass);
                    hREC_1D_in_PT_2D_bin[iPT2D]                     ->  Fill(iInvarMass);
                }
            //
            // >->-->-> Multiplicity
            //
                if ( evPhiCandidate.kHasMult && kDoMultiplicity )  {
                    hREC_1D_in_MT[iMult+1]                          ->  Fill(iInvarMass);
                    hREC_1D_in_MT_in_PT[iMult+1][iPT1D]             ->  Fill(iInvarMass);
                    hREC_1D_in_MT_in_PT_2D_bin[iMult+1][iPT2D]      ->  Fill(iInvarMass);
                    hREC_1D_in_MT[0]                                ->  Fill(iInvarMass);
                    hREC_1D_in_MT_in_PT[0][iPT1D]                   ->  Fill(iInvarMass);
                    hREC_1D_in_MT_in_PT_2D_bin[0][iPT2D]            ->  Fill(iInvarMass);
                }
            }
            //
            for ( Int_t jPhi = iPhi+1; jPhi < U_nAccept; jPhi++ )    {
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
                Float_t jRapidity           =   evPhiCandidate.Rap[U_AccCand[jPhi]];
                Int_t   jRap                =   evPhiCandidate.iRap[U_AccCand[jPhi]];
                        fHasRapidity        =   evPhiCandidate.kHasRap[U_AccCand[iPhi]] && evPhiCandidate.kHasRap[U_AccCand[jPhi]];
                Int_t   ijRap               =   fGetBinRap_(fabs(evPhiCandidate.Rap[U_AccCand[iPhi]]-evPhiCandidate.Rap[U_AccCand[jPhi]]));
                Bool_t  fHasDiffRap         =   ijRap != -1;
                //
                // >->-->-> Trigger
                //
                if ( !fCheckFill2 ) {
                    hTriggerEvt->Fill(2);
                    fCheckFill2 = true;
                }
                hTriggerEvt2D                                       ->  Fill(evPhiCandidate.pT[U_AccCand[iPhi]],evPhiCandidate.pT[U_AccCand[jPhi]],0.5);
                //
                if  ( fabs(iRapidity) < .5 && fabs(jRapidity) < .5 )    {
                    if  ( fHasDiffRap && kDoRapidity )     {
                //
                // >->-->-> Rapidity
                //
                        hREC_2D_in_RP[ijRap]                            ->  Fill(iInvarMass,jInvarMass,0.5);
                        hREC_2D_in_RP_in_PT[ijRap][iPT2D][jPT2D]        ->  Fill(iInvarMass,jInvarMass,0.5);
                    }
                //
                     if (kDoYield )    {
                //
                // >->-->-> Yield
                //
                         hREC_2D                                    ->  Fill(iInvarMass,jInvarMass);
                         hREC_2D_in_PT[iPT2D][jPT2D]                ->  Fill(iInvarMass,jInvarMass);
                     }
                //
                // >->-->-> Multiplicity
                //
                    if ( evPhiCandidate.kHasMult && kDoMultiplicity )  {
                        hREC_2D_in_MT[iMult+1]                      ->  Fill(iInvarMass,jInvarMass);
                        hREC_2D_in_MT_in_PT[iMult+1][iPT2D][jPT2D]  ->  Fill(iInvarMass,jInvarMass);
                        hREC_2D_in_MT[0]                            ->  Fill(iInvarMass,jInvarMass);
                        hREC_2D_in_MT_in_PT[0][iPT2D][jPT2D]        ->  Fill(iInvarMass,jInvarMass);
                    }
                }
                //
                if ( !kDoTrigger ) continue;
                //
                for ( Int_t kPhi = jPhi+1; kPhi < U_nAccept; kPhi++ )    {
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

                    for ( Int_t lPhi = kPhi+1; lPhi < U_nAccept; lPhi++ )    {
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

    if ( nEvents > 0 )  fStartTimer("Kaon Analysis");
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TKaonCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Kaon Analysis",iEvent,nEvents,kPrintIntervalPP);
        
        // Utilities
        TLorentzVector  LKaon1,    LKaon2,  LKaon3,    LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        U_nAccep2 = 0;
        
        for ( Int_t iKaon = 0; iKaon < evKaonCandidate.nKaon; iKaon++ )  {
            LKaon1.SetXYZM(evKaonCandidate.Px[iKaon],evKaonCandidate.Py[iKaon],evKaonCandidate.Pz[iKaon],kKaonMass);
            for ( Int_t jKaon = 0; jKaon < evKaonCandidate.nKaon; jKaon++ )  {
                LKaon2.SetXYZM(evKaonCandidate.Px[jKaon],evKaonCandidate.Py[jKaon],evKaonCandidate.Pz[jKaon],kKaonMass);
                for ( Int_t kKaon = 0; kKaon < evKaonCandidate.nKaon; kKaon++ )  {
                    LKaon3.SetXYZM(evKaonCandidate.Px[kKaon],evKaonCandidate.Py[kKaon],evKaonCandidate.Pz[kKaon],kKaonMass);
                    if ( evKaonCandidate.Charge[iKaon] == evKaonCandidate.Charge[jKaon] ) continue;
                    LPhi_candidate1 =   LKaon1+LKaon2;
                    if ( !fAcceptCandidate(LPhi_candidate1.Mag(),LPhi_candidate1.Pt()) ) continue;
                    if ( fGetBinRap_(LPhi_candidate1.Rapidity() - LKaon3.Rapidity()) == -1 ) continue;
                    // >> 1-Dimensional Analysis Fill
                    //
                    // >>-->> Utilities
                    //
                    // >>-->>-->> True Phis
                    //
                    Int_t   iPT1D               =   fGetBinPT1D(LPhi_candidate1.Pt());
                    //Int_t   iPT2D               =   evPhiCandidate.iPT2D[U_AccCand[iPhi]];
                    Float_t iInvarMass          =   LPhi_candidate1.Mag();
                    Float_t iRapidity           =   LPhi_candidate1.Rapidity()-LKaon3.Rapidity() ;
                    Int_t   iRap                =   fGetBinRap_(LPhi_candidate1.Rapidity()-LKaon3.Rapidity() );
                    
                    hREC_1D_in_KP[iRap]->               Fill(iInvarMass);
                    hREC_1D_in_KP_in_PT[iRap][iPT1D]->  Fill(iInvarMass);
                }
            }
        }
    }
    //
    if ( nEvents > 0 )  fStopTimer("Kaon Analysis");
    
    /*
    Struct_KaonCandidate    evKaonCandidateEvMix1;
    Struct_KaonCandidate    evKaonCandidateEvMix2;
    Struct_KaonCandidate    evKaonCandidateEvMix3;
    Struct_PhiCandidate     evPhiCandidate_EvMix1;
    Struct_PhiCandidate     evPhiCandidate_EvMix2;
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Recovering events
        TKaonCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Kaon Analysis",iEvent,nEvents,kPrintIntervalPP);

        //First events
        if ( iEvent == 0 )  {
            evKaonCandidateEvMix1 = evKaonCandidate;
            continue;
        }
        if ( iEvent == 1 )  {
            evKaonCandidateEvMix2 = evKaonCandidateEvMix1;
            evKaonCandidateEvMix1 = evKaonCandidate;
            continue;
        }
        if ( iEvent == 2 )  {
            evKaonCandidateEvMix3 = evKaonCandidateEvMix2;
            evKaonCandidateEvMix2 = evKaonCandidateEvMix1;
            evKaonCandidateEvMix1 = evKaonCandidate;
            continue;
        }
        // Utilities
        TLorentzVector  LKaon1,    LKaon2,  LKaon3,    LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        U_nAccep2 = 0;
        evPhiCandidate_EvMix1.nPhi = 0;
        evPhiCandidate_EvMix2.nPhi = 0;
        
        for ( Int_t iKaon = 0; iKaon < evKaonCandidate.nKaon; iKaon++ )  {
            LKaon1.SetXYZM(evKaonCandidate.Px[iKaon],evKaonCandidate.Py[iKaon],evKaonCandidate.Pz[iKaon],kKaonMass);
            for ( Int_t jKaon = 0; jKaon < evKaonCandidateEvMix1.nKaon; jKaon++ )  {
                if ( evKaonCandidate.Charge[iKaon] == evKaonCandidateEvMix1.Charge[jKaon] ) continue;
                LKaon2.SetXYZM(evKaonCandidateEvMix1.Px[jKaon],evKaonCandidateEvMix1.Py[jKaon],evKaonCandidateEvMix1.Pz[jKaon],kKaonMass);
                LPhi_candidate1 =   LKaon1  +   LKaon2;
                evPhiCandidate_EvMix1.Px[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Px();
                evPhiCandidate_EvMix1.Py[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Py();
                evPhiCandidate_EvMix1.Pz[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Pz();
                evPhiCandidate_EvMix1.pT[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Pt();
                evPhiCandidate_EvMix1.iKaon[evPhiCandidate_EvMix1.nPhi]   =   iKaon;
                evPhiCandidate_EvMix1.jKaon[evPhiCandidate_EvMix1.nPhi]   =   jKaon;
                evPhiCandidate_EvMix1.InvMass[evPhiCandidate_EvMix1.nPhi] =   LPhi_candidate1.Mag();
                evPhiCandidate_EvMix1.Rap[evPhiCandidate_EvMix1.nPhi]     =   LPhi_candidate1.Rapidity();
                evPhiCandidate_EvMix1.nPhi++;
            }
            for ( Int_t jKaon = 0; jKaon < evKaonCandidateEvMix2.nKaon; jKaon++ )  {
                if ( evKaonCandidate.Charge[iKaon] == evKaonCandidateEvMix2.Charge[jKaon] ) continue;
                LKaon2.SetXYZM(evKaonCandidateEvMix2.Px[jKaon],evKaonCandidateEvMix2.Py[jKaon],evKaonCandidateEvMix2.Pz[jKaon],kKaonMass);
                LPhi_candidate1 =   LKaon1  +   LKaon2;
                evPhiCandidate_EvMix1.Px[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Px();
                evPhiCandidate_EvMix1.Py[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Py();
                evPhiCandidate_EvMix1.Pz[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Pz();
                evPhiCandidate_EvMix1.pT[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Pt();
                evPhiCandidate_EvMix1.iKaon[evPhiCandidate_EvMix1.nPhi]   =   iKaon;
                evPhiCandidate_EvMix1.jKaon[evPhiCandidate_EvMix1.nPhi]   =   jKaon;
                evPhiCandidate_EvMix1.InvMass[evPhiCandidate_EvMix1.nPhi] =   LPhi_candidate1.Mag();
                evPhiCandidate_EvMix1.Rap[evPhiCandidate_EvMix1.nPhi]     =   LPhi_candidate1.Rapidity();
                evPhiCandidate_EvMix1.nPhi++;
            }
            for ( Int_t jKaon = 0; jKaon < evKaonCandidateEvMix3.nKaon; jKaon++ )  {
                if ( evKaonCandidate.Charge[iKaon] == evKaonCandidateEvMix3.Charge[jKaon] ) continue;
                LKaon2.SetXYZM(evKaonCandidateEvMix3.Px[jKaon],evKaonCandidateEvMix3.Py[jKaon],evKaonCandidateEvMix3.Pz[jKaon],kKaonMass);
                LPhi_candidate1 =   LKaon1  +   LKaon2;
                evPhiCandidate_EvMix1.Px[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Px();
                evPhiCandidate_EvMix1.Py[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Py();
                evPhiCandidate_EvMix1.Pz[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Pz();
                evPhiCandidate_EvMix1.pT[evPhiCandidate_EvMix1.nPhi]      =   LPhi_candidate1.Pt();
                evPhiCandidate_EvMix1.iKaon[evPhiCandidate_EvMix1.nPhi]   =   iKaon;
                evPhiCandidate_EvMix1.jKaon[evPhiCandidate_EvMix1.nPhi]   =   jKaon;
                evPhiCandidate_EvMix1.InvMass[evPhiCandidate_EvMix1.nPhi] =   LPhi_candidate1.Mag();
                evPhiCandidate_EvMix1.Rap[evPhiCandidate_EvMix1.nPhi]     =   LPhi_candidate1.Rapidity();
                evPhiCandidate_EvMix1.nPhi++;
            }
        }
        for ( Int_t iKaon = 0; iKaon < evKaonCandidateEvMix1.nKaon; iKaon++ )  {
            LKaon1.SetXYZM(evKaonCandidateEvMix1.Px[iKaon],evKaonCandidateEvMix1.Py[iKaon],evKaonCandidateEvMix1.Pz[iKaon],kKaonMass);
            for ( Int_t jKaon = 0; jKaon < evKaonCandidateEvMix1.nKaon; jKaon++ )  {
                if ( evKaonCandidateEvMix1.Charge[iKaon] == evKaonCandidateEvMix1.Charge[jKaon] ) continue;
                LKaon2.SetXYZM(evKaonCandidateEvMix1.Px[jKaon],evKaonCandidateEvMix1.Py[jKaon],evKaonCandidateEvMix1.Pz[jKaon],kKaonMass);
                LPhi_candidate1 =   LKaon1  +   LKaon2;
                evPhiCandidate_EvMix2.Px[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Px();
                evPhiCandidate_EvMix2.Py[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Py();
                evPhiCandidate_EvMix2.Pz[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Pz();
                evPhiCandidate_EvMix2.pT[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Pt();
                evPhiCandidate_EvMix2.iKaon[evPhiCandidate_EvMix2.nPhi]   =   iKaon;
                evPhiCandidate_EvMix2.jKaon[evPhiCandidate_EvMix2.nPhi]   =   jKaon;
                evPhiCandidate_EvMix2.InvMass[evPhiCandidate_EvMix2.nPhi] =   LPhi_candidate1.Mag();
                evPhiCandidate_EvMix2.Rap[evPhiCandidate_EvMix2.nPhi]     =   LPhi_candidate1.Rapidity();
                evPhiCandidate_EvMix2.nPhi++;
            }
        }
        for ( Int_t iKaon = 0; iKaon < evKaonCandidateEvMix2.nKaon; iKaon++ )  {
            LKaon1.SetXYZM(evKaonCandidateEvMix2.Px[iKaon],evKaonCandidateEvMix2.Py[iKaon],evKaonCandidateEvMix2.Pz[iKaon],kKaonMass);
            for ( Int_t jKaon = 0; jKaon < evKaonCandidateEvMix2.nKaon; jKaon++ )  {
                if ( evKaonCandidateEvMix2.Charge[iKaon] == evKaonCandidateEvMix2.Charge[jKaon] ) continue;
                LKaon2.SetXYZM(evKaonCandidateEvMix2.Px[jKaon],evKaonCandidateEvMix2.Py[jKaon],evKaonCandidateEvMix2.Pz[jKaon],kKaonMass);
                LPhi_candidate1 =   LKaon1  +   LKaon2;
                evPhiCandidate_EvMix2.Px[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Px();
                evPhiCandidate_EvMix2.Py[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Py();
                evPhiCandidate_EvMix2.Pz[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Pz();
                evPhiCandidate_EvMix2.pT[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Pt();
                evPhiCandidate_EvMix2.iKaon[evPhiCandidate_EvMix2.nPhi]   =   iKaon;
                evPhiCandidate_EvMix2.jKaon[evPhiCandidate_EvMix2.nPhi]   =   jKaon;
                evPhiCandidate_EvMix2.InvMass[evPhiCandidate_EvMix2.nPhi] =   LPhi_candidate1.Mag();
                evPhiCandidate_EvMix2.Rap[evPhiCandidate_EvMix2.nPhi]     =   LPhi_candidate1.Rapidity();
                evPhiCandidate_EvMix2.nPhi++;
            }
        }
        for ( Int_t iKaon = 0; iKaon < evKaonCandidateEvMix3.nKaon; iKaon++ )  {
            LKaon1.SetXYZM(evKaonCandidateEvMix3.Px[iKaon],evKaonCandidateEvMix3.Py[iKaon],evKaonCandidateEvMix3.Pz[iKaon],kKaonMass);
            for ( Int_t jKaon = 0; jKaon < evKaonCandidateEvMix3.nKaon; jKaon++ )  {
                if ( evKaonCandidateEvMix3.Charge[iKaon] == evKaonCandidateEvMix3.Charge[jKaon] ) continue;
                LKaon2.SetXYZM(evKaonCandidateEvMix3.Px[jKaon],evKaonCandidateEvMix3.Py[jKaon],evKaonCandidateEvMix3.Pz[jKaon],kKaonMass);
                LPhi_candidate1 =   LKaon1  +   LKaon2;
                evPhiCandidate_EvMix2.Px[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Px();
                evPhiCandidate_EvMix2.Py[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Py();
                evPhiCandidate_EvMix2.Pz[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Pz();
                evPhiCandidate_EvMix2.pT[evPhiCandidate_EvMix2.nPhi]      =   LPhi_candidate1.Pt();
                evPhiCandidate_EvMix2.iKaon[evPhiCandidate_EvMix2.nPhi]   =   iKaon;
                evPhiCandidate_EvMix2.jKaon[evPhiCandidate_EvMix2.nPhi]   =   jKaon;
                evPhiCandidate_EvMix2.InvMass[evPhiCandidate_EvMix2.nPhi] =   LPhi_candidate1.Mag();
                evPhiCandidate_EvMix2.Rap[evPhiCandidate_EvMix2.nPhi]     =   LPhi_candidate1.Rapidity();
                evPhiCandidate_EvMix2.nPhi++;
            }
        }
        evKaonCandidateEvMix3 = evKaonCandidateEvMix2;
        evKaonCandidateEvMix2 = evKaonCandidateEvMix1;
        evKaonCandidateEvMix1 = evKaonCandidate;
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate_EvMix1.nPhi; iPhi++ )  {
            if ( !fAcceptCandidate(evPhiCandidate_EvMix1.InvMass[iPhi],evPhiCandidate_EvMix1.pT[iPhi]) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate_EvMix2.nPhi; iPhi++ )  {
            if ( !fAcceptCandidate(evPhiCandidate_EvMix2.InvMass[iPhi],evPhiCandidate_EvMix2.pT[iPhi]) ) continue;
            U_AccCan2[U_nAccep2] = iPhi;
            U_nAccep2++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;
            //
            // >-> 1-Dimensional Analysis Fill   //
            //
            // >->-->-> Utilities
            //
            Float_t iRapidity   = evPhiCandidate_EvMix1.Rap[U_AccCand[iPhi]];
            //
            for ( Int_t jPhi = 0; jPhi < U_nAccep2; jPhi++ )    {
                // Must have at least 1 candidate
                if ( U_nAccep2 < 1 ) break;
                //
                // >-> 2-Dimensional Analysis Fill
                //
                // >->-->-> Utilities
                //
                Float_t jRapidity   = evPhiCandidate_EvMix2.Rap[U_AccCan2[jPhi]];
                //
                // >->-->-> Rapidity
                //
            }
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )    {
                // Must have at least 1 candidate
                if ( U_nAccept < 2 ) break;
                if ( iPhi == jPhi ) continue;
                //
                // >-> 1-Dimensional Analysis Fill   //
                //
                // >->-->-> Utilities
                //
                Float_t jRapidity   = evPhiCandidate_EvMix1.Rap[U_AccCand[jPhi]];
                //
                // >->-->-> Rapidity
                //
            }
        }
    }
    //
    if ( nEvents > 0 )  fStopTimer("Kaon Analysis");
    //
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
    if ( kDoTrigger )  {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Trigger")+kFolder).Data())));
        TFile *outFil1  =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Trigger")+kFolder).Data()),"recreate");
        //
        fHEventCount    ->Write();
        fHEvCountMlt    ->Write();
        hTriggerEvt     ->Write();
        hTriggerEvt1D   ->Write();
        hTriggerEvt2D   ->Write();
        //
        outFil1->Close();
    }
    //
    // >-> Yield Analysis
    //
    if ( kDoYield )  {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data())));
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/")).Data()));
        TFile *outFil2  =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Yield")+kFolder).Data()),"recreate");
        //
        fHEventCount    ->Write();
        fHEvCountMlt    ->Write();
        hREC_1D->Write();
        for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
        {
            hREC_1D_in_PT[iHisto]   ->Write();
        }
        //
        hREC_2D->Write();
        //
        for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)   {
            hREC_1D_in_PT_2D_bin[iHisto]   ->Write();
            for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)   {
                
                hREC_2D_in_PT[iHisto][jHisto]->Write();
                
                if ( hREC_2D_in_PT[iHisto][jHisto]->GetEntries() == 0 ) continue;
                
                gROOT->SetBatch(kTRUE);
                
                SetStyle();
                
                gStyle->SetPadTopMargin(0.2);
                gStyle->SetPadRightMargin(0.18);
                
                TCanvas    *cDrawHisto  =   new TCanvas("cDrawHisto","cDrawHisto",1000,1000);
                //
                //  X axis
                hREC_2D_in_PT[iHisto][jHisto]->GetXaxis()->SetLabelOffset(0.015);
                hREC_2D_in_PT[iHisto][jHisto]->GetXaxis()->SetTitleOffset(1.5);
                //
                //  Y axis
                hREC_2D_in_PT[iHisto][jHisto]->GetYaxis()->SetTitleOffset(1.75);
                //
                //  Z axis
                hREC_2D_in_PT[iHisto][jHisto]->GetZaxis()->SetTitle(Form("Counts/( %.1f MeV/#it{c}^{2} )",1000*kBinningPrecision2D));
                hREC_2D_in_PT[iHisto][jHisto]->GetZaxis()->SetTitleOffset(1.5);
                //
                hREC_2D_in_PT[iHisto][jHisto]->Draw("COLZ");
                
                uLatex->SetTextFont(60);
                uLatex->SetTextSize(0.05);
                uLatex->DrawLatexNDC(0.12, 0.95,"ALICE Performance");
                uLatex->SetTextFont(42);
                uLatex->SetTextSize(0.04);
                uLatex->DrawLatexNDC(0.12, 0.90,"pp #sqrt{#it{s}}= 7 TeV");
                uLatex->DrawLatexNDC(0.12, 0.85,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
                //
                uLatex->DrawLatexNDC(0.50, 0.90,Form("%.2f < #it{p}_{T,#phi_{1}} < %.2f GeV/#it{c}",fArrPT2D[iHisto],fArrPT2D[iHisto+1]));
                uLatex->DrawLatexNDC(0.50, 0.85,Form("%.2f < #it{p}_{T,#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[jHisto],fArrPT2D[jHisto+1]));
                
                cDrawHisto->SaveAs(Form("./result/Yield/PreProcessing/Plots/InvariantMass_%.1f_%.1f_%.1f_%.1f.pdf",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]));
                cDrawHisto->SaveAs(Form("./result/Yield/PreProcessing/Plots/InvariantMass_%.1f_%.1f_%.1f_%.1f.eps",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]));
                delete cDrawHisto;
                
                gROOT->SetBatch(kFALSE);
                
            }
        }
        //
        outFil2->Close();
    }
    //
    // >-> Multiplicity Analysis
    //
    if ( kDoMultiplicity )  {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data())));
        TFile *outFil3  =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        fHEventCount->Write();
        fHEvCountMlt->Write();
        fHEventCount_in_MT                          ->Write();
        fHEventCount_PhiCandidate_General_in_MT     ->Write();
        fHEventCount_PhiCandidate_Accepted_in_MT    ->Write();
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
        {
            hREC_1D_in_MT[iMult]->Write();
            hREC_2D_in_MT[iMult]->Write();

            for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
            {
                hREC_1D_in_MT_in_PT[iMult][jHisto]->Write();
            }
            for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
            {
                hREC_1D_in_MT_in_PT_2D_bin[iMult][jHisto]->Write();
                
                for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
                {
                    hREC_2D_in_MT_in_PT[iMult][jHisto][kHisto]->Write();
                }
            }
        }
        //
        outFil3->Close();
    }
    //
    // >-> Rapidity Analysis
    //
    if ( kDoRapidity )  {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Rapidity")+kFolder).Data())));
        TFile *outFil4  =   new TFile   (Form(kAnalysis_InvMassHist,(TString("Rapidity")+kFolder).Data()),"recreate");
        //
        fHEventCount->Write();
        fHEvCountMlt->Write();
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
        {
            hREC_1D_in_RP[iRap]->Write();
            hREC_1D_in_KP[iRap]->Write();
            hREC_2D_in_RP[iRap]->Write();

            for ( Int_t jHisto = 0; jHisto < nBinPT1D; jHisto++ )
            {
                hREC_1D_in_RP_in_PT[iRap][jHisto]->Write();
                hREC_1D_in_KP_in_PT[iRap][jHisto]->Write();
            }
            for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
            {
                hREC_1D_in_RP_in_PT_2D_bin[iRap][jHisto]->Write();
                
                for ( Int_t kHisto = 0; kHisto < nBinPT2D; kHisto++ )
                {
                    hREC_2D_in_RP_in_PT[iRap][jHisto][kHisto]->Write();
                }
            }
        }
        //
        outFil4->Close();
    }
    //
    fStopTimer("writing output objects to file(s)");
    // >-> Close input File
    //
    insFileDT->Close();
    //
}
