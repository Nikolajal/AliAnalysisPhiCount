#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: All Set!

void PreProcessing_MC ( string fFileName = "", TString fOption = "", Int_t nEventsCut = -1., TString kFolder = "" )   {
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
    TTree   *TPhiEfficiency =   (TTree*)insFileMC->Get(Form("%s%s",fPhiCandidateEff_Tree,"_name"));
    TTree   *TKaonCandidate =   nullptr;//(TTree*)insFileMC->Get(fKaonCandidateEff_Tree);
    TTree   *TPhiCandidate  =   (TTree*)insFileMC->Get(Form("%s%s",fPhiCandidate_Tree,"_name"));
    TTree   *TKaonEfficiency=   nullptr;//(TTree*)insFileMC->Get(fKaonCandidate_Tree);
    
    // Retrieving Event Count Histogram
    TList  *fQCOutputList   =   (TList*)insFileMC       ->Get("fQCOutputList_name");
    TH1D   *fHEventCount    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enum_FLL");
    TH1D   *fHEvCountMlt    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enum_V0M");
    
    // Define tree data structures
    Struct_PhiEfficiency    evPhiEfficiency;
    Struct_KaonEfficiency   evKaonEfficiency;
    Struct_PhiCandidate     evPhiCandidate;
    Struct_KaonCandidate    evKaonCandidate;

    // Setting the input Candidates in the Trees
    if ( !fSetCandidates(TPhiEfficiency,evPhiEfficiency,TKaonEfficiency,evKaonEfficiency) ) return;
    if ( !fSetCandidates(TPhiCandidate,evPhiCandidate,TKaonCandidate,evKaonCandidate) )     return;
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetAllBins();
    auto    fNBinning = (int)((fMaxPT1D-fMinPT1D)*10.);
    Float_t *fUniformBinning100MeV = new Float_t[fNBinning+1];
    fSetUniformBinning(fUniformBinning100MeV,fMinPT1D,fMaxPT1D,fNBinning);
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept,  U_nAccep2;
    
    // Creating the histograms-------------------------------------------------------------------------------

    // >> YIELD ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hProdDistrTRU;
    TH1F       *hProdDistrGEN;
    TH1F       *hProdDistrREC;
    TH1F       *hREC_1D;
    TH1F       *hGEN_1D;
    TH1F       *hREC_Rw_1D;
    TH1F       *hGEN_Rw_1D;
    TH1F       *hREC_IM_1D;
    TH1F       *hGEN_IM_1D;
    TH1F       *hTRU_1D;
    TH1F       *hEFF_1D;
    TH1F       *hEFF_IM_1D;
    TH1F       *hEFF_SL_1D;
    TH1F       *hGEN_INELVTX_1D;
    TH1F       *hGEN_INELFLL_1D;
    TH1F       *hREC_1D_in_2D_bin;
    TH1F       *hGEN_1D_in_2D_bin;
    TH1F       *hTRU_1D_in_2D_bin;
    TH1F       *hEFF_1D_in_2D_bin;
    TH1F       *hEFF_SL_1D_in_2D_bin;
    TH1F       *hGEN_INELVTX_1D_in_2D_bin;
    TH1F       *hGEN_INELFLL_1D_in_2D_bin;
    //
    //  Defining Efficiency and check utilities
    //
    hName       =   Form("hProdDistrTRU");
    hTitle      =   Form("hProdDistrTRU");
    hProdDistrTRU     =   new TH1F (hName,hTitle,10,-.5,9.5);
    SetAxis(hProdDistrTRU,"PT 1D");
    //
    hName       =   Form("hProdDistrGEN");
    hTitle      =   Form("hProdDistrGEN");
    hProdDistrGEN     =   new TH1F (hName,hTitle,10,-.5,9.5);
    SetAxis(hProdDistrGEN,"PT 1D");
    //
    hName       =   Form("hProdDistrREC");
    hTitle      =   Form("hProdDistrREC");
    hProdDistrREC     =   new TH1F (hName,hTitle,10,-.5,9.5);
    SetAxis(hProdDistrREC,"PT 1D");
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
    hName       =   Form("hREC_Rw_1D");
    hTitle      =   Form("hREC_Rw_1D");
    hREC_Rw_1D  =   new TH1F (hName,hTitle,fNBinning,fUniformBinning100MeV);
    SetAxis(hREC_Rw_1D,"PT 1D");
    //
    hName       =   Form("hGEN_Rw_1D");
    hTitle      =   Form("hGEN_Rw_1D");
    hGEN_Rw_1D  =   new TH1F (hName,hTitle,fNBinning,fUniformBinning100MeV);
    SetAxis(hGEN_Rw_1D,"PT 1D");
    //
    hName       =   Form("hREC_IM_1D");
    hTitle      =   Form("hREC_IM_1D");
    hREC_IM_1D  =   new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
    SetAxis(hREC_IM_1D,"PT 1D");
    //
    hName       =   Form("hGEN_IM_1D");
    hTitle      =   Form("hGEN_IM_1D");
    hGEN_IM_1D  =   new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
    SetAxis(hGEN_IM_1D,"PT 1D");
    //
    hName       =   Form("hTRU_1D");
    hTitle      =   Form("hTRU_1D");
    hTRU_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hTRU_1D,"PT 1D");
    //
    hName       =   Form("hGEN_INELVTX_1D");
    hTitle      =   Form("hGEN_INELVTX_1D");
    hGEN_INELVTX_1D =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hGEN_INELVTX_1D,"PT 1D");
    //
    hName       =   Form("hGEN_INELFLL_1D");
    hTitle      =   Form("hGEN_INELFLL_1D");
    hGEN_INELFLL_1D =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hGEN_INELFLL_1D,"PT 1D");
    //
    hName       =   Form("hEFF_1D");
    hTitle      =   Form("hEFF_1D");
    hEFF_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hEFF_1D,"PT 1D");
    //
    hName       =   Form("hEFF_IM_1D");
    hTitle      =   Form("hEFF_IM_1D");
    hEFF_IM_1D     =   new TH1F (hName,hTitle,nBinIM2D,fArrIM2D);
    SetAxis(hEFF_IM_1D,"PT 1D");
    //
    hName       =   Form("hEFF_SL_1D");
    hTitle      =   Form("hEFF_SL_1D");
    hEFF_SL_1D     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hEFF_SL_1D,"PT 1D");
    //
    hName       =   Form("hREC_1D_in_2D_bin");
    hTitle      =   Form("hREC_1D_in_2D_bin");
    hREC_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hREC_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hGEN_1D_in_2D_bin");
    hTitle      =   Form("hGEN_1D_in_2D_bin");
    hGEN_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hTRU_1D_in_2D_bin");
    hTitle      =   Form("hTRU_1D_in_2D_bin");
    hTRU_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hEFF_1D_in_2D_bin");
    hTitle      =   Form("hEFF_1D_in_2D_bin");
    hEFF_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hEFF_SL_1D_in_2D_bin");
    hTitle      =   Form("hEFF_SL_1D_in_2D_bin");
    hEFF_SL_1D_in_2D_bin     =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_SL_1D_in_2D_bin,"PT 1D");
    //
    hName       =   Form("hGEN_INELVTX_1D_in_2D_bin");
    hTitle      =   Form("hGEN_INELVTX_1D_in_2D_bin");
    hGEN_INELVTX_1D_in_2D_bin =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_INELVTX_1D_in_2D_bin,"PT 2D");
    //
    hName       =   Form("hGEN_INELFLL_1D_in_2D_bin");
    hTitle      =   Form("hGEN_INELFLL_1D_in_2D_bin");
    hGEN_INELFLL_1D_in_2D_bin =   new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_INELFLL_1D_in_2D_bin,"PT 2D");
    //
    // >>-->> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hREC_2D;
    TH2F       *hGEN_2D;
    TH2F       *hTRU_2D;
    TH2F       *hGEN_INELVTX_2D;
    TH2F       *hGEN_INELFLL_2D;
    TH2F       *hEFF_2D;
    TH2F       *hEFF_2D_fr_1D;
    TH2F       *hEFF_SL_2D;
    TH2F       *hEFF_SL_2D_fr_1D;
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
    hName       =   Form("hGEN_INELVTX_2D");
    hTitle      =   Form("hGEN_INELVTX_2D");
    hGEN_INELVTX_2D =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_INELVTX_2D,"PT 2D");
    //
    hName       =   Form("hGEN_INELFLL_2D");
    hTitle      =   Form("hGEN_INELFLL_2D");
    hGEN_INELFLL_2D =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_INELFLL_2D,"PT 2D");
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
    hName       =   Form("hEFF_SL_2D");
    hTitle      =   Form("hEFF_SL_2D");
    hEFF_SL_2D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_SL_2D,"PT 2D");
    //
    hName       =   Form("hEFF_SL_2D_fr_1D");
    hTitle      =   Form("hEFF_SL_2D_fr_1D");
    hEFF_SL_2D_fr_1D     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_SL_2D_fr_1D,"PT 2D");
    //

    // >> YIELD ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //>>    Declaring all histograms
    //
    // >>-->>-->>   Mass Resolution
    //
    TH1F  **hMRS_1D;
    TH1F  **hMRS_1D_in_2D_bin;
    TH2F ***hMRS_2D;
    //
    // >>-->>-->>   Mass Distribution
    //
    TH1F  **hMDS_1D;
    TH1F  **hMDS_1D_in_2D_bin;
    TH2F ***hMDS_2D;
    //
    // >>-->>-->>   True Mass Distribution
    //
    TH1F  **hTMD_1D;
    TH1F  **hTMD_1D_in_2D_bin;
    TH2F ***hTMD_2D;
    //
    //>>    Defining all histograms
    //
    hMRS_1D     =   new TH1F   *[nBinPT1D];
    hMDS_1D     =   new TH1F   *[nBinPT1D];
    hTMD_1D     =   new TH1F   *[nBinPT1D];
    for ( Int_t iTer = 0; iTer < nBinPT1D; iTer++ )  {
        hName   =   Form("hMRS_1D_%i",iTer);
        hTitle  =   Form("hMRS_1D_%i",iTer);
        hMRS_1D[iTer]   =   new TH1F(hName,hTitle,nBinIMRs,fArrIMRs);
        hName   =   Form("hMDS_1D_%i",iTer);
        hTitle  =   Form("hMDS_1D_%i",iTer);
        hMDS_1D[iTer]   =   new TH1F(hName,hTitle,nBinIM1D,fArrIM1D);
        hName   =   Form("hTMD_1D_%i",iTer);
        hTitle  =   Form("hTMD_1D_%i",iTer);
        hTMD_1D[iTer]   =   new TH1F(hName,hTitle,nBinIM1D,fArrIM1D);
    }
    //
    hMRS_2D             =   new TH2F  **[nBinPT2D];
    hMDS_2D             =   new TH2F  **[nBinPT2D];
    hTMD_2D             =   new TH2F  **[nBinPT2D];
    hMRS_1D_in_2D_bin   =   new TH1F   *[nBinPT2D];
    hMDS_1D_in_2D_bin   =   new TH1F   *[nBinPT2D];
    hTMD_1D_in_2D_bin   =   new TH1F   *[nBinPT2D];
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )  {
        hName   =   Form("hMRS_1D_in_2D_bin_%i",iTer);
        hTitle  =   Form("hMRS_1D_in_2D_bin_%i",iTer);
        hMRS_1D_in_2D_bin[iTer] =   new TH1F(hName,hTitle,nBinIMRs,fArrIMRs);
        hName   =   Form("hMDS_1D_in_2D_bin_%i",iTer);
        hTitle  =   Form("hMDS_1D_in_2D_bin_%i",iTer);
        hMDS_1D_in_2D_bin[iTer] =   new TH1F(hName,hTitle,nBinIM2D,fArrIM2D);
        hName   =   Form("hTMD_1D_in_2D_bin_%i",iTer);
        hTitle  =   Form("hTMD_1D_in_2D_bin_%i",iTer);
        hTMD_1D_in_2D_bin[iTer] =   new TH1F(hName,hTitle,nBinIM2D,fArrIM2D);
        //
        hMRS_2D[iTer]           =   new TH2F   *[nBinPT2D];
        hMDS_2D[iTer]           =   new TH2F   *[nBinPT2D];
        hTMD_2D[iTer]           =   new TH2F   *[nBinPT2D];
        for ( Int_t jTer = 0; jTer < nBinPT2D; jTer++ )  {
            hName   =   Form("hMRS_2D_%i_%i",iTer,jTer);
            hTitle  =   Form("hMRS_2D_%i_%i",iTer,jTer);
            hMRS_2D[iTer][jTer] =   new TH2F(hName,hTitle,nBinIMR2,fArrIMR2,nBinIMR2,fArrIMR2);
            hName   =   Form("hMDS_2D_%i_%i",iTer,jTer);
            hTitle  =   Form("hMDS_2D_%i_%i",iTer,jTer);
            hMDS_2D[iTer][jTer] =   new TH2F(hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            hName   =   Form("hTMD_2D_%i_%i",iTer,jTer);
            hTitle  =   Form("hTMD_2D_%i_%i",iTer,jTer);
            hTMD_2D[iTer][jTer] =   new TH2F(hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
        }
    }
    //
    
    // >> MULTIPLICITY ANALYSIS //

    // >>-->> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F      **hREC_1D_in_MT                   = new TH1F     *[nBinMult+1];
    TH1F      **hGEN_1D_in_MT                   = new TH1F     *[nBinMult+1];
    TH1F      **hGEN_INELVTX_1D_in_MT           = new TH1F     *[nBinMult+1];
    TH1F      **hTRU_1D_in_MT                   = new TH1F     *[nBinMult+1];
    TH1F      **hEFF_1D_in_MT                   = new TH1F     *[nBinMult+1];
    TH1F      **hEFF_SL_1D_in_MT                = new TH1F     *[nBinMult+1];
    TH1F      **hREC_1D_in_2D_bin_in_MT         = new TH1F     *[nBinMult+1];
    TH1F      **hGEN_1D_in_2D_bin_in_MT         = new TH1F     *[nBinMult+1];
    TH1F      **hGEN_INELVTX_1D_in_2D_bin_in_MT = new TH1F     *[nBinMult+1];
    TH1F      **hTRU_1D_in_2D_bin_in_MT         = new TH1F     *[nBinMult+1];
    TH1F      **hEFF_1D_in_2D_bin_in_MT         = new TH1F     *[nBinMult+1];
    TH1F      **hEFF_SL_1D_in_2D_bin_in_MT      = new TH1F     *[nBinMult+1];
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
    
    hName = Form("hEFF_SL_1D_in_MT_%i",0);
    hTitle= Form("hEFF_SL_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_SL_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hEFF_SL_1D_in_MT[0],"PT 1D");
    
    hName = Form("hGEN_INELVTX_1D_in_MT_%i",0);
    hTitle= Form("hGEN_INELVTX_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hGEN_INELVTX_1D_in_MT[0]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    SetAxis(hGEN_INELVTX_1D_in_MT[0],"PT 1D");
    
    hName = Form("hREC_1D_in_2D_bin_in_MT_%i",0);
    hTitle= Form("hREC_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hREC_1D_in_2D_bin_in_MT[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hREC_1D_in_2D_bin_in_MT[0],"PT 1D");
        
    hName = Form("hGEN_1D_in_2D_bin_in_MT_%i",0);
    hTitle= Form("hGEN_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hGEN_1D_in_2D_bin_in_MT[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_1D_in_2D_bin_in_MT[0],"PT 1D");
    
    hName = Form("hTRU_1D_in_2D_bin_in_MT_%i",0);
    hTitle= Form("hTRU_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hTRU_1D_in_2D_bin_in_MT[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hTRU_1D_in_2D_bin_in_MT[0],"PT 1D");
    
    hName = Form("hEFF_1D_in_2D_bin_in_MT_%i",0);
    hTitle= Form("hEFF_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_1D_in_2D_bin_in_MT[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_1D_in_2D_bin_in_MT[0],"PT 1D");
    
    hName = Form("hEFF_SL_1D_in_2D_bin_in_MT_%i",0);
    hTitle= Form("hEFF_SL_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_SL_1D_in_2D_bin_in_MT[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_SL_1D_in_2D_bin_in_MT[0],"PT 1D");
    
    hName = Form("hGEN_INELVTX_1D_in_2D_bin_in_MT_%i",0);
    hTitle= Form("hGEN_INELVTX_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hGEN_INELVTX_1D_in_2D_bin_in_MT[0]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
    SetAxis(hGEN_INELVTX_1D_in_2D_bin_in_MT[0],"PT 1D");
    
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
        
        hName = Form("hEFF_SL_1D_in_MT_%i",iMult);
        hTitle= Form("hEFF_SL_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_SL_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hEFF_SL_1D_in_MT[iMult],"PT 1D");
        
        hName = Form("hGEN_INELVTX_1D_in_MT_%i",iMult);
        hTitle= Form("hGEN_INELVTX_1D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hGEN_INELVTX_1D_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
        SetAxis(hGEN_INELVTX_1D_in_MT[iMult],"PT 1D");
        
        hName = Form("hREC_1D_in_2D_bin_in_MT_%i",iMult);
        hTitle= Form("hREC_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hREC_1D_in_2D_bin_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hREC_1D_in_2D_bin_in_MT[iMult],"PT 1D");
            
        hName = Form("hGEN_1D_in_2D_bin_in_MT_%i",iMult);
        hTitle= Form("hGEN_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hGEN_1D_in_2D_bin_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_1D_in_2D_bin_in_MT[iMult],"PT 1D");
        
        hName = Form("hTRU_1D_in_2D_bin_in_MT_%i",iMult);
        hTitle= Form("hTRU_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hTRU_1D_in_2D_bin_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_1D_in_2D_bin_in_MT[iMult],"PT 1D");
        
        hName = Form("hEFF_1D_in_2D_bin_in_MT_%i",iMult);
        hTitle= Form("hEFF_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_1D_in_2D_bin_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_1D_in_2D_bin_in_MT[iMult],"PT 1D");
        
        hName = Form("hEFF_SL_1D_in_2D_bin_in_MT_%i",iMult);
        hTitle= Form("hEFF_SL_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_SL_1D_in_2D_bin_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_SL_1D_in_2D_bin_in_MT[iMult],"PT 1D");
        
        hName = Form("hGEN_INELVTX_1D_in_2D_bin_in_MT_%i",iMult);
        hTitle= Form("hGEN_INELVTX_1D_in_2D_bin_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult],"PT 1D");
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
    TH2F      **hEFF_SL_2D_in_MT            = new TH2F     *[nBinMult+1];
    TH2F      **hEFF_SL_2D_in_MT_fr_1D      = new TH2F     *[nBinMult+1];
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
    
    hName = Form("hEFF_SL_2D_in_MT_%i",0);
    hTitle= Form("hEFF_SL_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_SL_2D_in_MT[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_SL_2D_in_MT[0],"PT 2D");
    
    hName = Form("hEFF_SL_2D_in_MT_fr_1D_%i",0);
    hTitle= Form("hEFF_SL_2D_in_MT_fr_1D Multiplicity [%.2f#;%.2f]",fArrMult[0],fArrMult[nBinMult]);
    hEFF_SL_2D_in_MT_fr_1D[0]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    SetAxis(hEFF_SL_2D_in_MT_fr_1D[0],"PT 2D");
    
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
        
        hName = Form("hEFF_SL_2D_in_MT_%i",iMult);
        hTitle= Form("hEFF_SL_2D_in_MT Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_SL_2D_in_MT[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_SL_2D_in_MT[iMult],"PT 2D");
        
        hName = Form("hEFF_SL_2D_in_MT_fr_1D_%i",iMult);
        hTitle= Form("hEFF_SL_2D_in_MT_fr_1D Multiplicity [%.2f#;%.2f]",fArrMult[iMult-1],fArrMult[iMult]);
        hEFF_SL_2D_in_MT_fr_1D[iMult]   = new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_SL_2D_in_MT_fr_1D[iMult],"PT 2D");
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
    TH1F      **hREC_1D_in_RP_in_2D_bin      = new TH1F     *[nBinRap_];
    TH1F      **hGEN_1D_in_RP_in_2D_bin      = new TH1F     *[nBinRap_];
    TH1F      **hTRU_1D_in_RP_in_2D_bin      = new TH1F     *[nBinRap_];
    TH1F      **hEFF_1D_in_RP_in_2D_bin      = new TH1F     *[nBinRap_];
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
        
        hName = Form("hREC_1D_in_RP_in_2D_bin_%i",iRap);
        hTitle= Form("hREC_1D_in_RP_in_2D_bin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hREC_1D_in_RP_in_2D_bin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hREC_1D_in_RP_in_2D_bin[iRap],"PT 1D");
            
        hName = Form("hGEN_1D_in_RP_in_2D_bin_%i",iRap);
        hTitle= Form("hGEN_1D_in_RP_in_2D_bin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hGEN_1D_in_RP_in_2D_bin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hGEN_1D_in_RP_in_2D_bin[iRap],"PT 1D");
        
        hName = Form("hTRU_1D_in_RP_in_2D_bin_%i",iRap);
        hTitle= Form("hTRU_1D_in_RP_in_2D_bin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hTRU_1D_in_RP_in_2D_bin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hTRU_1D_in_RP_in_2D_bin[iRap],"PT 1D");
        
        hName = Form("hEFF_1D_in_RP_in_2D_bin_%i",iRap);
        hTitle= Form("hEFF_1D_in_RP_in_2D_bin Rapidity [%.2f#;%.2f]",fArrRap_[iRap],fArrRap_[iRap+1]);
        hEFF_1D_in_RP_in_2D_bin[iRap]   = new TH1F (hName,hTitle,nBinPT2D,fArrPT2D);
        SetAxis(hEFF_1D_in_RP_in_2D_bin[iRap],"PT 1D");
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
    
    fStartTimer("Efficiency Utility Histograms Production");
    
    // Evaluating entries
    Int_t nEvents = (!TPhiEfficiency) ? 0 : ( nEventsCut == -1.? TPhiEfficiency->GetEntries() : min( nEventsCut, (int)TPhiEfficiency->GetEntries() ) );
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent+= 12 )    {
        // Recovering events
        TPhiEfficiency->GetEntry(iEvent);
        
        fPrintLoopTimer("Efficiency Utility Histograms Production",iEvent,nEvents,kPrintIntervalPP);
        
        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        //
        auto    fTru = 0;
        auto    fGen = 0;
        auto    fRec = 0;
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPhi; iPhi++ ) {
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],kPhiMesonMass_);
            
            if ( fabs(LPhi_candidate1.Rapidity()) < 0.5 )   {
                fTru++;
                if ( evPhiEfficiency.Selection[iPhi] >= 1 ) fGen++;
                if ( evPhiEfficiency.Selection[iPhi] >= 2 ) fRec++;
            }
            
            if ( !fAcceptCandidate(kPhiMesonMass_,LPhi_candidate1.Pt()) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        //
        hProdDistrTRU->Fill(fTru);
        hProdDistrGEN->Fill(fGen);
        hProdDistrREC->Fill(fRec);
        //
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
            Bool_t  fIsMBevent          =   ( evPhiEfficiency.TrueEventMask == 0 );
            Bool_t  fIsCUT_VTX          =   ( evPhiEfficiency.TrueEventMask == 0 || fCheckMask( evPhiEfficiency.TrueEventMask, 5 ) );
            //
            // >>-->>-->> Multiplicity
            //
            Int_t   iMult               =   fGetBinMult(evPhiEfficiency.Multiplicity);
            Bool_t  fHasMultiplicity    =   (iMult != -1);
            //
            // >>-->>-->> True Phis
            //
            Int_t   iSelection          =   (int)evPhiEfficiency.Selection[U_AccCand[iPhi]];
            Bool_t  iIsGen              =   (iSelection >= 1);
            Bool_t  iIsRec              =   (iSelection >= 2);
            Float_t iTransMom           =   LPhi_candidate1.Pt();
            Float_t iTrueIMass          =   evPhiCandidate.TrueInvMass[U_AccCand[iPhi]];
            Float_t iInvarMass          =   evPhiCandidate.InvMass[U_AccCand[iPhi]];
            Float_t iRapidity           =   LPhi_candidate1.Rapidity();
            Int_t   iRap                =   fGetBinRap_(LPhi_candidate1.Rapidity());
            Bool_t  fHasRapidity        =   fabs(LPhi_candidate1.Rapidity()) <0.5;
            //
            if ( fHasRapidity )  {
                //
                //  Mid-Rapidity Analyses
                if ( iIsGen )hGEN_INELFLL_1D                                ->  Fill(iTransMom);
                if ( iIsGen )hGEN_INELFLL_1D_in_2D_bin                      ->  Fill(iTransMom);
                //
                if ( fIsMBevent )   {
                    //
                    //  Minimum-Bias Analyses
                    //
                    // >>-->> Yield
                    //
                    hTRU_1D                                     ->  Fill(iTransMom);
                    hTRU_1D_in_2D_bin                           ->  Fill(iTransMom);
                    if ( iIsGen )   {
                        hGEN_1D                                 ->  Fill(iTransMom);
                        hGEN_Rw_1D                              ->  Fill(iTransMom);
                        hGEN_1D_in_2D_bin                       ->  Fill(iTransMom);
                    } if ( iIsRec )  {
                        hREC_1D                                 ->  Fill(iTransMom);
                        hREC_Rw_1D                              ->  Fill(iTransMom);
                        hREC_1D_in_2D_bin                       ->  Fill(iTransMom);
                    }
                    if ( fHasMultiplicity ) {
                        //
                        //  Multiplicity differentiation
                        hTRU_1D_in_MT[iMult+1]                  ->  Fill(iTransMom);
                        hTRU_1D_in_MT[0]                        ->  Fill(iTransMom);
                        hTRU_1D_in_2D_bin_in_MT[iMult+1]        ->  Fill(iTransMom);
                        hTRU_1D_in_2D_bin_in_MT[0]              ->  Fill(iTransMom);
                        if ( iIsGen )   {
                            hGEN_1D_in_MT[iMult+1]              ->  Fill(iTransMom);
                            hGEN_1D_in_MT[0]                    ->  Fill(iTransMom);
                            hGEN_1D_in_2D_bin_in_MT[iMult+1]    ->  Fill(iTransMom);
                            hGEN_1D_in_2D_bin_in_MT[0]          ->  Fill(iTransMom);
                        } if ( iIsRec )  {
                            hREC_1D_in_MT[iMult+1]              ->  Fill(iTransMom);
                            hREC_1D_in_MT[0]                    ->  Fill(iTransMom);
                            hREC_1D_in_2D_bin_in_MT[iMult+1]    ->  Fill(iTransMom);
                            hREC_1D_in_2D_bin_in_MT[0]          ->  Fill(iTransMom);
                        }
                    }
                }
                if ( fIsCUT_VTX && iIsGen )  {
                    //
                    //  Signal Loss Reference
                    hGEN_INELVTX_1D                             ->  Fill(iTransMom);
                    hGEN_INELVTX_1D_in_2D_bin                   ->  Fill(iTransMom);
                    if ( fHasMultiplicity ) {
                        //
                        //  Multiplicity differentiation
                        hGEN_INELVTX_1D_in_MT[0]                    ->  Fill(iTransMom);
                        hGEN_INELVTX_1D_in_MT[iMult+1]              ->  Fill(iTransMom);
                        hGEN_INELVTX_1D_in_2D_bin_in_MT[0]          ->  Fill(iTransMom);
                        hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult+1]    ->  Fill(iTransMom);
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
                Float_t jInvarMass          =   evPhiCandidate.InvMass[U_AccCand[jPhi]];
                Float_t jTrueIMass          =   evPhiCandidate.TrueInvMass[U_AccCand[jPhi]];
                Float_t jRapidity           =   LPhi_candidate2.Rapidity();
                Int_t   jRap                =   fGetBinRap_(LPhi_candidate2.Rapidity());
                        fHasRapidity        =   fabs(LPhi_candidate1.Rapidity())<0.5 && fabs(LPhi_candidate2.Rapidity())<0.5;
                Int_t   ijRap               =   fGetBinRap_(fabs(LPhi_candidate2.Rapidity()-LPhi_candidate1.Rapidity()));
                //
                if ( fHasRapidity )  {
                    if ( iIsGen && jIsGen ) hGEN_INELFLL_2D                                     ->  Fill(iTransMom,jTransMom,0.5);
                }
                if ( fHasRapidity && ( fIsCUT_VTX || fIsMBevent ) )  {
                    if ( iIsGen && jIsGen ) hGEN_INELVTX_2D                                     ->  Fill(iTransMom,jTransMom,0.5);
                }
                if ( fHasRapidity && fIsMBevent )   {
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
                        hTRU_1D_in_RP[ijRap]                             ->  Fill(iTransMom,0.5);
                        hTRU_1D_in_RP_in_2D_bin[ijRap]                    ->  Fill(iTransMom,0.5);
                        hTRU_1D_in_RP[ijRap]                             ->  Fill(jTransMom,0.5);
                        hTRU_1D_in_RP_in_2D_bin[ijRap]                    ->  Fill(jTransMom,0.5);
                        if ( iIsGen )   {
                            hGEN_1D_in_RP[ijRap]                         ->  Fill(iTransMom,0.5);
                            hGEN_1D_in_RP_in_2D_bin[ijRap]                ->  Fill(iTransMom,0.5);
                            hGEN_1D_in_RP[ijRap]                         ->  Fill(jTransMom,0.5);
                            hGEN_1D_in_RP_in_2D_bin[ijRap]                ->  Fill(jTransMom,0.5);
                        } if ( iIsRec )  {
                            hREC_1D_in_RP[ijRap]                         ->  Fill(iTransMom,0.5);
                            hREC_1D_in_RP_in_2D_bin[ijRap]                ->  Fill(iTransMom,0.5);
                            hREC_1D_in_RP[ijRap]                         ->  Fill(jTransMom,0.5);
                            hREC_1D_in_RP_in_2D_bin[ijRap]                ->  Fill(jTransMom,0.5);
                        }
                    }
                //
                // >>-->> Yield
                //
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
    fStopTimer("Efficiency Utility Histograms Production");
    //
    //>>    Evaluating entries and saving them for later
    nEvents = 0;//(!TPhiCandidate) ? 0 : ( nEventsCut == -1.? TPhiCandidate->GetEntries() : min( nEventsCut, TPhiCandidate->GetEntries() ) );
    //
    //>>    If there are actually entries start timer
    if ( nEvents > 0 )  fStartTimer("Resolution Utility Histograms Production");
    //
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        
        fPrintLoopTimer("Resolution Utility Histograms Production",iEvent,nEvents,kPrintIntervalPP);
        
        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        // Discarding Pile-up events
        //if ( kDoYield           &&  fCheckMask(evPhiCandidate.EventMask,1) ) continue;
        //if ( kDoMultiplicity    &&  fCheckMask(evPhiCandidate.TrueEventMask,1) ) continue;
        
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; iPhi++ )  {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fAcceptCandidate(evPhiCandidate.InvMass[iPhi],LPhi_candidate1.Pt()) ) continue;
            if ( !kDoRapidity && !fCutRapidity(LPhi_candidate1.Rapidity()) ) continue;
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
            Bool_t  fHasRapidity        =   fabs(evPhiCandidate.Rap[U_AccCand[iPhi]]) < 0.5;
            //
            // >->-->-> Yie
            //
            if ( fHasRapidity && iTrueIMass !=0 ) {
                hMRS_1D[iPT1D]           -> Fill(-iTrueIMass+iInvarMass);
                hMRS_1D_in_2D_bin[iPT2D] -> Fill(-iTrueIMass+iInvarMass);
                hMDS_1D[iPT1D]           -> Fill(iInvarMass);
                hMDS_1D_in_2D_bin[iPT2D] -> Fill(iInvarMass);
                hTMD_1D[iPT1D]           -> Fill(iTrueIMass);
                hTMD_1D_in_2D_bin[iPT2D] -> Fill(iTrueIMass);
                hGEN_IM_1D               -> Fill(iInvarMass);
                hREC_IM_1D               -> Fill(iInvarMass);
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
                if  ( fHasRapidity && iTrueIMass !=0  && jTrueIMass !=0 )    {
                //
                // >->-->-> Yield
                //
                    hMRS_2D[iPT2D][jPT2D]   ->Fill(-iTrueIMass+iInvarMass,-jTrueIMass+jInvarMass);
                    hMDS_2D[iPT2D][jPT2D]   ->Fill(iInvarMass,jInvarMass);
                    hTMD_2D[iPT2D][jPT2D]   ->Fill(iTrueIMass,jTrueIMass);
                //
                }
                //
            }
        }
    }
    //
    fStopTimer("Resolution Utility Histograms Production");
    //
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    //
    // >> YIELD ANALYSIS //
    //
    auto fNormEvent = fHEventCount->GetBinContent(kTrigger);
    auto fTotlEvent = hProdDistrTRU->GetEntries();
    /*
    for ( Int_t iCnt = 0; iCnt < fNormEvent - fTotlEvent; iCnt++ ) {
        hProdDistrTRU               ->Fill(0);
        hProdDistrGEN               ->Fill(0);
        hProdDistrREC               ->Fill(0);
    }
    */
    hEFF_1D                             ->Divide(hREC_1D,                   hGEN_1D,                1.,1.,"b");
    hEFF_IM_1D                          ->Divide(hREC_IM_1D,                hGEN_IM_1D,             1.,1.,"b");
    hEFF_1D_in_2D_bin                   ->Divide(hREC_1D_in_2D_bin,         hGEN_1D_in_2D_bin,      1.,1.,"b");
    hEFF_2D                             ->Divide(hREC_2D,                   hGEN_2D,                1.,1.,"b");
    hEFF_SL_1D                          ->Divide(hGEN_INELVTX_1D,           hGEN_1D,                1.,1.,"b");
    hEFF_SL_1D_in_2D_bin                ->Divide(hGEN_INELVTX_1D_in_2D_bin, hGEN_1D_in_2D_bin,      1.,1.,"b");
    hEFF_SL_2D                          ->Divide(hGEN_INELVTX_2D,           hGEN_2D,                1.,1.,"b");
    for ( Int_t iPT2D = 1; iPT2D <= nBinPT2D; iPT2D++ )  {
        for ( Int_t jPT2D = 1; jPT2D <= nBinPT2D; jPT2D++ )  {
            auto iEffVal    =   hEFF_1D_in_2D_bin->GetBinContent (iPT2D);
            auto jEffVal    =   hEFF_1D_in_2D_bin->GetBinContent (jPT2D);
            auto iEffErr    =   hEFF_1D_in_2D_bin->GetBinError   (iPT2D);
            auto jEffErr    =   hEFF_1D_in_2D_bin->GetBinError   (jPT2D);
            hEFF_2D_fr_1D->SetBinContent(iPT2D,jPT2D,iEffVal*jEffVal);
            hEFF_2D_fr_1D->SetBinError  (iPT2D,jPT2D,iEffVal*jEffVal*(iEffErr/iEffVal + jEffErr/jEffVal));
            //
            iEffVal    =   1./hEFF_SL_1D_in_2D_bin->GetBinContent (iPT2D);
            jEffVal    =   1./hEFF_SL_1D_in_2D_bin->GetBinContent (jPT2D);
            iEffErr    =   hEFF_SL_1D_in_2D_bin->GetBinError   (iPT2D) * iEffVal;
            jEffErr    =   hEFF_SL_1D_in_2D_bin->GetBinError   (jPT2D) * jEffVal;
            hEFF_SL_2D_fr_1D->SetBinContent(iPT2D,jPT2D,1./(iEffVal*jEffVal));
            hEFF_SL_2D_fr_1D->SetBinError  (iPT2D,jPT2D,(1./(iEffVal*jEffVal))*(iEffErr+jEffErr));
        }
    }
    //  Renormalisation
    //
    hREC_1D                     ->Scale(1.,"width");
    hREC_Rw_1D                  ->Scale(1.,"width");
    hREC_IM_1D                  ->Scale(1.,"width");
    hREC_1D_in_2D_bin           ->Scale(1.,"width");
    hREC_2D                     ->Scale(1.,"width");
    hREC_1D                     ->Scale(1./fNormEvent);
    hREC_Rw_1D                  ->Scale(1./fNormEvent);
    hREC_IM_1D                  ->Scale(1./fNormEvent);
    hREC_1D_in_2D_bin           ->Scale(1./fNormEvent);
    hREC_2D                     ->Scale(1./fNormEvent);
    hProdDistrTRU               ->Scale(1./fNormEvent);
    hProdDistrGEN               ->Scale(1./fNormEvent);
    hProdDistrREC               ->Scale(1./fNormEvent);
    //
    hGEN_1D                     ->Scale(1.,"width");
    hGEN_1D_in_2D_bin           ->Scale(1.,"width");
    hGEN_2D                     ->Scale(1.,"width");
    hGEN_Rw_1D                  ->Scale(1.,"width");
    hGEN_INELVTX_1D             ->Scale(1.,"width");
    hGEN_INELFLL_1D             ->Scale(1.,"width");
    hGEN_INELVTX_1D_in_2D_bin   ->Scale(1.,"width");
    hGEN_INELFLL_1D_in_2D_bin   ->Scale(1.,"width");
    hGEN_INELVTX_2D             ->Scale(1.,"width");
    hGEN_INELFLL_2D             ->Scale(1.,"width");
    hGEN_1D                     ->Scale(1./fNormEvent);
    hGEN_1D_in_2D_bin           ->Scale(1./fNormEvent);
    hGEN_2D                     ->Scale(1./fNormEvent);
    hGEN_Rw_1D                  ->Scale(1./fNormEvent);
    hGEN_INELVTX_1D             ->Scale(1./fNormEvent);
    hGEN_INELFLL_1D             ->Scale(1./fNormEvent);
    hGEN_INELVTX_1D_in_2D_bin   ->Scale(1./fNormEvent);
    hGEN_INELFLL_1D_in_2D_bin   ->Scale(1./fNormEvent);
    hGEN_INELVTX_2D             ->Scale(1./fNormEvent);
    hGEN_INELFLL_2D             ->Scale(1./fNormEvent);
    //
    hTRU_1D                     ->Scale(1.,"width");
    hTRU_1D_in_2D_bin           ->Scale(1.,"width");
    hTRU_2D                     ->Scale(1.,"width");
    hTRU_1D                     ->Scale(1./fNormEvent);
    hTRU_1D_in_2D_bin           ->Scale(1./fNormEvent);
    hTRU_2D                     ->Scale(1./fNormEvent);
    //
    
    // >> MULTIPLICITY ANALYSIS //
    //
    for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
    {
        hEFF_1D_in_MT[iMult]                ->Divide(hREC_1D_in_MT[iMult],          hGEN_1D_in_MT[iMult],           1.,1.,"b");
        hEFF_1D_in_2D_bin_in_MT[iMult]      ->Divide(hREC_1D_in_2D_bin_in_MT[iMult],hGEN_1D_in_2D_bin_in_MT[iMult], 1.,1.,"b");
        hEFF_2D_in_MT[iMult]                ->Divide(hREC_2D_in_MT[iMult],          hGEN_2D_in_MT[iMult],           1.,1.,"b");
        hEFF_SL_1D_in_MT[iMult]             ->Divide(hGEN_INELVTX_1D_in_MT[iMult],          hGEN_1D_in_MT[iMult],           1.,1.,"b");
        hEFF_SL_1D_in_2D_bin_in_MT[iMult]   ->Divide(hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult],hGEN_1D_in_2D_bin_in_MT[iMult], 1.,1.,"b");
        hREC_1D_in_MT[iMult]                ->Scale(1.,"width");
        hGEN_1D_in_MT[iMult]                ->Scale(1.,"width");
        hTRU_1D_in_MT[iMult]                ->Scale(1.,"width");
        hGEN_INELVTX_1D_in_MT[iMult]        ->Scale(1.,"width");
        hREC_1D_in_2D_bin_in_MT[iMult]      ->Scale(1.,"width");
        hGEN_1D_in_2D_bin_in_MT[iMult]      ->Scale(1.,"width");
        hTRU_1D_in_2D_bin_in_MT[iMult]      ->Scale(1.,"width");
        hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult]->Scale(1.,"width");
        hREC_2D_in_MT[iMult]                ->Scale(1.,"width");
        hGEN_2D_in_MT[iMult]                ->Scale(1.,"width");
        hTRU_2D_in_MT[iMult]                ->Scale(1.,"width");
        hREC_1D_in_MT[iMult]                ->Scale(1./fNormEvent);
        hGEN_1D_in_MT[iMult]                ->Scale(1./fNormEvent);
        hTRU_1D_in_MT[iMult]                ->Scale(1./fNormEvent);
        hGEN_INELVTX_1D_in_MT[iMult]        ->Scale(1./fNormEvent);
        hREC_1D_in_2D_bin_in_MT[iMult]      ->Scale(1./fNormEvent);
        hGEN_1D_in_2D_bin_in_MT[iMult]      ->Scale(1./fNormEvent);
        hTRU_1D_in_2D_bin_in_MT[iMult]      ->Scale(1./fNormEvent);
        hREC_2D_in_MT[iMult]                ->Scale(1./fNormEvent);
        hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult]->Scale(1./fNormEvent);
        hGEN_2D_in_MT[iMult]                ->Scale(1./fNormEvent);
        hTRU_2D_in_MT[iMult]                ->Scale(1./fNormEvent);
    }
    //
    
    // >> RAPIDITY ANALYSIS //
    //
    for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
    {
        hEFF_1D_in_RP[iRap]               ->Divide(hREC_1D_in_RP[iRap],         hGEN_1D_in_RP[iRap],          1.,1.,"b");
        hEFF_1D_in_RP_in_2D_bin[iRap]      ->Divide(hREC_1D_in_RP_in_2D_bin[iRap],hGEN_1D_in_RP_in_2D_bin[iRap], 1.,1.,"b");
        hEFF_2D_in_RP[iRap]               ->Divide(hREC_2D_in_RP[iRap],         hGEN_2D_in_RP[iRap],          1.,1.,"b");
        hREC_1D_in_RP[iRap]->Scale(1.,"width");
        hGEN_1D_in_RP[iRap]->Scale(1.,"width");
        hTRU_1D_in_RP[iRap]->Scale(1.,"width");
        hREC_1D_in_RP_in_2D_bin[iRap]->Scale(1.,"width");
        hGEN_1D_in_RP_in_2D_bin[iRap]->Scale(1.,"width");
        hTRU_1D_in_RP_in_2D_bin[iRap]->Scale(1.,"width");
        hREC_2D_in_RP[iRap]->Scale(1.,"width");
        hGEN_2D_in_RP[iRap]->Scale(1.,"width");
        hTRU_2D_in_RP[iRap]->Scale(1.,"width");
        hREC_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hGEN_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hTRU_1D_in_RP[iRap]->Scale(1./fNormEvent);
        hREC_1D_in_RP_in_2D_bin[iRap]->Scale(1./fNormEvent);
        hGEN_1D_in_RP_in_2D_bin[iRap]->Scale(1./fNormEvent);
        hTRU_1D_in_RP_in_2D_bin[iRap]->Scale(1./fNormEvent);
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
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Trigger")+kFolder).Data())));
        TFile *outFil1  =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Trigger")+kFolder).Data()),"recreate");
        //
        outFil1->Close();
    }
    //
    // >> Yield Analysis
    //
    if ( kDoYield ) {
        // >> All Analysis Utility
        //
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Dir_,(TString("Yield")+kFolder).Data())));
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s/Plots",Form(kMassResolution_Dir_,(TString("Yield")+kFolder).Data())));
        TFile *outFil0  =   new TFile   (Form(kMassResolution_Prod,(TString("Yield")+kFolder).Data()),"recreate");
        //
        for ( Int_t iTer = 0; iTer < nBinPT1D; iTer++ )  {
            hMRS_1D[iTer]->Write();
            hMDS_1D[iTer]->Write();
            hTMD_1D[iTer]->Write();
        }
        for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )  {
            hMRS_1D_in_2D_bin[iTer]->Write();
            hMDS_1D_in_2D_bin[iTer]->Write();
            hTMD_1D_in_2D_bin[iTer]->Write();
            for ( Int_t jTer = 0; jTer < nBinPT2D; jTer++ )  {
                hMRS_2D[iTer][jTer]->Write();
                hMDS_2D[iTer][jTer]->Write();
                hTMD_2D[iTer][jTer]->Write();
            }
        }
        //
        outFil0->Close();
        //
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data())));
        TFile *outFil2  =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()),"recreate");
        //
        hProdDistrTRU->Write();
        hProdDistrGEN->Write();
        hProdDistrREC->Write();
        hREC_1D->Write();
        hGEN_1D->Write();
        hREC_Rw_1D->Write();
        hGEN_Rw_1D->Write();
        hREC_IM_1D->Write();
        hGEN_IM_1D->Write();
        hTRU_1D->Write();
        hGEN_INELVTX_1D->Write();
        hGEN_INELFLL_1D->Write();
        hGEN_INELVTX_1D_in_2D_bin->Write();
        hGEN_INELFLL_1D_in_2D_bin->Write();
        hEFF_1D->Write();
        hEFF_IM_1D->Write();
        hEFF_SL_1D->Write();
        hREC_1D_in_2D_bin->Write();
        hGEN_1D_in_2D_bin->Write();
        hTRU_1D_in_2D_bin->Write();
        hEFF_1D_in_2D_bin->Write();
        hEFF_SL_1D_in_2D_bin->Write();
        hREC_2D->Write();
        hGEN_2D->Write();
        hTRU_2D->Write();
        hGEN_INELVTX_2D->Write();
        hGEN_INELFLL_2D->Write();
        hEFF_2D->Write();
        hEFF_SL_2D->Write();
        hEFF_2D_fr_1D->Write();
        hEFF_SL_2D_fr_1D->Write();
        //
        outFil2->Close();
        //
        //  -   //  Quality Control Plots
        //************************************************************
        SetStyle();
        //
        gROOT           ->  SetBatch();
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/")).Data()));
        //
        TCanvas        *cDrawEff    =   new TCanvas("","",1200,1200);
        gStyle          ->  SetOptStat(0);
        gPad            ->  SetLogx();
        gPad            ->  SetGridy();
        uSetHisto(hEFF_1D,"EFF 1D");
        hEFF_1D         ->  Draw("");
        
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.19, 0.77,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        
        cDrawEff        ->  SaveAs((TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/")+TString("hEFF_1D.pdf")).Data());
        delete cDrawEff;
        //
                        cDrawEff    =   new TCanvas("","",900,1200);
        gStyle          ->  SetOptStat(0);
        cDrawEff->Divide(3,4);
        for ( Int_t iPT2D = 1; iPT2D <= nBinPT2D; iPT2D++ )  {
            cDrawEff->cd(iPT2D);
            gPad            ->  SetLogx();
            gPad            ->  SetGridy();
            auto h12D = hEFF_2D_fr_1D   ->ProjectionX(Form("12D_%i",iPT2D),iPT2D,iPT2D);
            auto h_2D = hEFF_2D         ->ProjectionX(Form("_2D_%i",iPT2D),iPT2D,iPT2D);
            uSetHisto(h12D,"EFF 1D");
            uSetHisto(h_2D,"EFF2 1D");
            h12D    ->  Draw();
            h_2D    ->  Draw("SAME");
            
            if ( iPT2D == 1 )   {
                uLatex->SetTextFont(60);
                uLatex->SetTextSize(0.05);
                uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
                uLatex->SetTextFont(42);
                uLatex->SetTextSize(0.04);
                uLatex->DrawLatexNDC(0.19, 0.77,"pp #sqrt{#it{s}}= 7 TeV");
                uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
                uLatex->DrawLatexNDC(0.19, 0.65,Form("%.2f < #it{p}_{T}^{#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[iPT2D-1],fArrPT2D[iPT2D]));
            } else {
                uLatex->DrawLatexNDC(0.19, 0.83,Form("%.2f < #it{p}_{T}^{#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[iPT2D-1],fArrPT2D[iPT2D]));
            }
        }
        cDrawEff        ->  SaveAs((TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/")+TString("hEFF_12D.pdf")).Data());
        delete cDrawEff;
        //
                    cDrawEff    =   new TCanvas("","",1200,1200);
        gStyle          ->  SetOptStat(0);
        gPad            ->  SetLogx();
        gPad            ->  SetGridy();
        uSetHisto(hEFF_SL_1D,"EFF_SL 1D");
        hEFF_SL_1D      ->  SetMaximum(+4.);
        hEFF_SL_1D      ->  SetMinimum(-1.);
        hEFF_SL_1D      ->  Draw("");
        
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.19, 0.77,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        
        cDrawEff        ->  SaveAs((TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/")+TString("hEFF_SL_1D.pdf")).Data());
        delete cDrawEff;
        //
                        cDrawEff    =   new TCanvas("","",900,1200);
        gStyle          ->  SetOptStat(0);
        cDrawEff->Divide(3,4);
        for ( Int_t iPT2D = 1; iPT2D <= nBinPT2D; iPT2D++ )  {
            cDrawEff->cd(iPT2D);
            gPad            ->  SetLogx();
            gPad            ->  SetGridy();
            auto h12D_SL = hEFF_SL_2D_fr_1D   ->ProjectionX(Form("12D_SL_%i",iPT2D),iPT2D,iPT2D);
            auto h_2D_SL = hEFF_SL_2D         ->ProjectionX(Form("_2D_SL_%i",iPT2D),iPT2D,iPT2D);
            uSetHisto(h12D_SL,"EFF SL 1D");
            h12D_SL      ->  SetMaximum(+4.);
            h12D_SL      ->  SetMinimum(-1.);
            uSetHisto(h_2D_SL,"EFF2 SL 1D");
            h12D_SL    ->  Draw();
            h_2D_SL    ->  Draw("SAME");
            
            if ( iPT2D == 1 )   {
                uLatex->SetTextFont(60);
                uLatex->SetTextSize(0.05);
                uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
                uLatex->SetTextFont(42);
                uLatex->SetTextSize(0.04);
                uLatex->DrawLatexNDC(0.19, 0.77,"pp #sqrt{#it{s}}= 7 TeV");
                uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
                uLatex->DrawLatexNDC(0.19, 0.65,Form("%.2f < #it{p}_{T}^{#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[iPT2D-1],fArrPT2D[iPT2D]));
            } else {
                uLatex->DrawLatexNDC(0.19, 0.83,Form("%.2f < #it{p}_{T}^{#phi_{2}} < %.2f GeV/#it{c}",fArrPT2D[iPT2D-1],fArrPT2D[iPT2D]));
            }
        }
        cDrawEff        ->  SaveAs((TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/")+TString("hEFF_SL_12D.pdf")).Data());
        delete cDrawEff;
        //
        gROOT           ->  SetBatch(kFALSE);
        //************************************************************
    }
    //
    // >> Multiplicity Analysis
    //
    if ( kDoMultiplicity )  {
        // >> All Analysis Utility
        //
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Dir_,(TString("Multiplicity")+kFolder).Data())));
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s/Plots",Form(kMassResolution_Dir_,(TString("Multiplicity")+kFolder).Data())));
        TFile *outFil0  =   new TFile   (Form(kMassResolution_Prod,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        for ( Int_t iTer = 0; iTer < nBinPT1D; iTer++ )  {
            hMRS_1D[iTer]->Write();
            hMDS_1D[iTer]->Write();
            hTMD_1D[iTer]->Write();
        }
        for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )  {
            hMRS_1D_in_2D_bin[iTer]->Write();
            hMDS_1D_in_2D_bin[iTer]->Write();
            hTMD_1D_in_2D_bin[iTer]->Write();
            for ( Int_t jTer = 0; jTer < nBinPT2D; jTer++ )  {
                hMRS_2D[iTer][jTer]->Write();
                hMDS_2D[iTer][jTer]->Write();
                hTMD_2D[iTer][jTer]->Write();
            }
        }
        //
        outFil0->Close();
        //
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data())));
        TFile *outFil3  =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity")+kFolder).Data()),"recreate");
        //
        hProdDistrTRU->Write();
        hProdDistrGEN->Write();
        hProdDistrREC->Write();
        hREC_1D->Write();
        hGEN_1D->Write();
        hREC_Rw_1D->Write();
        hGEN_Rw_1D->Write();
        hREC_IM_1D->Write();
        hGEN_IM_1D->Write();
        hTRU_1D->Write();
        hGEN_INELVTX_1D->Write();
        hGEN_INELFLL_1D->Write();
        hGEN_INELVTX_1D_in_2D_bin->Write();
        hGEN_INELFLL_1D_in_2D_bin->Write();
        hEFF_1D->Write();
        hEFF_IM_1D->Write();
        hEFF_SL_1D->Write();
        hREC_1D_in_2D_bin->Write();
        hGEN_1D_in_2D_bin->Write();
        hTRU_1D_in_2D_bin->Write();
        hEFF_1D_in_2D_bin->Write();
        hEFF_SL_1D_in_2D_bin->Write();
        hREC_2D->Write();
        hGEN_2D->Write();
        hTRU_2D->Write();
        hGEN_INELVTX_2D->Write();
        hGEN_INELFLL_2D->Write();
        hEFF_2D->Write();
        hEFF_SL_2D->Write();
        hEFF_2D_fr_1D->Write();
        hEFF_SL_2D_fr_1D->Write();
        //
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )
        {
            hREC_1D_in_MT[iMult]->Write();
            hGEN_1D_in_MT[iMult]->Write();
            hGEN_INELVTX_1D_in_MT[iMult]->Write();
            hTRU_1D_in_MT[iMult]->Write();
            hEFF_1D_in_MT[iMult]->Write();
            hEFF_SL_1D_in_MT[iMult]->Write();
            hREC_1D_in_2D_bin_in_MT[iMult]->Write();
            hGEN_1D_in_2D_bin_in_MT[iMult]->Write();
            hGEN_INELVTX_1D_in_2D_bin_in_MT[iMult]->Write();
            hTRU_1D_in_2D_bin_in_MT[iMult]->Write();
            hEFF_1D_in_2D_bin_in_MT[iMult]->Write();
            hEFF_SL_1D_in_2D_bin_in_MT[iMult]->Write();
            hREC_2D_in_MT[iMult]->Write();
            hGEN_2D_in_MT[iMult]->Write();
            hTRU_2D_in_MT[iMult]->Write();
            hEFF_2D_in_MT[iMult]->Write();
            hEFF_SL_2D_in_MT[iMult]->Write();
        }
        //
        //  -   //  Quality Control Plots
        //************************************************************
        //
        //  Initia settings
        SetStyle();
        gROOT           ->  SetBatch();
        //
        //  Make output directory
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/")).Data()));
        //
        //  Efficiency comparison
        TCanvas*        cDrawEfficienciesComparison =   new TCanvas("","",1200,1200);
        gStyle          ->  SetOptStat(0);
        gPad            ->  SetLogx();
        gPad            ->  SetGridy();
        if ( !kDoYield ) uSetHisto       ( hEFF_1D, "EFF 1D" );
        //
        TLegend*        lEfficiencies   =   new TLegend(0.625,0.88,0.88,0.7);
        lEfficiencies   ->  SetNColumns(2);
        lEfficiencies   ->  SetFillColorAlpha(0.,0.);
        lEfficiencies   ->  SetLineColorAlpha(0.,0.);
        //
        lEfficiencies   ->  AddEntry(hEFF_1D,"INEL","EP");
        for ( Int_t iTer = 0; iTer <= nBinMult; iTer++ ) {
            uSetHisto   ( hEFF_1D_in_MT[iTer], "EFF 1D" );
            //
            if ( iTer == 0 )    hEFF_1D_in_MT[iTer]->SetMarkerStyle(kMarkers[5]);
            else                hEFF_1D_in_MT[iTer]->SetMarkerStyle(kMarkers[4]);
            if ( iTer >= 2 )    hEFF_1D_in_MT[iTer]->SetMarkerColor(kColors[iTer+1]);
            else                hEFF_1D_in_MT[iTer]->SetMarkerColor(kColors[iTer]);
            if ( iTer >= 2 )    hEFF_1D_in_MT[iTer]->SetLineColor(kColors[iTer+1]);
            else                hEFF_1D_in_MT[iTer]->SetLineColor(kColors[iTer]);
            //
            hEFF_1D_in_MT[iTer]    ->  Draw("SAME");
            if ( iTer == 0 )    lEfficiencies ->AddEntry(hEFF_1D_in_MT[iTer],Form("V0M [%3.0f;%3.0f]",fArrMult[0],fArrMult[nBinMult]),"EP");
            else                lEfficiencies ->AddEntry(hEFF_1D_in_MT[iTer],Form("V0M [%3.0f;%3.0f]",fArrMult[iTer-1],fArrMult[iTer]),"EP");
        }
        hEFF_1D         ->  Draw("SAME");
        lEfficiencies   ->  Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.19, 0.77,"pp #sqrt{#it{s}}= 5 TeV");
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        cDrawEfficienciesComparison ->  SaveAs((TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/")+TString("hEFF_in_Mult.pdf")).Data());
        delete cDrawEfficienciesComparison;
        //
        cDrawEfficienciesComparison =   new TCanvas("","",1200,1200);
        gStyle          ->  SetOptStat(0);
        gPad            ->  SetLogx();
        gPad            ->  SetGridy();
        if ( !kDoYield ) uSetHisto       ( hEFF_SL_1D, "EFF SL 1D" );
        //
        for ( Int_t iTer = 0; iTer <= nBinMult; iTer++ ) {
            uSetHisto   ( hEFF_SL_1D_in_MT[iTer], "EFF SL 1D" );
            //
            if ( iTer == 0 )    hEFF_SL_1D_in_MT[iTer]->SetMarkerStyle(kMarkers[5]);
            else                hEFF_SL_1D_in_MT[iTer]->SetMarkerStyle(kMarkers[4]);
            if ( iTer >= 2 )    hEFF_SL_1D_in_MT[iTer]->SetMarkerColor(kColors[iTer+1]);
            else                hEFF_SL_1D_in_MT[iTer]->SetMarkerColor(kColors[iTer]);
            if ( iTer >= 2 )    hEFF_SL_1D_in_MT[iTer]->SetLineColor(kColors[iTer+1]);
            else                hEFF_SL_1D_in_MT[iTer]->SetLineColor(kColors[iTer]);
            //
            hEFF_SL_1D_in_MT[iTer]    ->  Draw("SAME");
        }
        hEFF_SL_1D      ->  Draw("SAME");
        lEfficiencies   ->  Draw("SAME");
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.19, 0.77,"pp #sqrt{#it{s}}= 5 TeV");
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        cDrawEfficienciesComparison ->  SaveAs((TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/")+TString("hEFF_SL_in_Mult.pdf")).Data());
        delete cDrawEfficienciesComparison;
        gROOT           ->  SetBatch(kFALSE);
        //************************************************************
        //
        outFil3->Close();
    }
    //
    if ( kDoRapidity )  {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Rapidity")+kFolder).Data())));
        TFile *outFil4  =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Rapidity")+kFolder).Data()),"recreate");
        //
        for ( Int_t iRap = 0; iRap < nBinRap_; iRap++ )
        {
            hREC_1D_in_RP[iRap]->Write();
            hGEN_1D_in_RP[iRap]->Write();
            hTRU_1D_in_RP[iRap]->Write();
            hEFF_1D_in_RP[iRap]->Write();
            hREC_1D_in_RP_in_2D_bin[iRap]->Write();
            hGEN_1D_in_RP_in_2D_bin[iRap]->Write();
            hTRU_1D_in_RP_in_2D_bin[iRap]->Write();
            hEFF_1D_in_RP_in_2D_bin[iRap]->Write();
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
