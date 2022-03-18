#ifndef ALIANALYSISPHIPAIR_CONSTANTS_H
#define ALIANALYSISPHIPAIR_CONSTANTS_H

#include "./AliAnalysisPhiPair.h"

//------------------------------//
//      !TODO: CLEAN - UP       //
//------------------------------//
//
auto kCheckPileUp   = true;
auto is_pp_anl      = true;
auto is_pb_anl      = false;
auto kEnergy        = 7.00;
//
//>->->->->->   Branching Ratio
auto const  kBR                     =   0.492;
auto const  kBRerr                  =   0.005;
auto const  kSysHig_BR              =   kBRerr/kBR;
auto const  kSysLow_BR              =   kBRerr/kBR;
//>->->->->->   Branching Ratio
auto const  kTR                     =   1.;
auto const  kTRerr                  =   0.08;
auto const  kSysHig_TR              =   kTRerr/kTR;
auto const  kSysLow_TR              =   kTRerr/kTR;
//
//-//       TRIGGER EFFICIENCY
//
//  -   13 TeV
auto const  kTriggerEff13TeV        =   0.7448;
auto const  kTriggerEff13TeV_Hig    =   0.0190/0.7448;
auto const  kTriggerEff13TeV_Low    =   0.0190/0.7448;
//
//  -   7 TeV
auto const  kTriggerEff10bcdef      =   0.852;
auto const  kSysHig_TE              =   0.062/0.852;
auto const  kSysLow_TE              =   0.030/0.852;
//
//  -   5 TeV
auto const  kTriggerEff15n17pq      =   .7574;
auto const  kLHC15n_nEV             =   1.013736e8;
auto const  kLHC17pq_nEV            =   8.692481e8;
auto const  kLHC15n_fEV             =   (kLHC15n_nEV)/(kLHC15n_nEV+kLHC17pq_nEV);
auto const  kLHC17pq_fEV            =   (kLHC17pq_nEV)/(kLHC15n_nEV+kLHC17pq_nEV);
//                                       0-1    1-5     5-10    10-15   15-20   20-30   30-40   40-50   50-70   70-100
auto kTriggerEffMNBins = 10;
std::vector<Float_t> const   kTriggerEffMBins   =   {0.,    1.,     5.,     10.,    15.,    20.,    30.,    40.,    50.,    70.,     100.};
std::vector<Float_t> const   kTriggerEffM15n    =   {0.999, 0.998,  0.994,  0.990,  0.984,  0.972,  0.952,  0.930,  0.890,  0.706};
std::vector<Float_t> const   kTriggerEffM17pq   =   {0.999, 0.999,  0.999,  0.998,  0.998,  0.995,  0.989,  0.978,  0.947,  0.795};
std::vector<Float_t> const   kTriggerEffMult_p_p__05TeV =   {1.000, 1.000,  0.999,  0.998,  0.997,  0.995,  0.988,  0.975,  0.942,  0.795};
std::vector<Float_t> const   kTriggerEffMult_p_Pb_05TeV =   {0.999, 0.999,  0.999,  0.998,  0.999,  0.999,  0.999,  0.999,  0.999,  0.999};
std::vector<Float_t> const   kTriggerEffMult_p_p__07TeV =   {1.000, 1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000}; // Not used
std::vector<Float_t> const   kTriggerEffMult_p_p__13TeV =   {1.000, 1.000,  0.999,  0.998,  0.997,  0.995,  0.988,  0.975,  0.942,  0.795};

//-//   Analysis settings
auto        kDoMultiplicity         =   false;
auto        kDoYield                =   false;
auto        kDoCorrelation          =   false;
auto        kDoTrigger              =   false;
auto        kDoRapidity             =   false;

//-// InvMass range Pythia MC
auto const  fMinIMMC                =   0.75;
auto const  fMaxIMMC                =   1.25;



//  Systematics Options
//  >   Signal Extraction
std::vector<TString>    kSyst_SEX_1D_Options    =   {"RA","RB"};//,"RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4"  };
std::vector<TString>    kSyst_SEX_1D_Names      =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4"  };
std::vector<TString>    kSyst_SEX_1D_Legend     =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4"  };
//
std::vector<TString>    kSyst_SEX_2D_Options    =   {"RA","RB"};//,"RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4",  "BKG"  };
std::vector<TString>    kSyst_SEX_2D_Names      =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4",    "BKG"   };
std::vector<TString>    kSyst_SEX_2D_Legend     =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4",    "BKG"   };
//
//-//-//    PID
//
std::vector<TString>    kSyst_PID_XD_Options    =   {"PID1","PID2","PID3","PID4","PID5","PID6"};
std::vector<TString>    kSyst_PID_XD_Names      =   {"PID1","PID2","PID3","PID4","PID5","PID6"};
std::vector<TString>    kSyst_PID_XD_Legend     =   {"PID1","PID2","PID3","PID4","PID5","PID6"};
//
//-//-//    TRK
//
std::vector<TString>    kSyst_TRK_XD_Options    =   {"TRK1","TRK2","TRK3","TRK4","TRK5","TRK6","TRK7","TRK8","TRK9","TRK10","TRK11","TRK12"};//,"TRK13","TRK14"};
std::vector<TString>    kSyst_TRK_XD_Names      =   {"TRK1","TRK2","TRK3","TRK4","TRK5","TRK6","TRK7","TRK8","TRK9","TRK10","TRK11","TRK12"};//,"TRK13","TRK14"};
std::vector<TString>    kSyst_TRK_XD_Legend     =   {"TRK1","TRK2","TRK3","TRK4","TRK5","TRK6","TRK7","TRK8","TRK9","TRK10","TRK11","TRK12"};//,"TRK13","TRK14"};

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
//  --- File Names
//  --- --- Analysis
//  --- --- --- Pre Processing
TString const   kAnalysis_PreProc_Dir   =   TString("./result/%s/PreProcessing/");
TString const   kAnalysis_InvMassHist   =   kAnalysis_PreProc_Dir + TString("IM_Histograms.root");
TString const   kAnalysis_MCTruthHist   =   kAnalysis_PreProc_Dir + TString("IM_MonteCarloTruth.root");
TString const   kAnalysis_Plot_Direc    =   kAnalysis_PreProc_Dir + TString("Plots/");
//
//  --- --- --- Mass Resolution
TString const   kMassResolution_Dir_    =   TString("./result/%s/PreProcessing/");
TString const   kMassResolution_Prod    =   kMassResolution_Dir_ + TString("MR_Histograms.root");
TString const   kMassResolution_Anal    =   kMassResolution_Dir_ + TString("MR_Results.root");
TString const   kMassResolution_Plot    =   kMassResolution_Dir_ + TString("Plots/");
//
//  --- --- --- Signal Extraction
TString const   kAnalysis_SigExtr_Dir   =   TString("./result/%s/SignalExtraction/");
TString const   kAnalysis_SgExSys_Dir   =   TString("./result/%s/SignalExtraction/Systematics/");
TString const   kASigExtr_FitCheckPlt   =   kAnalysis_SigExtr_Dir + TString("FitCheck.root");
TString const   kASigExtr_FitChkPltSY   =   kAnalysis_SigExtr_Dir + TString("/%s/FitCheck.root");
TString const   kASigExtr_FitCheckRst   =   kAnalysis_SigExtr_Dir + TString("FitResults.root");
TString const   kASigExtr_FitChkRstSY   =   kAnalysis_SigExtr_Dir + TString("/%s/FitResults_%s.root");
TString const   kASigExtr_Plot_Direct   =   kAnalysis_SigExtr_Dir + TString("Plots/");
TString const   kASigExtrPlotDirectSY   =   kAnalysis_SigExtr_Dir + TString("/%s/Plots/");
//
//  --- --- --- Signal Extrapolation
TString const   kAnalysis_SigExtp_Dir   =   TString("./result/%s/SignalExtrapolation/");
TString const   kASigExtp_FitCheckPlt   =   kAnalysis_SigExtp_Dir + TString("FitCheck.root");
TString const   kASigExtp_FitCheckRst   =   kAnalysis_SigExtp_Dir + TString("FitResults.root");
TString const   kASigExtp_Plot_Direct   =   kAnalysis_SigExtp_Dir + TString("Plots/");
//
//  --- --- --- Systematics
TString const   kAnalysis_Systemt_Dir   =   TString("./result/%s/Systematics/");
TString const   kSystematicsPlot        =   kAnalysis_Systemt_Dir + TString("Plots/");
//
//  --- --- --- MonteCarlo Comparison
TString const   kProduction_MC_Dir      =   TString("./result/%s/MC_Production/");
TString const   kProduction_MC_Plot     =   kProduction_MC_Dir + TString("Plots/");
TString const   kProduction_MC_Ofl      =   kProduction_MC_Dir + TString("ComparePlots_%s.root");
//
//  --- --- --- Show Plots
TString const   kDIR_ShowPlots          =   TString("./result/%s/ShowPlots/");
//
//  --- --- --- Final Results
TString const   kDIR_FinalResults       =   TString("./result/FinalResults/%s/");
TString const   kDIR_FinalResultsPlots  =   kDIR_FinalResults + TString("/Plots/");
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Tree Names
//
auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaonCandidate_Tree     =   "KaonCandidate";
auto const  fKaonCandidateEff_Tree  =   "KaonEfficiency";
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- PDG Values
//
//-//-//    Particle Constants ( GEV )
//
auto const  kPhiMesonMass_          =   1.019455;
auto const  kPhiMesonMsErr          =   0.000020;
auto const  kPhiMesonWidth          =   0.004249;
auto const  kPhiMesonWdErr          =   0.000013;
auto const  kKaonMass               =   .493677;
auto const  kKaonMassUncert         =   .000013;
auto const  kTriggerEff             =   0.99;
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Binning
//
//-// InvMass bins 1D
const   Float_t     kBinningPrecision1D = 0.4 * MeV;     // From AliAnalysisUtility *MeV the value was in GeV is now in MeV, /MeV the value was in MeV and is now in GeV
const   Float_t     fMinIM1D  = 0.99;
const   Float_t     fMaxIM1D  = 1.08;
const   Int_t       nBinIM1D  = (int)((fMaxIM1D-fMinIM1D)/kBinningPrecision1D);
        Float_t*    fArrIM1D  = new Float_t [nBinIM1D+1];

//-// InvMass bins 2D
const   Float_t     kBinningPrecision2D = 0.8 * MeV;
const   Float_t     fMinIM2D    =   0.99;
const   Float_t     fMaxIM2D    =   1.08;
const   Int_t       nBinIM2D    =   (int)((fMaxIM2D-fMinIM2D)/kBinningPrecision2D);
        Float_t*    fArrIM2D    =   new Float_t [nBinIM2D+1];

//-// InvMass bins Resolution for REC
const   Float_t     kBinningPrecisionRC = 0.2 * MeV;
const   Float_t     fMinIMRC    =   1.00;
const   Float_t     fMaxIMRC    =   1.04;
const   Int_t       nBinIMRC    =   (int)((fMaxIMRC-fMinIMRC)/kBinningPrecisionRC);
        Float_t*    fArrIMRC    =   new Float_t [nBinIMRC+1];

//-// InvMass bins Resolution for TRUE
const   Float_t     kBinningPrecisionTR = 0.001 * MeV;
const   Float_t     fMinIMTR    =   1.00;
const   Float_t     fMaxIMTR    =   1.04;
const   Int_t       nBinIMTR    =   (int)((fMaxIMTR-fMinIMTR)/kBinningPrecisionTR);
        Float_t*    fArrIMTR    =   new Float_t [nBinIMTR+1];


//-// InvMass bins Resolution
const   Float_t     kBinningPrecisionRX = 0.2 * MeV;
const   Float_t     fMinIMRX    =   1.00;
const   Float_t     fMaxIMRX    =   1.04;
const   Int_t       nBinIMRX    =   (int)((fMaxIMRX-fMinIMRX)/kBinningPrecisionRX);
        Float_t*    fArrIMRX    =   new Float_t [nBinIMRX+1];

//-// pT bins 1D
        Int_t     nBinPT1D  =   -1;
        Float_t   fMinPT1D  =   -1;
        Float_t   fMaxPT1D  =   -1;
        Float_t  *fArrPT1D  =   new Float_t [1024];
        Float_t  *fArrPT1D_Comp  =   new Float_t [1024];

//-// pT bins 2D
        Int_t     nBinPT2D  =   -1;
        Float_t   fMinPT2D  =   -1;
        Float_t   fMaxPT2D  =   -1;
        Float_t  *fArrPT2D  =   new Float_t [1024];
        Float_t  *fArrPT2D_Comp  =   new Float_t [1024];

//-// Muliplicity bins
        Int_t     nBinMult  =   -1;
const   Float_t   fMinMult  =   0.0;
const   Float_t   fMaxMult  =   100.0;
        Float_t  *fArrMult  =   new Float_t [30];

//-// Muliplicity reference
        Int_t     nBinRMlt  =   -1;
        Float_t  *fArrRMlt  =   new Float_t [30];

//-// Phi-Correlation bins
        Int_t     nBinCrPh  =   5;
        Float_t   fMinCrPh  =   -180;
        Float_t   fMaxCrPh  =   180;
        Float_t  *fArrCrPh  =   new Float_t [1024];

//-// Utility Bins

//-// Rapidity bins
const   Int_t     nBinRap_  =   200;
const   Float_t   fMinRap_  =   -.5;
const   Float_t   fMaxRap_  =   0.5;
        Float_t  *fArrRap_  =   new Float_t [nBinRap_+1];

//-// N-Tuples bins
const   Int_t     nBinNTup  =   5;
const   Float_t   fMinNTup  =   -0.5;
const   Float_t   fMaxNTup  =   4.5;
        Float_t  *fArrNTup  =   new Float_t [nBinNTup+1];

//-// Systematics bins
const   Int_t     nBinSyst  =   120;
const   Float_t   fMinSyst  =   -.6;
const   Float_t   fMaxSyst  =   .6;
        Float_t  *fArrSyst  =   new Float_t [nBinSyst+1];

//-// Mass Resolution bins
const   Int_t     nBinIMRs  =   100;
const   Float_t   fMinIMRs  =   -.01;
const   Float_t   fMaxIMRs  =   .01;
        Float_t  *fArrIMRs  =   new Float_t [nBinIMRs+1];

//-// Mass Resolution bins
const   Int_t     nBinIMR2  =   40;
const   Float_t   fMinIMR2  =   -.01;
const   Float_t   fMaxIMR2  =   .01;
        Float_t  *fArrIMR2  =   new Float_t [nBinIMR2+1];



//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Dump Variables
//
Double_t                    fDumpVar_Double_t;
std::vector<TH1F*>          NULL_VECTOR;
RooFitResult   **           NULL_ROOFITPTR2  =   nullptr;
RooFitResult   ***          NULL_ROOFITPTR3  =   nullptr;

//------------------------------//
//       DATA STRUCTURES        //
//------------------------------//
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Trees
typedef struct  {
    //  --- Event Variables
    UChar_t     EventMask;
    Bool_t      kHasMult;
    Int_t       nPhi,       iMult;
    Float_t     Multiplicity, Spherocity, RTransverse;
    //  --- Pair Variables
    Bool_t      kHasRap [2048];
    Int_t       iPT1D   [2048], iPT2D   [2048], iRap    [2048], iKaon   [2048], jKaon   [2048];
    Float_t     Px      [2048], Py      [2048], Pz      [2048], pT  [2048], Rap [2048], InvMass [2048], TrueInvMass [2048], Phi [2048];
} Struct_PhiCandidate;
typedef struct  {
    //  --- Event Variables
    Bool_t      IsMB;
    Int_t       nPart;
    Float_t     Eta_10, Eta_08, Eta_05, V0M, V0A;
    //  --- Pair Variables
    Bool_t      kHasRap [2048];
    Int_t       iPT1D   [2048], iPT2D   [2048];
    Float_t     Px      [2048], Py      [2048], Pz      [2048], pT  [2048], Rap [2048], Mass [2048], Eta [2048], Phi [2048];
} Struct_MCParticle;
typedef struct  {
    UChar_t     EventMask,      Charge[2048];
    Char_t      SigmaTOF[2048], SigmaTPC[2048];
    Float_t     Multiplicity,   Px[2048],       Py[2048],       Pz[2048],   InvMass[2048],       Phi[2048];
    Int_t       nKaon;
    //UChar_t     nKaon;
} Struct_KaonCandidate;
typedef struct  {
    UChar_t     TrueEventMask,  EventMask,      Selection[2048];
    Float_t     Multiplicity,   Px[2048],       Py[2048],       Pz[2048],       InvMass[2048],  pT[2048],      Rap[2048],   Phi[2048];
    Bool_t      kHasRap[2048],  kIsGen[2048],   kIsRec[2048],   kHasMult;
    Bool_t      fTru,           fGen,           fRec,           IsMB,           IsVTX10;
    Int_t       iPT1D[2048],    iPT2D[2048],    iRap[2048],     iMult;
    Int_t       nPhi;
    //UChar_t     nPhi;
} Struct_PhiEfficiency;
typedef struct  {
    UChar_t     TrueEventMask,  EventMask,      Charge[2048],   Selection[2048];
    Float_t     Multiplicity,   Px[2048],       Py[2048],       Pz[2048],   InvMass[2048];
    Bool_t      nKaon,          ftru;
} Struct_KaonEfficiency;
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Event Utilities
enum class kEventMask {
    kVoid1, kPileUp, kPileUpMult, kINELgt0, kVoid5, kVoid6, kVoid7, kVoid8
};
enum kEventCount {
    kUnderFlow = 0, kALL = 1, kHasEvent = 2, kHasMCTracks = 3, kNoTrigger = 4, kTrigger = 5, kHasPID = 5, kIncmpDAQ = 6, kNoSPDVtx = 7, kVtxMismatch = 8, kVertex = 9, kVertex10 = 10, kPU_MB = 12, kPU_Mlt = 13
};
//
std::vector<std::pair<TString,TString>> kpp7TeVDataset = {std::pair<TString,TString>("LHC10b","LHC14j4b"),std::pair<TString,TString>("LHC10c","LHC14j4c"),std::pair<TString,TString>("LHC10d","LHC14j4d"),std::pair<TString,TString>("LHC10e","LHC14j4e")};
std::vector<std::pair<TString,TString>> kpp5TeVDataset = {std::pair<TString,TString>("LHC15TeV_DT","LHC15TeV_MC")};
//
//  Option choosing accepted strings
std::vector<TString>    kOptStrings_All             =   {"all","full","fll"};
std::vector<TString>    kOptStrings_Yield           =   {"yield","yld","std","standard"};
std::vector<TString>    kOptStrings_Multiplicity    =   {"multiplicity","mult","mlt"};
std::vector<TString>    kOptStrings_Correlation     =   {"correlation","corr","crph"};

//------------------------------//
//      BINNING UTILITIES       //
//------------------------------//
//
//-//   Setters
void
fSetBinIM1D
()      {
    for (int i = 0; i <= nBinIM1D; i++ )
    {
        fArrIM1D[i] = fMinIM1D+(i)*(fMaxIM1D - fMinIM1D)/(static_cast<Float_t>(nBinIM1D));
    }
}
void
fSetBinIM2D
()      {
    for (int i = 0; i <= nBinIM2D; i++ )
    {
        fArrIM2D[i] = fMinIM2D+(i)*(fMaxIM2D - fMinIM2D)/(static_cast<Float_t>(nBinIM2D));
    }
}
void
fSetBinIMRX
()      {
    for (int i = 0; i <= nBinIMRX; i++ )
    {
        fArrIMRX[i] = fMinIMRX+(i)*(fMaxIMRX - fMinIMRX)/(static_cast<Float_t>(nBinIMRX));
    }
}
void
fSetBinIMRC
()      {
    for (int i = 0; i <= nBinIMRC; i++ )
    {
        fArrIMRC[i] = fMinIMRC+(i)*(fMaxIMRC - fMinIMRC)/(static_cast<Float_t>(nBinIMRC));
    }
}
void
fSetBinIMTR
()      {
    for (int i = 0; i <= nBinIMTR; i++ )
    {
        fArrIMTR[i] = fMinIMTR+(i)*(fMaxIMTR - fMinIMTR)/(static_cast<Float_t>(nBinIMTR));
    }
}
void
fSetBinPT1D
()      {
    //  --- QM2022 Preliminary Binning pp@7TeV
    nBinPT1D        =   20;
    fArrPT1D[0]     =   0.4; //0.2
    fArrPT1D[1]     =   0.6; //0.1
    fArrPT1D[2]     =   0.7; //0.1
    fArrPT1D[3]     =   0.8; //0.1
    fArrPT1D[4]     =   0.9; //0.1
    fArrPT1D[5]     =   1.0; //0.2
    fArrPT1D[6]     =   1.1; //0.2
    fArrPT1D[7]     =   1.2; //0.2
    fArrPT1D[8]     =   1.3; //0.2
    fArrPT1D[9]     =   1.4; //0.2
    fArrPT1D[10]    =   1.6; //0.2
    fArrPT1D[11]    =   1.8; //0.2
    fArrPT1D[12]    =   2.0; //0.4
    fArrPT1D[13]    =   2.4; //0.4
    fArrPT1D[14]    =   2.8; //0.4
    fArrPT1D[15]    =   3.2; //0.4
    fArrPT1D[16]    =   3.6; //0.4
    fArrPT1D[17]    =   4.0; //1.0
    fArrPT1D[18]    =   5.0; //1.0
    fArrPT1D[19]    =   6.0; //2.0
    fArrPT1D[20]    =   8.0; //---
    //  --- Preliminary Binning pp@5TeV in mult ( Sushanta )
    nBinPT1D        =   17;
    fArrPT1D[0]     =   0.4; //0.2
    fArrPT1D[1]     =   0.6; //0.2
    fArrPT1D[2]     =   0.8; //0.1
    fArrPT1D[3]     =   1.0; //0.2
    fArrPT1D[4]     =   1.2; //0.2
    fArrPT1D[5]     =   1.4; //0.2
    fArrPT1D[6]     =   1.6; //0.2
    fArrPT1D[7]     =   1.8; //0.2
    fArrPT1D[8]     =   2.0; //0.4
    fArrPT1D[9]     =   2.4; //0.4
    fArrPT1D[10]    =   2.8; //0.4
    fArrPT1D[11]    =   3.2; //0.4
    fArrPT1D[12]    =   3.6; //0.4
    fArrPT1D[13]    =   4.0; //1.0
    fArrPT1D[14]    =   5.0; //1.0
    fArrPT1D[15]    =   6.0; //2.0
    fArrPT1D[16]    =   8.0; //2.0
    fArrPT1D[17]    =   10.; //---
    //  --- Maximum and Minimum assignment
    fMinPT1D        =   fArrPT1D[0];
    fMaxPT1D        =   fArrPT1D[nBinPT1D];
    fArrPT1D_Comp[0]=   0.00;
    for ( Int_t iPT1D = 0; iPT1D <= nBinPT1D; iPT1D++ ) fArrPT1D_Comp[iPT1D+1]  =   fArrPT1D[iPT1D];
}
void
fSetBinPT2D
()      {
    //  --- QM2022 Preliminary Binning pp@7TeV
    nBinPT2D        =   8;
    fArrPT2D[0]     =   0.40; //0.5
    fArrPT2D[1]     =   0.90; //0.3
    fArrPT2D[2]     =   1.20; //0.2
    fArrPT2D[3]     =   1.40; //0.3
    fArrPT2D[4]     =   1.70; //0.3
    fArrPT2D[5]     =   2.00; //0.5
    fArrPT2D[6]     =   2.50; //1.5
    fArrPT2D[7]     =   4.00; //4.0
    fArrPT2D[8]     =   10.0; //---
    //  --- Maximum and Minimum assignment
    fMinPT2D        =   fArrPT2D[0];
    fMaxPT2D        =   fArrPT2D[nBinPT2D];
    fArrPT2D_Comp[0]=   0.00;
    for ( Int_t iPT2D = 0; iPT2D <= nBinPT2D; iPT2D++ ) fArrPT2D_Comp[iPT2D+1]  =   fArrPT2D[iPT2D];
}
void
fSetBinMult
()      {
    if ( is_pp_anl )    {
        //  --- Standard Binning pp@5TeV
        nBinMult     =  5;
        fArrMult[0]  =  0.00;
        fArrMult[1]  =  5.00;
        fArrMult[2]  =  15.0;
        fArrMult[3]  =  30.0;
        fArrMult[4]  =  50.0;
        fArrMult[5]  =  100.0;
        //  --- Standard Binning pp@13TeV
        nBinMult     =  7;
        fArrMult[0]  =  0.00;
        fArrMult[1]  =  1.00;
        fArrMult[2]  =  5.00;
        fArrMult[3]  =  10.0;
        fArrMult[4]  =  15.0;
        fArrMult[5]  =  30.0;
        fArrMult[6]  =  50.0;
        fArrMult[7]  =  100.0;
        //  --- Preliminary Binning pp@5TeV in mult ( Sushanta )
        nBinMult     =  9;
        fArrMult[0]  =  0.00;
        fArrMult[1]  =  5.00;
        fArrMult[2]  =  10.0;
        fArrMult[3]  =  15.0;
        fArrMult[4]  =  20.0;
        fArrMult[5]  =  30.0;
        fArrMult[6]  =  40.0;
        fArrMult[7]  =  50.0;
        fArrMult[8]  =  70.0;
        fArrMult[9]  =  100.0;
    } else if ( is_pb_anl )    {
        nBinMult     =  3;
        fArrMult[0]  =  0.00;   //
        fArrMult[1]  =  15.0;   // 57.23
        fArrMult[2]  =  40.0;   // 37.08
        fArrMult[3]  =  100.;   // 15.95
    }
}
void
fSetBinRMlt
()      {
    if ( is_pp_anl )    {
        //  --- Standard Binning pp@5TeV
        nBinRMlt     =  5;
        fArrRMlt[0]  =  0;
        fArrRMlt[1]  =  15.48;
        fArrRMlt[2]  =  (12.07+10.40)/2.;
        fArrRMlt[3]  =  (9.17+7.76*2)/3.;
        fArrRMlt[4]  =  (6.30+5.16)/2.;
        fArrRMlt[5]  =  (3.90*2+2.38*3)/5.;
        //  --- Standard Binning pp@13TeV
        nBinRMlt     =  7;
        fArrRMlt[0]  =  0;
        fArrRMlt[1]  =  26.02;
        fArrRMlt[2]  =  20.02;
        fArrRMlt[3]  =  16.17;
        fArrRMlt[4]  =  13.77;
        fArrRMlt[5]  =  (12.04+10.02*2)/3.;
        fArrRMlt[6]  =  (7.95+6.32)/2.;
        fArrRMlt[7]  =  (4.50*2+2.55*3)/5.;
    } else if ( is_pb_anl )    {
        nBinRMlt     =  3;
        fArrRMlt[0]  =  0.00;   //
        fArrRMlt[1]  =  (44.96+36.05+30.33)/3.;  // 57.23
        fArrRMlt[2]  =  (30.33+23.12*4)/5.;  // 37.08
        fArrRMlt[3]  =  (15.89+9.63+4.13)/3.;  // 15.95
    }
}
void
fSetBinCrPh
()      {
    //  --- Standard Binning pp@13TeV
    nBinCrPh  =   8;
    fArrCrPh[0]  =  -105.;
    fArrCrPh[1]  =  -60.;
    fArrCrPh[2]  =  -15.;
    fArrCrPh[3]  =  +30.;
    fArrCrPh[4]  =  +75.;
    fArrCrPh[5]  =  +120.;
    fArrCrPh[6]  =  +165.;
    fArrCrPh[7]  =  +210.;
    fArrCrPh[8]  =  +255.;
    fMinCrPh  =   fArrCrPh[0];
    fMaxCrPh  =   fArrCrPh[nBinCrPh];
}
void
fSetBinRap_
()      {
    for (int i = 0; i <= nBinRap_; i++ )
    {
        fArrRap_[i] = fMinRap_+(i)*(fMaxRap_ - fMinRap_)/(static_cast<Float_t>(nBinRap_));
    }
    
    fArrRap_[0]  =  0.000;
    fArrRap_[1]  =  0.080;
    fArrRap_[2]  =  0.170;
    fArrRap_[3]  =  0.290;
    fArrRap_[4]  =  0.480;
    fArrRap_[5]  =  1.000;
    
}
void
fSetBinNTup
()      {
    for (int i = 0; i <= nBinNTup; i++ )
    {
        fArrNTup[i] = fMinNTup+(i)*(fMaxNTup - fMinNTup)/(static_cast<Float_t>(nBinNTup));
    }
}
void
fSetBinSyst
()      {
    for (int i = 0; i <= nBinSyst; i++ )
    {
        fArrSyst[i] = fMinSyst+(i)*(fMaxSyst - fMinSyst)/(static_cast<Float_t>(nBinSyst));
    }
}
void
fSetBinIMRs
()      {
    for (int i = 0; i <= nBinIMRs; i++ )
    {
        fArrIMRs[i] = fMinIMRs+(i)*(fMaxIMRs - fMinIMRs)/(static_cast<Float_t>(nBinIMRs));
    }
}
void
fSetBinIMR2
()      {
    for (int i = 0; i <= nBinIMR2; i++ )
    {
        fArrIMR2[i] = fMinIMR2+(i)*(fMaxIMR2 - fMinIMR2)/(static_cast<Float_t>(nBinIMR2));
    }
}
void
fSetAllBins
()      {
    fSetBinIM1D();
    fSetBinIM2D();
    fSetBinIMRX();
    fSetBinIMRC();
    fSetBinIMTR();
    fSetBinPT1D();
    fSetBinPT2D();
    fSetBinCrPh();
    fSetBinNTup();
    fSetBinMult();
    fSetBinRap_();
    fSetBinSyst();
    fSetBinIMRs();
    fSetBinIMR2();
    fSetBinRMlt();
}
//
//-//   Getters
Int_t
fGetBinIM1D
 (Float_t input_value )      {
    if ( input_value > fMaxIM1D ) return -1;
    for ( Int_t iBin = 0; iBin <= nBinIM1D; iBin++ )
    {
        if ( input_value <= fArrIM1D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinIM2D
 (Float_t input_value )      {
    if ( input_value > fMaxIM2D ) return -1;
    for ( Int_t iBin = 0; iBin <= nBinIM2D; iBin++ )
    {
        if ( input_value <= fArrIM2D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinPT1D
 (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinPT1D; iBin++ )
    {
        if ( input_value <= fArrPT1D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinPT2D
 (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinPT2D; iBin++ )
    {
        if ( input_value <= fArrPT2D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinMult
 (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinMult; iBin++ )
    {
        if ( input_value == fMinMult ) return 0;
        if ( input_value <= fArrMult[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinCrPh
 (Float_t input_value )      {
    auto    kDeltaPhi   =   ( input_value < fMinCrPh ? input_value + 360 : input_value > fMaxCrPh ? input_value -360 : input_value );
    //
    for ( Int_t iBin = 0; iBin <= nBinCrPh; iBin++ ) {
        if ( kDeltaPhi == fMinCrPh )          return 0;
        if ( kDeltaPhi <= fArrCrPh[iBin] )    return iBin -1;
    }
    return -1;
}
Float_t
fGetDltCrPh
 ( Float_t input_value )      {
    return ( input_value < fMinCrPh ? input_value + 360 : input_value > fMaxCrPh ? input_value -360 : input_value );
}
Int_t
fGetBinRap_
 (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinRap_; iBin++ )
    {
        if ( fabs(input_value) <= fArrRap_[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinNTup
 (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinNTup; iBin++ )
    {
        if ( input_value <= fArrNTup[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t
fGetBinMultEff
 (Float_t input_value )      {
    return -1;
    /*
    for ( Int_t iBin = 0; iBin <= nMltTrgECls; iBin++ )
    {
        if ( input_value <= kMltTrgECls[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
     */
}

//------------------------------//
//      CUTS UTILITIES          //
//------------------------------//
//
bool
fCutRapidity
 ( Double_t  dRapidity )         {
    if ( fabs(dRapidity) < 0.5 ) return true;
    return false;
}
bool
fCutInvariantMass
 ( Double_t  dInvariantMass )   {
    if ( dInvariantMass < fMinIM1D ) return false;
    if ( dInvariantMass > fMaxIM1D ) return false;
    return true;
}
bool
fCutTransverseMom
 ( Double_t  dTransverseMom )    {
    if ( dTransverseMom < fMinPT1D ) return false;
    if ( dTransverseMom > fMaxPT1D ) return false;
    return true;
}
bool
fCutMultiplicity
 ( Double_t  dMultiplicity )     {
    return true;
    if ( dMultiplicity < fMinMult ) return false;
    if ( dMultiplicity > fMaxMult ) return false;
    return true;
}
bool
fAcceptCandidate
 ( Double_t  dInvariantMass, Double_t dTransverseMom )  {
    if ( !fCutInvariantMass(dInvariantMass) ) return false;
    if ( !fCutTransverseMom(dTransverseMom) ) return false;
    return true;
}
bool
fCheckCoupleKaons
 ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi ) {
    if ( SPhiCandidates.iKaon[aAccCandidates[iPhi]] == SPhiCandidates.iKaon[aAccCandidates[jPhi]] ) return false;
    if ( SPhiCandidates.jKaon[aAccCandidates[iPhi]] == SPhiCandidates.jKaon[aAccCandidates[jPhi]] ) return false;
    if ( SPhiCandidates.iKaon[aAccCandidates[iPhi]] == SPhiCandidates.jKaon[aAccCandidates[jPhi]] ) return false;
    if ( SPhiCandidates.jKaon[aAccCandidates[iPhi]] == SPhiCandidates.iKaon[aAccCandidates[jPhi]] ) return false;
    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi )             {
    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi ) {
    // Non equal candidates
    if ( iPhi == jPhi ) return false;
                
    // Only non overlapping couples of Kaons
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,iPhi,jPhi) ) return false;

    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi, Int_t kPhi )             {
    // Non equal candidates
    if ( iPhi == kPhi ) return false;
    if ( jPhi == kPhi ) return false;
                
    // Only non overlapping couples of Kaons
    // >> iPhi vs kPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,iPhi,kPhi) ) return false;
    // >> jPhi vs kPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,jPhi,kPhi) ) return false;
    
    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi, Int_t kPhi, Int_t lPhi ) {
    // Non equal candidates
    if ( iPhi == lPhi ) return false;
    if ( jPhi == lPhi ) return false;
    if ( kPhi == lPhi ) return false;
                
    // Only non overlapping couples of Kaons
    // >> iPhi vs lPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,iPhi,lPhi) ) return false;
    // >> jPhi vs lPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,jPhi,lPhi) ) return false;
    // >> kPhi vs lPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,kPhi,lPhi) ) return false;
    
    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t iPhi )             {
    return true;
}
bool
fCheckCoupleKaons
 ( Struct_PhiCandidate SPhiCandidates, Int_t iPhi, Int_t jPhi ) {
    if ( SPhiCandidates.iKaon[iPhi] == SPhiCandidates.iKaon[jPhi] ) return false;
    if ( SPhiCandidates.jKaon[iPhi] == SPhiCandidates.jKaon[jPhi] ) return false;
    if ( SPhiCandidates.iKaon[iPhi] == SPhiCandidates.jKaon[jPhi] ) return false;
    if ( SPhiCandidates.jKaon[iPhi] == SPhiCandidates.iKaon[jPhi] ) return false;
    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t iPhi, Int_t jPhi ) {
    // Non equal candidates
    if ( iPhi == jPhi ) return false;
                
    // Only non overlapping couples of Kaons
    if ( !fCheckCoupleKaons(SPhiCandidates,iPhi,jPhi) ) return false;

    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t iPhi, Int_t jPhi, Int_t kPhi )             {
    // Non equal candidates
    if ( iPhi == kPhi ) return false;
    if ( jPhi == kPhi ) return false;
                
    // Only non overlapping couples of Kaons
    // >> iPhi vs kPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,iPhi,kPhi) ) return false;
    // >> jPhi vs kPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,jPhi,kPhi) ) return false;
    
    return true;
}
bool
fAcceptCandidate
 ( Struct_PhiCandidate SPhiCandidates, Int_t iPhi, Int_t jPhi, Int_t kPhi, Int_t lPhi ) {
    // Non equal candidates
    if ( iPhi == lPhi ) return false;
    if ( jPhi == lPhi ) return false;
    if ( kPhi == lPhi ) return false;
                
    // Only non overlapping couples of Kaons
    // >> iPhi vs lPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,iPhi,lPhi) ) return false;
    // >> jPhi vs lPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,jPhi,lPhi) ) return false;
    // >> kPhi vs lPhi
    if ( !fCheckCoupleKaons(SPhiCandidates,kPhi,lPhi) ) return false;
    
    return true;
}

#endif
