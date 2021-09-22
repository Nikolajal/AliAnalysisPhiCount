#ifndef ALIANALYSISPHIPAIR_CONSTANTS_H
#define ALIANALYSISPHIPAIR_CONSTANTS_H

#include "./AliAnalysisPhiPair.h"

//------------------------------//
//      !TODO: CLEAN - UP       //
//------------------------------//
//
//>->->->->->   Branching Ratio
auto const  kBR                     =   0.492;
auto const  kBRerr                  =   0.005;
auto const  kSysLow_BR              =   kBRerr/kBR;
auto const  kSysHig_BR              =   kBRerr/kBR;
//
//>->->->->->   Signal Extraction
auto const  kSysLow_1D_SE           =   0.06;
auto const  kSysHig_1D_SE           =   0.06;
auto const  kSysLow_2D_SE           =   0.15;
auto const  kSysHig_2D_SE           =   0.15;
//
//>->->->->->   Tracking Efficiency
auto const  kSysLow_TR              =   0.08;
auto const  kSysHig_TR              =   0.08;
//
//>->->->->->   PID Efficiency
auto const  kSysLow_PD              =   0.015;
auto const  kSysHig_PD              =   0.015;
//
auto const  kTriggerEff             =   0.85; // 5TeV 0.7574;// 7TeV 0.85;
auto const  nMltTrgECls             =   10;
Float_t const  kMultTrgEff      []  =   {.998822,0.995576,0.991524,0.986489,0.975743,0.9575743,0.937151,0.897753,0.696985};
Float_t const  kMltTrgECls      []  =   {0,5,10,15,20,30,40,50,70,100};
auto const  kSysLow_Trigger         =   0.01;
auto const  kSysHig_Trigger         =   0.01;
//
auto const  kTrackingEff            =   -1.0; //Unused
auto const  kSysLow_Tracking        =   0.01;
auto const  kSysHig_Tracking        =   0.01;

//
auto const  kSignalMiss1D           =   0.960419;
auto const  kSignalMiss2D           =   0.922404;

//-//   Analysis settings
auto        kDoMultiplicity         =   false;
auto        kDoYield                =   false;
auto        kDoTrigger              =   false;
auto        kDoRapidity             =   false;

//-// InvMass range Pythia MC
const   Float_t   fMinIMMC  =   0.75;
const   Float_t   fMaxIMMC  =   1.25;



//-// Systematics Options
//
std::vector<TString>    sOptions    =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4"  };
const   Int_t   nOptions            =   sOptions.size();
std::vector<TString>    sOption2    =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ", "RK","RL","RM","RN","WDT","RSH","RSL",/*"RSX",*/"DG2","DG4",  "BKG"   };
const   Int_t   nOption2            =   sOption2.size();
//
//-//-//    PID
//
const   Int_t   nPIDFiles           =   6;
std::vector<TString>    sOptiPID    =   {"#sigma_{TPC}^{vet} 3.3","#sigma_{TPC}^{vet} 2.7","#sigma_{TPC}^{aln} 5.5","#sigma_{TPC}^{aln} 4.5","#sigma_{TOF}^{vet} 3.3","#sigma_{TOF}^{vet} 2.7"};
const   TString sPID_DT_Name        =   "20210805_LHC10X%s.root";
const   TString sPID_MC_Name        =   "20210805_LHC14j4X%s.root";


//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- File Names
//
//-//-//    Analysis
//
//-//-//-//     Pre Processing
//
TString const   kAnalysis_PreProc_Dir   =   TString("./result/%s/PreProcessing/");
TString const   kAnalysis_InvMassHist   =   kAnalysis_PreProc_Dir + TString("IM_Histograms.root");
TString const   kAnalysis_MCTruthHist   =   kAnalysis_PreProc_Dir + TString("IM_MonteCarloTruth.root");
TString const   kAnalysis_Plot_Direc    =   kAnalysis_PreProc_Dir + TString("Plots/");
//
//-//-//-//     Mass Resolution
//
TString const   kMassResolution_Dir_    =   TString("./result/%s/MassResolution/");
TString const   kMassResolution_Prod    =   kMassResolution_Dir_ + TString("MR_Histograms.root");
TString const   kMassResolution_Anal    =   kMassResolution_Dir_ + TString("MR_Results.root");
TString const   kMassResolution_Plot    =   kMassResolution_Dir_ + TString("Plots/");
//
//-//-//-//     Signal Extraction
//
TString const   kAnalysis_SigExtr_Dir   =   TString("./result/%s/SignalExtraction/");
TString const   kAnalysis_SgExSys_Dir   =   TString("./result/%s/SignalExtraction/Systematics/");
TString const   kASigExtr_FitCheckPlt   =   kAnalysis_SigExtr_Dir + TString("FitCheck.root");
TString const   kASigExtr_FitCheckRst   =   kAnalysis_SigExtr_Dir + TString("FitResults.root");
TString const   kASigExtr_Plot_Direct   =   kAnalysis_SigExtr_Dir + TString("Plots/");
//
//-//-//-//     Signal Extrapolation
//
TString const   kAnalysis_SigExtp_Dir   =   TString("./result/%s/SignalExtrapolation/");
TString const   kASigExtp_FitCheckPlt   =   kAnalysis_SigExtp_Dir + TString("FitCheck.root");
TString const   kASigExtp_FitCheckRst   =   kAnalysis_SigExtp_Dir + TString("FitResults.root");
TString const   kASigExtp_Plot_Direct   =   kAnalysis_SigExtp_Dir + TString("Plots/");
//
//-//-//-//     Systematics
//
TString const   kAnalysis_Systemt_Dir   =   TString("./result/%s/Systematics/");
TString const   kSystematicsPlot        =   kAnalysis_Systemt_Dir + TString("Plots/");
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
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Binning
//
//-// InvMass bins 1D
const   Float_t     kBinningPrecision1D =   .5*MeV;     // From AliAnalysisUtility *MeV the value was in GeV is now in MeV, /MeV the value was in MeV and is now in GeV
const   Float_t     fMinIM1D  = 0.99;
const   Float_t     fMaxIM1D  = 1.08;
const   Int_t       nBinIM1D  = (int)((fMaxIM1D-fMinIM1D)/kBinningPrecision1D);
        Float_t    *fArrIM1D  = new Float_t [nBinIM1D+1];

//-// InvMass bins 2D
auto const  kBinningPrecision2D     =   1.*MeV;
const   Float_t   fMinIM2D  = 0.99;
const   Float_t   fMaxIM2D  = 1.08;
const   Int_t     nBinIM2D  = (int)((fMaxIM2D-fMinIM2D)/kBinningPrecision2D);
        Float_t * fArrIM2D  = new Float_t [nBinIM2D+1];

//-// pT bins 1D
const   Int_t     nBinPT1D  =   20;
const   Float_t   fMinPT1D  =   0.4;
const   Float_t   fMaxPT1D  =   10.;
        Float_t  *fArrPT1D  =   new Float_t [nBinPT1D+1];

//-// pT bins 2D
const   Int_t     nBinPT2D  =   10;
const   Float_t   fMinPT2D  =   0.4;
const   Float_t   fMaxPT2D  =   20.;
        Float_t  *fArrPT2D  =   new Float_t [nBinPT2D+1];

//-// Muliplicity bins
const   Int_t     nBinMult  =   10;
const   Float_t   fMinMult  =   0.0;
const   Float_t   fMaxMult  =   100.0;
        Float_t  *fArrMult  =   new Float_t [nBinMult+1];

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

//-// N-Tuples bins
const   Int_t     nBinSyst  =   120;
const   Float_t   fMinSyst  =   -.6;
const   Float_t   fMaxSyst  =   .6;
        Float_t  *fArrSyst  =   new Float_t [nBinSyst+1];

//-// N-Tuples bins
const   Int_t     nBinIMRs  =   100;
const   Float_t   fMinIMRs  =   -.01;
const   Float_t   fMaxIMRs  =   .01;
        Float_t  *fArrIMRs  =   new Float_t [nBinIMRs+1];

//-// N-Tuples bins
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
    UChar_t     nPhi,           EventMask,      iKaon[1024],    jKaon[1024],   Nature[1024];
    Float_t     Multiplicity,   Px[1024],       Py[1024],       Pz[1024],       pT[1024],      Rap[1024],       InvMass[1024],       TrueInvMass[1024];
    Bool_t      kHasRap[1024],  kHasMult;
    Int_t       iPT1D[1024],    iPT2D[1024],    iRap[1024],     iMult;
} Struct_PhiCandidate;
typedef struct  {
    UChar_t nKaon,          EventMask,      Charge[1024];
    Char_t  SigmaTOF[1024], SigmaTPC[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_KaonCandidate;
typedef struct  {
    UChar_t nPhi,           TrueEventMask,  EventMask,      Selection[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  fTru,           fGen,           fRec;
} Struct_PhiEfficiency;
typedef struct  {
    UChar_t nKaon,          TrueEventMask,  EventMask,      Charge[1024],   Selection[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  ftru;
} Struct_KaonEfficiency;
//
//-// -- -- -- -- -- -- -- -- -- -- -- -- -- Event Utilities
enum class kEventMask {
    kVoid1, kPileUp, kPileUpMult, kINELgt0, kVoid5, kVoid6, kVoid7, kVoid8
};
enum kEventCount {
    kUnderFlow = 0, kALL = 1, kHasEvent = 2, kHasMCTracks = 3, kNoTrigger = 4, kTrigger = 5, kHasPID = 5, kIncmpDAQ = 6, kNoSPDVtx = 7, kVtxMismatch = 8, kVertex = 9, kVertex10 = 10
};

//------------------------------//
//      BINNING UTILITIES       //
//------------------------------//
//
//-//   Setters
void    fSetBinIM1D ()      {
    for (int i = 0; i <= nBinIM1D; i++ )
    {
        fArrIM1D[i] = fMinIM1D+(i)*(fMaxIM1D - fMinIM1D)/(static_cast<Float_t>(nBinIM1D));
    }
}
void    fSetBinIM2D ()      {
    for (int i = 0; i <= nBinIM2D; i++ )
    {
        fArrIM2D[i] = fMinIM2D+(i)*(fMaxIM2D - fMinIM2D)/(static_cast<Float_t>(nBinIM2D));
    }
}
void    fSetBinNTup ()      {
    for (int i = 0; i <= nBinNTup; i++ )
    {
        fArrNTup[i] = fMinNTup+(i)*(fMaxNTup - fMinNTup)/(static_cast<Float_t>(nBinNTup));
    }
}
void    fSetBinPT1D ()      {
    /*
    fArrPT1D[0]     =   0.5; //0.1
    fArrPT1D[1]     =   0.7; //0.1
    fArrPT1D[2]     =   0.9; //0.1
    fArrPT1D[3]     =   1.2; //0.1
    fArrPT1D[4]     =   1.4; //0.1
    fArrPT1D[5]     =   1.6; //0.1
    fArrPT1D[6]     =   1.8; //0.2
    fArrPT1D[7]     =   2.0; //0.2
    fArrPT1D[8]     =   2.2; //0.2
    fArrPT1D[9]     =   2.6; //0.2
    fArrPT1D[10]    =   3.0; //0.2
    fArrPT1D[11]    =   3.5; //0.4
    fArrPT1D[12]    =   4.0; //0.4
    fArrPT1D[13]    =   5.0; //0.4
    fArrPT1D[14]    =   8.0; //0.4
    fArrPT1D[15]    =   12.; //0.4
    */
    /*
    fArrPT1D[0]     =   0.4; //0.1
    fArrPT1D[1]     =   0.6; //0.1
    fArrPT1D[2]     =   0.8; //0.1
    fArrPT1D[3]     =   1.0; //0.1
    fArrPT1D[4]     =   1.2; //0.1
    fArrPT1D[5]     =   1.4; //0.1
    fArrPT1D[6]     =   1.6; //0.2
    fArrPT1D[7]     =   1.8; //0.2
    fArrPT1D[8]     =   2.0; //0.2
    fArrPT1D[9]     =   2.5; //0.2
    fArrPT1D[10]    =   3.0; //0.2
    fArrPT1D[11]    =   3.5; //0.4
    fArrPT1D[12]    =   4.0; //0.4
    fArrPT1D[13]    =   4.5; //0.4
    fArrPT1D[14]    =   5.0; //0.4
    fArrPT1D[15]    =   6.0; //0.4
    fArrPT1D[16]    =   7.0; //0.4
    fArrPT1D[17]    =   8.0; //0.4
    fArrPT1D[18]    =   10.; //0.4
    fArrPT1D[19]    =   13.; //0.4
    fArrPT1D[20]    =   16.; //0.4
    fArrPT1D[21]    =   21.; //0.4
    */
    fArrPT1D[0]     =   0.4; //0.1
    fArrPT1D[1]     =   0.5; //0.1
    fArrPT1D[2]     =   0.6; //0.1
    fArrPT1D[3]     =   0.7; //0.1
    fArrPT1D[4]     =   0.8; //0.1
    fArrPT1D[5]     =   0.9; //0.1
    fArrPT1D[6]     =   1.0; //0.2
    fArrPT1D[7]     =   1.2; //0.2
    fArrPT1D[8]     =   1.4; //0.2
    fArrPT1D[9]     =   1.6; //0.2
    fArrPT1D[10]    =   1.8; //0.2
    fArrPT1D[11]    =   2.0; //0.4
    fArrPT1D[12]    =   2.4; //0.4
    fArrPT1D[13]    =   2.8; //0.4
    fArrPT1D[14]    =   3.2; //0.4
    fArrPT1D[15]    =   3.6; //0.4
    fArrPT1D[16]    =   4.0; //1.0
    fArrPT1D[17]    =   5.0; //1.0
    fArrPT1D[18]    =   6.0; //2.0
    fArrPT1D[19]    =   8.0; //2.0
    fArrPT1D[20]    =   10.;
}
void    fSetBinPT2D ()      {
    fArrPT2D[0]     =   0.40; //0.3
    fArrPT2D[1]     =   0.80; //0.2
    fArrPT2D[2]     =   1.00; //0.1
    fArrPT2D[3]     =   1.10; //0.1
    fArrPT2D[4]     =   1.20; //0.2
    fArrPT2D[5]     =   1.40; //0.2
    fArrPT2D[6]     =   1.60; //0.4
    fArrPT2D[7]     =   2.00; //0.8
    fArrPT2D[8]     =   2.80; //1.2
    fArrPT2D[9]     =   4.00; //6.0
    fArrPT2D[10]    =   10.0;
}
void    fSetBinMult ()      {
    fArrMult[0]  =  0.00;
    fArrMult[1]  =  1.00;
    fArrMult[2]  =  5.00;
    fArrMult[3]  =  10.0;
    fArrMult[4]  =  15.0;
    fArrMult[5]  =  20.0;
    fArrMult[6]  =  30.0;
    fArrMult[7]  =  40.0;
    fArrMult[8]  =  50.0;
    fArrMult[9]  =  70.0;
    fArrMult[10] =  100.;
}
void    fSetBinRap_ ()      {
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
void    fSetBinSyst ()      {
    for (int i = 0; i <= nBinSyst; i++ )
    {
        fArrSyst[i] = fMinSyst+(i)*(fMaxSyst - fMinSyst)/(static_cast<Float_t>(nBinSyst));
    }
}
void    fSetBinIMRs ()      {
    for (int i = 0; i <= nBinIMRs; i++ )
    {
        fArrIMRs[i] = fMinIMRs+(i)*(fMaxIMRs - fMinIMRs)/(static_cast<Float_t>(nBinIMRs));
    }
}
void    fSetBinIMR2 ()      {
    for (int i = 0; i <= nBinIMR2; i++ )
    {
        fArrIMR2[i] = fMinIMR2+(i)*(fMaxIMR2 - fMinIMR2)/(static_cast<Float_t>(nBinIMR2));
    }
}
void    fSetAllBins ()      {
    fSetBinIM1D();
    fSetBinIM2D();
    fSetBinPT1D();
    fSetBinPT2D();
    fSetBinNTup();
    fSetBinMult();
    fSetBinRap_();
    fSetBinSyst();
    fSetBinIMRs();
    fSetBinIMR2();
}
//
//-//   Getters
Int_t   fGetBinIM1D (Float_t input_value )      {
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
Int_t   fGetBinIM2D (Float_t input_value )      {
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
Int_t   fGetBinPT1D (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinPT1D; iBin++ )
    {
        if ( input_value <= fArrPT1D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t   fGetBinPT2D (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinPT2D; iBin++ )
    {
        if ( input_value <= fArrPT2D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t   fGetBinMult (Float_t input_value )      {
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
Int_t   fGetBinRap_ (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinRap_; iBin++ )
    {
        if ( fabs(input_value) <= fArrRap_[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t   fGetBinNTup (Float_t input_value )      {
    for ( Int_t iBin = 0; iBin <= nBinNTup; iBin++ )
    {
        if ( input_value <= fArrNTup[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}
Int_t   fGetBinMultEff (Float_t input_value )   {
    for ( Int_t iBin = 0; iBin <= nMltTrgECls; iBin++ )
    {
        if ( input_value <= kMltTrgECls[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
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

#endif
