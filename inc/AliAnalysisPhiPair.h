#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H

// Analysis Utility
#include "../AliAnalysisUtility/AliAnalysisUtility.h"

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//

enum            fFitResults1D
{
    Background, Signal, kMass, kSlop, kWidt
};

enum            fFitResults2D
{
    BackgBackg, BackgSignl, SignlBackg, SignlSignl
};

// Performance Values
auto const  kCPU_use                =   4;
auto const  kNCycle_                =   10;
auto const  kStatEvalCycles         =   200;
auto const  kBinningPrecision1D     =   0.25*MeV;
auto const  kBinningPrecision2D     =   0.5*MeV;
auto const  kPrintIntervalPP        =   1000000;

// Analysis Values
auto const  bPythiaTest             =   kFALSE;
auto const  kOnlyTrue               =   false;

//-// File Names                    // Re-Make the name with %s instead of the folder name

//-//-// Yield Analysis
    //PreProcessing
auto const  fYldPreProc             =   "./result/yield/InvariantMassHistograms.root";
auto const  fYldPrePrMC             =   "./result/yield/PPReference.root";
auto const  fYldSigExtr             =   "./result/yield/SEHistograms.root";
auto const  fYldSigChek             =   "./result/yield/FitCheckHisto.root";
auto const  fYldSigCorr             =   "./result/yield/SCHistograms.root";
auto const  fYldSigCh2k             =   "./result/yield/FitCheckHist2.root";

//-//-// Yield in Multiplicity Analysis
    //PreProcessing
auto const  fMltPreProc             =   "./result/multiplicity/InvariantMassHistograms.root";
auto const  fMltPrePrMC             =   "./result/multiplicity/PPReference.root";
auto const  fMltSigExtr             =   "./result/multiplicity/SEHistograms.root";
auto const  fMltSigChek             =   "./result/multiplicity/FitCheckHisto.root";
auto const  fMltSigCorr             =   "./result/multiplicity/SCHistograms.root";
auto const  fMltSigCh2k             =   "./result/multiplicity/FitCheckHist2.root";

//-//-// Trigger Analysis
    //PreProcessing
auto const  fTrgPreProc             =   "./result/trigger/DifferentialCandidateEventsHistograms.root";
auto const  fTrgPrePrMC             =   "./result/trigger/PPReference.root";

//-//-// Rapidity Analysis
    //PreProcessing
auto const  fRapPreProc             =   "./result/rapidity/PPRapidityDistributions.root";
auto const  fRapPrePrMC             =   "./result/rapidity/PPReference.root";
auto const  fRapSigExtr             =   "./result/rapidity/SEHistograms.root";
auto const  fRapSigChek             =   "./result/rapidity/FitCheckHisto.root";
auto const  fRapSigCorr             =   "./result/rapidity/SCHistograms.root";
auto const  fRapSigCh2k             =   "./result/rapidity/FitCheckHist2.root";

//-// Tree Names
auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaonCandidate_Tree     =   "KaonCandidate";
auto const  fKaonCandidateEff_Tree  =   "KaonEfficiency";

// Analysis Values
//
//>->   Analysis constants
//
//>->->->   Particles
//
auto const  kPhiMesonMass_          =   1.019455;   //  1.019455    +- 0.000020
auto const  kPhiMesonWidth          =   0.00426;    //  0.00426     +- 0.00004
auto const  kKaonMass               =   .493677;
auto const  kKaonMassUncert         =   .000013;
//
// TO BE GOT RID OF, Legacy compatibility
auto const  kPMas                   =   1.019455;   //  1.019455    +- 0.000020
auto const  kPWid                   =   0.00426;    //  0.00426     +- 0.00004
//
//>->->->   Detectors
//
auto const  kDetectorSlope          =   1.;
//
//>->->->   Systematics
//
//>->->->->->   Branching Ratio
auto const  kBR                     =   0.492;
auto const  kSysLow_BR              =   0.01;
auto const  kSysHig_BR              =   0.01;
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
auto const  kTriggerEff             =   .85; // 5TeV 0.7574;// 7TeV 0.85;
auto const  nMltTrgECls             =   9;
Float_t const  kMultTrgEff      []  =   {.998822,0.995576,0.991524,0.986489,0.975743,0.9575743,0.937151,0.897753,0.696985};
Float_t const  kMltTrgECls      []  =   {0,5,10,15,20,30,40,50,70,100};
auto const  kSysLow_Trigger         =   0.01;
auto const  kSysHig_Trigger         =   0.01;
//
auto const  kTrackingEff            =   -1.0; //Unused
auto const  kSysLow_Tracking        =   0.01;
auto const  kSysHig_Tracking        =   0.01;

//-//   Analysis settings
auto        kDoMultiplicity         =   false;
auto        kDoYield                =   false;
auto        kDoTrigger              =   false;
auto        kDoRapidity             =   false;

//-// InvMass range Pythia MC
const   Float_t   fMinIMMC  =   0.75;
const   Float_t   fMaxIMMC  =   1.25;

//-// InvMass bins 1D
const   Float_t   fMinIM1D  = 0.99;
const   Float_t   fMaxIM1D  = 1.08;
const   Int_t     nBinIM1D  = (int)((fMaxIM1D-fMinIM1D)/kBinningPrecision1D);
        Float_t * fArrIM1D  = new Float_t [nBinIM1D+1];

//-// InvMass bins 2D
const   Float_t   fMinIM2D  = 0.99;
const   Float_t   fMaxIM2D  = 1.08;
const   Int_t     nBinIM2D  = (int)((fMaxIM2D-fMinIM2D)/kBinningPrecision2D);
        Float_t * fArrIM2D  = new Float_t [nBinIM2D+1];

//-// pT bins 1D
const   Int_t     nBinPT1D  =   15;
const   Float_t   fMinPT1D  =   0.4;
const   Float_t   fMaxPT1D  =   12.;
        Float_t  *fArrPT1D  =   new Float_t [nBinPT1D+1];

//-// pT bins 2D
const   Int_t     nBinPT2D  =   7;
const   Float_t   fMinPT2D  =   0.4;
const   Float_t   fMaxPT2D  =   12.;
        Float_t  *fArrPT2D  =   new Float_t [nBinPT2D+1];

//-// Muliplicity bins
const   Int_t     nBinMult  =   5;
const   Float_t   fMinMult  =   0.0;
const   Float_t   fMaxMult  =   100.0;
        Float_t  *fArrMult  =   new Float_t [nBinMult+1];

//-// Rapidity bins
const   Int_t     nBinRap_  =   4;
const   Float_t   fMinRap_  =   0.;
const   Float_t   fMaxRap_  =   .5;
        Float_t  *fArrRap_  =   new Float_t [nBinRap_+1];

//-// N-Tuples bins
const   Int_t     nBinNTup  =   5;
const   Float_t   fMinNTup  =   -0.5;
const   Float_t   fMaxNTup  =   4.5;
        Float_t  *fArrNTup  =   new Float_t [nBinNTup+1];

//-// Systematics Options
//
//-//-//    SE
//
const   Bool_t  f1DOptio    =   true;
const   Bool_t  f2DOptio    =   true;
const   Int_t   nOptions    =   18;
const   string  sOptions[]  =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ","RK","RL","RI","RM","RN","W","CH3","CH5"};
const   Int_t   nOption2    =   0;
const   string  sOption2[]  =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ","RK","RL","RI","RM","RN","W","CH3","CH5","BKG"};
//
//-//-//    PID
//
const   Int_t   nPIDFiles       =   6;
const   TString sPID_DT_Name    =   "LHC10_kAnyINT_PID_%i.root";
const   TString sPID_MC_Name    =   "LHC14j4_kAnyINT_PID_%i.root";
//
//------------------------------//
//    DATA STRUCTURES           //
//------------------------------//

typedef struct
{
    UChar_t     nPhi,           EventMask,      iKaon[1024],    jKaon[1024],   Nature[1024];
    Float_t     Multiplicity,   Px[1024],       Py[1024],       Pz[1024],       pT[1024],      Rap[1024],       InvMass[1024];
    Bool_t      kHasRap[1024],  kHasMult;
    Int_t       iPT1D[1024],    iPT2D[1024],    iRap[1024],     iMult;
} Struct_PhiCandidate;

typedef struct
{
    UChar_t nKaon,          EventMask,      Charge[1024];
    Char_t  SigmaTOF[1024], SigmaTPC[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_KaonCandidate;

typedef struct
{
    UChar_t nPhi,           TrueEventMask,  EventMask,      Selection[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  fTru,           fGen,           fRec;
} Struct_PhiEfficiency;

typedef struct
{
    UChar_t nKaon,          TrueEventMask,  EventMask,      Charge[1024],   Selection[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  ftru;
} Struct_KaonEfficiency;

enum class kEventMask {
    kVoid1, kPileUp, kPileUpMult, kINELgt0, kVoid5, kVoid6, kVoid7, kVoid8
};

enum kEventCount {
    kUnderFlow = 0, kTrigger = 1, kHasEvent = 2, kHasMCTracks = 3, kHasPID = 4, kNoSPDVtx = 5, kVtxMismatch = 6, kVertex = 7, kVertex10 = 8
};

//------------------------------//
//    VARIABLES UTILITIES       //
//------------------------------//

void    fSetBinIM1D ()
{
    for (int i = 0; i <= nBinIM1D; i++ )
    {
        fArrIM1D[i] = fMinIM1D+(i)*(fMaxIM1D - fMinIM1D)/(static_cast<Float_t>(nBinIM1D));
    }
}

void    fSetBinIM2D ()
{
    for (int i = 0; i <= nBinIM2D; i++ )
    {
        fArrIM2D[i] = fMinIM2D+(i)*(fMaxIM2D - fMinIM2D)/(static_cast<Float_t>(nBinIM2D));
    }
}

void    fSetBinNTup ()
{
    for (int i = 0; i <= nBinNTup; i++ )
    {
        fArrNTup[i] = fMinNTup+(i)*(fMaxNTup - fMinNTup)/(static_cast<Float_t>(nBinNTup));
    }
}

void    fSetBinPT1D ()
{
    fArrPT1D[0]     =   0.4;
    fArrPT1D[1]     =   0.7;
    fArrPT1D[2]     =   0.9;
    fArrPT1D[3]     =   1.2;
    fArrPT1D[4]     =   1.4;
    fArrPT1D[5]     =   1.6;
    fArrPT1D[6]     =   1.8;
    fArrPT1D[7]     =   2.0;
    fArrPT1D[8]     =   2.2;
    fArrPT1D[9]     =   2.6;
    fArrPT1D[10]    =   3.0;
    fArrPT1D[11]    =   3.5;
    fArrPT1D[12]    =   4.0;
    fArrPT1D[13]    =   5.0;
    fArrPT1D[14]    =   8.0;
    fArrPT1D[15]    =   12.0;
}

void    fSetBinPT2D ()
{
    
    fArrPT2D[0]     =   0.40;
    fArrPT2D[1]     =   0.70;
    fArrPT2D[2]     =   0.90;
    fArrPT2D[3]     =   1.40;
    fArrPT2D[4]     =   2.00;
    fArrPT2D[5]     =   3.00;
    fArrPT2D[6]     =   5.00;
    fArrPT2D[7]     =   12.00;
}

void    fSetBinMult ()
{
    fArrMult[0]  =  0.00;
    fArrMult[1]  =  5.00;
    fArrMult[2]  =  15.0;
    fArrMult[3]  =  30.0;
    fArrMult[4]  =  50.0;
    fArrMult[5]  =  100.;
}

void    fSetBinRap_ ()
{
    fArrRap_[0]  =  0.0;
    fArrRap_[1]  =  0.08;
    fArrRap_[2]  =  0.16;
    fArrRap_[3]  =  0.32;
    fArrRap_[4]  =  0.5;
}

Int_t   fGetBinIM1D (Float_t input_value )
{
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

Int_t   fGetBinIM2D (Float_t input_value )
{
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

Int_t   fGetBinPT1D (Float_t input_value )
{
    for ( Int_t iBin = 0; iBin <= nBinPT1D; iBin++ )
    {
        if ( input_value <= fArrPT1D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}

Int_t   fGetBinPT2D (Float_t input_value )
{
    for ( Int_t iBin = 0; iBin <= nBinPT2D; iBin++ )
    {
        if ( input_value <= fArrPT2D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}

Int_t   fGetBinMult (Float_t input_value )
{
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

Int_t   fGetBinRap_ (Float_t input_value )
{
    for ( Int_t iBin = 0; iBin <= nBinRap_; iBin++ )
    {
        if ( fabs(input_value) <= fArrRap_[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}

Int_t   fGetBinNTup (Float_t input_value )
{
    for ( Int_t iBin = 0; iBin <= nBinNTup; iBin++ )
    {
        if ( input_value <= fArrNTup[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}

Int_t   fGetBinMultEff (Float_t input_value )
{
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

bool    fCutRapidity        ( Double_t  dRapidity )
{
    if ( fabs(dRapidity) < 0.5 ) return true;
    return false;
}

bool    fCutInvariantMass    ( Double_t  dInvariantMass )
{
    if ( dInvariantMass < fMinIM1D ) return false;
    if ( dInvariantMass > fMaxIM1D ) return false;
    return true;
}

bool    fCutTransverseMom   ( Double_t  dTransverseMom )
{
    if ( dTransverseMom < fMinPT1D ) return false;
    if ( dTransverseMom > fMaxPT1D ) return false;
    return true;
}

bool    fCutMultiplicity    ( Double_t  dMultiplicity )
{
    return true;
    if ( dMultiplicity < fMinMult ) return false;
    if ( dMultiplicity > fMaxMult ) return false;
    return true;
}

bool    fAcceptCandidate ( Double_t  dInvariantMass, Double_t dTransverseMom )
{
    if ( !fCutInvariantMass(dInvariantMass) ) return false;
    if ( !fCutTransverseMom(dTransverseMom) ) return false;
    return true;
}

bool    fCheckCoupleKaons( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi )
{
    if ( SPhiCandidates.iKaon[aAccCandidates[iPhi]] == SPhiCandidates.iKaon[aAccCandidates[jPhi]] ) return false;
    if ( SPhiCandidates.jKaon[aAccCandidates[iPhi]] == SPhiCandidates.jKaon[aAccCandidates[jPhi]] ) return false;
    if ( SPhiCandidates.iKaon[aAccCandidates[iPhi]] == SPhiCandidates.jKaon[aAccCandidates[jPhi]] ) return false;
    if ( SPhiCandidates.jKaon[aAccCandidates[iPhi]] == SPhiCandidates.iKaon[aAccCandidates[jPhi]] ) return false;
    return true;
}

bool    fAcceptCandidate ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi )
{
    return true;
}

bool    fAcceptCandidate ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi )
{
    // Non equal candidates
    if ( iPhi == jPhi ) return false;
                
    // Only non overlapping couples of Kaons
    if ( !fCheckCoupleKaons(SPhiCandidates,aAccCandidates,iPhi,jPhi) ) return false;

    return true;
}

bool    fAcceptCandidate ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi, Int_t kPhi )
{
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

bool    fAcceptCandidate ( Struct_PhiCandidate SPhiCandidates, Int_t aAccCandidates[], Int_t iPhi, Int_t jPhi, Int_t kPhi, Int_t lPhi )
{
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

//------------------------------//
//    HISTOGRAM UTILITIES       //
//------------------------------//
//
template < class Tclass >
void    SetAxis             ( Tclass * aTarget, string aOption = "" )   {
    if ( aOption.find("IM") != -1 )
    {
        if ( aOption.find("2D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("m_{K^{+}K^{-}} candidate #phi_{1} (GeV/c^{2})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("m_{K^{+}K^{-}} candidate #phi_{2} (GeV/c^{2})");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("1D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
        }
    }
    else if ( aOption.find("PT") != -1 )
    {
        if ( aOption.find("DD") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi#phi}}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("2D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
                
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
                
            // Z-Axis formatting
            aTarget->GetZaxis()->SetTitle("#frac{d^{3}N_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
            aTarget->GetZaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("1D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
    }
}
//
//--
//
template < class Tclass >
Double_t        HistIntegrals               ( Tclass* aTarget, string aOption = "", Double_t aLowBound = 0, Double_t aHigBound = 0 )    {
    
    Double_t    fResult = 0;
    Int_t       fnPowr;
    Bool_t      fWidth, fPower, fError;
    
    if ( aOption.find("W") != -1 )      fWidth  =   true;
    else                                fWidth  =   false;
    
    if ( aOption.find("x^") != -1 )   { fPower  =   true;   fnPowr  =   aOption.at(aOption.find("x^")+2)-'0'; }
    else                                fPower  =   false;
    
    if ( aOption.find("E") != -1 )      fError  =   true;
    else                                fError  =   false;
    
    // Starting iterations at 0
    for ( Int_t fiTer = 1; fiTer <= aTarget->GetNbinsX(); fiTer++ )
    {
        // Determining the pT point to calculate
        auto fPoint =   aTarget->GetBinCenter(fiTer);
        
        // Exiting the loop if out of bound
        if ( aTarget->GetBinLowEdge(fiTer)   <= aLowBound && aLowBound != aHigBound ) continue;
        if ( aTarget->GetBinLowEdge(fiTer+1) >= aHigBound && aLowBound != aHigBound ) break;
        
        Double_t        fCycleAdd   =   ( aTarget->GetBinContent(fiTer) );
        if ( fError )   fCycleAdd   =   ( aTarget->GetBinError(fiTer) );
        if ( fWidth )   fCycleAdd  *=   ( aTarget->GetBinWidth(fiTer) );
        if ( fPower )   fCycleAdd  *=   pow( ( aTarget->GetBinCenter(fiTer) ) , fnPowr );
        if ( fError )   fCycleAdd  *=   fCycleAdd;
        fResult += fCycleAdd;
    }
    
    if ( fError )   return sqrt(fResult);
    return fResult;
}
//
//------------------------------//
//    ANALYSISI SEPCIFIC Fncs   //
//------------------------------//
//
//_____________________________________________________________________________
//
Int_t kColor[] = {38,kBlue,kBlue+3,46,38};
Int_t kStyle[] = {26,9,10,25,22};
Int_t kWidth[] = {1,3,3,1,1};
//
//_____________________________________________________________________________
//
void                fChooseOption                   ( TString fOption ) {
    if ( fOption.IsNull() ) cout << "[INFO] No option chosen, standard all inclusive analysis enabled" <<endl;
    else    {
        kDoMultiplicity         =   false;
        kDoYield                =   false;
        kDoTrigger              =   false;
        kDoRapidity             =   false;
        if ( fOption.Contains("Multiplicity",TString::kIgnoreCase) )    { kDoMultiplicity = true;   cout << "[INFO] Mutliplicity option chosen" <<endl;}
        if ( fOption.Contains("Yield",TString::kIgnoreCase) )           { kDoYield = true;          cout << "[INFO] Yield option chosen" <<endl;}
        if ( fOption.Contains("Trigger",TString::kIgnoreCase) )         { kDoTrigger = true;        cout << "[INFO] Trigger option chosen" <<endl;}
        if ( fOption.Contains("Rapidity",TString::kIgnoreCase) )        { kDoRapidity = true;       cout << "[INFO] Rapidity option chosen" <<endl;}
    }
}
//
//_____________________________________________________________________________
//
Double_t            fGammaPhiValue                  ( Double_t fYieldPhi, Double_t fYieldPhiPhi )  {
    return  2*fYieldPhiPhi/fYieldPhi -fYieldPhi;
}
//
//_____________________________________________________________________________
//
Double_t            fGammaPhiError                  ( Double_t fYieldPhi, Double_t fYieldPhiPhi, Double_t fErrorPhi, Double_t fErrorPhiPhi)  {
    auto    fPar1   =   2*fErrorPhiPhi/fYieldPhi;
    auto    fPar2   =   (2*fYieldPhiPhi/(fYieldPhi*fYieldPhi)+1)*fErrorPhi;
    return  fPar1 + fPar2;
}
//
//_____________________________________________________________________________
//
void                fSetLevyTsalis                  ( Bool_t fIsConditional = false, Double_t fIntegral = 0.032 ) {
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas*0.9,kPMas*1.1);
    fLevyFit1D  ->  SetParameter(0,kPMas);  //Particle Mass?
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,5.,9.);
    fLevyFit1D  ->  SetParameter(1,7.18);    // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.2,.400);
    fLevyFit1D  ->  SetParameter(2,.300);   // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,1.e-9,1.);
    fLevyFit1D  ->  SetParameter(3,fIntegral);  //
    
    if ( false )   {
        // - // Setting up Fit parameters// - // Setting up Fit parameters
        
        // Mass
        fLevyFit1D  ->  SetParLimits(0,0.,100.);
        fLevyFit1D  ->  SetParameter(0,kPMas);
        
        // n-Parameter
        fLevyFit1D  ->  SetParLimits(1,2.0001,1000.);
        fLevyFit1D  ->  SetParameter(1,5.); // 6.7
        
        // T-Parameter
        fLevyFit1D  ->  SetParLimits(2,0.,100.);
        fLevyFit1D  ->  SetParameter(2,.3); // .272
        
        // dN/dy
        fLevyFit1D  ->  SetParLimits(3,0.5e-6,1.e-3);
        fLevyFit1D  ->  SetParameter(3,1.e-6);
    }
}
//
//_____________________________________________________________________________
//
void                fFitLevyTsalis                  ( TGraphAsymmErrors* gToBeFitted ) {
    return;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*  fSetSystErrors                  ( TGraphAsymmErrors* gStatistics ) {
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gStatistics);
    for ( Int_t iPnt = 0; iPnt < fResult->GetN(); iPnt++ ) {
        auto    fYValue =   fResult ->  GetPointY(iPnt);
        fResult ->  SetPointEYhigh  ( iPnt, fYValue*sqrt(kSysHig_BR*kSysHig_BR+kSysHig_TR*kSysHig_TR+kSysHig_PD*kSysHig_PD+kSysHig_1D_SE*kSysHig_1D_SE) );
        fResult ->  SetPointEYlow   ( iPnt, fYValue*sqrt(kSysLow_BR*kSysLow_BR+kSysLow_TR*kSysLow_TR+kSysLow_PD*kSysLow_PD+kSysLow_1D_SE*kSysLow_1D_SE) );
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
std::vector<TGraphAsymmErrors*>   fSetSystErrors    ( std::vector<TGraphAsymmErrors*>  gStatistics ) {
    std::vector<TGraphAsymmErrors*> fResult;
    for ( auto& iGraph : gStatistics )  {
        TGraphAsymmErrors  *fCurrentGraph =   new TGraphAsymmErrors(*iGraph);
        for ( Int_t iPnt = 0; iPnt < fCurrentGraph->GetN(); iPnt++ ) {
            auto    fYValue =   fCurrentGraph ->  GetPointY(iPnt);
            fCurrentGraph ->  SetPointEYhigh  ( iPnt, fYValue*sqrt(4*(kSysHig_BR*kSysHig_BR+kSysHig_TR*kSysHig_TR+kSysHig_PD*kSysHig_PD)+kSysHig_2D_SE*kSysHig_2D_SE) );
            fCurrentGraph ->  SetPointEYlow   ( iPnt, fYValue*sqrt(4*(kSysLow_BR*kSysLow_BR+kSysLow_TR*kSysLow_TR+kSysLow_PD*kSysLow_PD)+kSysLow_2D_SE*kSysLow_2D_SE) );
        }
        fResult.push_back(fCurrentGraph);
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t*           fExtrapolateModel               ( bool fIsConditional, TGraphAsymmErrors* gStatistics, TGraphAsymmErrors* gSystematics, Double_t fIntegral = 0.032, TString fName = "ExtrapolateSignal" )    {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    //  Starting the
    fSetLevyTsalis(fIsConditional,fIntegral);
    //
    //  Generating a Full Error Spectra to fit and extrapolating at low pT
    TGraphAsymmErrors      *gTotal      =   new TGraphAsymmErrors(*(fSumGraphErrors(gStatistics,gSystematics)));
    gTotal                              ->  Fit(fLevyFit1D,"IMREQ0SEX0","",fMinPT1D,fMaxPT1D);
    fResult[0]                          =   fLevyFit1D  ->Integral(0.,fMinPT1D);
    //
    //  !TODO: Scorporare la fuznione calcola mean pT
    //  Evaluating the Mean pT
    if ( !fIsConditional )  {
        fResult[5]                          =   (fMinPT1D)*fLevyFit1D->Moment(1,0.,fMinPT1D)*(fResult[0]);
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            fResult[5]                     +=   (fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D])*(fLevyFit1D->Moment(1,fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]))*(gTotal->GetPointY(iPT1D));
        }
        cout << fResult[5] << " -- " << endl;
    }   else    {
        fResult[5]                          =   (fMinPT2D)*fLevyFit1D->Moment(1,0.,fMinPT2D)*(fResult[0]);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            fResult[5]                     +=   (fArrPT2D[iPT2D+1] - fArrPT2D[iPT2D])*(fLevyFit1D->Moment(1,fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]))*(gTotal->GetPointY(iPT2D));
        }
    }
    //
    //  Save Fitresult for check later
    TCanvas                *cDrawFit    =   new TCanvas(Form("%s_%s",gStatistics->GetName(),fName.Data()),"");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gTotal->Draw();
    fLevyFit1D->Draw("same");
    cDrawFit->Write();
    cDrawFit->SaveAs(Form("result/tmp/%s_%s.pdf",gStatistics->GetName(),fName.Data()));
    delete cDrawFit;
    //
    //  Measuring Statistical error
    TCanvas *   cDrawFitStat    =   new TCanvas(Form("STAT_%s_%s",gStatistics->GetName(),fName.Data()),"");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gTotal->Draw();
    TH1D*   hStatIntegral   =   new TH1D(Form("hStatIntegral_%s_%s",gStatistics->GetName(),fName.Data()),"hStatIntegral",100000,0.0,.1);
    for ( Int_t iFit = 0; iFit < kStatEvalCycles; iFit++ )  {
        //  Set Standard Fit
        fSetLevyTsalis(fIsConditional,fIntegral);
        //
        //  Generating the Fit TGraph
        auto fSubject   =   fRandomizePoints(gSystematics,gStatistics);
        //
        fSubject    ->  Fit(fLevyFit1D,"IMREQ0SEX0","",fMinPT1D,fMaxPT1D);
        //
        hStatIntegral->Fill(fLevyFit1D  ->Integral(0.,fMinPT1D));
        auto fMemory = new TF1(*fLevyFit1D);
        fMemory->Draw("same");
    }
    cDrawFitStat->Write();
    cDrawFitStat->SaveAs(Form("result/tmp/STAT_%s_%s.pdf",gStatistics->GetName(),fName.Data()));
    delete cDrawFitStat;
    fResult[1]  =   hStatIntegral->GetRMS();
    fResult[2]  =   hStatIntegral->GetRMS();
    //
    //  Measuring Systematics error
    TCanvas *   cDrawFitSyst    =   new TCanvas(Form("SYST_%s_%s",gStatistics->GetName(),fName.Data()),"");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gTotal->Draw();
    TH1D*   hSystIntegral   =   new TH1D(Form("hSystIntegral_%s_%s",gStatistics->GetName(),fName.Data()),"hSystIntegral",100000,0.0,.1);
    for ( Int_t iFit = 0; iFit < kStatEvalCycles; iFit++ )  {
        //  Set Standard Fit
        fSetLevyTsalis(fIsConditional,fIntegral);
        //
        //  Generating the Fit TGraph
        auto fSubject   =   fRandomizePoints(gStatistics,gSystematics);
        //
        fSubject    ->  Fit(fLevyFit1D,"IMREQ0SEX0","",fMinPT1D,fMaxPT1D);
        //
        hSystIntegral->Fill(fLevyFit1D  ->Integral(0.,fMinPT1D));
        auto fMemory = new TF1(*fLevyFit1D);
        fMemory->Draw("same");
    }
    fResult[3]  =   hSystIntegral->GetRMS();
    fResult[4]  =   hSystIntegral->GetRMS();
    cDrawFitSyst->Write();
    cDrawFitSyst->SaveAs(Form("result/tmp/SYST_%s_%s.pdf",gStatistics->GetName(),fName.Data()));
    delete cDrawFitSyst;
    //
    fResult[6]  =   0;
    fResult[7]  =   0;
    fResult[8]  =   0;
    fResult[9]  =   0;
    hStatIntegral->Write();
    hSystIntegral->Write();
    //
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t*           fIntegrateModel                 ( TGraphAsymmErrors* gStatistics, TGraphAsymmErrors* gSystematics, TString fName = "IntegrateSignal" )      {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    Int_t   fNPoints =   gStatistics ->  GetN();
    if  ( fNPoints  != gSystematics ->  GetN() )
    {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult = new   Double_t    [5];
    //
    //  Calculate Integral and mean pT
    for ( Int_t iFill = 0; iFill < 5; iFill++ ) fResult[iFill]  =   0;
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fXBinCentr      =   ( gStatistics ->  GetPointX(iFit) );
        auto    fYValue         =   ( gStatistics ->  GetPointY(iFit) );
        auto    fXBinWidth      =   ( gStatistics ->  GetErrorXhigh(iFit)    +   gStatistics ->  GetErrorXlow(iFit) );
        auto    fYErrStatLow    =   ( gStatistics ->  GetErrorYlow(iFit) );
        auto    fYErrStatHigh   =   ( gStatistics ->  GetErrorYhigh(iFit) );
        auto    fYErrSystLow    =   ( gSystematics ->  GetErrorYlow(iFit) );
        auto    fYErrSystHigh   =   ( gSystematics ->  GetErrorYhigh(iFit) );
        //
        fResult[0]             +=   fXBinWidth*fYValue;
        fResult[1]             +=   fXBinWidth*fYErrStatLow*fXBinWidth*fYErrStatLow;
        fResult[2]             +=   fXBinWidth*fYErrStatHigh*fXBinWidth*fYErrStatHigh;
        fResult[3]             +=   fXBinWidth*fYErrSystLow;
        fResult[4]             +=   fXBinWidth*fYErrSystHigh;
    }
    for ( Int_t iTer = 1; iTer <= 2; iTer++ )    {
        auto fTemp      =   fResult[iTer];
        fResult[iTer]   =   sqrt(fTemp);
    }
    //
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t*           fMeasureFullYield               ( TGraphAsymmErrors* gStatistics, TGraphAsymmErrors* gSystematics, TString fName = "MeasureFullYield" )     {
    // Optimisation mode
    gROOT->SetBatch(true);
    //
    // Result format:  Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult             =   new Double_t        [10];
    //
    bool fIsConditional = false;
    if ( fName.First("2D") != -1  ) fIsConditional = true;
    //
    auto        fIntegralResults    =   fIntegrateModel     (gStatistics,gSystematics,fName);
    auto        fExtrapolResults    =   fExtrapolateModel   (fIsConditional,gStatistics,gSystematics,fIntegralResults[0],fName);
    //
    //  Mean Value of Result
    fResult[0]  =   fIntegralResults[0] +   fExtrapolResults[0];
    //
    //  Statistical Error of Result
    fResult[1]  =   sqrt(fIntegralResults[1]*fIntegralResults[1] +   fExtrapolResults[1]*fExtrapolResults[1]);
    fResult[2]  =   sqrt(fIntegralResults[2]*fIntegralResults[2] +   fExtrapolResults[2]*fExtrapolResults[2]);
    //
    //  Systematical Error of Result
    fResult[3]  =   fIntegralResults[3] +   fExtrapolResults[4];
    fResult[4]  =   fIntegralResults[3] +   fExtrapolResults[4];
    //
    //  Mean Value of pT
    fResult[5]  =   fExtrapolResults[5]/fResult[0];
    //
    //  Statistical Error of Result
    fResult[6]  =   fExtrapolResults[6];
    fResult[7]  =   fExtrapolResults[7];
    //
    //  Systematical Error of Result
    fResult[8]  =   fExtrapolResults[8];
    fResult[9]  =   fExtrapolResults[9];
    //
    // End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t            fEvaluateINELgt0                ( Int_t iMultBin, TH1  *hMultCounter)  {
    Double_t    fResult =   0;
    Int_t       fUppLimit   =   fArrMult[nBinMult];
    Int_t       fLowLimit   =   fArrMult[0];
    if  ( iMultBin > -1  && iMultBin < nBinMult)  {
        fUppLimit   =   fArrMult[iMultBin+1];
        fLowLimit   =   fArrMult[iMultBin];
    }   else    {
        cout << "[WARNING] Invalid index, returning full multiplicity count" << endl;
    }
    for ( Int_t iBin = 1; iBin <= hMultCounter->GetNbinsX(); iBin++ )   {
        auto fMultValue = hMultCounter->GetBinCenter(iBin);
        if ( fMultValue < fLowLimit )   continue;
        if ( fMultValue > fUppLimit )   break;
        fResult  +=  hMultCounter->GetBinContent(iBin)/kMultTrgEff[fGetBinMultEff(iBin-1.5)];
    }
    return      fResult;
}
//
//_____________________________________________________________________________
//
//------------------------------//
//    ANALYSISI LEGACY Fncs     //
//------------------------------//

TH1F*           fCheckPublishedData  ( TH1F* fMyResults, TH1F* hPublishedResults, TGraphAsymmErrors* gPublishedResults )    {
    TH1F   *fCheck  =   new TH1F(*hPublishedResults);
    fCheck->Divide(fMyResults,hPublishedResults);
    for ( int i = 0; i < fCheck->GetNbinsX(); i++ ) {
        auto pubError   =   gPublishedResults->GetErrorYhigh(i);
        auto pubValY   =   gPublishedResults->GetPointY(i);
        auto myError    =   fMyResults->GetBinError(i+1);
        auto myValY    =   fMyResults->GetBinContent(i+1);
        auto checkVal   =   fCheck->GetBinContent(i+1);
        fCheck->SetBinError (i+1,checkVal*sqrt((pubError*pubError/(pubValY*pubValY))+(myError*myError/(myValY*myValY))));
    }
    return fCheck;
}

TH1F*           fCheckPublishedData  ( TH1D* fMyResults, TH1F* hPublishedResults, TGraphAsymmErrors* gPublishedResults )    {
    TH1F   *fCheck  =   new TH1F(*hPublishedResults);
    fCheck->Divide(fMyResults,hPublishedResults);
    for ( int i = 0; i < fCheck->GetNbinsX(); i++ ) {
        auto pubError   =   gPublishedResults->GetErrorYhigh(i);
        auto pubValY   =   gPublishedResults->GetPointY(i);
        auto myError    =   fMyResults->GetBinError(i+1);
        auto myValY    =   fMyResults->GetBinContent(i+1);
        auto checkVal   =   fCheck->GetBinContent(i+1);
        fCheck->SetBinError (i+1,checkVal*sqrt((pubError*pubError/(pubValY*pubValY))+(myError*myError/(myValY*myValY))));
    }
    return fCheck;
}

TH1D           *fCheckPublishedData     ( TH1D* fMyResults, TH1* fPublishedResults )   {
    auto    fCheck  =   new TH1D(*fMyResults);
    fCheck->Divide(fMyResults,fPublishedResults);
    return  fCheck;
}

TH1D           *fCheckPublishedData     ( TH1D* fMyResults, TGraphAsymmErrors* fPublishedResults )   {
    auto    fCheck  =   new TH1D(*fMyResults);
    fCheck->Divide(fMyResults,(fPublishedResults->GetHistogram()));
    return  fCheck;
}

TH1F           *fCheckPublishedData     ( TGraphAsymmErrors* fMyResults, TH1* fPublishedResults )   {
    auto    fCheck  =   new TH1F(*(fMyResults->GetHistogram()));
    fCheck->Divide((fMyResults->GetHistogram()),fPublishedResults);
    return  fCheck;
}

TH1F           *fCheckPublishedData     ( TGraphAsymmErrors* fMyResults, TGraphAsymmErrors* fPublishedResults )   {
    auto    fCheck  =   new TH1F(*(fMyResults->GetHistogram()));
    fCheck->Divide((fMyResults->GetHistogram()),(fPublishedResults->GetHistogram()));
    return  fCheck;
}

int             fLegendSelect                   ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("RapTru") )       return 1;
    if ( !fOption.compare("Rap") )          return 3;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 2;
    else return -1;
}

void            fLegendMaker                    ( RooPlot * fRooPlot, const char * fSelect, TLegend * fLegend )
{
    switch (fLegendSelect(fSelect))
    {
        case 1:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        case 2:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBS"),   "Fit (Bkg #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSB"),   "Fit (Sig #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        case 3:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg Lin)",        "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooB2"),   "Fit (Bkg Gaus)",       "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

int             fAxisSelect                     ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 3;
    if ( !fOption.compare("Rap") )          return 4;
    if ( !fOption.compare("RapTru") )       return 4;
    else return -1;
}
 
void            fAxisMaker                      ( RooPlot * fRooPlot, const char * fSelect )
{
    switch (fAxisSelect(fSelect))
    {
        case 1:
            fRooPlot                    ->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
            break;
        case 2:
            fRooPlot                    ->GetXaxis()->SetTitle("m^{x}_{K^{+}K^{-}} (GeV/c^{2})");
            break;
        case 3:
            fRooPlot                    ->GetXaxis()->SetTitle("m^{y}_{K^{+}K^{-}} (GeV/c^{2})");
            break;
        case 4:
            fRooPlot                    ->GetXaxis()->SetTitle("|#Delta y| #phi_{1,2}");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

int             fPlotterSelect                  ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("RapTru") )       return 1;
    if ( !fOption.compare("Rap") )          return 3;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 2;
    else return -1;
}

void            fRooPlotPlotter                 ( RooPlot * fRooPlot, const char * fSelect, RooAddPdf fModel , RooDataHist * fData )
{
    switch (fPlotterSelect(fSelect))
    {
        case 1:
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(2),       Name("RooSS"));
            break;
        case 2:
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(38),    MarkerStyle(26) ,   Name("RooData"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  LineColor(4),       LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkg"), LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigSig"),LineColor(2),       Name("RooSS"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigBkg"),LineStyle(kDashed), LineColor(33),    Name("RooSB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkgSig"),LineStyle(kDashed), LineColor(36),    Name("RooBS"));
            break;
        case 3:
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fBk2"),             LineStyle(kDashed), LineColor(33),      Name("RooB2"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(2),       Name("RooSS"));
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

void            fRooPlotMaker                   ( RooPlot * fRooPlot, TLegend * fLegend, RooAddPdf fModel , RooDataHist * fData, const char * fSelect )
{
    fRooPlotPlotter(fRooPlot,fSelect,fModel,fData);
    fLegendMaker(fRooPlot,fSelect,fLegend);
    fAxisMaker(fRooPlot,fSelect);
}

void            SetBoundaries   ( string fOption, Double_t &aValMin, Double_t &aValMax )
{
    aValMin = 0.993;
    aValMax = 1.050;
    if ( fOption.find("RA") != -1 )
    {
        aValMin =   0.993;
        aValMax =   1.045;
    }
    if ( fOption.find("RB") != -1 )
    {
        aValMin =   0.993;
        aValMax =   1.055;
    }
    if ( fOption.find("RC") != -1 )
    {
        aValMin =   0.993;
        aValMax =   1.065;
    }
    if ( fOption.find("RD") != -1 )
    {
        aValMin =   0.993;
        aValMax =   1.075;
    }
    if ( fOption.find("RE") != -1 )
    {
        aValMin =   0.996;
        aValMax =   1.045;
    }
    if ( fOption.find("RF") != -1 )
    {
        aValMin =   0.996;
        aValMax =   1.050;
    }
    if ( fOption.find("RG") != -1 )
    {
        aValMin =   0.996;
        aValMax =   1.055;
    }
    if ( fOption.find("RH") != -1 )
    {
        aValMin =   0.996;
        aValMax =   1.065;
    }
    if ( fOption.find("RI") != -1 )
    {
        aValMin =   0.996;
        aValMax =   1.075;
    }
    if ( fOption.find("RJ") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.045;
    }
    if ( fOption.find("RK") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.050;
    }
    if ( fOption.find("RL") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.055;
    }
    if ( fOption.find("RM") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.065;
    }
    if ( fOption.find("RN") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.075;
    }
}

void            fCoreFitModelSetBoundaries      ( string fOption, Double_t &aValMin, Double_t &aValMax )
{
    SetBoundaries   ( fOption, aValMin, aValMax );
}

Double_t   *    fCoreFitModelOptionSelect       ( string fOption )
{
    Double_t   *fResult     =   new Double_t [7];
    fResult[0] = 0.;
    fResult[1] = 0.;
    fResult[2] = 1.;
    fResult[3] = 0.;
    fResult[4] = 0.;
    fResult[5] = 0.;
    fResult[6] = 1.;
    fCoreFitModelSetBoundaries ( fOption, fResult[0], fResult[1] );
    if  ( fOption.find("CH3") != -1 ) { fResult[2] = 0.; };
    if  ( fOption.find("CH5") != -1 ) { fResult[3] = 1.; };
    if  ( fOption.find("Wdt") != -1 ) { fResult[4] = kPhiMesonMass_*0.1; };
    if  ( fOption.find("Mss") != -1 ) { fResult[5] = kPhiMesonWidth*0.1; };
    if  ( bPythiaTest               ) { fResult[6] = 0.; };
    return fResult;
}

RooFitResult*   fCoreFitModel                   ( RooDataHist *fDataHist, RooRealVar fVariable, RooAddPdf &fModel_, Double_t kCh4Limits, Double_t kCh5Limits, Double_t kVmsLimits, Double_t kVwdLimits, Double_t kVslLimits )
{
    //------ Define what your model is made of -----//
    
    Int_t       nEntries        =   fDataHist->sumEntries();
    
    // Background PDF Coefficients
    RooRealVar  ChebychevPar0   =   RooRealVar      ("ch0","ch0",   0., -1, 1);
    RooRealVar  ChebychevPar1   =   RooRealVar      ("ch1","ch1",   0., -1, 1);
    RooRealVar  ChebychevPar2   =   RooRealVar      ("ch2","ch2",   0., -1, 1);
    RooRealVar  ChebychevPar3   =   RooRealVar      ("ch3","ch3",   0., -1, 1);
    RooRealVar  ChebychevPar4   =   RooRealVar      ("ch4","ch4",   0., 0. - kCh4Limits, 0. + kCh4Limits);
    RooRealVar  ChebychevPar5   =   RooRealVar      ("ch5","ch5",   0., 0. - kCh5Limits, 0. + kCh5Limits);

    //Signal
    RooRealVar  VoigtianMass_   =   RooRealVar      ("Vms","Vms",   kPhiMesonMass_, kPhiMesonMass_ - kVmsLimits, kPhiMesonMass_ + kVmsLimits);
    RooRealVar  VoigtianWidth   =   RooRealVar      ("Vwd","Vwd",   kPhiMesonWidth, kPhiMesonWidth - kVwdLimits, kPhiMesonWidth + kVwdLimits);
    RooRealVar  VoigtianSlope   =   RooRealVar      ("Vsl","Vsl",   kDetectorSlope, kDetectorSlope - kVslLimits, kDetectorSlope + kVslLimits);
    
    // Normalisation coefficients
    RooRealVar  SignalMagnit    =   RooRealVar      ("1SS","1SS",   0.5*nEntries, 0., nEntries);
    RooRealVar  BackgroundMg    =   RooRealVar      ("1BB","1BB",   0.5*nEntries, 0., nEntries);
    
    // PDFs
    RooVoigtian     fSignal     =   RooVoigtian     ("fSig","fSig", fVariable,  VoigtianMass_,  VoigtianWidth,  VoigtianSlope);
    RooChebychev    fBkgrnd     =   RooChebychev    ("fBkg","fBkg", fVariable,  RooArgSet( ChebychevPar0, ChebychevPar1, ChebychevPar2, ChebychevPar3, ChebychevPar4, ChebychevPar5 ));
                    fModel_     =   RooAddPdf       ("fMod","fMod", RooArgList( fSignal, fBkgrnd ),RooArgList( SignalMagnit, BackgroundMg ));
    
    return fModel_.fitTo(*fDataHist,Extended(kTRUE),SumW2Error(kTRUE),Save());
}

RooFitResult*   FitModel                        ( TH1F * _h_Data, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
   
    // Defining the Fit Options
    Double_t       *fFitOptions =   fCoreFitModelOptionSelect( fOption );
    
    //------ General Information on the Data Histogram -----//
    
    // Defining the RooFit Variable
    RooRealVar      fVariable   =       RooRealVar  ("fVariable","fVariable",fFitOptions[0],fFitOptions[1]);
    RooDataHist    *fDataHist   =   new RooDataHist ("","",fVariable,Import(*_h_Data));
    RooAddPdf       fModel_;
    
    // Fitting histogram and
    RooFitResult   *fFitResults     =   fCoreFitModel( fDataHist, fVariable, fModel_, fFitOptions[2], fFitOptions[3], fFitOptions[4], fFitOptions[5], fFitOptions[6] );
    
    // Saving to canvas on file if requested
    if ( fOption.find("S") != -1 )
    {
        Float_t fPTMax      =   fMaxPT1D;
        Float_t fPTMin      =   fMaxPT2D;
        Int_t   fPTindex    =   0;
        if ( fOption.find("PT=")    != -1 ) { fPTindex  =   10*(fOption.at(fOption.find("PT=")+3)-'0')+(fOption.at(fOption.find("PT=")+4)-'0'); }
        if ( fOption.find("1D")     != -1 ) { fPTMin    =   fArrPT1D[fPTindex];  fPTMax    =   fArrPT1D[fPTindex+1]; };
        if ( fOption.find("12D")    != -1 ) { fPTMin    =   fArrPT2D[fPTindex];  fPTMax    =   fArrPT2D[fPTindex+1]; };
        
        hName           =   "DT";
        hTitle          =   Form( "Invariant Mass of Kaons in pT %.2f-%.2f GeV/c", fPTMin, fPTMax );
        if ( bPythiaTest )  { Form( "Invariant Mass of Kaons in pT %.2f-%.2f GeV/c (MC)", fPTMin, fPTMax );   hName = "MC"; };
        
        // Canvas to plot
        TCanvas * fSaveToCanvas     =   new TCanvas(Form( "PT_%.1f_%.1f_1D_%s", fPTMin, fPTMax, hName ),
                                                    Form( "PT_%.1f_%.1f_1D_%s", fPTMin, fPTMax, hName ) );
        
        RooPlot * fSaveToFrame      =   fVariable.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           =   new TLegend   (0.12,0.60,0.30,0.85);
        
        fRooPlotMaker(fSaveToFrame,fLegend,fModel_,fDataHist,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return fFitResults;
}

RooFitResult*   FitModel                        ( TH1D * _h_Data, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
   
    // Defining the Fit Options
    Double_t       *fFitOptions =   fCoreFitModelOptionSelect( fOption );
    
    //------ General Information on the Data Histogram -----//
    
    // Defining the RooFit Variable
    RooRealVar      fVariable   =       RooRealVar  ("fVariable","fVariable",fFitOptions[0],fFitOptions[1]);
    RooDataHist    *fDataHist   =   new RooDataHist ("","",fVariable,Import(*_h_Data));
    RooAddPdf       fModel_;
    
    // Fitting histogram and
    RooFitResult   *fResult     =   fCoreFitModel( fDataHist, fVariable, fModel_, fFitOptions[2], fFitOptions[3], fFitOptions[4], fFitOptions[5], fFitOptions[6] );
    
    // Saving to canvas on file if requested
    if ( fOption.find("S") != -1 )
    {
        Float_t fPTMax      =   fMaxPT1D;
        Float_t fPTMin      =   fMaxPT2D;
        Int_t   fPTindex    =   0;
        if ( fOption.find("PT=")    != -1 ) { fPTindex  =   10*(fOption.at(fOption.find("PT=")+3)-'0')+(fOption.at(fOption.find("PT=")+4)-'0'); }
        if ( fOption.find("1D")     != -1 ) { fPTMin    =   fArrPT1D[fPTindex];  fPTMax    =   fArrPT1D[fPTindex+1]; };
        if ( fOption.find("12D")    != -1 ) { fPTMin    =   fArrPT2D[fPTindex];  fPTMax    =   fArrPT2D[fPTindex+1]; };
        
        hName           =   "DT";
        hTitle          =   Form( "Invariant Mass of Kaons in pT %.2f-%.2f GeV/c", fPTMin, fPTMax );
        if ( bPythiaTest )  { Form( "Invariant Mass of Kaons in pT %.2f-%.2f GeV/c (MC)", fPTMin, fPTMax );   hName = "MC"; };
        
        // Canvas to plot
        TCanvas * fSaveToCanvas     =   new TCanvas(Form( "PT_%.1f_%.1f_1D_%s", fPTMin, fPTMax, hName ),
                                                    Form( "PT_%.1f_%.1f_1D_%s", fPTMin, fPTMax, hName ) );
        
        RooPlot * fSaveToFrame      =   fVariable.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           =   new TLegend   (0.12,0.60,0.30,0.85);
        
        fRooPlotMaker(fSaveToFrame,fLegend,fModel_,fDataHist,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return fResult;
}

RooFitResult*   FitModel        ( TH1D * THdata, const char* fName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTDimension = 1, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    // Check there is a reasonable amount of entries
    if ( !fIsWorthFitting( THdata ) ) return nullptr;
    
    Bool_t  fCheb3  =   false;
    if  ( fOption.find("CH3") != -1 )   fCheb3  =   true;
    
    Bool_t  fCheb5  =   false;
    if  ( fOption.find("CH5") != -1 )   fCheb5  =   true;
    
    Bool_t  fWidth  =   false;
    if  ( fOption.find("W") != -1 )     fWidth  =   true;
    
    Bool_t  fMass_  =   false;
    if  ( fOption.find("M") != -1 )     fMass_  =   true;
    
    Double_t fInvMassValMax, fInvMassValMin;
    SetBoundaries(fOption,fInvMassValMin,fInvMassValMax);
    
    // Global Variables
    Int_t nEntries          = THdata->GetEntries();
    RooRealVar InvMass      = RooRealVar        ("InvMass","InvMass",fInvMassValMin,fInvMassValMax);
    RooDataHist* data       = new RooDataHist   (fName,fName,InvMass,Import(*THdata));
    RooDataHist* dataLoose  = new RooDataHist   (fName,fName,InvMass,Import(*fLooseErrors(THdata)));
    Int_t kNCycle           = kNCycle_;
    
    // Background PDF Coefficients
    RooRealVar ch0      = RooRealVar        ("ch0","ch0"      ,0.5,   -1, 1);//,0.5,-1,1);
    RooRealVar ch1      = RooRealVar        ("ch1","ch1"      ,-0.1, -1, 1);//,-0.1,-1,1);
    
    RooRealVar ch2, ch3;
    if ( fCheb3 && !fCheb5 )    ch2     = RooRealVar        ("ch2","ch2"        ,0.);
    else                        ch2     = RooRealVar        ("ch2","ch2"        ,0.,-1,1);
    
    if ( fCheb5 )               ch3     = RooRealVar        ("ch3","ch3"        ,0.,-1,1);
    else                        ch3     = RooRealVar        ("ch3","ch3"        ,0.);
    
    //Signal
    RooRealVar sMass, sWidt, sSlop;
    if ( fWidth )               sWidt   = RooRealVar        ("bWidt","bWidt"    ,kPWid);
    else                        sWidt   = RooRealVar        ("bWidt","bWidt"    ,kPWid,kPWid*0.9,kPWid*1.1);

    if ( fMass_ )               sMass   = RooRealVar        ("bMass","bMass"    ,kPMas);
    else                        sMass   = RooRealVar        ("bMass","bMass"    ,kPMas,kPMas*0.9,kPMas*1.1);
    
    if ( bPythiaTest )          sSlop   = RooRealVar        ("bSlop","bSlop"    ,0.);
    else                        sSlop   = RooRealVar        ("bSlop","bSlop"    ,0.001,0.,0.002);
    
    // Coefficients
    RooRealVar nSS      = RooRealVar        ("anSS","anSS"      ,0.5*nEntries,0.,nEntries);
    RooRealVar nBB      = RooRealVar        ("anBB","anBB"      ,0.5*nEntries,0.,nEntries);
    
    // PDFs
    RooVoigtian     fSig= RooVoigtian      ("fSig","fSig"      ,InvMass,sMass,sWidt,sSlop);
    RooChebychev    fBkg= RooChebychev     ("fBkg","fBkg"      ,InvMass,RooArgSet(ch0,ch1,ch2,ch3));
    RooAddPdf       fMod= RooAddPdf        ("fMod","fMod"      ,RooArgList(fBkg,fSig),RooArgList(nBB,nSS));
    
    
    RooFitResult* fFitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        if ( !kOnlyTrue )   {
            fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
            fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
            auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(Signal));
            if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
        }
        if ( kOnlyTrue )    {
            fFitResults      =   fSig.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
            fFitResults      =   fSig.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        }
    }
    
    if ( fSaveToFile )
    {
        hName                       = "InvMass";
        hTitle                      = "Invariant Mass of Kaons in pT 0-6 GeV";
        if ( PTindex != -1 && !bPythiaTest )
        {
            if ( PTDimension == 1 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV",fArrPT1D[PTindex],fArrPT1D[PTindex+1]);
            if ( PTDimension == 2 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1]);
        }
        if ( PTindex != -1 &&  bPythiaTest )
        {
            if ( PTDimension == 1 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV for MC",fArrPT1D[PTindex],fArrPT1D[PTindex+1]);
            if ( PTDimension == 2 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1]);
        }
        
        TCanvas * fSaveToCanvas;
        if ( PTDimension == 1 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName),
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName)
                                                );
        
        if ( PTDimension == 2 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName),
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName)
                                                );
        
        RooPlot * fSaveToFrame      = InvMass.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);
        
        fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        fSaveToCanvas               ->SaveAs(Form("result/SEFitCheck/PT_%.1f_%.1f_1D_%s.pdf",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName));
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return fFitResults;
}

RooFitResult*   FitModel        ( TH1F * THdata, TString fName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTDimension = 1, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    // Check there is a reasonable amount of entries
    if ( !fIsWorthFitting( THdata ) ) return nullptr;
    
    Bool_t  fCheb3  =   false;
    if  ( fOption.find("CH3") != -1 )   fCheb3  =   true;
    
    Bool_t  fCheb5  =   false;
    if  ( fOption.find("CH5") != -1 )   fCheb5  =   true;
    
    Bool_t  fWidth  =   false;
    if  ( fOption.find("W") != -1 )     fWidth  =   true;
    
    Bool_t  fMass_  =   false;
    if  ( fOption.find("M") != -1 )     fMass_  =   true;
    
    Double_t fInvMassValMax, fInvMassValMin;
    SetBoundaries(fOption,fInvMassValMin,fInvMassValMax);
    
    // Global Variables
    Int_t nEntries      = THdata->GetEntries();
    RooRealVar InvMass  = RooRealVar        ("InvMass","InvMass",fInvMassValMin,fInvMassValMax);
    RooDataHist* data   = new RooDataHist   (fName.Data(),fName.Data(),InvMass,Import(*THdata));
    RooDataHist* dataLoose  = new RooDataHist   (fName,fName,InvMass,Import(*fLooseErrors(THdata)));
    Int_t kNCycle       = kNCycle_;
    
    // Background PDF Coefficients
    RooRealVar ch0      = RooRealVar        ("ch0","ch0"      ,0.5,   -1, 1);//,0.5,-1,1);
    RooRealVar ch1      = RooRealVar        ("ch1","ch1"      ,-0.1, -1, 1);//,-0.1,-1,1);
    
    RooRealVar ch2, ch3;
    if ( fCheb3 && !fCheb5 )    ch2     = RooRealVar        ("ch2","ch2"        ,0.);
    else                        ch2     = RooRealVar        ("ch2","ch2"        ,0.,-1,1);
    
    if ( fCheb5 )               ch3     = RooRealVar        ("ch3","ch3"        ,0.,-1,1);
    else                        ch3     = RooRealVar        ("ch3","ch3"        ,0.);
    
    //Signal
    RooRealVar sMass, sWidt, sSlop;
    if ( fWidth )               sWidt   = RooRealVar        ("bWidt","bWidt"    ,kPWid);
    else                        sWidt   = RooRealVar        ("bWidt","bWidt"    ,kPWid,kPWid*0.9,kPWid*1.1);

    if ( fMass_ )               sMass   = RooRealVar        ("bMass","bMass"    ,kPMas);
    else                        sMass   = RooRealVar        ("bMass","bMass"    ,kPMas,kPMas*0.9,kPMas*1.1);
    
    if ( bPythiaTest )          sSlop   = RooRealVar        ("bSlop","bSlop"    ,0.);
    else                        sSlop   = RooRealVar        ("bSlop","bSlop"    ,0.001,0.,0.002);
    
    // Coefficients
    RooRealVar nSS      = RooRealVar        ("anSS","anSS"      ,0.5*nEntries,0.,nEntries);
    RooRealVar nBB      = RooRealVar        ("anBB","anBB"      ,0.5*nEntries,0.,nEntries);
    
    // PDFs
    RooVoigtian     fSig= RooVoigtian      ("fSig","fSig"      ,InvMass,sMass,sWidt,sSlop);
    RooChebychev    fBkg= RooChebychev     ("fBkg","fBkg"      ,InvMass,RooArgSet(ch0,ch1,ch2,ch3));
    RooAddPdf       fMod= RooAddPdf        ("fMod","fMod"      ,RooArgList(fBkg,fSig),RooArgList(nBB,nSS));
    
    RooFitResult* fFitResults;
    fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
    fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
    
    sSlop               .  setRange(0.,0.01);
    sSlop               .  setVal(0.00001);
    sWidt               .  setRange(kPWid*0.5,kPWid*1.5);
    sWidt               .  setVal(kPWid);
    sMass               .  setRange(kPMas*0.5,kPMas*1.5);
    sMass               .  setVal(kPMas);
    
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        if ( !kOnlyTrue )   {
            fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
            fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
            auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(Signal));
            if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
        }
        if ( kOnlyTrue )    {
            fFitResults      =   fSig.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
            fFitResults      =   fSig.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        }
    }
    
    if ( fSaveToFile )
    {
        hName                       = "InvMass";
        hTitle                      = "Invariant Mass of Kaons in pT 0-6 GeV";
        if ( PTindex != -1 && !bPythiaTest )
        {
            if ( PTDimension == 1 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV",fArrPT1D[PTindex],fArrPT1D[PTindex+1]);
            if ( PTDimension == 2 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1]);
        }
        if ( PTindex != -1 &&  bPythiaTest )
        {
            if ( PTDimension == 1 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV for MC",fArrPT1D[PTindex],fArrPT1D[PTindex+1]);
            if ( PTDimension == 2 ) hTitle = Form("Invariant Mass of Kaons in pT %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1]);
        }
        
        TCanvas * fSaveToCanvas;
        if ( PTDimension == 1 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()),
                                                Form("PT_%.1f_%.1f_1D_%s",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data())
                                                );
        
        if ( PTDimension == 2 )fSaveToCanvas   =   new TCanvas(
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data()),
                                                Form("PT_%.1f_%.1f_2D_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data())
                                                );
        
        RooPlot * fSaveToFrame      = InvMass.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);
        
        fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        if ( PTDimension == 1 )fSaveToCanvas               ->SaveAs(Form("result/SEFitCheck/PT_%.1f_%.1f_1D_%s.pdf",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()));
        if ( PTDimension == 2 )fSaveToCanvas               ->SaveAs(Form("result/SEFitCheck/PT_%.1f_%.1f_1D_%s.pdf",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data()));
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return fFitResults;
}

RooFitResult*   FitModel        ( TH2F * THdata, RooFitResult * fFitShapeX, RooFitResult * fFitShapeY, string fHistName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTjndex = -1, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    // Check there is a reasonable amount of entries
    if ( !fIsWorthFitting( THdata ) ) return nullptr;
    
    Bool_t  fBackg  =   false;
    if  ( fOption.find("BK") != -1 )     fBackg  =   true;
    
    Double_t fInvMassValMax, fInvMassValMin;
    fCoreFitModelSetBoundaries(fOption,fInvMassValMin,fInvMassValMax);
    
    // Global Variables
    Int_t nEntries      = THdata->GetEntries();
    RooRealVar varx     = RooRealVar        ("xInvMass2D","xInvMass2D",fInvMassValMin,fInvMassValMax);
    RooRealVar vary     = RooRealVar        ("yInvMass2D","yInvMass2D",fInvMassValMin,fInvMassValMax);
    RooDataHist* data   = new RooDataHist   (fHistName.c_str(),fHistName.c_str(),RooArgList(varx,vary),Import(*THdata));
    RooDataHist* dataLoose  = new RooDataHist   (fHistName.c_str(),fHistName.c_str(),RooArgList(varx,vary),Import(*fLooseErrors(THdata)));
    Int_t kNCycle       = kNCycle_;
    
    RooArgSet  *utilityx    =   new RooArgSet(fFitShapeX->floatParsFinal(),fFitShapeX->constPars());
    RooArgSet  *utilityy    =   new RooArgSet(fFitShapeY->floatParsFinal(),fFitShapeY->constPars());
    
    
    // Background
    RooRealVar ch0x, ch1x, ch2x, ch3x, ch0y, ch1y, ch2y, ch3y;
    
    if ( fBackg )
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0),utilityx->getRealValue("ch0",0)-0.1,utilityx->getRealValue("ch0",0)+0.1);
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0),utilityx->getRealValue("ch1",0)-0.1,utilityx->getRealValue("ch1",0)+0.1);
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0),utilityx->getRealValue("ch2",0)-0.1,utilityx->getRealValue("ch2",0)+0.1);
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0),utilityx->getRealValue("ch3",0)-0.1,utilityx->getRealValue("ch3",0)+0.1);
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0),utilityy->getRealValue("ch0",0)-0.1,utilityy->getRealValue("ch0",0)+0.1);
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0),utilityy->getRealValue("ch1",0)-0.1,utilityy->getRealValue("ch1",0)+0.1);
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0),utilityy->getRealValue("ch2",0)-0.1,utilityy->getRealValue("ch2",0)+0.1);
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0),utilityy->getRealValue("ch3",0)-0.1,utilityy->getRealValue("ch3",0)+0.1);
    }
    else
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0));
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0));
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0));
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0));
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0));
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0));
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0));
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0));
    }
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,utilityx->getRealValue("bMass",0));
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,utilityx->getRealValue("bWidt",0));
    RooRealVar pSlopex  = RooRealVar ("pSlopx","pSlopx" ,utilityx->getRealValue("bSlop",0));
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,utilityy->getRealValue("bMass",0));
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,utilityy->getRealValue("bWidt",0));
    RooRealVar pSlopey  = RooRealVar ("pSlopy","pSlopy" ,utilityy->getRealValue("sSlop",0));
    
    // Coefficients
    auto fS__x  =   utilityx->getRealValue("1nSS",0);
    auto fB__x  =   utilityx->getRealValue("1nBB",0);
    auto fS__y  =   utilityy->getRealValue("1nSS",0);
    auto fB__y  =   utilityy->getRealValue("1nBB",0);
    auto fTot_  =   fS__x*fS__y+fB__x*fB__y+fS__x*fB__y+fB__x*fS__y;
    RooRealVar n0       = RooRealVar ("anSS2D","anSS2D" ,nEntries*(fS__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n1       = RooRealVar ("anBB2D","anBB2D" ,nEntries*(fB__x*fB__y)/fTot_,0.,nEntries);
    RooRealVar n2       = RooRealVar ("anBS2D","anBS2D" ,nEntries*(fB__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n3       = RooRealVar ("anSB2D","anSB2D" ,nEntries*(fS__x*fB__y)/fTot_,0.,nEntries);
    
    // PDFs
    RooChebychev        fBkgx ("fBkgx","fBkgx"          ,varx,RooArgSet(ch0x,ch1x,ch2x,ch3x));
    RooVoigtian         fSigx ("fSigx","fSigx"          ,varx,pMassx,pWidthx,pSlopex);
    RooChebychev        fBkgy ("fBkgy","fBkgy"          ,vary,RooArgSet(ch0y,ch1y,ch2y,ch3y));
    RooVoigtian         fSigy ("fSigy","fSigy"          ,vary,pMassy,pWidthy,pSlopey);
    RooProdPdf          fBB   ("fBkg","fBkg"            ,fBkgx,fBkgy);
    RooProdPdf          fSB   ("fSigBkg","fSBWBkg"      ,fSigx,fBkgy);
    RooProdPdf          fBS   ("fBkgSig","fBkgSig"      ,fBkgx,fSigy);
    RooProdPdf          fSS   ("fSigSig","fSigSig"      ,fSigx,fSigy);
    RooAddPdf           fMod  ("fMod2D","fMod2D"        ,RooArgList(fBB,fSS,fSB,fBS),RooArgList(n1,n0,n3,n2));
    
    RooFitResult* fFitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        //fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        auto N_Raw      =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(Signal));
        if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError(),20.) ) break;
    }
    
    // Save to file
    if ( fSaveToFile )
    {
        int         nBinsPrint      =   3;
        double      dIncrement      =   (fInvMassValMax-fInvMassValMin)/nBinsPrint;
        TLatex*     latext          =   new TLatex();
        TCanvas*    cTotal          =   new TCanvas("","",0,45,1440,855);
                    cTotal          ->  SetTitle(Form("Slices of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
                    cTotal          ->  SetName(Form("PT_%.1f_%.1f__%.1f_%.1f_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()));
                    cTotal          ->  Divide(2,nBinsPrint);
        
                            varx.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                            vary.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()),
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame  =   vary.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);

                            varx.setRange("fDrawRange",fInvMassValMin+i*dIncrement,fInvMassValMin+(i+1)*dIncrement);
                            vary.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);

            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"yInvMass2D");
            
            cTotal->cd( i+1 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.3f < m^{x}_{K^{+}K^{-}} < %.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.3f < m^{x}_{K^{+}K^{-}} < %.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        varx.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        vary.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()),
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame      =   varx.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);
            
                                        varx.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        vary.setRange("fDrawRange",fInvMassValMin+i*dIncrement,fInvMassValMin+(i+1)*dIncrement);
                                                                            
            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"xInvMass2D");
            
            cTotal->cd( i+1 +3 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.2f < m^{y}_{K^{+}K^{-}} < %.2f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.2f < m^{y}_{K^{+}K^{-}} < %.2f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        varx.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        vary.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        cTotal ->Write();
        cTotal               ->SaveAs(Form("result/SEFitCheck/PT_%.1f_%.1f__%.1f_%.1f_%s.pdf",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()));
        delete cTotal;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    // Fit
    return fFitResults;
}

RooFitResult*   fExtrapolateModelLEGACYROOFIT               ( TH1F *HData, TString fName = "ExtrapolateSignal" ) {
    
    // Check is worth fitting
    if ( !fIsWorthFitting( HData ) ) return nullptr;
    
    auto nEntries = HData->GetEntries();
    
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    // Global Variables
    RooRealVar TransMom = RooRealVar        ("TransMom","TransMom",0.4,1.6);
    RooDataHist* RData  = new RooDataHist   (fName.Data(),fName.Data(),TransMom,Import(*HData));
    
    // Signal PDF Parameters
    RooRealVar          n_value ("n_value", "n_value",  6.7,    2.01,    20.);
    RooRealVar          exp_par ("exp_par", "exp_par",  .272,   0.1,    2.);
    RooRealVar          prt_mss ("prt_mss", "prt_mss",  kPhiMesonMass_);
    
    // Normalisation Coefficients
    RooRealVar          Sig_str ("Sig_str", "Sig_str",  0.033, 0.2, 0.4);
    
    // Formulas
    TString             LevyTsl ("((x[0]-1)*(x[0]-2))/(x[0]*x[1]*(x[0]*x[1]+x[2]*(x[0]-2)))*x[3]*(TMath::Power(1+(sqrt(x[2]*x[2]+x[3]*x[3])-x[2])/(x[0]*x[1]),(-x[0])))");

    // PDFs
    RooGenericPdf       fModel  ("fModel",  "fModel",   LevyTsl.Data(), RooArgList( n_value,    exp_par,    prt_mss,  TransMom));
    
    RooFitResult *result;
    
    //TransMom.setRange("FitRange",0.4,1.6);
    result = fModel.chi2FitTo(*RData,SumW2Error(kTRUE),Save(),Minos(true),NumCPU(kCPU_use));
    result = fModel.fitTo(*RData,SumW2Error(kTRUE),Save(),Minos(true),NumCPU(kCPU_use));
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    
    auto fSaveToCanvas   =   new TCanvas();
    gPad->SetLogy();
    
    RooPlot * fSaveToFrame      = TransMom.frame(Name(fName.Data()),Title(fName.Data()));
    
    RData                           ->plotOn(fSaveToFrame,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
    fModel                          .plotOn (fSaveToFrame,      LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
    
    fSaveToFrame                ->Draw("same");
    fSaveToCanvas               ->Write();
    fSaveToCanvas               ->SaveAs("check.pdf");
    
    return result;
}

//
template < class Tclass >
Double_t*       fExtrapolateModel               ( Tclass *THdata, TString fName = "ExtrapolateSignal" ) {
    // Optimisation mode
    gROOT->SetBatch(true);
    
    // Result format: Integral, Stat err, Syst err, Mean pT, Stat err, Syst err
    Double_t   *fResult = new   Double_t    [6];
    
    // Setting -1. for the default
    for ( Int_t iFill = 0; iFill < 6; iFill++ ) fResult[iFill]  =   -1.;
    
    // Set Standard Fit
    fSetLevyTsalis();
    
    // Fit the Spectra
    THdata->Fit(fLevyFit1D,"IMREQ0S","",fMinPT1D,10.);
    
    // Save to further checks
    TCanvas * fCheckFit = new TCanvas();
    gPad->SetLogy();
    THdata      ->Draw("same");
    fLevyFit1D  ->Draw("same");
    fCheckFit   ->Write();
    fCheckFit   ->SaveAs(Form("tmp/%s.pdf",fName.Data()));
    
    fResult[0]  =   fLevyFit1D->Integral(0.,fMinPT1D);
    fResult[1]  =   fLevyFit1D->IntegralError(0.,fMinPT1D);
    
    return fResult;
    
    // End Optimisation mode
    gROOT->SetBatch(false);
}
//
//_____________________________________________________________________________
//
template < class Tclass >
Double_t*           fIntegrateModel                 ( Tclass *THdata, TString fName = "IntegrateSignal" )    {
    // Optimisation mode
    gROOT->SetBatch(true);
    
    // Result format: Integral, Stat err, Syst err, Mean pT, Stat err, Syst err
    Double_t   *fResult = new   Double_t    [6];
    
    // Setting -1. for the default
    for ( Int_t iFill = 0; iFill < 6; iFill++ ) fResult[iFill]  =   -1.;
    
    fResult[0]  =   THdata  ->IntegralAndError(-1.,100,fResult[1],"width");
    //fSetSystErr(THdata)     ->IntegralAndError(-1.,100,fResult[2],"width");
    
    return fResult;
    
    // End Optimisation mode
    gROOT->SetBatch(false);
}
//
//_____________________________________________________________________________
//
template < class Tclass >
Double_t*           fMeasureFullYield               ( Tclass *THdata, TString fName = "MeasureFullYield" ) {
    // Optimisation mode
    gROOT->SetBatch(true);
    
    // Result format: Integral, Stat err, Syst err, Mean pT, Stat err, Syst err
    Double_t   *fResult             =   new Double_t        [6];
    
    // Setting -1. for the default
    for ( Int_t iFill = 0; iFill < 6; iFill++ ) fResult[iFill]  =   -1.;
    
    auto        fIntegralResults    =   fIntegrateModel     (THdata,fName);
    auto        fExtrapolResults    =   fExtrapolateModel   (THdata,fName);
    
    fResult[0]  =   fIntegralResults[0] +   fExtrapolResults[0];
    fResult[1]  =   fIntegralResults[1] +   fExtrapolResults[1];
    
    return fResult;
    
    // End Optimisation mode
    gROOT->SetBatch(false);
}
//
//_____________________________________________________________________________
//
RooFitResult*   FitModelRap                         ( TH1D * THdata, TH1D * THbkg1,  TH1D * THbkg2,  const char* fName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTDimension = 1, string fOption = "" )
{
    //FitModelRap(      hBKG_SIG_Rap,       hBKG_BKG_Rap,       hBKG_BKG_BKG_SIG_Rap);
    
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    // Check there is a reasonable amount of entries
    if ( !fIsWorthFitting( THdata ) ) return nullptr;
    
    // Global Variables
    Int_t nDataEntries          = THdata->GetEntries();
    Int_t nBkg1Entries          = THbkg1->GetEntries();
    Int_t nBkg2Entries          = THbkg2->GetEntries();
    RooRealVar fRap_        = RooRealVar        ("Rapidity","Rapidity",-1.,1.);
    RooDataHist* data       = new RooDataHist   (fName,fName,fRap_,Import(*THdata));
    RooDataHist* dataLoose  = new RooDataHist   (fName,fName,fRap_,Import(*fLooseErrors(THdata)));
    RooDataHist* bkg1       = new RooDataHist   (fName,fName,fRap_,Import(*THbkg1));
    RooDataHist* bkg1Loose  = new RooDataHist   (fName,fName,fRap_,Import(*fLooseErrors(THbkg1)));
    RooDataHist* bkg2       = new RooDataHist   (fName,fName,fRap_,Import(*THbkg2));
    RooDataHist* bkg2Loose  = new RooDataHist   (fName,fName,fRap_,Import(*fLooseErrors(THbkg2)));
    Int_t kNCycle           = kNCycle_;
    
    // Combinatorial Background
    //
    //>->-> Variables
    //
    RooRealVar         *rBkg1_Cf00      =   new RooRealVar      ("rBkg1_Cf00","rBkg1_Cf00",     .05, 0., 1.);
    RooRealVar         *rBkg1_Cf01      =   new RooRealVar      ("rBkg1_Cf01","rBkg1_Cf01",     .5, -100., 100);
    RooRealVar         *rBkg1_Cf02      =   new RooRealVar      ("rBkg1_Cf02","rBkg1_Cf02",     .005, 0., 1);
    //
    RooRealVar         *rBkg1_Mean      =   new RooRealVar      ("rBkg1_Mean","rBkg1_Mean",     0.);
    RooRealVar         *rBkg1_Widt      =   new RooRealVar      ("rBkg1_Widt","rBkg1_Widt",     .005, 0., 1);
    // Formulas
    TString             fBkg1_Expo  ("TMath::Sqrt( +x[3]*x[3]*(x[0]-1)*(x[0]-1)*(x[0]+1)*(x[0]+1) + TMath::Power( 1 + TMath::Sqrt( (x[0]*x[0] + x[1]*x[1] ) ) , -2*x[2]) )");
    //
    // Coefficients
    RooRealVar         *rBkg1_Expo__Norm=   new RooRealVar     ("rBkg1_Expo__Norm",  "rBkg1_Expo__Norm",  0.5*nBkg1Entries,   0., nBkg1Entries);
    RooRealVar         *rBkg1_Gauss_Norm=   new RooRealVar     ("rBkg1_Gauss_Norm",  "rBkg1_Gauss_Norm",  0.5*nBkg1Entries,   0., nBkg1Entries);
    //
    // PDFs
    RooGenericPdf      *fBkg1_Expo_     =   new RooGenericPdf   ("fBkg1_Cstm_",     "fBkg1_Cstm_",   fBkg1_Expo.Data(), RooArgList( fRap_,    *rBkg1_Cf00, *rBkg1_Cf01, *rBkg1_Cf02));
    RooGaussian        *fBkg1_Gauss     =   new RooGaussian     ("fBkg1_Gauss",     "fBkg1_Gauss",   fRap_,    *rBkg1_Mean,    *rBkg1_Widt);
    RooAddPdf          *fBkg1_Model     =   new RooAddPdf       ("fBkg1_Model",     "fBkg1_Model",  RooArgList( *fBkg1_Expo_, *fBkg1_Gauss ), RooArgList(*rBkg1_Expo__Norm, *rBkg1_Gauss_Norm) );
    //
    // Combinatorial Background Fit
    RooFitResult* fFitRsltBkg1;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )    {
        fFitRsltBkg1    =   fBkg1_Model->fitTo(*bkg1Loose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        fFitRsltBkg1    =   fBkg1_Model->fitTo(*bkg1,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        //auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(0));
        //if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
    }
    //
    // Combinatorial Background + Combinatorial Phis
    //
    //>->-> Recovering Previous Shape
    //
    RooRealVar         *rBkg2_Cf00      =   new RooRealVar      ("rBkg2_Cf00",  "rBkg2_Cf00",   rBkg1_Cf00->getValV());
    RooRealVar         *rBkg2_Cf01      =   new RooRealVar      ("rBkg2_Cf01",  "rBkg2_Cf01",   rBkg1_Cf01->getValV());
    RooRealVar         *rBkg2_Cf02      =   new RooRealVar      ("rBkg2_Cf02",  "rBkg2_Cf02",   rBkg1_Cf02->getValV());
    //
    RooRealVar         *rBkg2_Mean      =   new RooRealVar      ("rBkg2_Mean",  "rBkg2_Mean",   rBkg1_Mean->getValV());
    RooRealVar         *rBkg2_Widt      =   new RooRealVar      ("rBkg2_Widt",  "rBkg2_Widt",   rBkg1_Widt->getValV());
    //
    RooGenericPdf      *fBkg2_Expo_     =   new RooGenericPdf   ("fBkg2_Expo_", "fBkg2_Expo_",  fBkg1_Expo.Data(), RooArgList( fRap_,    *rBkg2_Cf00, *rBkg2_Cf01, *rBkg2_Cf02));
    RooGaussian        *fBkg2_Gauss     =   new RooGaussian     ("fBkg2_Gauss", "fBkg2_Gauss",  fRap_,    *rBkg2_Mean,    *rBkg2_Widt);
    //
    RooRealVar         *rBkg2_Expo__Norm=   new RooRealVar     ("rBkg2_Expo__Norm",  "rBkg2_Expo__Norm",  rBkg1_Expo__Norm->getValV()/nBkg1Entries);
    RooRealVar         *rBkg2_Gauss_Norm=   new RooRealVar     ("rBkg2_Gauss_Norm",  "rBkg2_Gauss_Norm",  rBkg1_Gauss_Norm->getValV()/nBkg1Entries);
    //
    RooAddPdf          *fBkg2_MdBk1     =   new RooAddPdf       ("fBkg2_MdBk1", "fBkg2_MdBk1",  RooArgList( *fBkg2_Expo_, *fBkg2_Gauss ), RooArgList(*rBkg2_Expo__Norm, *rBkg2_Gauss_Norm) );
    //
    //>->-> Variables
    //
    RooRealVar         *rBkg2_Slop      =   new RooRealVar      ("rBkg1_Cf00",  "rBkg1_Cf00",     .5,    0.1,     1.);
    //
    // Formulas
    TString             fBkg2_Lin       ("TMath::Max(0.,1.*-TMath::Abs(x[0])*(x[1]) + x[1])");
    //
    // Coefficients
    RooRealVar         *rBkg2_Bkg1__Norm=   new RooRealVar     ("rBkg2_Bkg1__Norm",  "rBkg2_Bkg1__Norm",  0.5*nBkg2Entries,   0., nBkg2Entries);
    RooRealVar         *rBkg2_Tria__Norm=   new RooRealVar     ("rBkg2_Tria__Norm",  "rBkg2_Tria__Norm",  0.5*nBkg2Entries,   0., nBkg2Entries);
    //
    // PDFs
    RooGenericPdf      *fBkg2_Tria     =   new RooGenericPdf   ("fBkg2_Tria",       "fBkg2_Tria",   fBkg2_Lin.Data(), RooArgList( fRap_,    *rBkg2_Slop));
    RooAddPdf          *fBkg2_Model    =   new RooAddPdf       ("fBkg2_Model",      "fBkg2_Model",  RooArgList( *fBkg2_Tria, *fBkg2_MdBk1 ), RooArgList(*rBkg2_Tria__Norm, *rBkg2_Bkg1__Norm) );
    //
    RooFitResult* fFitRsltBkg2;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )    {
        fFitRsltBkg2    =   fBkg2_Model->fitTo(*bkg2Loose,  Extended(kTRUE),    SumW2Error(kTRUE),  Save(), NumCPU(kCPU_use));
        fFitRsltBkg2    =   fBkg2_Model->fitTo(*bkg2,       Extended(kTRUE),    SumW2Error(kTRUE),  Save(), NumCPU(kCPU_use));
        //auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(0));
        //if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
    }
    //
    // Signal Extraction
    //
    //>->-> Recovering Previous Shape
    //
    RooRealVar         *rSign_BkSl      =   new RooRealVar      ("rSign_BkSl",      "rSign_BkSl",       rBkg2_Slop->getValV());
    //
    RooRealVar         *rSign_Bkg1__Norm=   new RooRealVar      ("rSign_Bkg1__Norm", "rSign_Bkg1__Norm",rBkg2_Bkg1__Norm->getValV()/nBkg2Entries);
    RooRealVar         *rSign_BkTr__Norm=   new RooRealVar      ("rSign_BkTr__Norm", "rSign_BkTr__Norm",rBkg2_Tria__Norm->getValV()/nBkg2Entries);
    //
    RooGenericPdf      *fSign_BkTr     =   new RooGenericPdf    ("fSign_BkTr",       "fSign_BkTr",      fBkg2_Lin.Data(), RooArgList( fRap_,    *rSign_BkSl));
    RooAddPdf          *fSign_MdBk2    =   new RooAddPdf        ("fSign_MdBk2",      "fSign_MdBk2",     RooArgList( *fSign_BkTr, *fBkg2_MdBk1 ), RooArgList(*rSign_BkTr__Norm, *rSign_Bkg1__Norm) );
    //
    //
    //>->-> Variables
    //
    RooRealVar         *rSign_Mean      =   new RooRealVar      ("rSign_Mean","rSign_Mean",     0.);
    RooRealVar         *rSign_Widt      =   new RooRealVar      ("rSign_Widt","rSign_Widt",     .5,    0.1,     1.);
    //
    RooRealVar         *rSign_Slop      =   new RooRealVar      ("rSign_Slop",  "rSign_Slop",     .5,    0.1,     1.);
    //
    // Coefficients
    RooRealVar         *rSign_FBkg_Norm=   new RooRealVar       ("rSign_FBkg_Norm",  "rSign_FBkg_Norm",  0.5*nDataEntries,     0., nDataEntries);
    RooRealVar         *rSign_Tria_Norm=   new RooRealVar       ("rSign_Tria_Norm",  "rSign_Tria_Norm",  0.5*nDataEntries,     0., nDataEntries);
    RooRealVar         *rSign_Sign_Norm=   new RooRealVar       ("rSign_Sign_Norm",  "rSign_Sign_Norm",  0.5*nDataEntries,     0., nDataEntries);
    //
    // PDFs
    RooGenericPdf      *fSign_Tria      =   new RooGenericPdf   ("fSign_Tria",       "fSign_Tria",  fBkg2_Lin.Data(), RooArgList( fRap_,    *rSign_Slop));
    RooGaussian        *fSign_Sign      =   new RooGaussian     ("fSign_Sign",      "fSign_Sign",   fRap_,  *rSign_Mean,   *rSign_Widt);
    RooAddPdf          *fSign_Model     =   new RooAddPdf       ("fSign_Model",     "fSign_Model",  RooArgList( *fSign_Sign,    *fSign_Tria,  *fSign_MdBk2 ), RooArgList(*rSign_Sign_Norm, *rSign_Tria_Norm, *rSign_FBkg_Norm) );
    //
    RooFitResult* fFitRsltSign;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )    {
        fFitRsltSign    =   fSign_Model->fitTo(*dataLoose,  Extended(kTRUE),    SumW2Error(kTRUE),  Save(), NumCPU(kCPU_use));
        fFitRsltSign    =   fSign_Model->fitTo(*data,       Extended(kTRUE),    SumW2Error(kTRUE),  Save(), NumCPU(kCPU_use));
        //auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(0));
        //if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
    }
    //
    
    /*
     
     //>->-> Variables
     //
     //
     // Formulas
     TString             fBkg2_Lin       ("TMath::Max(0.,1.*-TMath::Abs(x[0])*(x[1]) + x[1])");
     //
     // Coefficients
     RooRealVar         *rBkg2_Bkg1__Norm=   new RooRealVar     ("rBkg2_Bkg1__Norm",  "rBkg2_Bkg1__Norm",  0.5*nBkg2Entries,   0., nBkg2Entries);
     RooRealVar         *rBkg2_Tria__Norm=   new RooRealVar     ("rBkg2_Tria__Norm",  "rBkg2_Tria__Norm",  0.5*nBkg2Entries,   0., nBkg2Entries);
     //
     // PDFs
     RooGenericPdf      *fBkg2_Tria     =   new RooGenericPdf   ("fBkg2_Tria",       "fBkg2_Tria",   fBkg2_Lin.Data(), RooArgList( fRap_,    *rBkg2_Slop));
     */
    if ( true )
    {
        TCanvas * fSaveToCanvas     = new TCanvas();
        RooPlot * fSaveToFrame      = fRap_.frame();
        TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);
        
        bkg1                        ->plotOn(fSaveToFrame,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
        fBkg1_Model                 ->plotOn (fSaveToFrame,      LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        fSaveToCanvas               ->SaveAs(Form("result/SEFitCheck/RaP_bkg1_%s.pdf",fName));
        fSaveToCanvas               ->SetLogy();
        fSaveToCanvas               ->SaveAs(Form("result/SEFitCheck/RaP_bkg1_logy_%s.pdf",fName));
        
        TCanvas * fSaveToCanva2     = new TCanvas();
        RooPlot * fSaveToFram2      = fRap_.frame();
        TLegend * fLegen2           = new TLegend   (0.12,0.60,0.30,0.85);
            
        bkg2                        ->plotOn    (fSaveToFram2,  MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
        fBkg2_Model                 ->plotOn    (fSaveToFram2,  LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
        fBkg2_Model                 ->plotOn    (fSaveToFram2,  Components("fBkg2_Tria"),       LineStyle(kDashed), LineColor(38),      Name("RooBB"));
        fBkg2_Model                 ->plotOn    (fSaveToFram2,  Components("fBkg2_MdBk1"),      LineStyle(kDashed), LineColor(33),      Name("RooBS"));

        fLegen2                     ->SetFillColor(kWhite);
        fLegen2                     ->SetLineColor(kWhite);
        fLegen2                     ->AddEntry(fSaveToFram2->findObject("RooData"), "Data",                 "EP");
        fLegen2                     ->AddEntry(fSaveToFram2->findObject("RooMod"),  "Full Model",           "L");
        fLegen2                     ->AddEntry(fSaveToFram2->findObject("RooBB"),   "Combinatorial Phis",   "L");
        fLegen2                     ->AddEntry(fSaveToFram2->findObject("RooBS"),   "Combinatorial Kaons",  "L");
    
        fSaveToFram2                ->Draw("same");
        fLegen2                     ->Draw("same");
        fSaveToCanva2               ->Write ();
        fSaveToCanva2               ->SaveAs(Form("result/SEFitCheck/RaP_bkg2_%s.pdf",fName));
        fSaveToCanva2               ->SetLogy();
        fSaveToCanva2               ->SaveAs(Form("result/SEFitCheck/RaP_bkg2_logy_%s.pdf",fName));
        
        TCanvas * fSaveToCanva3     = new TCanvas();
        RooPlot * fSaveToFram3      = fRap_.frame();
        TLegend * fLegen3           = new TLegend   (0.12,0.60,0.30,0.85);
            
        data                        ->plotOn    (fSaveToFram3,  MarkerColor(38),            MarkerStyle(26),    Name("RooData"));
        fSign_Model                 ->plotOn    (fSaveToFram3,  LineColor(4),               LineStyle(kDashed), Name("RooMod"));
        fSign_Model                 ->plotOn    (fSaveToFram3,  Components("fBkg2_MdBk1"),  LineStyle(kDashed), LineColor(38),      Name("RooBB"));
        fSign_Model                 ->plotOn    (fSaveToFram3,  Components("fSign_Tria"),   LineStyle(kDashed), LineColor(2),       Name("RooBS"));
        fSign_Model                 ->plotOn    (fSaveToFram3,  Components("fSign_Sign"),   LineStyle(kSolid),  LineColor(2),       Name("RooSS"));
        
        fLegen3                     ->SetFillColor(kWhite);
        fLegen3                     ->SetLineColor(kWhite);
        fLegen3                     ->AddEntry(fSaveToFram3->findObject("RooData"), "Data",                 "EP");
        fLegen3                     ->AddEntry(fSaveToFram3->findObject("RooMod"),  "Full Model",           "L");
        fLegen3                     ->AddEntry(fSaveToFram3->findObject("RooBB"),   "Kaon Combinatorial",   "L");
        fLegen3                     ->AddEntry(fSaveToFram3->findObject("RooBS"),   "Secondary Signal",  "L");
        fLegen3                     ->AddEntry(fSaveToFram3->findObject("RooSS"),   "Signal",  "L");
            
        fSaveToFram3                ->Draw("same");
        fLegen3                     ->Draw("same");
        fSaveToCanva3               ->Write ();
        fSaveToCanva3               ->SaveAs(Form("result/SEFitCheck/RaP_data_%s.pdf",fName));
        fSaveToCanva3               ->SetLogy();
        fSaveToCanva3               ->SaveAs(Form("result/SEFitCheck/RaP_data_logy_%s.pdf",fName));
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return nullptr; //fFitResults;
    
    /*
    // Combinatorial Background
    RooRealVar          Intercpt ("Intercpt", "Intercpt",  1.,    .5,    10.);
    
    // Formulas
    TString             AbsLin ("-TMath::Abs(x[0])*(x[1])+x[1]");//((x[0]-1)*(x[0]-2))/(x[0]*x[1]*(x[0]*x[1]+x[2]*(x[0]-2)))*x[3]*(TMath::Power(1+(sqrt(x[2]*x[2]+x[3]*x[3])-x[2])/(x[0]*x[1]),(-x[0])))");
     
    RooRealVar ch0      = RooRealVar        ("ch0","ch0"      ,0.5,   -1, 1);//,0.5,-1,1);
    RooRealVar ch1      = RooRealVar        ("ch1","ch1"      ,-0.1, -1, 1);//,-0.1,-1,1);
    RooRealVar ch2      = RooRealVar        ("ch2","ch2"      ,0.01,-1, 1);//,0.01,-1,1);
    RooRealVar ch3      = RooRealVar        ("ch3","ch3"      ,-0.05,-1, 1);//,-0.05,-1,1);
    
    RooRealVar ch4, ch5;
    if ( fCheb3 && !fCheb5 )    ch4     = RooRealVar        ("ch4","ch4"        ,0.);
    else                        ch4     = RooRealVar        ("ch4","ch4"        ,0.,-1,1);
    
    if ( fCheb5 )               ch5     = RooRealVar        ("ch5","ch5"        ,0.,-1,1);
    else                        ch5     = RooRealVar        ("ch5","ch5"        ,0.);
    
    //Signal
    RooRealVar sMass, sWidt, sSlop;
    if ( fWidth )               sWidt   = RooRealVar        ("sWidt","sWidt"    ,kPWid);
    else                        sWidt   = RooRealVar        ("sWidt","sWidt"    ,kPWid,kPWid*0.9,kPWid*1.1);
    
    if ( fMass_ )               sMass   = RooRealVar        ("sMass","sMass"    ,kPMas);
    else                        sMass   = RooRealVar        ("sMass","sMass"    ,kPMas,kPMas*0.9,kPMas*1.1);
    
    if ( bPythiaTest )          sSlop   = RooRealVar        ("sSlop","sSlop"    ,0.);
    else                        sSlop   = RooRealVar        ("sSlop","sSlop"    ,0.001,0.,0.002);
    
    // Coefficients
    RooRealVar nSS      = RooRealVar        ("1nSS","1nSS"      ,0.5*nEntries,0.,nEntries);
    RooRealVar nBB      = RooRealVar        ("1nBB","1nBB"      ,0.5*nEntries,0.,nEntries);
    
    // PDFs
    RooVoigtian     fSig= RooVoigtian      ("fSig","fSig"      ,fRap_,sMass,sWidt,sSlop);
    RooChebychev    fBkg= RooChebychev     ("fBkg","fBkg"      ,fRap_,RooArgSet(ch0,ch1,ch2,ch3,ch4,ch5));
    RooAddPdf       fMod= RooAddPdf        ("fMod","fMod"      ,RooArgList(fBkg,fSig),RooArgList(nBB,nSS));
    
    
    RooFitResult* fFitRsltBkg1;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        fFitResults      =   fModel.fitTo(*bkg2Loose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        fFitResults      =   fModel.fitTo(*bkg2,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        //auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(0));
        //if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
    }
    */
}
//
//_____________________________________________________________________________







#endif
