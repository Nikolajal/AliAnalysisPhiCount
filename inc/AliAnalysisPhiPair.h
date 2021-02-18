#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H

// Analysis Utility
#include "./AliAnalysisUtility/AliAnalysisUtility.h"

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//

enum            fFitResults1D
{
    Background, Signal
};

enum            fFitResults2D
{

    BackgBackg, BackgSignl, SignlBackg, SignlSignl
};

// Performance Values
auto const  kCPU_use                =   3;
auto const  kNCycle_                =   5;
auto const  kStatEvalCycles         =   10;

// Analysis Values
auto const  bPythiaTest             =   kFALSE;

//-// File Names

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

//-//-// Others
auto const  fInvMasHist             =   "./result/InvariantMassHistograms.root";
auto const  fEfficiHist             =   "./result/Efficiencies_MCTruth.root";
auto const  fFitResHist             =   "./result/InvariantMassFitResultsPlots.root";
auto const  fFitResults             =   "./result/InvariantMassFitResults.root";
auto const  fAnlResHist             =   "./result/AnalysisResultsPlots.root";
auto const  fAnlResults             =   "./result/AnalysisResults.root";
auto const  fSystError_             =   "./result/Syst_SigExt.root";

//-// Tree Names
auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaonCandidate_Tree     =   "KaonCandidate";
auto const  fKaonCandidateEff_Tree  =   "KaonEfficiency";

// Analysis Values
//-// Analysis constants
auto const  kPhiMesonMass_          =   1.019455;   //  1.019455    +- 0.000020
auto const  kPhiMesonWidth          =   0.00426;    //  0.00426     +- 0.00004
auto const  kKaonMass               =   .493677;
auto const  kKaonMassUncert         =   .000013;
auto const  kDetectorSlope          =   1.;
auto const  kBranchingRtio          =   0.492;
auto const  kVertexEfficnc          =   1;
auto const  kTriggerEfficnc         =   1;
auto const  kRapidityIntvl          =   1;

//-// Analysis systematics (%)
auto const  kBRSystematics          =   0.01;

auto const  kPMas                   =   1.019455;   //  1.019455    +- 0.000020
auto const  kPWid                   =   0.00426;    //  0.00426     +- 0.00004

//-//   Analysis settings
auto        kDoMultiplicity         =   true;
auto        kDoYield                =   true;
auto        kDoTrigger              =   true;

//-// InvMass range Pythia MC
const   Float_t   fMinIMMC  =   0.75;
const   Float_t   fMaxIMMC  =   1.25;

//-// InvMass bins 1D
const   Int_t     nBinIM1D  =   135;
const   Float_t   fMinIM1D  =   0.99;
const   Float_t   fMaxIM1D  =   1.08;
        Float_t * fArrIM1D  =   new Float_t [nBinIM1D+1];

//-// InvMass bins 2D
const   Int_t     nBinIM2D  =   135;
const   Float_t   fMinIM2D  =   0.99;
const   Float_t   fMaxIM2D  =   1.08;
        Float_t * fArrIM2D  =   new Float_t [nBinIM2D+1];

//-// pT bins 1D
const   Int_t     nBinPT1D  =   21;
const   Float_t   fMinPT1D  =   0.40;
const   Float_t   fMaxPT1D  =   21.0;
        Float_t  *fArrPT1D  =   new Float_t [nBinPT1D+1];

//-// pT bins 2D
const   Int_t     nBinPT2D  =   10;
const   Float_t   fMinPT2D  =   0.40;
const   Float_t   fMaxPT2D  =   21.0;
        Float_t  *fArrPT2D  =   new Float_t [nBinPT2D+1];

//-// Muliplicity bins
const   Int_t     nBinMult  =   5;
const   Float_t   fMinMult  =   0.0;
const   Float_t   fMaxMult  =   100.0;
        Float_t  *fArrMult  =   new Float_t [nBinMult+1];

//-// Rapidity bins
const   Int_t     nBinRap_  =   20;
const   Float_t   fMinRap_  =   0.0;
const   Float_t   fMaxRap_  =   1.0;
        Float_t  *fArrRap_  =   new Float_t [nBinRap_+1];

//-// N-Tuples bins
const   Int_t     nBinNTup  =   5;
const   Float_t   fMinNTup  =   -0.5;
const   Float_t   fMaxNTup  =   4.5;
        Float_t  *fArrNTup  =   new Float_t [nBinNTup+1];

//------------------------------//
//    DATA STRUCTURES           //
//------------------------------//

typedef struct
{
    UChar_t nPhi,           Multiplicity,   iKaon[1024],    jKaon[1024],   Nature[1024];
    Float_t Px[1024],       Py[1024],       Pz[1024],       pT[1024],      Rap[1024],        InvMass[1024];
} Struct_PhiCandidate;

typedef struct
{
    UChar_t nKaon,          Multiplicity,   Charge[1024];
    Char_t  SigmaTOF[1024], SigmaTPC[1024];
    Float_t Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_KaonCandidate;

typedef struct
{
    UChar_t nPhi,           Multiplicity,   Selection[1024];
    Float_t Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  fTru,           fGen,           fRec;
} Struct_PhiEfficiency;

typedef struct
{
    UChar_t nKaon,          Multiplicity,   Charge[1024],   Selection[1024];
    Float_t Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  ftru;
} Struct_KaonEfficiency;

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
    fArrPT1D[0] =   0.4;
    fArrPT1D[1] =   0.6;
    fArrPT1D[2] =   0.8;
    fArrPT1D[3] =   1.0;
    fArrPT1D[4] =   1.2;
    fArrPT1D[5] =   1.4;
    fArrPT1D[6] =   1.6;
    fArrPT1D[7] =   1.8;
    fArrPT1D[8] =   2.0;
    fArrPT1D[9]=   2.5;
    fArrPT1D[10]=   3.0;
    fArrPT1D[11]=   3.5;
    fArrPT1D[12]=   4.0;
    fArrPT1D[13]=   4.5;
    fArrPT1D[14]=   5.0;
    fArrPT1D[15]=   6.0;
    fArrPT1D[16]=   7.0;
    fArrPT1D[17]=   8.0;
    fArrPT1D[18]=   10.0;
    fArrPT1D[19]=   13.0;
    fArrPT1D[20]=   16.0;
    fArrPT1D[21]=   21.0;
}

void    fSetBinPT2D ()
{
    fArrPT2D[0] =   0.40;
    fArrPT2D[1] =   0.68;
    fArrPT2D[2] =   0.82;
    fArrPT2D[3] =   0.95;
    fArrPT2D[4] =   1.1;
    fArrPT2D[5] =   1.3;
    fArrPT2D[6] =   1.6;
    fArrPT2D[7] =   2.3;
    fArrPT2D[8] =  3.0;
    fArrPT2D[9] =  5.0;
    fArrPT2D[10] =  21.;
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
    fArrRap_[1]  =  0.025;
    fArrRap_[2]  =  0.05;
    fArrRap_[3]  =  0.075;
    fArrRap_[4]  =  0.1;
    fArrRap_[5]  =  0.125;
    fArrRap_[6]  =  0.150;
    fArrRap_[7]  =  0.175;
    fArrRap_[8]  =  0.2;
    fArrRap_[9]  =  0.25;
    fArrRap_[10]  =  0.3;
    fArrRap_[11]  =  0.35;
    fArrRap_[12]  =  0.4;
    fArrRap_[13]  =  0.45;
    fArrRap_[14]  =  0.5;
    fArrRap_[15]  =  0.55;
    fArrRap_[16]  =  0.6;
    fArrRap_[17]  =  0.7;
    fArrRap_[18]  =  0.8;
    fArrRap_[19]  =  0.9;
    fArrRap_[20]  =  1.0;
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
        if ( input_value <= fArrRap_[iBin] )
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

//------------------------------//
//      CUTS UTILITIES          //
//------------------------------//

bool    fCutRapidity        ( Double_t  dRapidity )
{
    if ( fabs(dRapidity) <= 0.5 ) return true;
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
    if ( dTransverseMom <= 0.4 ) return false;
    if ( dTransverseMom >= fMaxPT1D ) return false;
    return true;
}

bool    fCutMultiplicity    ( Double_t  dMultiplicity )
{
    return true;
    if ( dMultiplicity < fMinMult ) return false;
    if ( dMultiplicity > fMaxMult ) return false;
    return true;
}

bool    fAcceptCandidate ( Double_t  dRapidity, Double_t  dInvariantMass, Double_t dTransverseMom, Double_t  dMultiplicity )
{
    if ( !fCutRapidity(dRapidity)           ) return false;
    if ( !fCutInvariantMass(dInvariantMass) ) return false;
    if ( !fCutTransverseMom(dTransverseMom) ) return false;
    if ( !fCutMultiplicity(dMultiplicity)   ) return false;
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
        if ( fOption.Contains("Multiplicity",TString::kIgnoreCase) )    { kDoMultiplicity = true;   cout << "[INFO] Mutliplicity option chosen" <<endl;}
        if ( fOption.Contains("Yield",TString::kIgnoreCase) )           { kDoYield = true;          cout << "[INFO] Yield option chosen" <<endl;}
        if ( fOption.Contains("Trigger",TString::kIgnoreCase) )         { kDoTrigger = true;        cout << "[INFO] Trigger option chosen" <<endl;}
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
    fLevyFit1D  ->  SetParLimits(0,0.,10.);
    fLevyFit1D  ->  SetParameter(0,kPMas);  //Particle Mass?
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,2.1,100.);
    fLevyFit1D  ->  SetParameter(1,6.7);    // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.21,10.);
    fLevyFit1D  ->  SetParameter(2,.272);   // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,1.e-9,1.);
    fLevyFit1D  ->  SetParameter(3,fIntegral);  //
    
    if ( fIsConditional )   {
        // - // Setting up Fit parameters
        
        // Mass
        fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
        fLevyFit1D  ->  SetParameter(0,kPMas);
        
        // n-Parameter
        fLevyFit1D  ->  SetParLimits(1,3.5,7.5);
        fLevyFit1D  ->  SetParameter(1,4.5); // 6.7
        
        // T-Parameter
        fLevyFit1D  ->  SetParLimits(2,.18,.45);
        fLevyFit1D  ->  SetParameter(2,.272); // .272
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
        fResult ->  SetPointEYhigh  ( iPnt, fYValue*sqrt(kBRSystematics*kBRSystematics) );
        fResult ->  SetPointEYlow   ( iPnt, fYValue*sqrt(kBRSystematics*kBRSystematics) );
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
            fCurrentGraph ->  SetPointEYhigh  ( iPnt, fYValue*sqrt(4*kBRSystematics*kBRSystematics) );
            fCurrentGraph ->  SetPointEYlow   ( iPnt, fYValue*sqrt(4*kBRSystematics*kBRSystematics) );
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
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err, Syst err
    Double_t   *fResult     =   new Double_t    [10];
    //
    //  Measuring mean value
    fSetLevyTsalis(fIsConditional,fIntegral);
    TGraphAsymmErrors   *   gTotal  =   new TGraphAsymmErrors(*(fSumGraphErrors(gStatistics,gSystematics)));
    gTotal->  Fit(fLevyFit1D,"IMRE0SEX0","",0.4,10.);
    fResult[0]  =   fLevyFit1D  ->Integral(0.,0.4);
    //
    TCanvas *   cDrawFit    =   new TCanvas(Form("gTotal_%s",fName.Data()),Form("gTotal_%s",fName.Data()));
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gTotal->Draw();
    fLevyFit1D->Draw("same");
    cDrawFit->Write();
    cDrawFit->SaveAs(Form("result/tmp/gTotal_%s.pdf",fName.Data()));
    delete cDrawFit;
    //
    //  Measuring Statistical error
    TH1D*   hStatIntegral   =   new TH1D(Form("hStatIntegral_%s",fName.Data()),"hStatIntegral",100000,0.0,.1);
    for ( Int_t iFit = 0; iFit < kStatEvalCycles; iFit++ )  {
        //  Set Standard Fit
        fSetLevyTsalis(fIsConditional,fIntegral);
        //
        //  Generating the Fit TGraph
        auto fSubject   =   fRandomizePoints(gStatistics,gSystematics);
        //
        fSubject    ->  Fit(fLevyFit1D,"IMREQ0SEX0","",0.4,10.);
        //
        hStatIntegral->Fill(fLevyFit1D  ->Integral(0.,0.4));
    }
    fResult[1]  =   hStatIntegral->GetRMS();
    //
    //  Measuring Systematics error
    TH1D*   hSystIntegral   =   new TH1D(Form("hSystIntegral_%s",fName.Data()),"hSystIntegral",100000,0.0,.1);
    for ( Int_t iFit = 0; iFit < kStatEvalCycles; iFit++ )  {
        //  Set Standard Fit
        fSetLevyTsalis(fIsConditional,fIntegral);
        //
        //  Generating the Fit TGraph
        auto fSubject   =   fRandomizePoints(gSystematics,gStatistics);
        //
        fSubject    ->  Fit(fLevyFit1D,"IMREQ0SEX0","",0.4,10.);
        //
        hSystIntegral->Fill(fLevyFit1D  ->Integral(0.,0.4));
    }
    fResult[2]  =   hSystIntegral->GetRMS();
    //
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
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err, Syst err
    Double_t   *fResult = new   Double_t    [10];
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
        fResult[1]             +=   fXBinWidth*fYErrStatLow;
        fResult[2]             +=   fXBinWidth*fYErrStatHigh;
        fResult[3]             +=   fXBinWidth*fYErrSystLow;
        fResult[4]             +=   fXBinWidth*fYErrSystHigh;
        /*
        fResult[5]             +=   fXBinWidth*fYValue;
        fResult[6]             +=   fXBinWidth*fYErrStatLow;
        fResult[7]             +=   fXBinWidth*fYErrStatHigh;
        fResult[8]             +=   fXBinWidth*fYErrSystLow;
        fResult[9]             +=   fXBinWidth*fYErrSystHigh;
        */
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
    // Result format: Integral, Stat err, Syst err, Mean pT, Stat err, Syst err
    Double_t   *fResult             =   new Double_t        [6];
    //
    bool fIsConditional = false;
    if ( fName.First("2D") != -1  ) fIsConditional = true;
    //
    auto        fIntegralResults    =   fIntegrateModel     (gStatistics,gSystematics,fName);
    auto        fExtrapolResults    =   fExtrapolateModel   (fIsConditional,gStatistics,gSystematics,fIntegralResults[0],fName);
    //
    fResult[0]  =   fIntegralResults[0] +   fExtrapolResults[0];
    fResult[1]  =   fIntegralResults[1] +   fExtrapolResults[1];
    fResult[2]  =   fIntegralResults[2] +   fExtrapolResults[2];
    //
    // !TODO: Revise the combination method
    fResult[3]  =   fIntegralResults[3] +   fExtrapolResults[3];
    fResult[4]  =   fIntegralResults[4] +   fExtrapolResults[4];
    fResult[5]  =   fIntegralResults[5] +   fExtrapolResults[5];
    //
    // End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
//
//------------------------------//
//    ANALYSISI LEGACY Fncs     //
//------------------------------//

TH1F*            fCheckPublishedResults( TH1F* fMyResults, TH1F* hPublishedResults, TGraphAsymmErrors* gPublishedResults )    {
    TH1F   *fCheck  =   new TH1F(*hPublishedResults);
    fCheck->Divide(fMyResults,hPublishedResults);
    for ( int i = 0; i < fCheck->GetNbinsX(); i++ ) {
        auto pubError   =   gPublishedResults->GetErrorXhigh(i);
        auto myError    =   fMyResults->GetBinError(i+1);
        fCheck->SetBinError (i,sqrt(pubError*pubError+myError*myError));
    }
    return fCheck;
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

void            fCoreFitModelSetBoundaries      ( string fOption, Double_t &aValMin, Double_t &aValMax )
{
    aValMin = 0.99;
    aValMax = 1.05;
    if ( fOption.find("RA") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.040;
    }
    if ( fOption.find("RB") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.060;
    }
    if ( fOption.find("RC") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.070;
    }
    if ( fOption.find("RD") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.040;
    }
    if ( fOption.find("RE") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.050;
    }
    if ( fOption.find("RF") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.060;
    }
    if ( fOption.find("RG") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.070;
    }
    if ( fOption.find("RH") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.040;
    }
    if ( fOption.find("RI") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.050;
    }
    if ( fOption.find("RJ") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.060;
    }
    if ( fOption.find("RK") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.070;
    }
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

void            SetBoundaries   ( string fOption, Double_t &aValMin, Double_t &aValMax )
{
    aValMin = 0.99;
    aValMax = 1.05;
    if ( fOption.find("RA") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.040;
    }
    if ( fOption.find("RB") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.060;
    }
    if ( fOption.find("RC") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.070;
    }
    if ( fOption.find("RD") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.040;
    }
    if ( fOption.find("RE") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.050;
    }
    if ( fOption.find("RF") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.060;
    }
    if ( fOption.find("RG") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.070;
    }
    if ( fOption.find("RH") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.040;
    }
    if ( fOption.find("RI") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.050;
    }
    if ( fOption.find("RJ") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.060;
    }
    if ( fOption.find("RK") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.070;
    }
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
    RooRealVar ch0      = RooRealVar        ("ch0","ch0"      ,1,   -1, 1);//,0.5,-1,1);
    RooRealVar ch1      = RooRealVar        ("ch1","ch1"      ,0.1, -1, 1);//,-0.1,-1,1);
    RooRealVar ch2      = RooRealVar        ("ch2","ch2"      ,0.01,-1, 1);//,0.01,-1,1);
    RooRealVar ch3      = RooRealVar        ("ch3","ch3"      ,0.05,-1, 1);//,-0.05,-1,1);
    
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
    RooVoigtian     fSig= RooVoigtian      ("fSig","fSig"      ,InvMass,sMass,sWidt,sSlop);
    RooChebychev    fBkg= RooChebychev     ("fBkg","fBkg"      ,InvMass,RooArgSet(ch0,ch1,ch2,ch3,ch4,ch5));
    RooAddPdf       fMod= RooAddPdf        ("fMod","fMod"      ,RooArgList(fBkg,fSig),RooArgList(nBB,nSS));
    
    
    RooFitResult* fFitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(Signal));
        if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
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
    Int_t kNCycle       = 1;
    
    // Background PDF Coefficients
    RooRealVar ch0      = RooRealVar        ("ch0","ch0"      ,0.5,-1,1);//,0.5,-1,1);
    RooRealVar ch1      = RooRealVar        ("ch1","ch1"      ,-0.1,-1,1);//,-0.1,-1,1);
    RooRealVar ch2      = RooRealVar        ("ch2","ch2"      ,0.,-1,1);//,0.01,-1,1);
    RooRealVar ch3      = RooRealVar        ("ch3","ch3"      ,0.,-1,1);//,-0.05,-1,1);
    
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
    RooVoigtian     fSig= RooVoigtian      ("fSig","fSig"      ,InvMass,sMass,sWidt,sSlop);
    RooChebychev    fBkg= RooChebychev     ("fBkg","fBkg"      ,InvMass,RooArgSet(ch0,ch1,ch2,ch3,ch4,ch5));
    RooAddPdf       fMod= RooAddPdf        ("fMod","fMod"      ,RooArgList(fBkg,fSig),RooArgList(nBB,nSS));
    
    RooFitResult* fFitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(Signal));
        if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
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
        fSaveToCanvas               ->SaveAs(Form("result/SEFitCheck/PT_%.1f_%.1f_1D_%s.pdf",fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()));
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
    Int_t kNCycle       = 1;
    
    RooArgSet  *utilityx    =   new RooArgSet(fFitShapeX->floatParsFinal(),fFitShapeX->constPars());
    RooArgSet  *utilityy    =   new RooArgSet(fFitShapeY->floatParsFinal(),fFitShapeY->constPars());
    
    
    // Background
    RooRealVar ch0x, ch1x, ch2x, ch3x, ch4x, ch5x, ch0y, ch1y, ch2y, ch3y, ch4y, ch5y;
    
    if ( fBackg )
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0),utilityx->getRealValue("ch0",0)*0.75,utilityx->getRealValue("ch0",0)*1.25);
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0),utilityx->getRealValue("ch1",0)*0.75,utilityx->getRealValue("ch1",0)*1.25);
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0),utilityx->getRealValue("ch2",0)*0.75,utilityx->getRealValue("ch2",0)*1.25);
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0),utilityx->getRealValue("ch3",0)*0.75,utilityx->getRealValue("ch3",0)*1.25);
        ch4x     = RooRealVar ("ch4x","ch4x"     ,utilityx->getRealValue("ch4",0),utilityx->getRealValue("ch4",0)*0.75,utilityx->getRealValue("ch4",0)*1.25);
        ch5x     = RooRealVar ("ch5x","ch5x"     ,utilityx->getRealValue("ch5",0),utilityx->getRealValue("ch5",0)*0.75,utilityx->getRealValue("ch5",0)*1.25);
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0),utilityy->getRealValue("ch0",0)*0.75,utilityy->getRealValue("ch0",0)*1.25);
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0),utilityy->getRealValue("ch1",0)*0.75,utilityy->getRealValue("ch1",0)*1.25);
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0),utilityy->getRealValue("ch2",0)*0.75,utilityy->getRealValue("ch2",0)*1.25);
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0),utilityy->getRealValue("ch3",0)*0.75,utilityy->getRealValue("ch3",0)*1.25);
        ch4y     = RooRealVar ("ch4y","ch4y"     ,utilityy->getRealValue("ch4",0),utilityy->getRealValue("ch4",0)*0.75,utilityy->getRealValue("ch4",0)*1.25);
        ch5y     = RooRealVar ("ch5y","ch5y"     ,utilityy->getRealValue("ch5",0),utilityy->getRealValue("ch5",0)*0.75,utilityy->getRealValue("ch5",0)*1.25);
    }
    else
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0));
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0));
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0));
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0));
        ch4x     = RooRealVar ("ch4x","ch4x"     ,utilityx->getRealValue("ch4",0));
        ch5x     = RooRealVar ("ch5x","ch5x"     ,utilityx->getRealValue("ch5",0));
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0));
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0));
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0));
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0));
        ch4y     = RooRealVar ("ch4y","ch4y"     ,utilityy->getRealValue("ch4",0));
        ch5y     = RooRealVar ("ch5y","ch5y"     ,utilityy->getRealValue("ch5",0));
    }
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,utilityx->getRealValue("sMass",0));
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,utilityx->getRealValue("sWidt",0));
    RooRealVar pSlopex  = RooRealVar ("pSlopx","pSlopx" ,utilityx->getRealValue("sSlop",0));
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,utilityy->getRealValue("sMass",0));
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,utilityy->getRealValue("sWidt",0));
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
    RooChebychev        fBkgx ("fBkgx","fBkgx"          ,varx,RooArgSet(ch0x,ch1x,ch2x,ch3x,ch4x,ch5x));
    RooVoigtian         fSigx ("fSigx","fSigx"          ,varx,pMassx,pWidthx,pSlopex);
    RooChebychev        fBkgy ("fBkgy","fBkgy"          ,vary,RooArgSet(ch0y,ch1y,ch2y,ch3y,ch4y,ch5x));
    RooVoigtian         fSigy ("fSigy","fSigy"          ,vary,pMassy,pWidthy,pSlopey);
    RooProdPdf          fBB   ("fBkg","fBkg"            ,fBkgx,fBkgy);
    RooProdPdf          fSB   ("fSigBkg","fSBWBkg"      ,fSigx,fBkgy);
    RooProdPdf          fBS   ("fBkgSig","fBkgSig"      ,fBkgx,fSigy);
    RooProdPdf          fSS   ("fSigSig","fSigSig"      ,fSigx,fSigy);
    RooAddPdf           fMod  ("fMod2D","fMod2D"        ,RooArgList(fBB,fSS,fSB,fBS),RooArgList(n1,n0,n3,n2));
    
    RooFitResult* fFitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        fFitResults      =   fMod.fitTo(*dataLoose,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save(),NumCPU(kCPU_use));
        auto N_Raw      =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().at(Signal));
        if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
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
Double_t*           fExtrapolateModel               ( Tclass *THdata, TString fName = "ExtrapolateSignal" ) {
    // Optimisation mode
    gROOT->SetBatch(true);
    
    // Result format: Integral, Stat err, Syst err, Mean pT, Stat err, Syst err
    Double_t   *fResult = new   Double_t    [6];
    
    // Setting -1. for the default
    for ( Int_t iFill = 0; iFill < 6; iFill++ ) fResult[iFill]  =   -1.;
    
    // Set Standard Fit
    fSetLevyTsalis();
    
    // Fit the Spectra
    THdata->Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
    
    // Save to further checks
    TCanvas * fCheckFit = new TCanvas();
    gPad->SetLogy();
    THdata      ->Draw("same");
    fLevyFit1D  ->Draw("same");
    fCheckFit   ->Write();
    fCheckFit   ->SaveAs(Form("tmp/%s.pdf",fName.Data()));
    
    fResult[0]  =   fLevyFit1D->Integral(0.,0.4);
    fResult[1]  =   fLevyFit1D->IntegralError(0.,0.4);
    
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

#endif
