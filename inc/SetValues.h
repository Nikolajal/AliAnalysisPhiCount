#ifndef SETVALUES_H
#define SETVALUES_H

// C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>

// ROOT
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"

// RooFit
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"

// RooFitFunction
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"

using namespace std;
using namespace RooFit;
using namespace std::chrono;

auto const oFileMonCar  = "./result/outGeneratorMC.root";
auto const oFilePreP1D  = "./result/outPreProcessing1D.root";
auto const oFilePreP2D  = "./result/outPreProcessing2D.root";
auto const oFileDataFt  = "./result/outDataFormat.root";
auto const oFileAnalys  = "./result/outAnalysis.root";
auto const oFileEffici  = "./result/outEfficiency.root";
auto const oFileHistog  = "./result/outHistogram.root";
auto const iFileMCEffi  = "./result/outGeneratorMC_Efficiency.root";
auto const hPhiEff      = "hPhiEff";
auto const PTreeNameK2  = "PythiaTreeK2";
auto const PTreeNamePhi = "PythiaTreePhi";
auto const PEvents      = 1e7;

/*Change_____________*/
// InvMass range
auto const nBins            = 100;       //nBinIM0D
auto const minBound         = 0.99;
auto const maxBound         = 1.09;

// pT cuts
const Int_t     nBin_pT     = 40;
const Float_t   nMin_pT     = 0.;
const Float_t   nMax_pT     = 4.;
//___________________

// InvMass range
const Int_t     nBinIM1D    = 90;
const Float_t   fMinIM1D    = 0.99;
const Float_t   fMaxIM1D    = 1.09;

// InvMass range
const Int_t     nBinIM2D    = 60;
const Float_t   fMinIM2D    = 0.99;
const Float_t   fMaxIM2D    = 1.09;

// pT cuts
const Int_t     nBinPT1D    = 40;
const Float_t   fMinPT1D    = 0.;
const Float_t   fMaxPT1D    = 4.;

// pT cuts
const Int_t     nBinPT2D    = 4;
const Float_t   fMinPT2D    = 0.;
const Float_t   fMaxPT2D    = 4.;

// pT cuts
const Int_t     nBinPhi1D   = 40;
const Float_t   fMinPhi1D   = 0.;
const Float_t   fMaxPhi1D   = 4.;

// pT cuts
const Int_t     nBinPhi2D   = 4;
const Float_t   fMinPhi2D   = 0.;
const Float_t   fMaxPhi2D   = 4.;

typedef struct
{
    Int_t nKaon, particleID[1024], mother1[1024], mother2[1024], motherID[1024];
    Float_t px[1024], py[1024], pz[1024], pT[1024], e[1024];
} EVKAON;

typedef struct
{
    Int_t   nKaonCouple, iKaon[1024], jKaon[1024];
    Bool_t  bPhi[1024], bRec[1024];
    Float_t InvMass[1024], pT[1024];
} EVKAONCOUPLE;

typedef struct
{
    Int_t   nPhi;
    Float_t pT[1024];
} EVPHI;

Float_t fBound_pT (Int_t index)
{
    return (index)*(nMax_pT - nMin_pT)/(static_cast<Float_t>(nBin_pT));
}

Float_t fBound2D_pT (Int_t index)
{
    return (index)*(fMaxPT2D - fMinPT2D)/(static_cast<Float_t>(nBinPT2D));
}

#endif
