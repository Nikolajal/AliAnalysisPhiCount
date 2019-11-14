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

auto const  oFileMonCar = "./result/outGeneratorMC.root";
auto const  oFilePreP1D = "./result/outPreProcessing1D.root";
auto const  oFilePreP2D = "./result/outPreProcessing2D.root";
auto const  oFilePreBKG = "./result/outPreProcessingBKG.root";
auto const  oFilePrBKG2 = "./result/outPreProcessingBKGModel.root";
auto const  oFileDataFt = "./result/outDataFormat.root";
auto const  oFileAnalys = "./result/outAnalysis.root";
auto const  oFileAnal1D = "./result/outAnalysis1D.root";
auto const  oFileAnal2D = "./result/outAnalysis2D.root";
auto const  oFileEffici = "./result/outEfficiency.root";
auto const  oFileHist1D = "./result/outHistogram1D.root";
auto const  oFileHist2D = "./result/outHistogram2D.root";
auto const  iFileMCEffi = "./result/outGeneratorMC_Efficiency.root";
auto const  hPhiEff     = "hPhiEff";
auto const  PTreeNameK2 = "PythiaTreeK2";
auto const  PTreeNameKS = "PythiaTreeKS";
auto const  PTreeNameKD = "PythiaTreeKD";
auto const  PTreeNamePhi= "PythiaTreePhi";
int  const  PEvents     = 1e5;
bool        BKG2        = false;

// InvMass range 1D
const   Int_t     nBinIM1D  = 100;
const   Float_t   fMinIM1D  = 0.99;
const   Float_t   fMaxIM1D  = 1.09;
        Float_t * fArrIM1D  = new Float_t [nBinIM1D+1];

// InvMass range 2D
const   Int_t     nBinIM2D  = 100;
const   Float_t   fMinIM2D  = 0.99;
const   Float_t   fMaxIM2D  = 1.09;
        Float_t * fArrIM2D  = new Float_t [nBinIM2D+1];

// pT cuts 1D
        Int_t     nBinPT1D  = 1;
const   Float_t   fMinPT1D  = 0.;
const   Float_t   fMaxPT1D  = 4.;
        Float_t * fArrPT1D  = new Float_t [nBinPT1D+1];

// pT cuts 2D
        Int_t     nBinPT2D  = 1;
const   Float_t   fMinPT2D  = 0.;
const   Float_t   fMaxPT2D  = 4.;
        Float_t * fArrPT2D  = new Float_t [nBinPT2D+1];

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

void    vSetBinsIM1D ()
{
    for (int i = 0; i <= nBinIM1D; i++ )
    {
        fArrIM1D[i] = fMinIM1D+(i)*(fMaxIM1D - fMinIM1D)/(static_cast<Float_t>(nBinIM1D));
    }
}

void    vSetBinsIM2D ()
{
    for (int i = 0; i <= nBinIM2D; i++ )
    {
        fArrIM2D[i] = fMinIM2D+(i)*(fMaxIM2D - fMinIM2D)/(static_cast<Float_t>(nBinIM2D));
    }
}

void    vSetBinsPT1D ()
{
    /*
    fArrPT1D[0] =   0.;
    fArrPT1D[1] =   0.5;
    fArrPT1D[2] =   1.;
    fArrPT1D[3] =   0.75;
    fArrPT1D[4] =   1.;
    fArrPT1D[5] =   1.25;
    fArrPT1D[6] =   1.5;
    fArrPT1D[7] =   2.;
    fArrPT1D[8] =   2.5;
    fArrPT1D[9] =   3.;
    fArrPT1D[10]=   4.;
    /*/
    for (int i = 0; i <= nBinPT1D; i++ )
    {
        fArrPT1D[i] = fMinPT1D+(i)*(fMaxPT1D - fMinPT1D)/(static_cast<Float_t>(nBinPT1D));
    }
    
}

void    vSetBinsPT2D ()
{
    /*
    fArrPT2D[0] =   0.;
    fArrPT2D[1] =   0.25;
    fArrPT2D[2] =   0.5;
    fArrPT2D[3] =   0.75;
    fArrPT2D[4] =   1.;
    fArrPT2D[5] =   1.25;
    fArrPT2D[6] =   1.5;
    fArrPT2D[7] =   2.;
    fArrPT2D[8] =   2.5;
    fArrPT2D[9] =   3.;
    fArrPT2D[10]=   4.;
    /*/
    for (int i = 0; i <= nBinPT2D; i++ )
    {
        fArrPT2D[i] = fMinPT2D+(i)*(fMaxPT2D - fMinPT2D)/(static_cast<Float_t>(nBinPT2D));
    }
    
}

#endif
