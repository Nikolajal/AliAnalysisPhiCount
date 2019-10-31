#ifndef SETFUNCTIONS_H
#define SETFUNCTIONS_H
#include "SetValues.h"

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

typedef struct
{
    RooFitResult    * FitRes;
    RooAddPdf       * fMod;
} FIT1D_RESULT;

typedef struct
{
    FIT1D_RESULT    Array[1024];
} FIT1D_RESULT_ARRAY;

enum fitresults1D
{
    ich0,
    ich1,
    ich2,
    ich3,
    ich4,
    inBkg,
    inSig,
    ipMass,
    ipWidth
};
enum fitresults2D
{
    inBB,
    inBS,
    inBkg2D,
    inSB,
    inSS,
};

void ModelFit1D (RooDataHist * data, RooRealVar var, FIT1D_RESULT & result, Int_t index)
{
    /* Global variables */
    auto nEntries = data->sumEntries();
    
    /* Background */
    RooRealVar * ch0    = new RooRealVar ("ch0","ch0"       ,0.,-1.,1.);
    RooRealVar * ch1    = new RooRealVar ("ch1","ch1"       ,0.,-1.,1.);
    RooRealVar * ch2    = new RooRealVar ("ch2","ch2"       ,0.,-1.,1.);
    RooRealVar * ch3    = new RooRealVar ("ch3","ch3"       ,0.,-1.,1.);
    RooRealVar * ch4    = new RooRealVar ("ch4","ch4"       ,0.,-1.,1.);
    RooRealVar * nBkg   = new RooRealVar ("nBkg","nBkg"     ,0.5*nEntries,0.,nEntries);
    // PDF
    RooChebychev * fBkg = new RooChebychev("fBkg","fBkg"    ,var,RooArgSet(*ch0,*ch1,*ch2,*ch3,*ch4));
    
    /* Signal */
    RooRealVar * pMass  = new RooRealVar ("pMass","pMass"   ,1.020,1.010,1.030);
    RooRealVar * pWidth = new RooRealVar ("pWidth","pWidth" ,0.005,0.0005,0.01);
    RooRealVar * nSig   = new RooRealVar ("nSig","nSig"     ,0.5*nEntries,0.,nEntries);
    // PDF
    RooBreitWigner * fSig = new RooBreitWigner ("fSig","fSig",var,*pMass,*pWidth);
    
    /* Model */
    result.fMod         = new RooAddPdf ("fMod","fMod",RooArgList(*fBkg,*fSig),RooArgList(*nBkg,*nSig));
    result.FitRes       = result.fMod->fitTo(*data,Save());
    
    /* Save to File */
    TH1* Mod1D = result.fMod->createHistogram(Form("Mod1D_%f_%f",fBound_pT(index),fBound_pT(index+1)),var,Binning(20));
    Mod1D->Write();
}

void ModelFit2D_Preprocess (RooDataHist ** data, RooRealVar var, FIT1D_RESULT_ARRAY & result)
{
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        ModelFit1D(data[iHisto],var,result.Array[iHisto],iHisto);
    }
}

void ModelFit2D (RooDataHist * data, RooRealVar varx, RooRealVar vary, FIT1D_RESULT & result, FIT1D_RESULT_ARRAY & preprocess, Int_t index, Int_t jndex)
{
    /* Global variables */
    RooRealVar * nBB    = new RooRealVar ("nBB","nBB"           ,0.5,0.,20000);
    RooRealVar * nSB    = new RooRealVar ("nSB","nSB"           ,0.5,0.,20000);
    RooRealVar * nBS    = new RooRealVar ("nBS","nBS"           ,0.5,0.,20000);
    RooRealVar * nSS    = new RooRealVar ("nSS","nSS"           ,0.5,0.,20000);
    RooRealVar * nBkg   = new RooRealVar ("nBkg","nBkg"         ,0.5,0.,20000);
    
    
    /* Background */
    // XBkg
    RooRealVar * ch0x   = new RooRealVar ("ch0x","ch0x"         ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ich0))->getVal()));
    RooRealVar * ch1x   = new RooRealVar ("ch1x","ch1x"         ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ich1))->getVal()));
    RooRealVar * ch2x   = new RooRealVar ("ch2x","ch2x"         ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ich2))->getVal()));
    RooRealVar * ch3x   = new RooRealVar ("ch3x","ch3x"         ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ich3))->getVal()));
    RooRealVar * ch4x   = new RooRealVar ("ch4x","ch4x"         ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ich4))->getVal()));
    // YBkg
    RooRealVar * ch0y   = new RooRealVar ("ch0y","ch0y"         ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ich0))->getVal()));
    RooRealVar * ch1y   = new RooRealVar ("ch1y","ch1y"         ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ich1))->getVal()));
    RooRealVar * ch2y   = new RooRealVar ("ch2y","ch2y"         ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ich2))->getVal()));
    RooRealVar * ch3y   = new RooRealVar ("ch3y","ch3y"         ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ich3))->getVal()));
    RooRealVar * ch4y   = new RooRealVar ("ch4y","ch4y"         ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ich4))->getVal()));
    
    // PDF
    RooChebychev * fBkgx= new RooChebychev("fBkgx","fBkgx"      ,varx,RooArgSet(*ch0x,*ch1x,*ch2x,*ch3x,*ch4x));
    RooChebychev * fBkgy= new RooChebychev("fBkgy","fBkgy"      ,vary,RooArgSet(*ch0y,*ch1y,*ch2y,*ch3y,*ch4y));
    
    /* Signal */
    RooRealVar * pMassx = new RooRealVar ("pMassx","pMassx"     ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ipMass))->getVal()));
    RooRealVar * pWidthx= new RooRealVar ("pWidthx","pWidthx"   ,(static_cast<RooRealVar*>(preprocess.Array[index].FitRes->floatParsFinal().at(ipWidth))->getVal()));
    RooRealVar * pMassy = new RooRealVar ("pMassy","pMassy"     ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ipMass))->getVal()));
    RooRealVar * pWidthy= new RooRealVar ("pWidthy","pWidthy"   ,(static_cast<RooRealVar*>(preprocess.Array[jndex].FitRes->floatParsFinal().at(ipWidth))->getVal()));
                        
    // PDF
    RooBreitWigner * fSigx= new RooBreitWigner ("fSigx","fSigx" ,varx,*pMassx,*pWidthx);
    RooBreitWigner * fSigy= new RooBreitWigner ("fSigy","fSigy" ,vary,*pMassy,*pWidthy);

    /* Model */
    
    //Backgrounds
    RooProdPdf* fBB     = new RooProdPdf("fBB","fBB"            ,*fBkgx,*fBkgy);
    RooProdPdf* fSB     = new RooProdPdf("fSB","fSB"            ,*fSigx,*fBkgy);
    RooProdPdf* fBS     = new RooProdPdf("fBS","fBS"            ,*fBkgx,*fSigy);
    RooAddPdf * fBkg    = new RooAddPdf ("fBkg","fBkg"          ,RooArgList(*fBB,*fSB,*fBS),RooArgList(*nBB,*nSB,*nBS));
    
    //Signal
    RooProdPdf* fSS     = new RooProdPdf("fBB","fBB"            ,*fBkgx,*fBkgy);
    
    //Total
    result.fMod         = new RooAddPdf ("fMod","fMod"          ,RooArgList(*fBkg,*fSS),RooArgList(*nBkg,*nSS));
    result.FitRes       = result.fMod->fitTo(*data,Save());
    
    /* Save to File */
    TH1* Mod2D = result.fMod->createHistogram(Form("Mod1D_%f_%f_%f_%f",fBound2D_pT(index),fBound2D_pT(index+1),fBound2D_pT(jndex),fBound2D_pT(jndex+1)),varx,Binning(20),YVar(vary,Binning(20)));
    Mod2D->Write();
}


#endif
