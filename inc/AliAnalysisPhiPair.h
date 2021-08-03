// Global file
//
// !TODO: 1.
//
#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H
//
//>>    Analysis Utility
#include "../AliAnalysisUtility/AliAnalysisUtility.h"
//
//>>    Analysis constants and binning utilities
#include "./AliAnalysisPhiPair_Constants.h"
//
//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
//>>    Performance Regulation Values
auto const  kCPU_use                =   1;
auto const  kNCycle_                =   3;
auto const  kStatEvalCycles         =   4;
auto const  kPrintIntervalPP        =   1000000;
auto const  k2DErrorLimit           =   10.;
auto const  kCPUStrategy            =   1;
auto const  kFitOffset              =   true;
auto const  kFitInitHesse           =   true;
auto const  kFitMinos               =   true;
auto const  kFitMinuitStrategy      =   2;
auto const  kSaveToFile             =   true;
//
//_____________________________________________________________________________
//
//>>    Analysis Values
auto const  bPythiaTest             =   kFALSE;
auto const  kOnlyTrue               =   false;
//
//_____________________________________________________________________________
//
//*****                 !TODO: CLEAN
std::vector<std::pair<TF1*,std::vector<float>>> fSystFitFunctions  =   {  {fMTExponential,{1.2,1.4,1.6,2.0}}, {fBoseEinstein,{1.2,1.4,1.6,2.0}}, {fBoltzmann,{1.2,1.4,1.6,2.0}}, {fPowerLaw,{2.0,2.8,4.0}}, {fLevyTsallis,{1.2,1.4,1.6,2.0,2.8,4.0}} };//, {fBGBlastWave,{1.2,1.4,1.6,2.0,2.8,4.0}} };
//
//*****

Float_t             fMeasureMeanPT                  ( TF1 * fLowFit, TGraphAsymmErrors * gTotal, Bool_t fReFit = false );
Float_t             fMeasureMeanPT                  ( TF1 * fLowFit, TH1F * gTotal, Bool_t fReFit = false );
//
//-------------------------------------//
//      Analysis Fit Functions         //
//-------------------------------------//
//
//  --  --  General Analysis Functions  --  --  //
//
//_____________________________________________________________________________
//
void                        SetBoundaries                   ( TString fOption, Double_t &aValMin, Double_t &aValMax )   {
    aValMin = 0.998;
    aValMax = 1.065;
    //----
    if ( fOption.Contains("RA") )
    {
        aValMin =   0.996;
        aValMax =   1.015;
    }
    if ( fOption.Contains("RB") )
    {
        aValMin =   0.996;
        aValMax =   1.062;
    }
    if ( fOption.Contains("RC") )
    {
        aValMin =   0.996;
        aValMax =   1.065;
    }
    if ( fOption.Contains("RD") )
    {
        aValMin =   0.996;
        aValMax =   1.068;
    }
    if ( fOption.Contains("RE") )
    {
        aValMin =   0.996;
        aValMax =   1.071;
    }
    //----
    if ( fOption.Contains("RF") )
    {
        aValMin =   0.998;
        aValMax =   1.059;
    }
    if ( fOption.Contains("RG") )
    {
        aValMin =   0.998;
        aValMax =   1.062;
    }
    if ( fOption.Contains("RH") )
    {
        aValMin =   0.998;
        aValMax =   1.068;
    }
    if ( fOption.Contains("RI") )
    {
        aValMin =   0.998;
        aValMax =   1.071;
    }
    //----
    if ( fOption.Contains("RJ") )
    {
        aValMin =   1.000;
        aValMax =   1.059;
    }
    if ( fOption.Contains("RK") )
    {
        aValMin =   1.000;
        aValMax =   1.062;
    }
    if ( fOption.Contains("RL") )
    {
        aValMin =   1.000;
        aValMax =   1.065;
    }
    if ( fOption.Contains("RM") )
    {
        aValMin =   1.000;
        aValMax =   1.068;
    }
    if ( fOption.Contains("RN") )
    {
        aValMin =   1.000;
        aValMax =   1.071;
    }
}
//
//_____________________________________________________________________________
//
//
//  --  --  General Customisation Functions  --  --  //
//
//_____________________________________________________________________________
//
int                         fLegendSelect                   ( string fOption )  {
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("RapTru") )       return 1;
    if ( !fOption.compare("Rap") )          return 3;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 2;
    else return -1;
}
//
//_____________________________________________________________________________
//
void                        fLegendMaker                    ( RooPlot * fRooPlot, const char * fSelect, TLegend * fLegend )
{
    fLegend                     ->SetFillColorAlpha(kWhite,0.);
    fLegend                     ->SetLineColorAlpha(kWhite,0.);
    switch (fLegendSelect(fSelect))
    {
        case 1:
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        case 2:
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBS"),   "Fit (Bkg #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSB"),   "Fit (Sig #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        case 3:
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
//
//_____________________________________________________________________________
//
int                         fAxisSelect                     ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 3;
    if ( !fOption.compare("Rap") )          return 4;
    if ( !fOption.compare("RapTru") )       return 4;
    else return -1;
}
//
//_____________________________________________________________________________
//
void                        fAxisMaker                      ( RooPlot * fRooPlot, const char * fSelect )
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
//
//_____________________________________________________________________________
//
int                         fPlotterSelect                  ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )    return 1;
    if ( !fOption.compare("RapTru") )       return 1;
    if ( !fOption.compare("Rap") )          return 3;
    if ( !fOption.compare("xInvMass2D") )   return 2;
    if ( !fOption.compare("yInvMass2D") )   return 2;
    else return -1;
}
//
//_____________________________________________________________________________
//
void                        fRooPlotPlotter                 ( RooPlot * fRooPlot, const char * fSelect, RooAddPdf fModel , RooDataHist * fData )
{
    switch (fPlotterSelect(fSelect))
    {
        case 1:
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(4),                   LineStyle(kSolid), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(2),       Name("RooSS"));
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            break;
        case 2:
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(38),    MarkerStyle(26) ,   Name("RooData"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  LineColor(4),       LineStyle(kSolid), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkg"), LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigSig"),LineColor(2),       Name("RooSS"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigBkg"),LineStyle(kDashed), LineColor(33),    Name("RooSB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkgSig"),LineStyle(kDashed), LineColor(36),    Name("RooBS"));
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(38),    MarkerStyle(26) ,   Name("RooData"));
            break;
        case 3:
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(4),                   LineStyle(kSolid), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fBk2"),             LineStyle(kDashed), LineColor(33),      Name("RooB2"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(2),       Name("RooSS"));
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}
//
//_____________________________________________________________________________
//
void                        fRooPlotMaker                   ( RooPlot * fRooPlot, TLegend * fLegend, RooAddPdf fModel , RooDataHist * fData, const char * fSelect )
{
    fRooPlotPlotter(fRooPlot,fSelect,fModel,fData);
    fLegendMaker(fRooPlot,fSelect,fLegend);
    fAxisMaker(fRooPlot,fSelect);
}
//
//_____________________________________________________________________________
//
//
//  --  --  1D Analysis Functions  --  --  //
//
//_____________________________________________________________________________
//
RooFitResult               *fFitCoreModel                   ( TH1 * THdata, TH1F* hSlopReference, TString fName = "", TString fOption = "", Int_t PTindex = -1, Int_t PTDimension = 1, TString fPathToSave = "./result/SEFitCheck" )    {
    //
    //>>    Silencing TCanvas Pop-Up and speeding Fit procedures
    gROOT->SetBatch(kTRUE);
    //
    //>>    Selecting all Options
    Bool_t  fDegree2,   fDegree4,   fLosWidt,   fUsePoly,   fAlwaysFit, fSaveToFile,    fUseHighRes,    fUseLowRes,    fUseFreeRes;
    //
    fDegree2    =   fOption.Contains("DG2", TString::kIgnoreCase);
    fDegree4    =   fOption.Contains("DG4", TString::kIgnoreCase);
    fLosWidt    =   fOption.Contains("WDT", TString::kIgnoreCase);
    fUsePoly    =   fOption.Contains("POL", TString::kIgnoreCase);
    fUseHighRes =   fOption.Contains("RSH", TString::kIgnoreCase);
    fUseLowRes  =   fOption.Contains("RSL", TString::kIgnoreCase);
    fUseFreeRes =   fOption.Contains("RSX", TString::kIgnoreCase);
    fAlwaysFit  =   fOption.Contains("FIT", TString::kIgnoreCase);
    fSaveToFile =   fOption.Contains("SAVE",TString::kIgnoreCase);
    //
    //>>    Check there is a reasonable amount of entries form the general method, or if it has been overrided
    if ( !fIsWorthFitting( THdata ) || fAlwaysFit ) return nullptr;
    //
    //>>    Global Variables
    Double_t                        fInvMassValMax, fInvMassValMin;
    SetBoundaries(  fOption,        fInvMassValMin, fInvMassValMax  );
    Int_t           nEntries    =   THdata->GetEntries();
    RooRealVar      InvMass     =   RooRealVar        ("InvMass",   "InvMass",      fInvMassValMin, fInvMassValMax                  );
    RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         InvMass,        Import(*THdata)                 );
    RooDataHist*    dataLoose   =   new RooDataHist   ("DataLoose", "DataLoose",    InvMass,        Import(*fLooseErrors(THdata))   );
    Int_t           kNCycle     =   kNCycle_;
    Float_t         fRescaleRes =   1.;
    if ( fUseHighRes )  fRescaleRes = 1.1;
    if ( fUseLowRes  )  fRescaleRes = 0.9;
    //
    //>>    Background PDF Coefficients
    RooRealVar ch1, ch2, ch3, ch4;
                                ch1     =   RooRealVar      ("ch1",     "ch1",      0.5,    -1,         1   );
                                ch2     =   RooRealVar      ("ch2",     "ch2",      -.06,   -1,         1   );
    if ( fDegree2 && !fDegree4 )ch3     =   RooRealVar      ("ch3",     "ch3",      0.  );
    else                        ch3     =   RooRealVar      ("ch3",     "ch3",      0.01,   -1,         1   );
    if ( fDegree4 )             ch4     =   RooRealVar      ("ch4",     "ch4",      0.,     -1,         1   );
    else                        ch4     =   RooRealVar      ("ch4",     "ch4",      0.  );
    //
    //>>    Signal PDF Coefficients
    RooRealVar sMass, sWidt, sSlop;
                                sMass   =   RooRealVar      ("bMass",   "bMass",    kPhiMesonMass_,  kPhiMesonMass_*0.5,  kPhiMesonMass_*1.5);
    if ( fLosWidt )             sWidt   =   RooRealVar      ("bWidt",   "bWidt",    kPhiMesonWidth,  kPhiMesonWidth*0.5,  kPhiMesonWidth*1.5);
    else                        sWidt   =   RooRealVar      ("bWidt",   "bWidt",    kPhiMesonWidth);
    if ( bPythiaTest )          sSlop   =   RooRealVar      ("bSlop",   "bSlop",    0.);
    else if ( !fUseFreeRes )    sSlop   =   RooRealVar      ("bSlop",   "bSlop",    fRescaleRes*hSlopReference->GetBinContent(PTindex+1)/1000.);
    else                        sSlop   =   RooRealVar      ("bSlop",   "bSlop",    fRescaleRes*hSlopReference->GetBinContent(PTindex+1)/1000.);
    //
    //>>    Normalisation Coefficients
    RooRealVar  nSS,    nBB;
                                nSS     =   RooRealVar      ("anSS",    "anSS",     .5*nEntries,    0., nEntries);
                                nBB     =   RooRealVar      ("anBB",    "anBB",     .5*nEntries,    0., nEntries);
    //
    //>>    Building the PDFs
    RooVoigtian                 fSig    =   RooVoigtian     ("fSig",    "fSig",     InvMass,    sMass,  sWidt,  sSlop);
    RooChebychev                fBkg    =   RooChebychev    ("fBkg",    "fBkg",     InvMass,    RooArgSet(ch1,ch2,ch3,ch4));
    RooPolynomial               fBkgPol =   RooPolynomial   ("fBkg",    "fBkg",     InvMass,    RooArgSet(ch1,ch2,ch3,ch4));
    RooAddPdf                  *fMod;
    if ( fUsePoly )   {
        fMod                            =   new RooAddPdf   ("fMod",    "fMod",     RooArgList(fBkgPol,fSig),   RooArgList(nBB,nSS));
        ch1.removeRange();
        ch2.removeRange();
        ch3.removeRange();
        ch4.removeRange();
    }   else    {
        fMod                            =   new RooAddPdf   ("fMod",    "fMod",     RooArgList(fBkg,fSig),      RooArgList(nBB,nSS));
    }
    //
    RooFitResult* fFitResults;
    fFitResults      =   fMod->fitTo(*dataLoose,Save(),NumCPU(kCPU_use,kCPUStrategy));
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )    {
        if ( !kOnlyTrue )   {
            fFitResults      =   fMod->fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
            auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("anSS"));
            if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError()) ) break;
        }
        if ( kOnlyTrue )    {
            fFitResults      =   fSig.fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
        }
    }
    //
    if ( fSaveToFile || kSaveToFile )
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
        
        fRooPlotMaker(fSaveToFrame,fLegend,*fMod,data,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        if ( PTDimension == 1 )fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f_1D_%s.pdf",fPathToSave.Data(),fArrPT1D[PTindex],fArrPT1D[PTindex+1],fName.Data()));
        if ( PTDimension == 2 )fSaveToCanvas               ->SaveAs(Form("%s/PT_%.1f_%.1f_1D_%s.pdf",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fName.Data()));
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(kFALSE);
    
    return fFitResults;
}
//
//_____________________________________________________________________________
//
std::vector<TH1F*>          FitModel                        ( TH1F **hTarget, TH1F* hSlopReference, RooFitResult** &fFitresultsStore, Int_t fDimension, TString fOption = "", TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "" )  {
    std::vector<TH1F*>  fResults;
    if ( !fFitresultsStore )  {  fFitresultsStore = new RooFitResult*[((fDimension == 1)? nBinPT1D : nBinPT2D)]; }
    for ( Int_t iFit = 0; iFit < ((fDimension == 1)? nBinPT1D : nBinPT2D); iFit++ )
    {
        //>>    Fit to model
        fFitresultsStore[iFit]      =   fFitCoreModel(hTarget[iFit],hSlopReference,fNameFile,fOption,iFit,fDimension,fTargetPath);
        //
        //>>    SegFault Protection
        if ( !fFitresultsStore[iFit] ) continue;
        //
        //>>    Building Raw Count histograms
        if ( iFit == 0 )   {
            Int_t   iTer = 0;
            for ( auto fCoeff : fFitresultsStore[iFit]->floatParsFinal() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                if  ( fDimension == 1)  fResults.push_back( new TH1F(Form("%s%s_%s","h1D_",    N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h1D_",     N_Raw->GetName()),nBinPT1D,fArrPT1D) );
                else                    fResults.push_back( new TH1F(Form("%s%s_%s","h2Dbin_", N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h2Dbin_",  N_Raw->GetName()),nBinPT2D,fArrPT2D) );
                if ( strncmp(N_Raw->GetName(),"anSS",4) == 0 )   {
                    if  ( fDimension == 1)  fResults.at(iTer)->SetName("hRAW_1D");
                    else                    fResults.at(iTer)->SetName("hRAW_1D_in_2D_bin");
                }
                iTer++;
            }
            for ( auto fCoeff : fFitresultsStore[iFit]->constPars() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                if  ( fDimension == 1)  fResults.push_back( new TH1F(Form("%s%s_%s","h1D_",    N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h1D_",     N_Raw->GetName()),nBinPT1D,fArrPT1D) );
                else                    fResults.push_back( new TH1F(Form("%s%s_%s","h2Dbin_", N_Raw->GetName(), fNameFile.Data()),  Form("%s%s","h2Dbin_",  N_Raw->GetName()),nBinPT2D,fArrPT2D) );
            }
        }
        //
        //>>    Filling Raw Count Histograms
        Int_t   iTer = 0;
        for ( auto fCoeff : fFitresultsStore[iFit]->floatParsFinal() )   {
            auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
            fResults.at(iTer)->SetBinContent          (iFit+1,N_Raw->getVal());
            fResults.at(iTer)->SetBinError            (iFit+1,N_Raw->getError());
            iTer++;
        }
        for ( auto fCoeff : fFitresultsStore[iFit]->constPars() )   {
            auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
            fResults.at(iTer)->SetBinContent          (iFit+1,N_Raw->getVal());
            fResults.at(iTer)->SetBinError            (iFit+1,N_Raw->getError());
            iTer++;
        }
    }
    //
    return  fResults;
}
std::vector<TH1F*>          FitModel                        ( TH1F **hTarget, TH1F* hSlopReference, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return FitModel(hTarget,hSlopReference,NULL_ROOFITPTR2,1,fOption,fTargetPath,fNameFile);
}
std::vector<TH1F*>          FitModel                        ( TH1F **hTarget, TH1F* hSlopReference, RooFitResult** &fFitresultsStore, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return FitModel(hTarget,hSlopReference,fFitresultsStore,2,fOption,fTargetPath,fNameFile);
}
//
//_____________________________________________________________________________
//
//
//  --  --  2D Analysis Functions  --  --  //
//
//_____________________________________________________________________________
//
RooFitResult               *fFitCoreModel                   ( TH2F * THdata, TH1F* hSlopReference, RooFitResult * fFitShapeX, RooFitResult * fFitShapeY, TString fName = "", TString fOption = "", Int_t PTindex = -1, Int_t PTjndex = -1, TString fPathToSave = "./result/SEFitCheck" )    {
    //
    //>>    Silencing TCanvas Pop-Up and speeding Fit procedures
    gROOT->SetBatch(kTRUE);
    //
    //>>    Selecting all Options
    Bool_t  fFreeBkg,   fUsePoly,   fAlwaysFit, fSaveToFile;
    //
    fFreeBkg    =   fOption.Contains("BKG", TString::kIgnoreCase);
    fUsePoly    =   fOption.Contains("POL", TString::kIgnoreCase);
    fAlwaysFit  =   fOption.Contains("FIT", TString::kIgnoreCase);
    fSaveToFile =   fOption.Contains("SAVE",TString::kIgnoreCase);
    //
    //>>    Check there is a reasonable amount of entries form the general method, or if it has been overrided
    if ( !fIsWorthFitting( THdata ) || fAlwaysFit ) return nullptr;
    //
    //>>    Global Variables
    Double_t                        fInvMassValMax, fInvMassValMin;
    SetBoundaries(  fOption,        fInvMassValMin, fInvMassValMax  );
    Int_t           nEntries    =   THdata->GetEntries();
    RooRealVar      xInvMass    =   RooRealVar        ("xInvMass",  "xInvMass",     fInvMassValMin, fInvMassValMax                                  );
    RooRealVar      yInvMass    =   RooRealVar        ("yInvMass",  "yInvMass",     fInvMassValMin, fInvMassValMax                                  );
    RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         RooArgList(xInvMass,yInvMass),  Import(*THdata)                 );
    RooDataHist*    dataLoose   =   new RooDataHist   ("DataLoose", "DataLoose",    RooArgList(xInvMass,yInvMass),  Import(*fLooseErrors(THdata))   );
    Int_t kNCycle               =   kNCycle_;
    //
    //>>    Recovering the previous fits
    RooArgSet      *fXShapes    =   new RooArgSet(fFitShapeX->floatParsFinal(),fFitShapeX->constPars());
    RooArgSet      *fYShapes    =   new RooArgSet(fFitShapeY->floatParsFinal(),fFitShapeY->constPars());
    //
    //>>    Background
    RooRealVar  ch1x, ch2x, ch3x, ch4x, ch1y, ch2y, ch3y, ch4y;
    //
    if ( fFreeBkg )
    {
        ch1x     = RooRealVar ("ch1x","ch1x"     ,fXShapes->getRealValue("ch1",0),fXShapes->getRealValue("ch1",0)-0.1,fXShapes->getRealValue("ch1",0)+0.1);
        ch2x     = RooRealVar ("ch2x","ch2x"     ,fXShapes->getRealValue("ch2",0),fXShapes->getRealValue("ch2",0)-0.1,fXShapes->getRealValue("ch2",0)+0.1);
        ch3x     = RooRealVar ("ch3x","ch3x"     ,fXShapes->getRealValue("ch3",0),fXShapes->getRealValue("ch3",0)-0.1,fXShapes->getRealValue("ch3",0)+0.1);
        ch4x     = RooRealVar ("ch4x","ch4x"     ,fXShapes->getRealValue("ch4",0),fXShapes->getRealValue("ch4",0)-0.1,fXShapes->getRealValue("ch4",0)+0.1);
        ch1y     = RooRealVar ("ch1y","ch1y"     ,fYShapes->getRealValue("ch1",0),fYShapes->getRealValue("ch1",0)-0.1,fYShapes->getRealValue("ch1",0)+0.1);
        ch2y     = RooRealVar ("ch2y","ch2y"     ,fYShapes->getRealValue("ch2",0),fYShapes->getRealValue("ch2",0)-0.1,fYShapes->getRealValue("ch2",0)+0.1);
        ch3y     = RooRealVar ("ch3y","ch3y"     ,fYShapes->getRealValue("ch3",0),fYShapes->getRealValue("ch3",0)-0.1,fYShapes->getRealValue("ch3",0)+0.1);
        ch4y     = RooRealVar ("ch4y","ch4y"     ,fYShapes->getRealValue("ch4",0),fYShapes->getRealValue("ch4",0)-0.1,fYShapes->getRealValue("ch4",0)+0.1);
    }
    else
    {
        ch1x     = RooRealVar ("ch1x","ch1x"     ,fXShapes->getRealValue("ch1",0));
        ch2x     = RooRealVar ("ch2x","ch2x"     ,fXShapes->getRealValue("ch2",0));
        ch3x     = RooRealVar ("ch3x","ch3x"     ,fXShapes->getRealValue("ch3",0));
        ch4x     = RooRealVar ("ch4x","ch4x"     ,fXShapes->getRealValue("ch4",0));
        ch1y     = RooRealVar ("ch1y","ch1y"     ,fYShapes->getRealValue("ch1",0));
        ch2y     = RooRealVar ("ch2y","ch2y"     ,fYShapes->getRealValue("ch2",0));
        ch3y     = RooRealVar ("ch3y","ch3y"     ,fYShapes->getRealValue("ch3",0));
        ch4y     = RooRealVar ("ch4y","ch4y"     ,fYShapes->getRealValue("ch4",0));
    }
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,fXShapes->getRealValue("bMass",0));
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,fXShapes->getRealValue("bWidt",0));
    RooRealVar pSlopex  = RooRealVar ("pSlopx","pSlopx" ,fXShapes->getRealValue("bSlop",0));
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,fYShapes->getRealValue("bMass",0));
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,fYShapes->getRealValue("bWidt",0));
    RooRealVar pSlopey  = RooRealVar ("pSlopy","pSlopy" ,fYShapes->getRealValue("bSlop",0));
    
    // Coefficients
    auto fS__x  =   fXShapes->getRealValue("anSS",0);
    auto fB__x  =   fXShapes->getRealValue("anBB",0);
    auto fS__y  =   fYShapes->getRealValue("anSS",0);
    auto fB__y  =   fYShapes->getRealValue("anBB",0);
    auto fTot_  =   fS__x*fS__y+fB__x*fB__y+fS__x*fB__y+fB__x*fS__y;
    
    RooRealVar n0       = RooRealVar ("anSS2D","anSS2D" ,nEntries*(fS__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n1       = RooRealVar ("anBB2D","anBB2D" ,nEntries*(fB__x*fB__y)/fTot_,0.,nEntries);
    RooRealVar n2       = RooRealVar ("anBS2D","anBS2D" ,nEntries*(fB__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n3       = RooRealVar ("anSB2D","anSB2D" ,nEntries*(fS__x*fB__y)/fTot_,0.,nEntries);
    
    // PDFs
    RooChebychev        fBkgx   ("fBkgx","fBkgx"        ,xInvMass,RooArgSet(ch1x,ch2x,ch3x,ch4x));
    RooPolynomial       fBk2x   ("fBkgx","fBkgx"        ,xInvMass,RooArgSet(ch1x,ch2x,ch3x,ch4x));
    RooVoigtian         fSigx   ("fSigx","fSigx"        ,xInvMass,pMassx,pWidthx,pSlopex);
    RooChebychev        fBkgy   ("fBkgy","fBkgy"        ,yInvMass,RooArgSet(ch1y,ch2y,ch3y,ch4y));
    RooPolynomial       fBk2y   ("fBkgy","fBkgy"        ,yInvMass,RooArgSet(ch1y,ch2y,ch3y,ch4y));
    RooVoigtian         fSigy   ("fSigy","fSigy"        ,yInvMass,pMassy,pWidthy,pSlopey);
    RooProdPdf         *fBB;
    if ( fUsePoly )   {
        fBB= new RooProdPdf         ("fBkg","fBkg"          ,fBk2x,fBk2y);
    }   else    {
        fBB= new RooProdPdf         ("fBkg","fBkg"          ,fBkgx,fBkgy);
    }
    RooProdPdf         *fSB;
    if ( fUsePoly )   {
        fSB= new RooProdPdf         ("fSigBkg","fSigBkg"    ,fSigx,fBk2y);
    }   else    {
        fSB= new RooProdPdf         ("fSigBkg","fSigBkg"    ,fSigx,fBkgy);
    }
    RooProdPdf         *fBS;
    if ( fUsePoly )   {
        fBS= new RooProdPdf         ("fBkgSig","fBkgSig"          ,fBk2x,fSigy);
    }   else    {
        fBS= new RooProdPdf         ("fBkgSig","fBkgSig"          ,fBkgx,fSigy);
    }
    RooProdPdf          fSS     ("fSigSig","fSigSig"    ,fSigx,fSigy);
    RooAddPdf           fMod    ("fMod2D","fMod2D"      ,RooArgList(*fBB,fSS,*fSB,*fBS),RooArgList(n1,n0,n3,n2));
    
    RooFitResult* fFitResults;
    fFitResults      =   fMod.fitTo(*dataLoose,Save(),NumCPU(kCPU_use,kCPUStrategy));
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
        fFitResults      =   fMod.fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
        auto N_Raw      =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("anSS2D"));
        if ( fIsResultAcceptable(N_Raw->getVal(),N_Raw->getError(),k2DErrorLimit) ) break;
    }
    // Save to file
    if ( fSaveToFile || kSaveToFile )
    {
        int         nBinsPrint      =   4;
        double      dIncrement      =   (fInvMassValMax-fInvMassValMin)/nBinsPrint;
        TLatex*     latext          =   new TLatex();
        TCanvas*    cTotal          =   new TCanvas("","",0,45,1440,855);
                    cTotal          ->  SetTitle(Form("Slices of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
                    cTotal          ->  SetName(Form("PT_%.1f_%.1f__%.1f_%.1f_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()));
                    cTotal          ->  Divide(2,nBinsPrint);
        
                            xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                            yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()),
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame  =   yInvMass.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);

                            xInvMass.setRange("fDrawRange",fInvMassValMin+i*dIncrement,fInvMassValMin+(i+1)*dIncrement);
                            yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);

            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"yInvMass2D");
            
            cTotal->cd( 2*i+1 );
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
                                        xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()),
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame      =   xInvMass.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);
            
                                        xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        yInvMass.setRange("fDrawRange",fInvMassValMin+i*dIncrement,fInvMassValMin+(i+1)*dIncrement);
                                                                            
            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"xInvMass2D");
            
            cTotal->cd( 2*i+2 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.3f < m^{y}_{K^{+}K^{-}} < %.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.3f < m^{y}_{K^{+}K^{-}} < %.3f",fInvMassValMin+dIncrement*i,fInvMassValMin+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        xInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
                                        yInvMass.setRange("fDrawRange",fInvMassValMin,fInvMassValMax);
        cTotal ->Write();
        cTotal               ->SaveAs(Form("%s/PT_%.1f_%.1f__%.1f_%.1f_%s.pdf",fPathToSave.Data(),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fName.Data()));
        delete cTotal;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    // Fit
    return fFitResults;
}
//
//_____________________________________________________________________________
//
std::vector<TH2F*>          FitModel                        ( TH1F  **hShapeFit, TH1F* hSlopReference, TH2F***hTarget, RooFitResult*** &fFitresultsStore = NULL_ROOFITPTR3, TString fOption = "", TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", std::vector<TH1F*> &f1Din2DbinCheck = NULL_VECTOR )  {
    std::vector<TH2F*>  fResults;
    RooFitResult      **fShapeStore = nullptr;
    if ( !fFitresultsStore )    {
        fFitresultsStore = new RooFitResult**[nBinPT2D];
        for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ ) {
            fFitresultsStore[iTer]  =    new RooFitResult    *[nBinPT2D];
        }
    }
    //
    //>>    Recover 1D Shapes
    f1Din2DbinCheck     =   FitModel(hShapeFit,hSlopReference,fShapeStore,2,fOption,fTargetPath,fNameFile);
    //
    //>>    Fit 2D Histograms
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        for ( Int_t jFit = iFit; jFit < nBinPT2D; jFit++ ) {
            //
            //>>    Protection Against SegFault
            if ( !fShapeStore[iFit] ) continue;
            if ( !fShapeStore[jFit] ) continue;
            //
            //>>    Fit
            fFitresultsStore[iFit][jFit]       =  fFitCoreModel(hTarget[iFit][jFit],hSlopReference,fShapeStore[iFit],fShapeStore[jFit],fNameFile,fOption,iFit,jFit,fTargetPath);
            //
            //>>    Protection Against SegFault
            if ( !fFitresultsStore[iFit][jFit] ) continue;
            //
            //>>    Building Raw Count histograms
            if ( iFit == 0 && jFit == 0 )   {
                Int_t   iTer = 0;
                for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                    auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                    fResults.push_back( new TH2F(Form("%s_%s",N_Raw->GetName(),fNameFile.Data()),N_Raw->GetName(),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D) );
                    if ( strncmp(N_Raw->GetName(),"anSS2D",6) == 0 ) fResults.at(iTer)->SetName("hRAW_2D");
                    iTer++;
                }
                for ( auto fCoeff : fFitresultsStore[iFit][jFit]->constPars() )   {
                    auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                    fResults.push_back( new TH2F(Form("%s_%s",N_Raw->GetName(),fNameFile.Data()),N_Raw->GetName(),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D) );
                }
            }
            //
            //>>    Filling Raw Count Histograms
            Int_t   iTer = 0;
            for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,N_Raw->getVal());
                fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,N_Raw->getError());
                if ( iFit == jFit ) {
                    fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,2.*N_Raw->getVal());
                    fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,2.*N_Raw->getError());
                }
                iTer++;
            }
            for ( auto fCoeff : fFitresultsStore[iFit][jFit]->constPars() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,N_Raw->getVal());
                fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,N_Raw->getError());
                if ( iFit == jFit ) {
                    fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,2.*N_Raw->getVal());
                    fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,2.*N_Raw->getError());
                }
                iTer++;
            }
        }
    }
    return fResults;
}
std::vector<TH2F*>          FitModel                        ( TH2F***hTarget, TH1F* hSlopReference,  RooFitResult**  fShapeStore, RooFitResult*** &fFitresultsStore = NULL_ROOFITPTR3, TString fOption = "", TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", std::vector<TH1F*> &f1Din2DbinCheck = NULL_VECTOR )  {
    std::vector<TH2F*>  fResults;
    if ( !fFitresultsStore )    {
        fFitresultsStore = new RooFitResult**[nBinPT2D];
        for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ ) {
            fFitresultsStore[iTer]  =    new RooFitResult    *[nBinPT2D];
        }
    }
    //
    //>>    Fit 2D Histograms
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        for ( Int_t jFit = 0; jFit < nBinPT2D; jFit++ ) {
            //
            //>>    Protection Against SegFault
            if ( !fShapeStore[iFit] ) continue;
            if ( !fShapeStore[jFit] ) continue;
            //
            //>>    Fit
            fFitresultsStore[iFit][jFit]       =  fFitCoreModel(hTarget[iFit][jFit],hSlopReference,fShapeStore[iFit],fShapeStore[jFit],fNameFile,fOption,iFit,jFit,fTargetPath);
            //
            //>>    Protection Against SegFault
            if ( !fFitresultsStore[iFit][jFit] ) continue;
            //
            //>>    Building Raw Count histograms
            if ( iFit == 0 && jFit == 0 )   {
                Int_t   iTer = 0;
                for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                    auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                    fResults.push_back( new TH2F(N_Raw->GetName(),N_Raw->GetName(),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D) );
                    if ( strncmp(N_Raw->GetName(),"anSS2D",6) == 0 ) fResults.at(iTer)->SetName("hRAW_2D");
                    iTer++;
                }
            }
            //
            //>>    Filling Raw Count Histograms
            Int_t   iTer = 0;
            for ( auto fCoeff : fFitresultsStore[iFit][jFit]->floatParsFinal() )   {
                auto N_Raw      = static_cast<RooRealVar*>(fCoeff);
                fResults.at(iTer)->SetBinContent          (iFit+1,jFit+1,N_Raw->getVal());
                fResults.at(iTer)->SetBinError            (iFit+1,jFit+1,N_Raw->getError());
                iTer++;
            }
        }
    }
    return fResults;
}
std::vector<TH2F*>          FitModel                        ( TH1F  **hShapeFit, TH1F* hSlopReference, TH2F***hTarget, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "",  std::vector<TH1F*> &f1Din2DbinCheck = NULL_VECTOR )  {
    return  FitModel(hShapeFit,hSlopReference,hTarget,NULL_ROOFITPTR3,fOption,fTargetPath,fNameFile);
}
std::vector<TH2F*>          FitModel                        ( TH1F  **hShapeFit, TH1F* hSlopReference, TH2F***hTarget, std::vector<TH1F*> &f1Din2DbinCheck, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return  FitModel(hShapeFit,hSlopReference,hTarget,NULL_ROOFITPTR3,fOption,fTargetPath,fNameFile,f1Din2DbinCheck);
}
std::vector<TH2F*>          FitModel                        ( TH2F***hTarget, TH1F* hSlopReference, RooFitResult**  fShapeStore, std::vector<TH1F*> &f1Din2DbinCheck, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "" )  {
    return  FitModel(hTarget,hSlopReference,fShapeStore,NULL_ROOFITPTR3,"",fTargetPath,fNameFile,f1Din2DbinCheck);
}
std::vector<TH2F*>          FitModel                        ( TH2F***hTarget, TH1F* hSlopReference, RooFitResult**  fShapeStore, std::vector<TH1F*> &f1Din2DbinCheck, TString fTargetPath = "./result/SEFitCheck", TString fNameFile = "", TString fOption = "" )  {
    return  FitModel(hTarget,hSlopReference,fShapeStore,NULL_ROOFITPTR3,fOption,fTargetPath,fNameFile,f1Din2DbinCheck);
}
//
//_____________________________________________________________________________

//------------------------------//
//    !TODO: Clean the mess below      //
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
            //aTarget->GetZaxis()->SetTitle("#frac{d^{3}N_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
            //aTarget->GetZaxis()->SetTitleOffset(1.15);
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
        if ( fOption.Contains("Multiplicity",TString::kIgnoreCase) )    { kDoMultiplicity = true;   cout << "[INFO] Multiplicity option chosen" <<endl;}
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
void                fSetFunction                    ( TF1* fFitFunction = fLevyTsallis, Double_t fIntegral = 0.032 ) {
    fSetAllFunctions();
    //
    //____________________________________________Mass
    auto    nMass   =   fFitFunction    ->GetParNumber("Mass");
    if ( nMass != -1 )   {
        fFitFunction    ->  FixParameter (nMass,kPhiMesonMass_);
    }
    //
    //____________________________________________dN/dy
    auto    ndNdy   =   fFitFunction    ->GetParNumber("dN_dy");
    if ( ndNdy != -1 )   {
        fFitFunction    ->  SetParLimits(ndNdy,0.,1.);
        fFitFunction    ->  SetParameter(ndNdy,fIntegral);
    }
    //
    //____________________________________________dN/dy
    auto    nEnne   =   fFitFunction    ->GetParNumber("n");
    if ( nEnne != -1 )   {
        fFitFunction    ->  SetParLimits(nEnne,0.,100.);
        fFitFunction    ->  SetParameter(nEnne,7.18);
        if ( strncmp(fFitFunction->GetName(),"LevyTsallis",11) == 0 )   {
            fFitFunction  ->  SetParLimits(nEnne,2.1,15.);
            fFitFunction  ->  SetParameter(nEnne,7.18);
        }
        if ( strncmp(fFitFunction->GetName(),"PowerLaw",13) == 0 )          fFitFunction  ->  SetParameter(nEnne,13);
    }
    //
    //____________________________________________dN/dy
    auto    nSlop   =   fFitFunction    ->GetParNumber("T");
    if ( nSlop != -1 )   {
        fFitFunction    ->  SetParLimits(nSlop,0.05,30);
        fFitFunction    ->  SetParameter(nSlop,.330);
        if ( strncmp(fFitFunction->GetName(),"LevyTsallis",11) == 0 )       fFitFunction  ->  SetParameter(nSlop,.300);
        if ( strncmp(fFitFunction->GetName(),"MTExponential",13) == 0 )     fFitFunction  ->  SetParameter(nSlop,.560);
        if ( strncmp(fFitFunction->GetName(),"PTExponential",13) == 0 )     fFitFunction  ->  SetParameter(nSlop,.560);
        if ( strncmp(fFitFunction->GetName(),"FermiDirac",13) == 0 )        fFitFunction  ->  SetParameter(nSlop,.330);
        if ( strncmp(fFitFunction->GetName(),"Boltzmann",13) == 0 )         fFitFunction  ->  SetParameter(nSlop,.350);
        if ( strncmp(fFitFunction->GetName(),"BoseEinstein",13) == 0 )      fFitFunction  ->  SetParameter(nSlop,.350);
        if ( strncmp(fFitFunction->GetName(),"PowerLaw",13) == 0 )          fFitFunction  ->  SetParameter(nSlop,.550);
    }
}
void                fSetFunction                    ( TF2* fFitFunction , Double_t fIntegral = 3.e3 ) {
    fSetAllFunctions();
    //
    //____________________________________________Mass
    auto    xMass   =   fFitFunction    ->GetParNumber("XMass");
    auto    yMass   =   fFitFunction    ->GetParNumber("YMass");
    fFitFunction    ->  FixParameter (xMass,kPhiMesonMass_);
    fFitFunction    ->  FixParameter (yMass,kPhiMesonMass_);
    //
    //____________________________________________dN/dy
    auto    XdNdy   =   fFitFunction    ->GetParNumber("XdN_dy");
    auto    YdNdy   =   fFitFunction    ->GetParNumber("YdN_dy");
    auto    _dNdy   =   fFitFunction    ->GetParNumber("dN_dy");
    fFitFunction    ->  SetParLimits(XdNdy,0.,1.e6);
    fFitFunction    ->  SetParameter(XdNdy,fIntegral);
    fFitFunction    ->  SetParLimits(YdNdy,0.,1.e6);
    fFitFunction    ->  SetParameter(YdNdy,fIntegral);
    fFitFunction    ->  SetParLimits(_dNdy,0.,1.e6);
    fFitFunction    ->  SetParameter(_dNdy,fIntegral);
    //
    //____________________________________________dN/dy
    auto    XEnne   =   fFitFunction    ->GetParNumber("Xn");
    auto    YEnne   =   fFitFunction    ->GetParNumber("Yn");
    fFitFunction    ->  SetParLimits(XEnne,0.,100.);
    fFitFunction    ->  SetParameter(XEnne,7.18);
    fFitFunction    ->  SetParLimits(YEnne,0.,100.);
    fFitFunction    ->  SetParameter(YEnne,7.18);
    //
    //____________________________________________dN/dy
    auto    xSlop   =   fFitFunction    ->GetParNumber("XT");
    auto    ySlop   =   fFitFunction    ->GetParNumber("YT");
    fFitFunction    ->  SetParLimits(xSlop,0.05,30);
    fFitFunction    ->  SetParameter(xSlop,.330);
    fFitFunction    ->  SetParLimits(ySlop,0.05,30);
    fFitFunction    ->  SetParameter(ySlop,.330);
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
TGraphAsymmErrors*
fSetSystErrors                                      ( TGraphAsymmErrors*                gStatistics ) {
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gStatistics);
    for ( Int_t iPnt = 0; iPnt < fResult->GetN(); iPnt++ ) {
        auto    fYValue =   fResult ->  GetPointY(iPnt);
        fResult ->  SetPointEYhigh  ( iPnt, fYValue*sqrt(kSysHig_BR*kSysHig_BR+kSysHig_TR*kSysHig_TR+kSysHig_PD*kSysHig_PD+kSysHig_1D_SE*kSysHig_1D_SE) );
        fResult ->  SetPointEYlow   ( iPnt, fYValue*sqrt(kSysLow_BR*kSysLow_BR+kSysLow_TR*kSysLow_TR+kSysLow_PD*kSysLow_PD+kSysLow_1D_SE*kSysLow_1D_SE) );
    }
    return fResult;
}
std::vector<TGraphAsymmErrors*>
fSetSystErrors                                      ( std::vector<TGraphAsymmErrors*>   gStatistics ) {
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
TH1F*
fSetSystErrors                                      ( TH1F*                             hStatistics ) {
    TH1F               *fResult =   new TH1F            (*hStatistics);
    for ( Int_t iBin = 0; iBin < fResult->GetNbinsX(); iBin++ ) {
        auto    fBinContent     =   fResult->GetBinContent( iBin );
        fResult->SetBinError( iBin, fBinContent*sqrt(kSysHig_BR*kSysHig_BR+kSysHig_TR*kSysHig_TR+kSysHig_PD*kSysHig_PD+kSysHig_1D_SE*kSysHig_1D_SE) );
    }
    return  fResult;
}
std::vector<TH1F*>
fSetSystErrors                                      ( std::vector<TH1F*>                hStatistics ) {
    std::vector<TH1F*> fResult;
    for ( auto& iHisto : hStatistics )  {
        auto    fCurrentHisto   =   new TH1F(*iHisto);
        for ( Int_t iBin = 0; iBin < iHisto->GetNbinsX(); iBin++ ) {
            auto    fBinContent     =   fCurrentHisto->GetBinContent( iBin );
            fCurrentHisto->SetBinError( iBin, fBinContent*sqrt(4*(kSysHig_BR*kSysHig_BR+kSysHig_TR*kSysHig_TR+kSysHig_PD*kSysHig_PD)+kSysHig_2D_SE*kSysHig_2D_SE) );
        }
        fResult.push_back(fCurrentHisto);
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
std::vector<float>  fEvaluateError                  ( bool fIsConditional, TGraphAsymmErrors* gStatic, TGraphAsymmErrors* gMoveable, Float_t fIntegral, TF1* fFitFunc, Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fName = "", TString fFolder = "")  {
    auto fExtrap    =   fFitFunc->Integral(0.,fMinPT1D);
    auto fMeanPT    =   fMeasureMeanPT(fFitFunc,gStatic);
    TGraphAsymmErrors      *gTotal      =   new TGraphAsymmErrors(*(fSumErrors(gStatic,gMoveable)));
    std::vector<float> fResult;
    //
    //  TCanvas w/ options
    auto    cDrawResult     =   new TCanvas("","",1500,500);
    auto    fText           =   new TLatex();
    gStyle->SetOptStat(0);
    cDrawResult             ->  Divide(3,1);
    cDrawResult             ->  cd(1);
    gPad->SetLogy();
    //
    //  Draw Spectra
    gTotal->Draw();
    gTotal->SetMarkerStyle(20);
    gTotal->SetMarkerColor(kBlue);
    gTotal->SetLineColor(kBlue);
    //
    //  Generate utility histogram to evaluate statistical error
    std::vector<float> fStatVariation;
    std::vector<float> fStatVariationMPT;
    for ( Int_t iFit = 0; iFit < kStatEvalCycles; iFit++ )  {
        //  Set Standard Fit
        fSetFunction(fFitFunc,fIntegral);
        //
        //  Generating the Fit TGraph
        auto fSubject   =   uRandomizePoints(gStatic,gMoveable);
        //
        fSubject    ->  Fit(fFitFunc,   "IMESQ","", fMinimumFitRange,fMaximumFitRange);
        //
        fStatVariation.push_back(   fFitFunc  ->Integral(0.,fMinPT1D));
        //
        fStatVariationMPT.push_back(fMeasureMeanPT(fFitFunc,fSubject));
        //
        //  Draw Function
        fFitFunc->DrawCopy("same");
    }
    cDrawResult             ->  cd(2);
    uCleanOutsiders         ( fStatVariation );
    auto    hStatIntegral   =   uBuildTH1F( fStatVariation, 50 );
    hStatIntegral           ->  SetNameTitle( Form("hIntegral_%s",fName.Data()), "dN/dy variations" );
    hStatIntegral           ->  Fit("gaus","IMESQ","");
    auto StatErrorMean      =   hStatIntegral->GetFunction("gaus")->GetParameter(1);
    auto StatError          =   hStatIntegral->GetFunction("gaus")->GetParameter(2);
    fResult                 .push_back( StatError );
    fResult                 .push_back( StatError );
    hStatIntegral           ->  Draw();
    fText   ->  DrawLatexNDC(0.45,0.850,Form("Mean: %.3e %s",100*(StatErrorMean/fExtrap -1),"%"));
    fText   ->  DrawLatexNDC(0.45,0.800,Form("#sigma_{Dev}: %.3e %s",100*(StatError)/(fExtrap),"%"));
    //
    cDrawResult             ->  cd(3);
    uCleanOutsiders         ( fStatVariationMPT );
    auto    hStatIntegralMPT   =   uBuildTH1F( fStatVariationMPT, 50 );
    hStatIntegralMPT           ->  SetNameTitle( Form("hIntegralMPT_%s",fName.Data()), "Mean p_{T} variations" );
    hStatIntegralMPT           ->  Fit("gaus","IMESQ","");
    auto StatErMPTMean      =   hStatIntegralMPT->GetFunction("gaus")->GetParameter(1);
    auto StatErMPT          =   hStatIntegralMPT->GetFunction("gaus")->GetParameter(2);
    //auto StatErMPTMean      =   hStatIntegralMPT->GetMean();
    //auto StatErMPT          =   hStatIntegralMPT->GetRMS();
    auto ReferenceMPT       =   (fMeanPT)/(fExtrap+fIntegral);
    fResult                 .push_back( StatErMPT );
    fResult                 .push_back( StatErMPT );
    hStatIntegralMPT           ->  Draw();
    fText   ->  DrawLatexNDC(0.45,0.850,Form("Mean: %.3e %s",100*(StatErMPTMean/ReferenceMPT -1),"%"));
    fText   ->  DrawLatexNDC(0.45,0.800,Form("#sigma_{Dev}: %.3e %s",100*(StatErMPT)/(ReferenceMPT),"%"));
    //
    cDrawResult             ->  cd(1);
    gTotal->Draw("same EP");
    //
    //  Write info on the Extrapolation:
    fText   ->  DrawLatexNDC(0.45,0.825,Form("#frac{dN}{dy}: %.3e",fResult[0]));
    fText   ->  DrawLatexNDC(0.45,0.650,Form("#sigma_{high}: %.3e %s",100*(fResult[2])/(fResult[0]),"%"));
    //
    cDrawResult->Write();
    if ( fIsConditional )   cDrawResult->SaveAs(Form("result/yield/ExtrapolateCheck/2D/ErrorFits_%s.pdf",fName.Data()));
    else                    cDrawResult->SaveAs(Form("result/yield/ExtrapolateCheck/1D/ErrorFits_%s.pdf",fName.Data()));
    delete cDrawResult;
    //
    return fResult;
}
std::vector<float>  fEvaluateError                  ( bool fIsConditional, TH1F* gStatic, TH1F* gMoveable, Float_t fIntegral, Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fName = "", TString fFolder = "")  {
    //
    auto fExtrap    =   fLevyTsallis->Integral(0.,fMinPT1D);
    auto fMeanPT    =   fMeasureMeanPT(fLevyTsallis,gStatic);
    TH1F      *gTotal      =    fSumErrors(gStatic,gMoveable);
    gTotal->SetName("gTotal");
    std::vector<float> fResult;
    //
    //  TCanvas w/ options
    auto    cDrawResult     =   new TCanvas("cDrawResult","cDrawResult",1800,600);
    auto    fText           =   new TLatex();
    gStyle->SetOptStat(0);
    cDrawResult             ->  Divide(3,1);
    cDrawResult             ->  cd(1);
    gPad->SetLogy();
    //
    //  Draw Spectra
    gTotal->SetMarkerStyle(26);
    gTotal->SetMarkerSize(2);
    gTotal->SetMarkerColor(38);
    gTotal->SetLineColor(kBlack);
    gTotal->Draw();
    //
    //  Generate utility histogram to evaluate statistical error
    std::vector<float> fStatVariation;
    std::vector<float> fStatVariationMPT;
    for ( Int_t iFit = 0; iFit < kStatEvalCycles; iFit++ )  {
        //  Set Standard Fit
        fSetFunction(fLevyTsallis,fIntegral);
        //
        auto dump = new TCanvas();
        //  Generating the Fit TGraph
        auto fSubject   =   uRandomizePoints(gStatic,gMoveable);
        //
        fSubject    ->  Fit(fLevyTsallis,   "EMRQS","", fMinimumFitRange,fMaximumFitRange);
        fSubject    ->  Fit(fLevyTsallis,   "EMRQSI","", fMinimumFitRange,fMaximumFitRange);
        //
        delete dump;
        //
        fStatVariation.push_back(   fLevyTsallis  ->Integral(0.,fMinPT1D));
        //
        fStatVariationMPT.push_back(fMeasureMeanPT(fLevyTsallis,fSubject));
        //
        //  Draw Function
        cDrawResult->cd(1);
        fSubject->GetFunction("LevyTsallis")->DrawCopy("SAME");
    }
    cDrawResult             ->  cd(2);
    uCleanOutsiders         ( fStatVariation );
    auto    hStatIntegral   =   uBuildTH1F( fStatVariation, 50 );
    hStatIntegral           ->  SetNameTitle( Form("hIntegral_%s",fName.Data()), "dN/dy variations" );
    hStatIntegral           ->  Fit("gaus","I EM R Q S","");
    auto StatErrorMean      =   hStatIntegral->GetFunction("gaus")->GetParameter(1);
    auto StatError          =   hStatIntegral->GetFunction("gaus")->GetParameter(2);
    fResult                 .push_back( StatError );
    fResult                 .push_back( StatError );
    hStatIntegral           ->  Draw();
    fText   ->  DrawLatexNDC(0.45,0.850,Form("Mean: %.3e %s",100*(StatErrorMean/fExtrap -1),"%"));
    fText   ->  DrawLatexNDC(0.45,0.800,Form("#sigma_{Dev}: %.3e %s",100*(StatError)/(fExtrap),"%"));
    //
    cDrawResult             ->  cd(3);
    uCleanOutsiders         ( fStatVariationMPT );
    auto    hStatIntegralMPT   =   uBuildTH1F( fStatVariationMPT, 50 );
    hStatIntegralMPT           ->  SetNameTitle( Form("hIntegralMPT_%s",fName.Data()), "Mean p_{T} variations" );
    hStatIntegralMPT           ->  Fit("gaus","I EM R Q S","");
    auto StatErMPTMean      =   hStatIntegralMPT->GetFunction("gaus")->GetParameter(1);
    auto StatErMPT          =   hStatIntegralMPT->GetFunction("gaus")->GetParameter(2);
    auto ReferenceMPT       =   (fMeanPT)/(fExtrap+fIntegral);
    fResult                 .push_back( StatErMPT );
    fResult                 .push_back( StatErMPT );
    hStatIntegralMPT           ->  Draw();
    fText   ->  DrawLatexNDC(0.45,0.850,Form("Mean: %.3e %s",100*(StatErMPTMean/fMeanPT -1),"%"));
    fText   ->  DrawLatexNDC(0.45,0.800,Form("#sigma_{Dev}: %.3e %s",100*(StatErMPT/fMeanPT),"%"));
    //
    cDrawResult             ->  cd(1);
    gTotal->Draw("same EP");
    //
    cDrawResult->Write();
    if ( fIsConditional )   cDrawResult->SaveAs(Form("%s/2D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResult->SaveAs(Form("%s/1D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResult;
    //
    return fResult;
}
std::vector<float>  fEvaluateError                  ( bool fIsConditional, TH1F* gStatic, TH1F* gMoveable, std::vector<std::pair<TF1*,std::vector<float>>> fSystFunc, Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fName = "", TString fFolder = "")  {
    //
    auto fExtrap    =   fLevyTsallis->Integral(0.,fMinPT1D);
    auto fIntegral  =   fLevyTsallis->Integral(fMinPT1D,fMaxPT1D);
    auto fMeanPT    =   fMeasureMeanPT(fLevyTsallis,gStatic);
    TH1F      *gTotal      =    fSumErrors(gStatic,gMoveable);
    gTotal->SetName("gTotal");
    std::vector<float> fResult;
    //
    //  TCanvas w/ options
    auto    cDrawResult     =   new TCanvas("cDrawResult","cDrawResult",1800,600);
    auto    fText           =   new TLatex();
    gStyle->SetOptStat(0);
    cDrawResult             ->  Divide(3,1);
    cDrawResult             ->  cd(1);
    gPad->SetLogy();
    //
    //  Draw Spectra
    gTotal->SetMarkerStyle(26);
    gTotal->SetMarkerSize(2);
    gTotal->SetMarkerColor(kRed);
    gTotal->SetLineColor(kBlack);
    gTotal->GetXaxis()->SetRangeUser(0.,3.0);
    gTotal->Draw();
    //
    //  Generate utility histogram to evaluate statistical error
    std::vector<float> fStatVariation;
    std::vector<float> fStatVariationMPT;
    TLegend    *lAll    =   new TLegend(0.35,0.15,0.6,0.3);
    lAll->SetLineColorAlpha(kWhite,0.);
    lAll->SetFillColorAlpha(kWhite,0.);
    auto iTer = 0;
    auto jTer = 0;
    for ( auto iFuncRange : fSystFunc )  {
        //
        //  Prepping Fit Function
        fSetAllFunctions();
        auto    iFunc   =   iFuncRange.first;
        iFunc->SetLineColor(fGetRainbowColor(iTer,true));
        //
        //  Prepping Associated Ranges
        for ( auto iRange : iFuncRange.second )  {
            fSetFunction(iFunc,fIntegral);
            //
            gTotal      ->  Fit(iFunc,          "EMRQS","", fMinimumFitRange,iRange);
            gTotal      ->  Fit(iFunc,          "EMRQSI","",fMinimumFitRange,iRange);
            //
            fStatVariation.push_back(iFunc  ->Integral(0.,fMinPT1D));
            fStatVariationMPT.push_back(iFunc  ->Integral(0.,fMinPT1D));
            //
            //  Draw Function
            iFunc->SetRange(0.,iRange);
            cDrawResult             ->  cd(1);
            iFunc->DrawCopy("same");
            cDrawResult             ->  cd(2);
            iFunc->DrawCopy("same");
            jTer++;
        }
        //
        //  TLegend
        lAll    ->  AddEntry(iFunc,iFunc->GetName(),"L");
        iTer++;
    }
    cDrawResult             ->  cd(2);
    uCleanOutsiders         ( fStatVariation );
    auto    hStatIntegral   =   uBuildTH1F( fStatVariation );
    hStatIntegral           ->  SetNameTitle( Form("hIntegral_%s",fName.Data()), "dN/dy variations" );
    //hStatIntegral           ->  Fit("gaus","I EM R Q S","");
    auto StatErrorMean      =   hStatIntegral->GetMean();//GetFunction("gaus")->GetParameter(1);
    auto StatError          =   hStatIntegral->GetRMS();//->GetFunction("gaus")->GetParameter(2);
    fResult                 .push_back( StatError );
    fResult                 .push_back( StatError );
    hStatIntegral           ->  Draw();
    fText   ->  DrawLatexNDC(0.45,0.850,Form("Mean: %.3e %s",100*(StatErrorMean/fExtrap -1),"%"));
    fText   ->  DrawLatexNDC(0.45,0.800,Form("#sigma_{Dev}: %.3e %s",100*(StatError)/(fExtrap),"%"));
    //
    cDrawResult             ->  cd(3);
    uCleanOutsiders         ( fStatVariationMPT );
    auto    hStatIntegralMPT   =   uBuildTH1F( fStatVariationMPT );
    hStatIntegralMPT           ->  SetNameTitle( Form("hIntegralMPT_%s",fName.Data()), "Mean p_{T} variations" );
    //hStatIntegralMPT           ->  Fit("gaus","I EM R Q S","");
    auto StatErMPTMean      =   hStatIntegralMPT->GetMean();//->GetFunction("gaus")->GetParameter(1);
    auto StatErMPT          =   hStatIntegralMPT->GetRMS();//->GetFunction("gaus")->GetParameter(2);
    auto ReferenceMPT       =   (fMeanPT)/(fExtrap+fIntegral);
    fResult                 .push_back( StatErMPT );
    fResult                 .push_back( StatErMPT );
    hStatIntegralMPT           ->  Draw();
    fText   ->  DrawLatexNDC(0.45,0.850,Form("Mean: %.3e %s",100*(StatErMPTMean/fMeanPT -1),"%"));
    fText   ->  DrawLatexNDC(0.45,0.800,Form("#sigma_{Dev}: %.3e %s",100*(StatErMPT/fMeanPT),"%"));
    //
    cDrawResult             ->  cd(1);
    lAll->Draw("same");
    gTotal->Draw("same EP");
    //
    cDrawResult->Write();
    if ( fIsConditional )   cDrawResult->SaveAs(Form("%s/2D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResult->SaveAs(Form("%s/1D/ErrorFits_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResult;
    //
    return fResult;
}
//
//_____________________________________________________________________________
//
Double_t*           fExtrapolateModel               ( bool fIsConditional, TGraphAsymmErrors* gStatistics, TGraphAsymmErrors* gSystematics, Double_t fIntegral = 0.032, TString fName = "ExtrapolateSignal", TF1* fFitFunc = fLevyTsallis, Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fFolder = "" )    {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    //  Initialising the Fit Function
    fSetFunction(fFitFunc,fIntegral);
    fFitFunc->SetLineColor(kRed);
    fFitFunc->SetRange(0.,fMaxPT1D);
    //
    //  Setting the Fit Range
    Double_t    fMaxFitInt  =   TMath::Min( fMaximumFitRange, (double)fMaxPT1D );
    Double_t    fMinFitInt  =   TMath::Max( fMinimumFitRange, (double)fMinPT1D );
    //
    //  Generating a Total Error Spectra to fit and extrapolating at low pT
    TGraphAsymmErrors      *gTotal      =   new TGraphAsymmErrors(*(fSumErrors(gStatistics,gSystematics)));
    //
    //  Fitting a first time to evaluate integral in non-measured region
    gTotal                              ->  Fit(fFitFunc,"BIMRSQEX","",fMinFitInt,fMaxFitInt);
    fResult[0]                          =   fFitFunc  ->Integral(0.,fMinPT1D);
    //
    fResult[5]                          =   fMeasureMeanPT(fFitFunc,gTotal);
    //
    //  TCanvas w/ options
    TCanvas                *cDrawResult =   new TCanvas(Form("%s_%s",gStatistics->GetName(),fName.Data()));
    gStyle  ->SetOptStat(0);
    gPad    ->SetLogy();
    //
    //  Draw Spectra w/ function
    gTotal  ->Draw();
    fFitFunc->Draw("same");
    //
    //  Write info on the Extrapolation:
    TLatex  *fText   =   new TLatex();
    fText   ->  DrawLatexNDC(0.5,0.825,Form("#frac{dN}{dy} in [%.1f;%.1f]: %.7f",0.,fMinPT1D,fResult[0]));
    fText   ->  DrawLatexNDC(0.5,0.750,Form("Function: %s",fFitFunc->GetName()));
    fText   ->  DrawLatexNDC(0.5,0.700,Form("Fit range: [%.1f;%.1f]",fMinFitInt,fMaxFitInt));
    //  Save To File
    cDrawResult->Write();
    //  Graphical Check
    if ( fIsConditional )   cDrawResult->   SaveAs(Form("%s/2D/EvaluationFit_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResult->   SaveAs(Form("%s/1D/EvaluationFit_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResult;
    //
    auto fStatResults                   =   fEvaluateError(fIsConditional,gSystematics,gStatistics,fIntegral,fFitFunc,fMinFitInt,fMaxFitInt,TString("Stat_")+fName,fFolder);
    auto fSystResults                   =   fEvaluateError(fIsConditional,gStatistics,gSystematics,fIntegral,fFitFunc,fMinFitInt,fMaxFitInt,TString("Syst_")+fName,fFolder);
    //
    fResult[1] = fResult[2] = fStatResults.at(0);
    fResult[3] = fResult[4] = fSystResults.at(0);
    fResult[6] = fResult[7] = fStatResults.at(2);
    fResult[8] = fResult[9] = fSystResults.at(2);
    /*
    // !TODO: Scorporare la multifucntion
    //  Generate utility histogram to evaluate statistical error
    std::vector<float> fMultiFitVariation;
    TLegend    *lAll    =   new TLegend();
    auto iTer = 0;
    for ( auto iFuncRange : fSystFitFunctions )  {
        //
        //  Prepping Fit Function
        auto    iFunc   =   iFuncRange.first;
        iFunc->SetLineColor(fGetRainbowColor(iTer,true));
        //
        //  Prepping Associated Ranges
        for ( auto iRange : iFuncRange.second )  {
            fSetFunction(iFunc,fIntegral);
            gTotal      ->  Fit(iFunc,"IMES","",fMinFitInt,iRange);
            //
            fMultiFitVariation.push_back(iFunc  ->Integral(0.,fMinPT1D));
            //
            //  Draw Function
            iFunc->SetRange(0.,iRange);
            cDrawResult             ->  cd(1);
            iFunc->DrawCopy("same");
            cDrawResult             ->  cd(2);
            iFunc->DrawCopy("same");
        }
        //
        //  TLegend
        lAll    ->  AddEntry(iFunc,iFunc->GetName(),"L");
        iTer++;
    }
    cDrawResult             ->  cd(3);
    //
    uCleanOutsiders         ( fMultiFitVariation );
    auto    hMultiFitHisto  =   uBuildTH1F( fMultiFitVariation );
    hMultiFitHisto          ->  SetNameTitle( Form("hMultiFit_%s",fName.Data()), "hMultiFit" );
    fResult[3]              =   hMultiFitHisto->GetRMS();
    fResult[4]              =   hMultiFitHisto->GetRMS();
    hMultiFitHisto          ->  Draw();
    //
    cDrawResult             ->  cd(1);
    //
    lAll->Draw("same");
    //
    //  Write info on the Extrapolation:
    fText   ->  DrawLatexNDC(0.5,0.825,Form("#frac{dN}{dy}: %.3e",fResult[0]));
    fText   ->  DrawLatexNDC(0.5,0.750,Form("Mean: %.3e",hMultiFitHisto->GetMean()));
    fText   ->  DrawLatexNDC(0.5,0.700,Form("#sigma_{low}: %.3e",fResult[3]));
    fText   ->  DrawLatexNDC(0.5,0.650,Form("#sigma_{high}: %.3e",fResult[4]));
    //
    cDrawResult->Write();
    if ( fIsConditional )   cDrawResult->SaveAs(Form("result/yield/ExtrapolateCheck/2D/MultiFitsError_%s.pdf",fName.Data()));
    else                    cDrawResult->SaveAs(Form("result/yield/ExtrapolateCheck/1D/MultiFitsError_%s.pdf",fName.Data()));
    delete cDrawResult;
    //
     */
    //_____________________________________
    //
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
Double_t*           fExtrapolateModel               ( bool fIsConditional, TH1F* gStatistics, TH1F* gSystematics, Double_t fIntegral = 0.032, TString fName = "ExtrapolateSignal", Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fFolder = ""  )    {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    //  Initialising the Fit Function
    fSetFunction(fLevyTsallis,fIntegral);
    fLevyTsallis->SetLineColor(kRed);
    fLevyTsallis->SetRange(0.,fMaxPT1D);
    //
    //  Setting the Fit Range
    Double_t    fMaxFitInt  =   TMath::Min( fMaximumFitRange, (double)fMaxPT1D );
    Double_t    fMinFitInt  =   TMath::Max( fMinimumFitRange, (double)fMinPT1D );
    //
    //  Generating a Total Error Spectra to fit and extrapolating at low pT
    TH1F      *gTotal      =   new TH1F(*(fSumErrors(gStatistics,gSystematics)));
    //
    //  Fitting a first time to evaluate integral in non-measured region
    gTotal                              ->  Fit(fLevyTsallis,"EMRQS","",fMinFitInt,fMaxFitInt);
    gTotal                              ->  Fit(fLevyTsallis,"EMRQSI","",fMinFitInt,fMaxFitInt);
    fResult[0]                          =   fLevyTsallis  ->Integral(0.,fMinPT1D);
    //
    fResult[5]                          =   fMeasureMeanPT(fLevyTsallis,gTotal);
    //
    //  TCanvas w/ options
    TCanvas                *cDrawResult =   new TCanvas(Form("%s_%s",gStatistics->GetName(),fName.Data()),"",1600,800);
    cDrawResult->Divide(2,1);
    gStyle  ->SetOptStat(0);
    gTotal->SetMarkerStyle(26);
    gTotal->SetMarkerColor(38);
    gTotal->SetLineColor(kBlack);
    //
    //  Draw Spectra w/ function
    cDrawResult->cd(1);
    gPad    ->SetLogy();
    gTotal->Draw();
    fLevyTsallis->Draw("same");
    //
    //  Write info on the Extrapolation:
    TLatex  *fText   =   new TLatex();
    fText   ->  DrawLatexNDC(0.5,0.825,Form("#frac{dN}{dy} in [%.1f;%.1f]: %.7f",0.,fMinPT1D,fResult[0]));
    fText   ->  DrawLatexNDC(0.5,0.750,Form("Function: %s",fLevyTsallis->GetName()));
    fText   ->  DrawLatexNDC(0.5,0.700,Form("Fit range: [%.1f;%.1f]",fMinFitInt,fMaxFitInt));
    //
    cDrawResult->cd(2);
    gPad    ->SetLogy();
    auto gTotalCloseUp = (TH1F*)gTotal->Clone();
    gTotalCloseUp->GetXaxis()->SetRangeUser(0.01,2.5);
    gTotalCloseUp->Draw();
    fLevyTsallis->Draw("same");
    //
    //  Save To File
    cDrawResult->Write();
    //
    //  Graphical Check
    if ( fIsConditional )   cDrawResult->   SaveAs(Form("%s/2D/EvaluationFit_%s.pdf",fFolder.Data(),fName.Data()));
    else                    cDrawResult->   SaveAs(Form("%s/1D/EvaluationFit_%s.pdf",fFolder.Data(),fName.Data()));
    delete cDrawResult;
    //
    fRandomGen->SetSeed(1);
    auto fStatResults                   =   fEvaluateError(fIsConditional,gSystematics,gStatistics,fIntegral,fMinFitInt,fMaxFitInt,TString("Stat_")+fName,fFolder);
    fRandomGen->SetSeed(1);
    auto fSystResults                   =   fEvaluateError(fIsConditional,gStatistics,gSystematics,fIntegral,fMinFitInt,fMaxFitInt,TString("Syst_")+fName,fFolder);
    //
    auto fSystResuFit                   =   fEvaluateError(fIsConditional,gStatistics,gSystematics,fSystFitFunctions,fMinFitInt,fMaxFitInt,TString("SFit_")+fName,fFolder);
    //
    fResult[1] = fResult[2] = fStatResults.at(0);
    fResult[3] = fResult[4] = SquareSum( { fSystResults.at(0), fSystResuFit.at(0) } );
    fResult[6] = fResult[7] = fStatResults.at(2);
    fResult[8] = fResult[9] = SquareSum( { fSystResults.at(2), fSystResuFit.at(2) } );
    //_____________________________________
    //
    //  End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
Double_t*           fExtrapolateModel               ( std::vector<TH1F*>  gStatistics, std::vector<TH1F*>  gSystematics, Double_t fIntegral = 0.032, TString fName = "ExtrapolateSignal", Double_t fMaximumFitRange = fMaxPT1D, Double_t fMinimumFitRange = fMinPT1D, TString fFolder = ""  )    {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    for ( int i(0); i < 10; i++ ) fResult[i] = 0;
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
Double_t*           fIntegrateModel                 ( std::vector<TH1F*>  gStatistics, std::vector<TH1F*>  gSystematics, TString fName = "IntegrateSignal", TString fFolder = ""  )     {
    //  Optimisation mode
    gROOT->SetBatch(true);
    //
    //  Result format: Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult     =   new Double_t    [10];
    //
    for ( int i(0); i < 10; i++ ) fResult[i] = 0;
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
    bool fIsConditional = fName.Contains("2D");
    //
    auto        fIntegralResults    =   fIntegrateModel     (gStatistics,gSystematics,fName);
    auto        fExtrapolResults    =   fExtrapolateModel   (fIsConditional,gStatistics,gSystematics,fIntegralResults[0],fName);
    //
    //  Mean Value of Result
    fResult[0]  =   fIntegralResults[0] +   fExtrapolResults[0];
    //
    //  Statistical Error of Result
    fResult[1]  =   fIntegralResults[1] +   fExtrapolResults[1];
    fResult[2]  =   fIntegralResults[2] +   fExtrapolResults[2];
    //
    //  Systematical Error of Result
    fResult[3]  =   fIntegralResults[3] +   fExtrapolResults[4];
    fResult[4]  =   fIntegralResults[3] +   fExtrapolResults[4];
    //
    //  Mean Value of pT
    fResult[5]  =   fExtrapolResults[5]/fResult[0];
    //
    //  Statistical Error of Result
    fResult[6]  =   fResult[5]*( fExtrapolResults[6]/fExtrapolResults[5] + fResult[1]/fResult[0] );
    fResult[7]  =   fResult[5]*( fExtrapolResults[7]/fExtrapolResults[5] + fResult[2]/fResult[0] );
    //
    //  Systematical Error of Result
    fResult[8]  =   fResult[5]*( fExtrapolResults[8]/fExtrapolResults[5] + fResult[3]/fResult[0] );
    fResult[9]  =   fResult[5]*( fExtrapolResults[9]/fExtrapolResults[5] + fResult[4]/fResult[0] );
    //
    // End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
Double_t*           fMeasureFullYield               ( TH1F* gStatistics, TH1F* gSystematics, TString fName = "MeasureFullYield", TString fFolder = ""  )     {
    // Optimisation mode
    gROOT->SetBatch(true);
    //
    // Result format:  Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult             =   new Double_t        [10];
    //
    bool fIsConditional = fName.Contains("2D");
    //
    auto        fIntegralStat       =   0.;
    auto        fIntegralSyst       =   0.;
    auto        fIntegral           =   gStatistics->IntegralAndError(-1,1000,fIntegralStat,"width");
                                        gSystematics->IntegralAndError(-1,1000,fIntegralSyst,"width");
    auto        fExtrapolResults    =   fExtrapolateModel   (fIsConditional,gStatistics,gSystematics,fIntegral,fName,fMinPT1D,fMaxPT1D,fFolder);
    //
    //  Mean Value of Result
    fResult[0]  =   fIntegral +   fExtrapolResults[0];
    //
    //  Statistical Error of Result
    fResult[1]  =   fIntegralStat +   fExtrapolResults[1];
    fResult[2]  =   fIntegralStat +   fExtrapolResults[2];
    //
    //  Systematical Error of Result
    fResult[3]  =   fIntegralSyst  +   fExtrapolResults[4];
    fResult[4]  =   fIntegralSyst  +   fExtrapolResults[4];
    //
    //  Mean Value of pT
    fResult[5]  =   fExtrapolResults[5]/fResult[0];
    //
    //  Statistical Error of Result
    fResult[6]  =   fResult[5]*( fExtrapolResults[6]/fExtrapolResults[5] + fResult[1]/fResult[0] );
    fResult[7]  =   fResult[5]*( fExtrapolResults[7]/fExtrapolResults[5] + fResult[2]/fResult[0] );
    //
    //  Systematical Error of Result
    fResult[8]  =   fResult[5]*( fExtrapolResults[8]/fExtrapolResults[5] + fResult[3]/fResult[0] );
    fResult[9]  =   fResult[5]*( fExtrapolResults[9]/fExtrapolResults[5] + fResult[4]/fResult[0] );
    //
    // End Optimisation mode
    gROOT->SetBatch(false);
    //
    return fResult;
}
Double_t*           fMeasureFullYield               ( std::vector<TH1F*>  gStatistics, std::vector<TH1F*>  gSystematics, TString fName = "MeasureFullYield", TString fFolder = "" )     {
    // Optimisation mode
    gROOT->SetBatch(true);
    //
    // Result format:  Integral, Stat err low, Stat err high, Syst err low, syst err high, Mean pT, Stat err low, Stat err high, Syst err low, syst err high
    Double_t   *fResult             =   new Double_t        [10];
    //
    auto        fIntegralResults    =   fIntegrateModel     (gStatistics,gSystematics,fName,fFolder);
    auto        fExtrapolResults    =   fExtrapolateModel   (gStatistics,gSystematics,fIntegralResults[0],fName,fMinPT1D,fMaxPT1D,fFolder);
    //
    //  Mean Value of Result
    fResult[0]  =   fIntegralResults[0] +   fExtrapolResults[0];
    //
    //  Statistical Error of Result
    fResult[1]  =   fIntegralResults[1] +   fExtrapolResults[1];
    fResult[2]  =   fIntegralResults[2] +   fExtrapolResults[2];
    //
    //  Systematical Error of Result
    fResult[3]  =   fIntegralResults[3] +   fExtrapolResults[4];
    fResult[4]  =   fIntegralResults[3] +   fExtrapolResults[4];
    //
    //  Mean Value of pT
    fResult[5]  =   fExtrapolResults[5]/fResult[0];
    //
    //  Statistical Error of Result
    fResult[6]  =   fResult[5]*( fExtrapolResults[6]/fExtrapolResults[5] + fResult[1]/fResult[0] );
    fResult[7]  =   fResult[5]*( fExtrapolResults[7]/fExtrapolResults[5] + fResult[2]/fResult[0] );
    //
    //  Systematical Error of Result
    fResult[8]  =   fResult[5]*( fExtrapolResults[8]/fExtrapolResults[5] + fResult[3]/fResult[0] );
    fResult[9]  =   fResult[5]*( fExtrapolResults[9]/fExtrapolResults[5] + fResult[4]/fResult[0] );
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
//__________________________________________________________________________
//
void        fSetPhiCandidate                    ( TTree* TPhiCandidate, Struct_PhiEfficiency &evPhiEfficiency )    {
    TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiEfficiency.EventMask);
    TPhiCandidate-> SetBranchAddress    ("TrueEventMask",   &evPhiEfficiency.TrueEventMask);
    TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiEfficiency.Multiplicity);
    TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiEfficiency.nPhi);
    TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiEfficiency.Px);
    TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiEfficiency.Py);
    TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiEfficiency.Pz);
    TPhiCandidate-> SetBranchAddress    ("Selection",       &evPhiEfficiency.Selection);
}
//
//_____________________________________________________________________________
//
void        fSetKaonCandidate                    ( TTree* TKaonCandidate, Struct_KaonEfficiency &evKaonEfficiency )    {
    TKaonCandidate->SetBranchAddress    ("EventMask",       &evKaonEfficiency.EventMask);
    TKaonCandidate->SetBranchAddress    ("TrueEventMask",   &evKaonEfficiency.TrueEventMask);
    TKaonCandidate->SetBranchAddress    ("Multiplicity",    &evKaonEfficiency.Multiplicity);
    TKaonCandidate->SetBranchAddress    ("nKaon",           &evKaonEfficiency.nKaon);
    TKaonCandidate->SetBranchAddress    ("Px",              &evKaonEfficiency.Px);
    TKaonCandidate->SetBranchAddress    ("Py",              &evKaonEfficiency.Py);
    TKaonCandidate->SetBranchAddress    ("Pz",              &evKaonEfficiency.Pz);
    TKaonCandidate->SetBranchAddress    ("Selection",       &evKaonEfficiency.Selection);
}
//
//_____________________________________________________________________________
//
Bool_t      fSetCandidates                      ( TTree* TPhiCnd, Struct_PhiEfficiency &fTargetPhi, TTree* TKaonCn, Struct_KaonEfficiency &fTargetKaon )    {
    if ( !TPhiCnd && !TKaonCn ) return false;
    if ( TPhiCnd )  fSetPhiCandidate(TPhiCnd,fTargetPhi);
    if ( TKaonCn )  fSetKaonCandidate(TKaonCn,fTargetKaon);
    return true;
}
//
//_____________________________________________________________________________
//
void        fSetPhiCandidate                    ( TTree* TPhiCandidate, Struct_PhiCandidate &evPhiCandidate )    {
    TPhiCandidate-> SetBranchAddress    ("EventMask",       &evPhiCandidate.EventMask);
    TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
    TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
    TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
    TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
    TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
    TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
    TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
    TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
    TPhiCandidate-> SetBranchAddress    ("TrueInvMass",     &evPhiCandidate.TrueInvMass);
}
//
//_____________________________________________________________________________
//
void        fSetKaonCandidate                    ( TTree* TKaonCandidate, Struct_KaonCandidate &evKaonCandidate )    {
    TKaonCandidate-> SetBranchAddress   ("EventMask",       &evKaonCandidate.EventMask);
    TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
    TKaonCandidate-> SetBranchAddress   ("nKaon",           &evKaonCandidate.nKaon);
    TKaonCandidate-> SetBranchAddress   ("Px",              &evKaonCandidate.Px);
    TKaonCandidate-> SetBranchAddress   ("Py",              &evKaonCandidate.Py);
    TKaonCandidate-> SetBranchAddress   ("Pz",              &evKaonCandidate.Pz);
    TKaonCandidate-> SetBranchAddress   ("Charge",          &evKaonCandidate.Charge);
    TKaonCandidate-> SetBranchAddress   ("TOFSigma",        &evKaonCandidate.SigmaTOF);
    TKaonCandidate-> SetBranchAddress   ("TPCSigma",        &evKaonCandidate.SigmaTPC);
}
//
//_____________________________________________________________________________
//
Bool_t      fSetCandidates                      ( TTree* TPhiCnd, Struct_PhiCandidate &fTargetPhi, TTree* TKaonCn, Struct_KaonCandidate &fTargetKaon )    {
    if ( !TPhiCnd && !TKaonCn ) return false;
    if ( TPhiCnd )  fSetPhiCandidate(TPhiCnd,fTargetPhi);
    if ( TKaonCn )  fSetKaonCandidate(TKaonCn,fTargetKaon);
    return true;
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
    result = fModel.chi2FitTo(*RData,Save(),Minos(true),NumCPU(kCPU_use));
    result = fModel.fitTo(*RData,Save(),Minos(true),NumCPU(kCPU_use));
    
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
Double_t*       fExtrapolateModel               ( Tclass *THdata, TString fName = "ExtrapolateSignal", TF1* fFitFunc = fLevyTsallis) {
    // Optimisation mode
    gROOT->SetBatch(true);
    
    // Result format: Integral, Stat err, Syst err, Mean pT, Stat err, Syst err
    Double_t   *fResult = new   Double_t    [6];
    
    // Setting -1. for the default
    for ( Int_t iFill = 0; iFill < 6; iFill++ ) fResult[iFill]  =   -1.;
    
    // Set Standard Fit
    fSetFunction(fFitFunc);
    
    // Fit the Spectra
    THdata->Fit(fFitFunc,"IMREQ0S","",fMinPT1D,10.);
    
    // Save to further checks
    TCanvas * fCheckFit = new TCanvas();
    gPad->SetLogy();
    THdata      ->Draw("same");
    fFitFunc  ->Draw("same");
    fCheckFit   ->Write();
    fCheckFit   ->SaveAs(Form("./yield/ExtrapolateCheck/%s.pdf",fName.Data()));
    
    fResult[0]  =   fFitFunc->Integral(0.,fMinPT1D);
    fResult[1]  =   fFitFunc->IntegralError(0.,fMinPT1D);
    
    // End Optimisation mode
    gROOT->SetBatch(false);
    
    return fResult;
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
        fFitRsltBkg1    =   fBkg1_Model->fitTo(*bkg1Loose,Extended(kTRUE),Save(),NumCPU(kCPU_use));
        fFitRsltBkg1    =   fBkg1_Model->fitTo(*bkg1,Extended(kTRUE),Save(),NumCPU(kCPU_use));
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
        fFitRsltBkg2    =   fBkg2_Model->fitTo(*bkg2Loose,  Extended(kTRUE),      Save(), NumCPU(kCPU_use));
        fFitRsltBkg2    =   fBkg2_Model->fitTo(*bkg2,       Extended(kTRUE),      Save(), NumCPU(kCPU_use));
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
        fFitRsltSign    =   fSign_Model->fitTo(*dataLoose,  Extended(kTRUE),      Save(), NumCPU(kCPU_use));
        fFitRsltSign    =   fSign_Model->fitTo(*data,       Extended(kTRUE),      Save(), NumCPU(kCPU_use));
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
}


Float_t             fMeasureMeanPT                  ( TF1 * fLowFit, TGraphAsymmErrors * gTotal, Bool_t fReFit = false )   {
    Float_t fResult = 0.;
    if ( fReFit )   {
        gTotal->Fit(fLowFit);
    }
    std::vector<float>  fIntegral;
    std::vector<float>  fMeanPT;
    std::vector<float>  fWidth;
    auto    nEntries        =   gTotal->GetN();
    fIntegral.push_back(fLowFit->Integral(0.,gTotal->GetPointX(0) - gTotal->GetErrorXlow(0)));
    fMeanPT.push_back(  fLowFit->Moment(1,0.,gTotal->GetPointX(0) - gTotal->GetErrorXlow(0)));
    fWidth.push_back(   1.);
    for ( Int_t iPnt = 0; iPnt < nEntries; iPnt++ )   {
        auto    fXLow   =   gTotal->GetPointX(iPnt) - gTotal->GetErrorXlow(iPnt);
        auto    fXHig   =   gTotal->GetPointX(iPnt) + gTotal->GetErrorXhigh(iPnt);
        fWidth.push_back(   fXHig - fXLow );
        fIntegral.push_back(gTotal->GetPointY(iPnt));
        fMeanPT.push_back(  fLowFit->Moment(1,fXLow,fXHig));
    }
    fIntegral.push_back(fLowFit->Integral(gTotal->GetPointX(nEntries-1) + gTotal->GetErrorXhigh(nEntries-1),1.e2));
    fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetPointX(nEntries-1) + gTotal->GetErrorXhigh(nEntries-1),1.e2));
    fWidth.push_back(   1.);
    for ( Int_t iTer = 0; iTer < fIntegral.size(); iTer++ ) {
        fResult +=  fIntegral.at(iTer)*fMeanPT.at(iTer)*fWidth.at(iTer);
    }
    return fResult;
    /*
    
    //--------------------
    //  !TODO: Scorporare la funzione che calcola mean pT
    //  Evaluating the Mean pT
    if ( !fIsConditional )  {
        fResult[5]                          =   (fMinPT1D)*fFitFunc->Moment(1,0.,fMinPT1D)*(fResult[0]);
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            fResult[5]                     +=   (fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D])*(fFitFunc->Moment(1,fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]))*(gTotal->GetPointY(iPT1D));
        }
    }   else    {
        fResult[5]                          =   (fMinPT2D)*fFitFunc->Moment(1,0.,fMinPT2D)*(fResult[0]);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            fResult[5]                     +=   (fArrPT2D[iPT2D+1] - fArrPT2D[iPT2D])*(fFitFunc->Moment(1,fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]))*(gTotal->GetPointY(iPT2D));
        }
    }
    //--------------------
     */
}
Float_t             fMeasureMeanPT                  ( TF1 * fLowFit, TH1F * gTotal, Bool_t fReFit = false )   {
    Float_t fResult = 0.;
    if ( fReFit )   {
        gTotal->Fit(fLowFit);
    }
    std::vector<float>  fIntegral;
    std::vector<float>  fMeanPT;
    std::vector<float>  fWidth;
    auto    nEntries        =   gTotal->GetNbinsX();
    fIntegral.push_back(fLowFit->Integral(0.,gTotal->GetBinLowEdge(1)));
    fMeanPT.push_back(  fLowFit->Moment(1,0.,gTotal->GetBinLowEdge(1)));
    fWidth.push_back(   1.);
    for ( Int_t iBin = 1; iBin <= nEntries; iBin++ )   {
        fWidth.push_back(   gTotal->GetBinWidth(iBin) );
        fIntegral.push_back(gTotal->GetBinContent(iBin));
        fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetBinLowEdge(iBin),gTotal->GetBinLowEdge(iBin+1)));
    }
    fIntegral.push_back(fLowFit->Integral(gTotal->GetBinLowEdge(nEntries+1),1.e2));
    fMeanPT.push_back(  fLowFit->Moment(1,gTotal->GetBinLowEdge(nEntries+1),1.e2));
    fWidth.push_back(   1.);
    for ( Int_t iTer = 0; iTer < fIntegral.size(); iTer++ ) {
        fResult +=  fIntegral.at(iTer)*fMeanPT.at(iTer)*fWidth.at(iTer);
    }
    return fResult;
    /*
    
    //--------------------
    //  !TODO: Scorporare la funzione che calcola mean pT
    //  Evaluating the Mean pT
    if ( !fIsConditional )  {
        fResult[5]                          =   (fMinPT1D)*fFitFunc->Moment(1,0.,fMinPT1D)*(fResult[0]);
        for ( Int_t iPT1D = 0; iPT1D < nBinPT1D; iPT1D++ )  {
            fResult[5]                     +=   (fArrPT1D[iPT1D+1] - fArrPT1D[iPT1D])*(fFitFunc->Moment(1,fArrPT1D[iPT1D],fArrPT1D[iPT1D+1]))*(gTotal->GetPointY(iPT1D));
        }
    }   else    {
        fResult[5]                          =   (fMinPT2D)*fFitFunc->Moment(1,0.,fMinPT2D)*(fResult[0]);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
            fResult[5]                     +=   (fArrPT2D[iPT2D+1] - fArrPT2D[iPT2D])*(fFitFunc->Moment(1,fArrPT2D[iPT2D],fArrPT2D[iPT2D+1]))*(gTotal->GetPointY(iPT2D));
        }
    }
    //--------------------
     */
}

Double_t            fIntegralOverHalf               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    for ( Int_t iTer = 0; iTer < hTarget.size(); iTer++ )   {
        for ( Int_t jTer = iTer+1; jTer < hTarget.size(); jTer++ )   {
            fResult             +=  hTarget.at(iTer)->GetBinContent(jTer);
            fError              +=  hTarget.at(iTer)->GetBinError(jTer)*hTarget.at(iTer)->GetBinError(jTer);
        }
    }
    sqrt(fError);
    return fResult;
}
Double_t            fIntegralDiagonal               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    for ( Int_t iTer = 0; iTer < hTarget.size(); iTer++ )   {
        fResult             +=  hTarget.at(iTer)->GetBinContent(iTer);
        fError              +=  hTarget.at(iTer)->GetBinError(iTer)*hTarget.at(iTer)->GetBinError(iTer);
    }
    sqrt(fError);
    return fResult;
}
Double_t            fIntegralDiagHalf               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    fResult +=  fIntegralOverHalf( hTarget, fDumpVar_Double_t );
    fError  +=  fDumpVar_Double_t*fDumpVar_Double_t;
    fResult +=  fIntegralDiagonal( hTarget, fDumpVar_Double_t );
    fError  +=  fDumpVar_Double_t*fDumpVar_Double_t;
    sqrt(fError);
    return  fResult;
}
Double_t            fIntegralCorrHalf               ( std::vector<TH1F*>  hTarget, Double_t &fError = fDumpVar_Double_t  ) {
    Double_t    fResult =   0.;
    fResult +=  2*fIntegralOverHalf( hTarget, fDumpVar_Double_t );
    fError  +=  4*fDumpVar_Double_t*fDumpVar_Double_t;
    fResult +=  fIntegralDiagonal( hTarget, fDumpVar_Double_t );
    fError  +=  fDumpVar_Double_t*fDumpVar_Double_t;
    sqrt(fError);
    return  fResult;
}

std::pair<std::vector<float>,TH1F*>     uIntegralError                  ( std::vector<TH1F*>  hTargetStat, std::vector<TH1F*>  hTargetSyst )  {
    //
    
    // auto fError = ; 
    
    //
    //  Result Container
    std::pair<std::vector<float>,TH1F*> fResult;
    //
    //  Generate utility histogram to evaluate error
    std::vector<float>      uIntegralVar;
    std::vector<float>*     uIntegralVarCond = new std::vector<float> [hTargetStat.size()];
    std::vector<float>      uExtrapolVar;
    std::vector<float>      uExtrapolVa2;
    std::vector<float>*     uExtrapolVarCond = new std::vector<float> [hTargetStat.size()];
    //
    //  Evaluate Central Value
    TH1F   *hConditionalYieldStat   =   new TH1F("hConditionalYieldStat","hConditionalYieldStat",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldSyst   =   new TH1F("hConditionalYieldSyst","hConditionalYieldSyst",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldINT_   =   new TH1F("hConditionalYieldINT_","hConditionalYieldINT_",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldEXT_   =   new TH1F("hConditionalYieldEXT_","hConditionalYieldEXT_",nBinPT2D,fArrPT2D);
    TH1F   *hConditionalYieldTest   =   new TH1F("hConditionalYieldTest","hConditionalYieldTest",nBinPT2D,fArrPT2D);
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        fSetAllFunctions();
        fSetFunction(fLevyTsallis);
        auto    fNewTarget  =   fSumErrors(hTargetStat.at(iPT2D),hTargetSyst.at(iPT2D));
        fNewTarget          ->  Fit(fLevyTsallis,   "EMRQS","",     fMinPT1D,fMaxPT1D);
        fNewTarget          ->  Fit(fLevyTsallis,   "EMRQSI","",    fMinPT1D,fMaxPT1D);
        //
        TCanvas*cTest = new TCanvas();
        fNewTarget->DrawCopy();
        fLevyTsallis->DrawCopy("same");
        cTest->SaveAs(Form("./result/Yield/ExtrapolateCheck/ttt_%i.pdf",iPT2D));
        delete cTest;
        //
        //  Integral
        auto    fIntegral   =   fNewTarget->Integral("width");
        //
        //  Extrapolation
        auto    fExtraplt   =   fLevyTsallis->Integral(0.,fMinPT1D);
        auto    fExtraplE   =   fLevyTsallis->IntegralError(0.,fMinPT1D);
        //
        //  Set Values
        hConditionalYieldStat->SetBinContent(iPT2D+1,fIntegral+fExtraplt);
        hConditionalYieldSyst->SetBinContent(iPT2D+1,fIntegral+fExtraplt);
        hConditionalYieldINT_->SetBinContent(iPT2D+1,fIntegral);
        hConditionalYieldEXT_->SetBinContent(iPT2D+1,fExtraplt);
        hConditionalYieldEXT_->SetBinError  (iPT2D+1,fExtraplE);
        //
        //  Set test
        hConditionalYieldTest->SetBinContent(iPT2D+1,fNewTarget->GetBinContent  (iPT2D+1));
        hConditionalYieldTest->SetBinError  (iPT2D+1,fNewTarget->GetBinError    (iPT2D+1));
    }
    //
    //  Save to File
    hConditionalYieldStat->Write();
    hConditionalYieldSyst->Write();
    hConditionalYieldINT_->Write();
    hConditionalYieldEXT_->Write();
    hConditionalYieldTest->Write();
    //
    //  Evaluate Variations
    auto cCheckMultiFit = new TCanvas();
    gStyle->SetOptStat(0);
    cCheckMultiFit->Divide(4,3);
    TH1F   *hConditionalUtil   =   new TH1F("hConditionalUtil","hConditionalUtil",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 0; iTer < kStatEvalCycles; iTer++ )  {
        //
        //  Randomise points Symmetrically
        auto fNewTarget = uRandomizePointsSymm(hTargetStat,hTargetSyst);
        //
        //  Loop over Randomised Conditional Spectra
        auto        jTer        =   0;
        Double_t    fIntegral   =   0.;
        Double_t    fExtraplt   =   0.;
        for ( auto hHistogram : fNewTarget )    {
            //
            //  Recover Integral of current Spectra & save it
            auto fCurrentIntegral = hHistogram->Integral("width");
            uIntegralVarCond[jTer].push_back(fCurrentIntegral);
            //
            //  Cumulative Integral
            fIntegral   +=  fCurrentIntegral*(fArrPT2D[jTer+1]-fArrPT2D[jTer]);
            //
            //  Prepping the Fit
            fSetAllFunctions();
            fSetFunction(fLevyTsallis);
            hHistogram->Fit(fLevyTsallis,"IQMRE0S");
            auto    fCurrentExtraplt    =   fLevyTsallis->Integral(0.,fMinPT1D);
            uExtrapolVarCond[jTer].push_back(fCurrentExtraplt);
            //
            //  Cumulative Extrapolation
            fExtraplt   +=  fCurrentExtraplt;
            //
            //  Make the Conditional Yield Histogram
            hConditionalUtil->SetBinContent (jTer+1, fCurrentIntegral+fCurrentExtraplt);
            //hConditionalUtil->SetBinError   (jTer+1, fError[jTer]);
            hConditionalYieldTest->SetBinContent( jTer+1, hHistogram->GetBinContent(jTer));
            //
            //  Graphics
            cCheckMultiFit->cd( jTer + 1 );
            gPad->SetLogy();
            if ( iTer == 0 )                    fSumErrors(hTargetStat.at(jTer),hTargetSyst.at(jTer))->Draw("SAME");
            fLevyTsallis->DrawCopy("same");
            if ( iTer == kStatEvalCycles-1 )    fSumErrors(hTargetStat.at(jTer),hTargetSyst.at(jTer))->Draw("SAME");
            //
            jTer++;
        }
        //
        uIntegralVar.push_back( fIntegral );
        //
        fSetAllFunctions();
        fSetFunction(fLevyTsallis);
        hConditionalUtil->Fit(fLevyTsallis,"IQMRE0S");
        cCheckMultiFit->cd( 11 );
        gPad->SetLogy();
        if ( iTer == 0 )                    hConditionalUtil->Draw("SAME");
        fLevyTsallis->DrawCopy("same");
        if ( iTer == kStatEvalCycles-1 )    hConditionalUtil->Draw("SAME");
        //
        uExtrapolVar.push_back( fLevyTsallis->Integral(0.,fMinPT1D) );
        //
        fSetAllFunctions();
        fSetFunction(fLevyTsallis);
        hConditionalYieldTest->Fit(fLevyTsallis,"IQMRE0S","",fMinPT1D,2.);
        cCheckMultiFit->cd( 12 );
        gPad->SetLogy();
        if ( iTer == 0 )                    hConditionalYieldTest->Draw("SAME");
        fLevyTsallis->DrawCopy("same");
        if ( iTer == kStatEvalCycles-1 )    hConditionalYieldTest->Draw("SAME");
        //
        uExtrapolVa2.push_back( fLevyTsallis->Integral(0.,fMinPT1D) );
    }
    cCheckMultiFit->SaveAs("./result/Yield/ExtrapolateCheck/cCheckMultiFit.pdf");
    delete cCheckMultiFit;
    //
    //  Analysis of Integral Variations
    auto cCheckMultiIntegral    = new TCanvas();
    gStyle->SetOptStat(0);
    cCheckMultiIntegral->Divide(4,3);
    //  First-Differential Spectrum
    cCheckMultiIntegral->cd(1);
    //  Generating the Utility Histogram
    auto    hIntegral       =   uBuildTH1F( uIntegralVar, 50, -hConditionalYieldINT_->Integral("width") );
    hIntegral               ->  SetNameTitle( "hIntegral", "dN/dy variations" );
    fSetAllFunctions();
    fSetFunction(fGauss);
    fGauss->SetParameter(0,hIntegral->Integral());
    fGauss->FixParameter(1,0);
    fGauss->SetParameter(2,hIntegral->GetRMS());
    hIntegral               ->  Fit(fGauss,"IMRE0S");
    
    hIntegral               ->  Draw();
    hIntegral               ->  Write();
    fGauss                  ->  DrawCopy("SAME");
    fGauss                  ->  Write();
    
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cCheckMultiIntegral->cd(iTer+2);
        
        auto    hIntegralCond       =   uBuildTH1F( uIntegralVarCond[iTer], 50, -hConditionalYieldINT_->GetBinContent(iTer+1) );
        hIntegralCond               ->  SetNameTitle( Form("hIntegral_%i",iTer), "dN/dy variations" );
        hIntegralCond               ->  Fit(fGauss,"IMRE0S");
        
        hIntegralCond               ->  Draw();
        hIntegralCond               ->  Write();
        fGauss                      ->  DrawCopy("SAME");
        fGauss                      ->  Write();
    }
    //
    cCheckMultiIntegral->SaveAs(Form("result/yield/Integral.pdf"));
    delete cCheckMultiIntegral;
    //
    //
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    hConditionalYieldStat->Fit(fLevyTsallis,"IQMRE0S");
    //
    //  Analysis of Extrapolation Variations
    auto cCheckMultiExtrap    = new TCanvas();
    gStyle->SetOptStat(0);
    cCheckMultiExtrap->Divide(4,3);
    //  Second-Differential Spectrum
    cCheckMultiExtrap->cd(1);
    //  Generating the Utility Histogram
    auto    hExtrap         =   uBuildTH1F( uExtrapolVar, 50, -fLevyTsallis->Integral(0.,fMinPT1D) );
    hExtrap                 ->  SetNameTitle( "hExtrap", "dN/dy variations" );
    fSetAllFunctions();
    fSetFunction(fGauss);
    fGauss->SetParameter(0,hExtrap->Integral());
    fGauss->FixParameter(1,0);
    fGauss->SetParameter(2,hExtrap->GetRMS());
    hExtrap                 ->  Fit(fGauss,"IMRE0S");
    
    hExtrap                 ->  Draw();
    hExtrap                 ->  Write();
    fGauss                  ->  DrawCopy("SAME");
    fGauss                  ->  Write();
    
    cCheckMultiExtrap->cd(12);
    //  Generating the Utility Histogram
    auto    hExtra2         =   uBuildTH1F( uExtrapolVa2, 50, -fLevyTsallis->Integral(0.,fMinPT1D) );
    hExtra2                 ->  SetNameTitle( "hExtrap", "dN/dy variations" );
    fSetAllFunctions();
    fSetFunction(fGauss);
    fGauss->SetParameter(0,hExtra2->Integral());
    fGauss->FixParameter(1,0);
    fGauss->SetParameter(2,hExtra2->GetRMS());
    hExtrap                 ->  Fit(fGauss,"IMRE0S");
    
    hExtra2                 ->  Draw();
    hExtra2                 ->  Write();
    fGauss                  ->  DrawCopy("SAME");
    fGauss                  ->  Write();
    
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cCheckMultiExtrap->cd(iTer+2);
        
        auto    hExtrapCond         =   uBuildTH1F( uExtrapolVarCond[iTer], 50, -hConditionalYieldEXT_->GetBinContent(iTer+1) );
        hExtrapCond                 ->  SetNameTitle( Form("hExtrap_%i",iTer), "dN/dy variations" );
        hExtrapCond                 ->  Fit(fGauss,"IMRE0S");
        
        hExtrapCond                 ->  Draw();
        hExtrapCond                 ->  Write();
        fGauss                      ->  DrawCopy("SAME");
        fGauss                      ->  Write();
    }
    //
    cCheckMultiExtrap->SaveAs(Form("result/yield/Extrap.pdf"));
    delete cCheckMultiExtrap;
    
    
    
    /*
    // RESULT CHECK
    
    auto cDrawResult = new TCanvas();
    cDrawResult->Divide(4,3);

    cDrawResult->cd(1);
    auto    hStatIntegral   =   uBuildTH1F( uIntegralVar, 50, -hConditionalYieldINT_->Integral("width") );
    hStatIntegral           ->  SetNameTitle( "  ssdd", "dN/dy variations" );
    auto fErr = 0.;
    auto    jTer    =   0;
    for ( auto hIntegrand : hTargetStat )    {
        auto gfs = 0.;
        hIntegrand->IntegralAndError(-1,10000,gfs,"width");
        fErr += gfs*gfs;
    }
    hStatIntegral->Fit("gaus","IMREQ0S");
    auto fResultsss = hStatIntegral->GetFunction("gaus");
    auto siggg = fResultsss->GetParameter(2);
    hStatIntegral->Draw();
    fResultsss->Draw("same");
    //Latx->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*fIntegralStat));
    //Latx->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*sqrt(fErr)));
    Latx->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
    Latx->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hStatIntegral->GetRMS()));
    
    hStatIntegral->SetName("NEW_-1");
    hStatIntegral->Write();
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cDrawResult->cd(iTer+2);
        uCleanOutsiders(  uIntegralVarCond[iTer] );
        auto    hConditional   =   uBuildTH1F( uIntegralVarCond[iTer], 50, -hConditionalYieldINT_->GetBinContent(iTer+1) );
        hConditional           ->  SetNameTitle( Form("NEW_%i",iTer), "dN/dy variations" );
        hConditional->Fit("gaus","IMREQ0S");
        hConditional->Draw();
        fResultsss = hConditional->GetFunction("gaus");
        siggg = fResultsss->GetParameter(2);
        fResultsss->Draw("same");
        //Latx->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*hRES_2D_Cond2_Stat_INT->GetBinError(iTer+1)));
        Latx->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
        Latx->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hConditional->GetRMS()));
        hConditional->Write();
    }
    
    cDrawResult->SaveAs(Form("result/yield/Integral.pdf"));
    delete cDrawResult;
    
    cDrawResult = new TCanvas();
    cDrawResult->Divide(4,3);
    cDrawResult->cd(1);
    
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    hConditionalYieldStat->Fit(fLevyTsallis,"IQMRE0S");
    
    auto hStatExtrapol   =   uBuildTH1F( uExtrapolVar, 50, -fLevyTsallis->Integral(0.,fMinPT1D) );
    hStatExtrapol           ->  SetNameTitle( "  hStatExtrapol", "dN/dy variations" );
    hStatExtrapol->Fit("gaus","IMRQE0S");
    hStatExtrapol->Draw();
    auto fResulFits = hStatExtrapol->GetFunction("gaus");
    siggg = fResulFits->GetParameter(2);
    fResulFits->DrawCopy("same");
    //Latx->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*fExtrapolResults[1]));
    Latx->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
    Latx->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hStatIntegral->GetRMS()));
    
    hStatIntegral->SetName("ENEW_-1");
    hStatIntegral->Write();
    for ( Int_t iTer = 0; iTer < hTargetStat.size(); iTer++ )    {
        cDrawResult->cd(iTer+2);
        uCleanOutsiders(  uExtrapolVarCond[iTer] );
        auto    hConditional   =   uBuildTH1F( uExtrapolVarCond[iTer], 50, -hConditionalYieldEXT_->GetBinContent(iTer+1) );
        hConditional           ->  SetNameTitle( Form("ENEW_%i",iTer), "dN/dy variations" );
        hConditional->Fit("gaus","IMRE0QS");
        hConditional->Draw();
        fResultsss = hConditional->GetFunction("gaus");
        siggg = fResultsss->GetParameter(2);
        fResultsss->Draw("same");
        //Latx->DrawLatexNDC(0.60,0.80,Form("OLD : %.4f",1.e6*hRES_2D_Cond2_Stat_EXT->GetBinError(iTer+1)));
        Latx->DrawLatexNDC(0.60,0.75,Form("GSS : %.4f",1.e6*siggg));
        Latx->DrawLatexNDC(0.60,0.70,Form("RMS : %.4f",1.e6*hConditional->GetRMS()));
        
        hConditional->Write();
    }
    
    cDrawResult->SaveAs(Form("result/yield/Extrapol.pdf"));
    //delete cDrawResult;
    */
    
    return fResult;
}

#endif
