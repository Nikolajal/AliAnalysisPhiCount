#include "../inc/AliAnalysisPhiPair.h"
#include "./Analysis/SignalExtraction.C"
#include "./PreProcessing.C"
// !TODO: All Set!

void TestMacro   ()  {
    
    Double_t fMin_, fMax_;
    for ( Int_t i = 0; i < nOptions; i++ )  {
        SetBoundaries(sOptions[i],fMin_,fMax_);
        cout << sOptions[i] << " min:" << fMin_ << " max:" << fMax_ << endl;
    }
    
    return;
    
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    auto nEntries = 100000;
    TH1F* THdata = new TH1F("htest","htest",1000,-20,20);
    
    for ( int i = 0; i < nEntries; i++ )  {
        THdata->Fill( uRandomGen->Gaus(0,5) );
        THdata->Fill( uRandomGen->Uniform(-20,20) );
    }
    
    RooRealVar      InvMass     =   RooRealVar        ("InvMass",   "InvMass",      -20, 20                  );
    RooRealVar      InvMassLoose=   RooRealVar        ("InvMassLoose",   "InvMassLoose",      -5, 5                  );
    
    RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         InvMass,        Import(*THdata)                 );
    RooDataHist*    dataLoose   =   new RooDataHist   ("DataLoose", "DataLoose",    InvMassLoose,        Import(*THdata)   );
    
    
    
    RooRealVar ch1, ch2;
                                ch1     =   RooRealVar      ("ch1",     "ch1",      0.,     -10,         10   );
                                ch2     =   RooRealVar      ("ch2",     "ch2",      0.,     -10,         10   );
    
    RooChebychev                fBkg1   =   RooChebychev    ("fBkg1",    "fBkg1",     InvMass,    RooArgSet(ch1));
    RooChebychev                fBkg2   =   RooChebychev    ("fBkg2",    "fBkg2",     InvMassLoose,    RooArgSet(ch2));
    
    RooRealVar mn1, st1;
    mn1     =   RooRealVar      ("mn1",     "mn1",      0.,     -10,         10   );
    st1     =   RooRealVar      ("st1",     "st1",      0.,     0,         100   );
    
    RooRealVar mn2, st2;
    mn2     =   RooRealVar      ("mn2",     "mn2",      0.,     -10,         10   );
    st2     =   RooRealVar      ("st2",     "st2",      0.,     0,         100   );
    
    RooGaussian                fSig1   =   RooGaussian    ("fSig1",    "fSig1",     InvMass,    mn1,st1);
    RooGaussian                fSig2   =   RooGaussian    ("fSig2",    "fSig2",     InvMassLoose,    mn1,st2);
    
    RooRealVar  nSS1,    nBB1, nSS2,    nBB2;
    nSS1     =   RooRealVar      ("anSS1",    "anSS1",     .5*nEntries,    0., 2*nEntries);
    nBB1     =   RooRealVar      ("anBB1",    "anBB1",     .5*nEntries,    0., 2*nEntries);
    nSS2     =   RooRealVar      ("anSS2",    "anSS2",     .5*nEntries,    0., 2*nEntries);
    nBB2     =   RooRealVar      ("anBB2",    "anBB2",     .5*nEntries,    0., 2*nEntries);
    
    
    RooAddPdf                  *fMod1 =   new RooAddPdf   ("fMod1",    "fMod1",     RooArgList(fBkg1,fSig1),   RooArgList(nBB1,nSS1));
    RooAddPdf                  *fMod2 =   new RooAddPdf   ("fMod2",    "fMod2",     RooArgList(fBkg2,fSig2),   RooArgList(nBB2,nSS2));
    
    fMod1->fitTo(*data,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
    fMod2->fitTo(*dataLoose,Extended(kTRUE),Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
    
    RooRealVar      vTest     =   RooRealVar        ("vTest",   "vTest",      -1000, 1000                  );
    RooRealVar      vMean     =   RooRealVar        ("vMean",   "vMean",      kPhiMesonMass_);
    RooRealVar      vStdv     =   RooRealVar        ("vStdv",   "vStdv",      kPhiMesonWidth);
    RooRealVar      vMea2     =   RooRealVar        ("vMea2",   "vMea2",      kPhiMesonMass_);
    RooRealVar      vStd2     =   RooRealVar        ("vStd2",   "vStd2",      kPhiMesonWidth);
    RooRealVar      vSlop     =   RooRealVar        ("vSlop",   "vSlop",      0.0);
    
    RooBreitWigner hTest = RooBreitWigner("hTest","hTest",vTest,vMean,vStdv);
    RooBreitWigner hTes2 = RooBreitWigner("hTes2","hTes2",vTest,vMea2,vStd2);
    
    vTest.setRange("Full",0,10);
    vTest.setRange("Meas",0.99,1.06);
    
    cout << hTest.analyticalIntegral(1,"Full") << endl;
    cout << hTest.analyticalIntegral(1,"Meas") << endl;
    cout << hTest.analyticalIntegral(1,"Meas")/hTest.analyticalIntegral(1,"Full") << endl;
    
    RooProdPdf hTes3 = RooProdPdf( "hTes3","hTes3",hTest,hTes2);
    
    auto hTestTF = hTes3.asTF(vTest,RooArgList(vMean,vStdv,vMea2,vStd2));
    
    cout << hTestTF->Integral(0,10) << endl;
    cout << hTestTF->Integral(0.99,1.06) << endl;
    cout << hTestTF->Integral(0.99,1.06)/(hTestTF->Integral(0,10)) << endl;
    
    TCanvas *c1 = new TCanvas();
    THdata->Draw();
    
    
    /*
    
    
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
    
    */
    
}
    
    
    /*
    auto fPhi1 = 0.;
    auto fPhi2 = 0.;
    auto fPhi11= 0.;
    auto fPhi22= 0.;
    auto nEvents_ = 1.e8;
    for ( int i = 0; i < nEvents_; i++ )    {
        auto fPhi = fRandomGen->Poisson(0.3);
        auto fPhi_= fPhi + ( fRandomGen->Uniform(0,1) > 0.95 && fPhi > 0);
        if ( fPhi > 0 ) fPhi1 += TMath::Binomial(fPhi,1);
        if ( fPhi > 1 ) fPhi2 += TMath::Binomial(fPhi,2);
        if ( fPhi_ > 0 ) fPhi11+= TMath::Binomial(fPhi_,1);
        if ( fPhi_ > 1 ) fPhi22+= TMath::Binomial(fPhi_,2);
    }
    fPhi1 /= nEvents_;
    fPhi2 /= nEvents_;
    fPhi11/= nEvents_;
    fPhi22/= nEvents_;
    cout << " 1: " << fPhi1 <<  " 2: " << fPhi2 << endl;
    cout << " 2/1 " << fPhi2/fPhi1 <<  " 2/(1*1): " << fPhi2/(fPhi1*fPhi1) << endl;
    cout << " \gamma " << 2*fPhi2/(fPhi1) - fPhi1 << endl;
    cout << " 1: " << fPhi11 <<  " 2: " << fPhi22 << endl;
    cout << " 2/1 " << fPhi22/fPhi11 <<  " 2/(1*1): " << fPhi22/(fPhi11*fPhi11) << endl;
    cout << " \gamma " << 2*fPhi22/(fPhi11) - fPhi11 << endl;
    return;
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (Form(kASigExtp_FitCheckRst,"Yield"));
    TFile*  insFile_DT_Yiel2            =   new TFile   (Form(kASigExtr_FitCheckRst,"Yield"));
    TFile*  insFile_EF_Yield            =   new TFile   (Form(kAnalysis_MCTruthHist,"Yield"));
    
    auto    h2D =   (TH2F*)(insFile_DT_Yiel2->Get("hRAW_2D"));
    auto    hEf =   (TH2F*)(insFile_EF_Yield->Get("hEFF_2D_fr_1D"));
    auto    h1D =   (TH2F*)(insFile_DT_Yiel2->Get("hRAW_1D"));
    auto    hE1 =   (TH2F*)(insFile_EF_Yield->Get("hEFF_1D"));
    
    auto    hRt =   (TH2F*)h2D->Clone();
    auto    hR1 =   (TH1F*)h1D->Clone();
    hRt->Divide(h2D,hEf);
    hR1->Divide(h1D,hE1);
    
    // >-> GENERAL ANALYSIS //
    //
    TH1D       *hEvntEff;
    TH1D       *hEvntMlt;
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yiel2->Get(hName));
    //
    hName       =   "fQC_Event_Enum_Mult";
    hEvntMlt    =   (TH1D*)(insFile_DT_Yiel2->Get(hName));
    
    // Scaling for efficiencies
    //
    auto        kN_Trg          =   (hEvntEff->GetBinContent(kEventCount::kTrigger));
    auto        kN_Vtx          =   (hEvntEff->GetBinContent(kEventCount::kVertex));
    auto        kN_MB           =   (hEvntEff->GetBinContent(kEventCount::kVertex10));
    Double_t    f1DCorrection   =   (1./kBR)        *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    
    hRt->Scale(f2DCorrection);
    hR1->Scale(f1DCorrection);
    
    TFile*  fout            =   new TFile   (Form(kAnalysis_SigExtr_Dir,"Yield")+TString("/dick.root"),"recreate");
    
    hRt->Write();
    
    fSetAllFunctions();
    fSetFunction();
    fSetFunction(fLevyTsallis2D);
    hR1->Fit(fLevyTsallis,"IM");
    fLevyTsallis2D->SetParameter(1,fLevyTsallis->GetParameter(1));
    fLevyTsallis2D->SetParameter(4,fLevyTsallis->GetParameter(1));
    fLevyTsallis2D->SetParameter(2,fLevyTsallis->GetParameter(2));
    fLevyTsallis2D->SetParameter(5,fLevyTsallis->GetParameter(2));
    fLevyTsallis2D->SetParameter(6,TMath::Power(fLevyTsallis->GetParameter(3),2));
    hRt->Fit(fLevyTsallis2D,"IM");
    
    fBreitWigner->FixParameter(0,kPhiMesonMass_);
    fBreitWigner->FixParameter(1,kPhiMesonWidth);
    fBreitWigner2D->FixParameter(0,kPhiMesonMass_);
    fBreitWigner2D->FixParameter(1,kPhiMesonWidth);
    cout <<   fBreitWigner->GetParameter(0) << endl;
    cout <<   fBreitWigner->GetParameter(1) << endl;
    cout <<   fBreitWigner2D->GetParameter(0) << endl;
    cout <<   fBreitWigner2D->GetParameter(1) << endl;
    
    auto fTotal     =   fBreitWigner->Integral(2*kKaonMass,100.);
    auto fPartial   =   fBreitWigner->Integral(0.996,1.065);
    
    auto fTotal2D   =   fBreitWigner2D->Integral(2*kKaonMass,100.,2*kKaonMass,100.);
    auto fPartial2D =   fBreitWigner2D->Integral(0.996,1.065,0.996,1.065);
    
    auto fStripe    = fLevyTsallis2D->Integral  (0.,    0.4,    0.4,    100.) + fLevyTsallis2D->Integral  (0.4,    10.,    10.,    100.);
    auto fBlock     = fLevyTsallis2D->Integral  (0.,    0.4,    0.,     0.4) + fLevyTsallis2D->Integral  (10.,    100.,    10.,     100.);
    auto f1_  =fBlock+2.*fStripe;
    auto f2_ = hRt->Integral(-1,1000,-1,1000,"width");
    cout << "Out of Measure:    " << f1_ << endl;
    cout << "In Measure:        " << f2_ << endl;
    cout << "FRAC:              " << fPartial/fTotal << endl;
    cout << "FRAC:              " << fPartial2D/fTotal2D << endl;
    
    fout -> Close();

    /*
    auto    hMeasured   =   (TH1F*)(insFile_DT_Yield->Get("hRES_1D_Stat"));
    auto    hMeasure2   =   (TH1F*)(insFile_DT_Yield->Get("hRES_1D_Syst"));
    auto    hREC        =   (TH1F*)(insFile_EF_Yield->Get("hREC_Rw_1D"));
    auto    hGEN        =   (TH1F*)(insFile_EF_Yield->Get("hGEN_Rw_1D"));
    auto    hTotal      =   fSumErrors(hMeasured,hMeasure2);
    
    
    TFile*  fout            =   new TFile   (Form(kAnalysis_SigExtr_Dir,"Yield")+TString("/dick.root"),"recreate");
    
    fSetAllFunctions();
    fSetFunction();
    
    TCanvas* c2 = new TCanvas();
    gPad->SetLogx();
    gStyle->SetOptStat(0);
    
    hTotal->Fit(fLevyTsallis,"IM");
    
    hTotal->Draw();
    fLevyTsallis->Draw("same");
    c2->SaveAs(Form(kAnalysis_SigExtr_Dir,"Yield")+TString("/dick2.pdf"));
    
    ReweightEfficiency(hTotal,fLevyTsallis,hGEN,hREC,fout,0,2,1.e-6);
    
    fout->Close();
    
    
    TFile*  fin            =   new TFile   (Form(kAnalysis_SigExtr_Dir,"Yield")+TString("/dick.root"));
    
    auto fREC0 = (TH1F*)(fin->Get("hREC_Rw_1D_i0"));
    auto fGEN0 = (TH1F*)(fin->Get("hGEN_Rw_1D_i0"));
    auto fRECL = (TH1F*)(fin->Get("hREC_Rw_1D_i2"));
    auto fGENL = (TH1F*)(fin->Get("hGEN_Rw_1D_i2"));
    
    fREC0->Divide(fGEN0);
    fRECL->Divide(fGENL);
    
    TCanvas* c1 = new TCanvas();
    c1->Divide(2,1);
    c1->cd(1);
    gPad->SetLogx();
    gStyle->SetOptStat(0);
    fREC0->DrawCopy();
    fRECL->DrawCopy();
    c1->cd(2);
    gPad->SetLogx();
    
    fREC0->Divide(fRECL);
    fREC0->DrawCopy("HIST");
    
    c1->SaveAs(Form(kAnalysis_SigExtr_Dir,"Yield")+TString("/dick.pdf"));
    
    fin->Close();
    
    */
    /*
    auto fTotal     =   0.;
    auto fPartial   =   0.;
    auto fTota2     =   0.;
    auto fPartia2   =   0.;
    auto step       =   1.e-4;
    auto start      =   2*.493677;
    auto stop       =   100;
    
    for ( int i = start*(1./step); i < stop*(1./step); i++) {
        auto fValue =   step*TMath::BreitWigner(i*step,1.019455,0.004249);
        fTotal  +=  fValue;
        if ( i*step <= 1.065 && i*step >= 0.996 ) fPartial  +=  fValue;
        for ( int j = start*(1./step); j < stop*(1./step); j++) {
            auto fValu2 =   fValue*step*TMath::BreitWigner(j*step,1.019455,0.004249);
            fTota2  +=  fValu2;
            if ( i*step <= 1.065 && i*step >= 0.996 && j*step <= 1.065 && j*step >= 0.996 ) fPartia2  +=  fValu2;
        }
    }
    
    cout << fTotal << " - " << fPartial << " - " << fPartial/fTotal << endl;
    cout << fTota2 << " - " << fPartia2 << " - " << fPartia2/fTota2 << endl;
    */


    
    /*
    TFile  *fAsym = new TFile("/Users/nikolajal/alidock/AliAnalysisPhiCount/result/Yield/SignalExtraction/FitResults_Asymm.root");
    TFile  *fSymm = new TFile("/Users/nikolajal/alidock/AliAnalysisPhiCount/result/Yield/SignalExtraction/FitResults_Symm.root");
    
    TH2F* hAsym = (TH2F*)(fAsym->Get("hRAW_2D"));
    TH2F* hSymm = (TH2F*)(fSymm->Get("hRAW_2D"));
    TH2F* hAsymErr = (TH2F*)(fAsym->Get("hRAW_2D")->Clone());
    TH2F* hSymmErr = (TH2F*)(fSymm->Get("hRAW_2D")->Clone());
    
    fSetAllBins ();
    
    TH1F* hAsym1D = new TH1F("Asymm","Asymm",nBinPT2D,fArrPT2D);
    TH1F* hSymm1D = new TH1F("Symm","Symm",nBinPT2D,fArrPT2D);
    
    for ( Int_t iBin = 0; iBin < hSymm->GetNbinsX(); iBin++ ) {
        hAsym1D->SetBinContent  ( iBin+1, hAsym->GetBinError  (iBin+1,iBin+1) );
        hAsym1D->SetBinError    ( iBin+1, 0/*hAsym->GetBinError    (iBin+1,iBin+1)*/// );
        //hSymm1D->SetBinContent  ( iBin+1, hSymm->GetBinError  (iBin+1,iBin+1) );
        //hSymm1D->SetBinError    ( iBin+1, 0/*hSymm->GetBinError    (iBin+1,iBin+1)*/ );
    /*
    hAsymErr->SetBinContent  ( iBin+1,iBin+1,hAsym->GetBinError  (iBin+1,iBin+1) );
        hSymmErr->SetBinContent  ( iBin+1,iBin+1,hSymm->GetBinError  (iBin+1,iBin+1) );
    }
    auto hRatio = (TH1F*)hAsym1D->Clone();
    hRatio->Divide(hSymm1D,hAsym1D);
    hRatio->SetMarkerColor(kGreen-1);
    hRatio->SetLineColor(kGreen-1);
    hRatio->SetMarkerStyle(26);
    hRatio->SetTitle("Ratio Symm/Asymm");
    
    TCanvas    *c1  =   new TCanvas("","",1400,700);
    c1->Divide(2,1);
    c1->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    hAsym1D->SetMarkerColor(kRed);
    hAsym1D->SetLineColor(kRed);
    hAsym1D->SetMarkerStyle(23);
    hAsym1D->Draw();
    hSymm1D->SetMarkerColor(kBlue);
    hSymm1D->SetLineColor(kBlue);
    hSymm1D->SetMarkerStyle(25);
    hSymm1D->SetFillColorAlpha(kBlue,0.15);
    hSymm1D->Draw("same E2");
    gPad->BuildLegend();
    c1->cd(2);
    gStyle->SetOptStat(0);
    hRatio->Draw();
    delete c1;
    
    TCanvas    *c2  =   new TCanvas("","",1400,1400);
    auto hRatio2D = (TH1F*)hAsymErr->Clone();
    hRatio2D->Divide(hSymmErr,hAsymErr);
    hRatio2D->SetTitle("Errors Symm/Asymm");
    gPad->SetLogx();
    gPad->SetLogy();
    hRatio2D->Draw("colz text");
    delete c2;
    
    fSetLevyTsallis2D();
    fSetFunction(fLevyTsallis2D);
    hAsym->Fit(fLevyTsallis2D,"");
}
    /*
     
    Float_t *fArrPT2D_Tst = new Float_t[12];
    
    fArrPT2D_Tst[0]     =   0.00; //0.4
    fArrPT2D_Tst[1]     =   0.40; //0.3
    fArrPT2D_Tst[2]     =   0.70; //0.2
    fArrPT2D_Tst[3]     =   0.90; //0.1
    fArrPT2D_Tst[4]     =   1.00; //0.2
    fArrPT2D_Tst[5]     =   1.20; //0.2
    fArrPT2D_Tst[6]     =   1.40; //0.2
    fArrPT2D_Tst[7]     =   1.60; //0.4
    fArrPT2D_Tst[8]     =   2.00; //0.8
    fArrPT2D_Tst[9]     =   2.80; //1.2
    fArrPT2D_Tst[10]    =   4.00; //6.0
    fArrPT2D_Tst[11]    =   10.0;
    
    /*
    //Retrieving Event data
    TFile *insFileMC        =   new TFile   ("/Volumes/[HD][Nikolajal]_Toshiba/Dataset/2010/Sim/2021_05_12/LHC14j4_STD.root");
    
    //Retrieving Event data TTree
    TTree   *TPhiEfficiency =   (TTree*)insFileMC->Get(fPhiCandidateEff_Tree);
    TTree   *TKaonCandidate =   nullptr;//(TTree*)insFileMC->Get(fKaonCandidateEff_Tree);
    TTree   *TPhiCandidate  =   (TTree*)insFileMC->Get(fPhiCandidate_Tree);
    TTree   *TKaonEfficiency=   nullptr;//(TTree*)insFileMC->Get(fKaonCandidate_Tree);
    
    // Retrieving Event Count Histogram
    TList  *fQCOutputList   =   (TList*)insFileMC       ->Get("fQCOutputList");
    TH1D   *fHEventCount    =   (TH1D*) fQCOutputList   ->FindObject("fQC_Event_Enumerate");
    
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
    fSetAllBins ();
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
    TH2F       *hTST_2D;
    //
    //  Defining Efficiency and check utilities
     //
     hName       =   Form("hTST_2D");
     hTitle      =   Form("hTST_2D");
     hTST_2D     =   new TH2F (hName,hTitle,11,fArrPT2D_Tst,11,fArrPT2D_Tst);
     SetAxis(hTST_2D,"PT 2D");
     //
     hName       =   Form("hTS2_2D");
     hTitle      =   Form("hTS2_2D");
     hTS2_2D     =   new TH2F (hName,hTitle,1000,0.,10.,1000,0.,10.);
     SetAxis(hTS2_2D,"PT 2D");
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    fStartTimer("Efficiency Utility Histograms Production");
    
    auto nEventsCut = -1.;
    
    // Evaluating entries
    Int_t nEvents = (!TPhiEfficiency) ? 0 : ( nEventsCut == -1.? TPhiEfficiency->GetEntries() : nEventsCut);
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEvents; iEvent++ )    {
        // Recovering events
        TPhiEfficiency->GetEntry(iEvent);
        
        fPrintLoopTimer("Efficiency Utility Histograms Production",iEvent,nEvents,kPrintIntervalPP);

        // Utilities
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        for ( Int_t iPhi = 0; iPhi < evPhiEfficiency.nPhi; iPhi++ ) {
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[iPhi],evPhiEfficiency.Py[iPhi],evPhiEfficiency.Pz[iPhi],1.019455);
            if ( !fAcceptCandidate(1.019455,LPhi_candidate1.Pt()) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )    {
            // Must have at least 1 candidate
            if ( U_nAccept < 1 ) break;

            // Building First Candidate
            LPhi_candidate1.SetXYZM(evPhiEfficiency.Px[U_AccCand[iPhi]],evPhiEfficiency.Py[U_AccCand[iPhi]],evPhiEfficiency.Pz[U_AccCand[iPhi]],1.019455);

            // >> 1-Dimensional Analysis Fill
            //
            // >>-->> Utilities
            //
            // >>-->>-->> Event
            //
            Bool_t  fNoVtxRec           =   (fCheckMask(evPhiEfficiency.TrueEventMask,0) || fCheckMask(evPhiEfficiency.TrueEventMask,1));
            Bool_t  fNoVtxCut           =   (fCheckMask(evPhiEfficiency.TrueEventMask,2) );
            Bool_t  fIsMBevent          =   (evPhiEfficiency.TrueEventMask == 0);
            //
            // >>-->>-->> Multiplicity
            //
            Int_t   iMult               =   fGetBinMult(evPhiEfficiency.Multiplicity);
            Bool_t  fHasMultiplicity    =   iMult != -1;
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
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )    {
                // Must have at least 2 candidates
                if ( U_nAccept < 2 ) break;
                
                // Protection against auto-correlation
                if ( iPhi == jPhi ) continue;

                // Building Second Candidate
                LPhi_candidate2.SetXYZM(evPhiEfficiency.Px[U_AccCand[jPhi]],evPhiEfficiency.Py[U_AccCand[jPhi]],evPhiEfficiency.Pz[U_AccCand[jPhi]],1.019455);

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
                if ( fHasRapidity )   {
                //
                // >>-->> Yield
                //
                hTST_2D                                     ->  Fill(iTransMom,jTransMom,0.5);
                hTS2_2D                                     ->  Fill(iTransMom,jTransMom,0.5);
                //
                }
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
    auto fNormEvent = fHEventCount->GetBinContent(1);
     hTST_2D->Scale(1.,"width");
     hTS2_2D->Scale(1.,"width");
     hTST_2D->Scale(1./fNormEvent);
     hTS2_2D->Scale(1./fNormEvent);
    //
    TFile *outFil0  =   new TFile   ("./result/TEST/test.root","recreate");
    //
     hTST_2D->Write();
     hTS2_2D->Write();
    //
    outFil0->Close();
    
    
    TFile  *inFile  =   new TFile("./result/TEST/test.root");
    
    TH2F * hTST_2D = (TH2F*)(inFile->Get("hTST_2D"));
    TH1F * hDiag    =   new TH1F("hDiag",   "hDiag",    11, fArrPT2D_Tst);
    TH1F * hSlice   =   new TH1F("hSlice",  "hSlice",   11, fArrPT2D_Tst);
    
    
    fSetAllFunctions();
    fSetFunction();
    
    TH1F * hCheck = new TH1F("hCheck","hCheck", 11, fArrPT2D_Tst);
    for ( Int_t iTer = 1; iTer <= 14; iTer++ ) {
        hDiag->SetBinContent( iTer, hTST_2D->GetBinContent( iTer, iTer ) );
        hDiag->SetBinError  ( iTer, hTST_2D->GetBinError  ( iTer, iTer ) );
        for ( Int_t jTer = 1; jTer <= 14; jTer++ ) {
            hSlice->SetBinContent( jTer, hTST_2D->GetBinContent( 1, jTer ) );
            hSlice->SetBinError  ( jTer, hTST_2D->GetBinError  ( 1, jTer ) );
        }
        hTST_2D->ProjectionY(Form("_%i",iTer),iTer,iTer)->Fit(fLevyTsallis,"","IMREQ0S",0.4,2.0);
        hCheck->SetBinContent   ( iTer, fLevyTsallis->Integral(0.,0.4) );
        hCheck->SetBinError     ( iTer, fLevyTsallis->IntegralError(0.,0.4) );
    }
    hCheck->Scale(1./(0.4),"");
    
    TLatex *Latx = new TLatex();
    
    
    TCanvas *c3 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    
    hCheck->SetLineColor(kRed);
    hCheck->SetMarkerColor(kRed);
    hCheck->SetMarkerStyle(24);
    
    hDiag->SetLineColor(kBlue);
    hDiag->SetMarkerColor(kBlue);
    hDiag->SetMarkerStyle(22);
    
    fLevyTsallis2->FixParameter(0,1.019455);
    fLevyTsallis2->SetParLimits(1,2.01,10);
    fLevyTsallis2->SetParLimits(2,0.,10);
    fLevyTsallis2->SetParLimits(3,0.,1.);
    
    fSetAllFunctions();
    fSetFunction(fBoseEinstein);
    hDiag->Fit(fLevyTsallis2,"IMRE0S","",0.4,12.);
    auto fFucnt = hDiag->GetFunction("LevyTsallis2");
    fFucnt->SetRange(0.,12.);
    
    TCanvas *c1 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hDiag->Draw();
    fFucnt->DrawCopy("same");
    Latx->DrawLatexNDC(0.6,0.85,Form("TRGe6: %.3f",1.e6*hDiag->GetBinContent( 1. )));
    Latx->DrawLatexNDC(0.6,0.80,Form("MSRe6: %.3f",1.e6*(1./0.4)*fFucnt->Integral( 0., .4 )));
    Latx->DrawLatexNDC(0.6,0.75,Form("MSRe6: %.3f",(1.e6*hDiag->GetBinContent( 1. ))/(1.e6*(1./0.4)*fFucnt->Integral( 0., .4 ))));
    c1->SaveAs("./result/TEST/c1.pdf");
    delete c1;
    
    c3->cd();
    hDiag->DrawCopy();
    fFucnt->DrawCopy("same");
    
    hSlice->SetLineColor(kGreen-1);
    hSlice->SetMarkerColor(kGreen-1);
    hSlice->SetMarkerStyle(22);
    
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    hSlice->Fit(fLevyTsallis,"IMRE0S","",0.4,12.);
    
    TCanvas *c2 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hSlice->Draw();
    fLevyTsallis->DrawCopy("same");
    hCheck->Draw("same");
    Latx->DrawLatexNDC(0.6,0.85,Form("TRGe6: %.3f",1.e6*hSlice->GetBinContent( 1. )));
    Latx->DrawLatexNDC(0.6,0.80,Form("MSRe6: %.3f",1.e6*(1./0.4)*fLevyTsallis->Integral( 0., .4 )));
    Latx->DrawLatexNDC(0.6,0.75,Form("TRGe6: %.3f",1.e6*hSlice->GetBinContent( 2. )));
    Latx->DrawLatexNDC(0.6,0.70,Form("MSRe6: %.3f",1.e6*(1./0.3)*fLevyTsallis->Integral( 0.4, .7 )));
    c2->SaveAs("./result/TEST/c2.pdf");
    delete c2;
    
    c3->cd();
    hSlice->DrawCopy("same");
    fLevyTsallis->DrawCopy("same");
    (gPad->BuildLegend(0.1,0.1,0.4,0.3,""))->SetName("gggg");
    
    hDiag->Fit(fLevyTsallis2,"IMRE0QS","",0.4,4.);
    auto fValue_Diag = fLevyTsallis2 ->Integral(0.,0.4);
    auto fError_Diag = fLevyTsallis2 ->IntegralError(0.,0.4);
    Latx->DrawLatexNDC(0.4,0.85,Form("INT Diag (x10^{6}): %.3f #pm %.3f",1.e6*(1./0.4)*fValue_Diag,1.e6*(1./0.4)*fError_Diag));
    hSlice->Fit(fLevyTsallis,"IMRE0QS","",0.4,1.5);
    Latx->DrawLatexNDC(0.4,0.80,Form("INT Slic (x10^{6}): %.3f #pm %.3f",1.e6*(1./0.4)*fLevyTsallis  ->Integral(0.,0.4),1.e6*(1./0.4)*fLevyTsallis  ->IntegralError(0.,0.4)));
    Latx->DrawLatexNDC(0.4,0.75,Form("INT Ratio Diag/Slic   : %.3f",fValue_Diag/(fLevyTsallis  ->Integral(0.,0.4))));
    Latx->DrawLatexNDC(0.4,0.70,Form("INT Ratio Diag/Truth   : %.3f",(fValue_Diag)/(0.4*hTST_2D->GetBinContent(1,1))));
    Latx->DrawLatexNDC(0.4,0.65,Form("INT Ratio Slic/Truth   : %.3f",(fLevyTsallis  ->Integral(0.,0.4))/(0.4*hTST_2D->GetBinContent(1,1))));
    
    c3->SaveAs("./result/TEST/c3.pdf");
    
    ((TH1F*)c3->GetPrimitive("hDiag_copy"))->GetXaxis()->SetRangeUser(0.,1.5);
    auto tl = ((TLegend*)c3->GetPrimitive("gggg"));
    tl->SetX1NDC(0.3);
    tl->SetX2NDC(0.7);
    tl->SetY1NDC(0.1);
    tl->SetY2NDC(0.3);
    
    c3->SaveAs("./result/TEST/c3_CU.pdf");
    delete c3;
    
    TFile  *inFil2  =   new TFile("./result/TEST/Yield/SignalExtraction/FitResults.root");
    
    TH1F * hSliceTrue = (TH1F*)(inFil2->Get("hConditionalYieldEXT_"));
    hSliceTrue->Scale(1./(0.4),"");
    hSliceTrue->SetTitle("hSlic");
    TH1F * hDiagnTrue = (TH1F*)(inFil2->Get("hConditionalYieldTest"));
    hDiagnTrue->SetTitle("hDiag");
    
    c3 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    
    hDiagnTrue->SetLineColor(kBlue);
    hDiagnTrue->SetMarkerColor(kBlue);
    hDiagnTrue->SetMarkerStyle(22);
    
    fLevyTsallis2->FixParameter(0,1.019455);
    fLevyTsallis2->SetParLimits(1,2.01,10);
    fLevyTsallis2->SetParLimits(2,0.,10);
    
    fSetAllFunctions();
    fSetFunction(fBoseEinstein);
    hDiagnTrue->Fit(fLevyTsallis2,"IMRE0S","",0.4,4.);
    fFucnt = hDiagnTrue->GetFunction("LevyTsallis2");
    fFucnt->SetRange(0.,12.);
    
    c1 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hDiagnTrue->Draw();
    fFucnt->DrawCopy("same");
    c1->SaveAs("./result/TEST/c1_2.pdf");
    delete c1;
    
    c3->cd();
    fFucnt->SetRange(0.,10.);
    fFucnt->DrawCopy();
    hDiagnTrue->DrawCopy("same");
    
    hSliceTrue->SetLineColor(kGreen-1);
    hSliceTrue->SetMarkerColor(kGreen-1);
    hSliceTrue->SetMarkerStyle(22);
    
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    hSliceTrue->Fit(fLevyTsallis,"IMRE0S","",0.4,4.);
    
    c2 = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hSliceTrue->Draw();
    fLevyTsallis->DrawCopy("same");
    c2->SaveAs("./result/TEST/c2_2.pdf");
    delete c2;
    
    c3->cd();
    fLevyTsallis->DrawCopy("same");
    hSliceTrue->DrawCopy("same");
    (gPad->BuildLegend())->SetName("gggg");
    
    hDiagnTrue->Fit(fLevyTsallis2,"IMRE0QS","",0.4,4.0);
    fValue_Diag = fLevyTsallis2 ->Integral(0.,0.4);
    fError_Diag = fLevyTsallis2 ->IntegralError(0.,0.4);
    Latx->DrawLatexNDC(0.4,0.85,Form("INT Diag (x10^{6}): %.3f #pm %.3f",1.e6*(1./0.4)*fValue_Diag,1.e6*(1./0.4)*fError_Diag));
    hSliceTrue->Fit(fLevyTsallis,"IMRE0QS","",0.4,1.5);
    Latx->DrawLatexNDC(0.4,0.80,Form("INT Slic (x10^{6}): %.3f #pm %.3f",1.e6*(1./0.4)*fLevyTsallis  ->Integral(0.,0.4),1.e6*(1./0.4)*fLevyTsallis  ->IntegralError(0.,0.4)));
    Latx->DrawLatexNDC(0.4,0.75,Form("INT Ratio Diag/Slic   : %.3f",fValue_Diag/(fLevyTsallis  ->Integral(0.,0.4))));
    
    c3->SaveAs("./result/TEST/c3_2.pdf");
    
    ((TF1*)c3->GetPrimitive("LevyTsallis2"))->SetRange(0.,2.0);
    ((TF1*)c3->GetPrimitive("LevyTsallis2"))->SetMaximum(5.e-3);
    tl = ((TLegend*)c3->GetPrimitive("gggg"));
    tl->SetX1NDC(0.3);
    tl->SetX2NDC(0.7);
    tl->SetY1NDC(0.1);
    tl->SetY2NDC(0.3);
    
    c3->SaveAs("./result/TEST/c3_2_CU.pdf");
    delete c3;
    
    auto kNBin  = 2;
    auto kMin   = 0.;
    auto kMax   = 2.;
    
    TH2D* Check     =   new TH2D("Check","Check",kNBin, kMin, kMax,kNBin, kMin, kMax);
    TH2D* Chec2     =   new TH2D("Chec2","Chec2",kNBin, kMin, kMax,kNBin, kMin, kMax);
    TH2D* Chec3     =   new TH2D("Chec3","Chec3",kNBin, kMin, kMax,kNBin, kMin, kMax);
    TH2D* Chec4     =   new TH2D("Chec4","Chec4",kNBin, kMin, kMax,kNBin, kMin, kMax);
    TH2D* Chec5     =   new TH2D("Chec5","Chec5",kNBin, kMin, kMax,kNBin, kMin, kMax);
    
    TH1D* Slic1     =   new TH1D("Slic1","Slic1",kNBin, kMin, kMax);
    TH1D* Slic2     =   new TH1D("Slic2","Slic2",kNBin, kMin, kMax);
    float * fPTArray = new Float_t[100];
    
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    
    for ( Int_t iEvt = 0; iEvt < 1.e6; iEvt++ ) {
        auto fNPhiMesons = 2;
        auto hchoose = fRandomGen->Uniform(0,1);
        for ( Int_t iPhi = 0; iPhi < fNPhiMesons; iPhi++ ) {
            
            fPTArray[0] = fabs(fRandomGen->Uniform(0.,2.));
            fPTArray[1] = fabs(fRandomGen->Uniform(0.,2.));
            
            /*
            auto fPT1 = fabs(fRandomGen->Uniform(0.,10.));
            auto fPT2 = fabs(fRandomGen->Uniform(0.,12.));
            
            auto minmaxchoose = fRandomGen->Uniform(0,.1);
            
            if ( minmaxchoose < 0.5 )   {
                fPTArray[0] = min(fPT1,fPT2);
                fPTArray[1] = max(fPT1,fPT2);
            }   else    {
                fPTArray[1] = min(fPT1,fPT2);
                fPTArray[0] = max(fPT1,fPT2);
            }
        }
        for ( Int_t iPhi = 0; iPhi < fNPhiMesons; iPhi++ ) {
            for ( Int_t jPhi = 0; jPhi < fNPhiMesons; jPhi++ ) {
                if ( iPhi == jPhi ) continue;
                Check->Fill( fPTArray[iPhi], fPTArray[jPhi],0.5);
                if ( fPTArray[iPhi] < 1. && fPTArray[iPhi] >= 0. ) Slice->Fill( fPTArray[jPhi] );
                if ( fPTArray[iPhi] < 2. && fPTArray[iPhi] >= 1. ) Slice->Fill( fPTArray[jPhi] );
                if ( iPhi >= jPhi ) continue;
                Chec2->Fill( fPTArray[iPhi], fPTArray[jPhi]);
                Chec3->Fill( min(fPTArray[iPhi],fPTArray[jPhi]), max(fPTArray[iPhi],fPTArray[jPhi]));
            }
        }
    }
    for ( Int_t iPhi = 1; iPhi <= 12; iPhi++ ) {
        for ( Int_t jPhi = 1; jPhi <= 12; jPhi++ ) {
            Chec4->SetBinContent( iPhi, jPhi, 0.5*Chec3->GetBinContent(min(iPhi,jPhi),max(iPhi,jPhi)) );
            if ( iPhi == jPhi ) Chec4->SetBinContent( iPhi, jPhi, Chec3->GetBinContent(min(iPhi,jPhi),max(iPhi,jPhi)) );
        }
    }
    
    TH1D* TestDiag  =   new TH1D("TestDiag","TestDiag",kNBin, kMin, kMax);
    TH1D* TestDia2  =   new TH1D("TestDia2","TestDia2",kNBin, kMin, kMax);
    TH1D* TestDia3  =   new TH1D("TestDia3","TestDia3",kNBin, kMin, kMax);
    TH1D* TestDia4  =   new TH1D("TestDia4","TestDia4",kNBin, kMin, kMax);
    
    for ( Int_t iTer = 0; iTer < kNBin; iTer++ )   {
        auto fCurrentBin = Check->GetBinContent(iTer+1,iTer+1);
        auto fCurrentBiE = Check->GetBinError(iTer+1,iTer+1);
        TestDiag->SetBinContent (iTer+1,fCurrentBin);
        TestDiag->SetBinError   (iTer+1,fCurrentBiE);
        
        fCurrentBin = Chec2->GetBinContent(iTer+1,iTer+1);
        fCurrentBiE = Chec2->GetBinError(iTer+1,iTer+1);
        TestDia2->SetBinContent (iTer+1,fCurrentBin);
        TestDia2->SetBinError   (iTer+1,fCurrentBiE);
        
        fCurrentBin = Chec3->GetBinContent(iTer+1,iTer+1);
        fCurrentBiE = Chec3->GetBinError(iTer+1,iTer+1);
        TestDia3->SetBinContent (iTer+1,fCurrentBin);
        TestDia3->SetBinError   (iTer+1,fCurrentBiE);
        
        fCurrentBin = Chec4->GetBinContent(iTer+1,iTer+1);
        fCurrentBiE = Chec4->GetBinError(iTer+1,iTer+1);
        TestDia4->SetBinContent (iTer+1,fCurrentBin);
        TestDia4->SetBinError   (iTer+1,fCurrentBiE);
    }
    
    
    cout << " Target: " << Check->GetBinContent(1,1) << endl;
    
    
    TFile * f1 = new TFile("test.root","recreate");
    Slice->Write();
    Check->Write();
    Chec2->Write();
    Chec3->Write();
    Chec4->Write();
    TestDiag->Write();
    TestDia2->Write();
    TestDia3->Write();
    TestDia4->Write();
    
    return;
    /*
    c1 = new TCanvas();
    gPad->SetLogz();
    Check->Draw("colz");
    c1->SaveAs("./result/test2.pdf");
    gPad->SetLogy();
    TestSide->SetMinimum(1.e-8);
    TestSide->Draw();
    //TestSide->Fit("landau","IMRE0S","",0.4,3.);
    auto ffunction = TestSide->GetFunction("landau");
    ffunction->SetRange(0.,12.);
    ffunction->DrawCopy("SAME");

    //cout << " Side: " << ffunction->Integral(0.,0.4) << " .-. " << TestSide->Integral(1,1,"width") << endl;
    
    TestDiag->Draw("SAME");
    TestDiag->Fit("landau","IMRE0S","",0.4,3.);
    ffunction = TestDiag->GetFunction("landau");
    ffunction->SetRange(0.,12.);
    ffunction->DrawCopy("SAME");
    
    cout << " Diag: " << ffunction->Integral(0.,0.4) << " .-. " << TestDiag->Integral(1,1,"width") << endl;
    
    c1->SaveAs("./result/test.pdf");
    delete c1;
    
    TH1D * hFirstConditional    =   new TH1D    ("hFirstConditional",   "hFirstConditional",    10,0,10);
    TH1D * hSecondConditional   =   new TH1D    ("hSecondConditional",  "hSecondConditional",   10,0,10);
    
    TH2D * hMeasureHistogram    =   new TH2D    ("hMeasureHistogram",   "hMeasureHistogram",    10,0,10, 10,0,10);
    TH2D * hLegend              =   new TH2D    ("hLegend",             "hLegend",              10,0,10, 10,0,10);
    TH1D * hSecondMeasured      =   new TH1D    ("hSecondMeasured",     "hSecondMeasured",      10,0,10);
    
    int kNEvents = int(1.e7);
    Float_t*fMemoryMomentum = new Float_t[10];
    Double_t fPairCount = 0.;
    Double_t fPairCoun2 = 0.;
    for ( Int_t iTer = 0; iTer < kNEvents; iTer++ )   {
        auto    fPhiProduced    =   (int)(min(fabs(fRandomGen->Gaus(1,2)),10.));
        
        if ( fPhiProduced < 2 ) continue;
        
        for ( Int_t iPhi = 0; iPhi < fPhiProduced; iPhi++ )   {
            fMemoryMomentum[iPhi]   =   min(fabs(fRandomGen->Gaus(1,4)),9.99999);
        }
        
        Bool_t  fFill   =   false;
        for ( Int_t iPhi = 0; iPhi < fPhiProduced; iPhi++ )   {
            if ( fMemoryMomentum[iPhi] < 1 ) fFill   =   true;
            for ( Int_t jPhi = 0; jPhi < fPhiProduced; jPhi++ )   {
                if ( iPhi == jPhi ) continue;
                hSecondConditional  ->  Fill( fMemoryMomentum[iPhi] );
                hMeasureHistogram   ->  Fill( fMemoryMomentum[iPhi], fMemoryMomentum[jPhi],   0.5 );
                fPairCoun2 += 0.5;
                if ( fFill ) hFirstConditional->Fill( fMemoryMomentum[jPhi] );
                if ( iPhi < jPhi ) continue;
                fPairCount++;
            }
            fFill   =   false;
        }
    }
    fPairCount /= kNEvents;
    fPairCoun2 /= kNEvents;
    hFirstConditional   ->Scale(1./kNEvents);
    hSecondConditional  ->Scale(1./kNEvents);
    hMeasureHistogram   ->Scale(1./kNEvents);
    
    for ( int iBin = 0; iBin < 10; iBin++ ) {
        auto fTarget = hMeasureHistogram->ProjectionX(Form("_%i",iBin),iBin+1,iBin+1);
        hSecondMeasured->SetBinContent( iBin+1, fTarget->Integral() );
        for ( int jBin = 10; jBin > iBin; jBin-- ) {
            hLegend->SetBinContent( iBin, jBin, iBin );
        }
    }
    
    TCanvas*cTest_ = new TCanvas("","",1400,800);
    gStyle->SetOptStat(0);
    hFirstConditional->SetMarkerStyle(25);
    hFirstConditional->SetMarkerColor(kRed);
    hFirstConditional->SetLineColor(kRed);
    
    hFirstConditional->SetMinimum(0);
    hFirstConditional->SetMaximum(0.1);
    hFirstConditional->Draw("HIST");
    auto hMeasuredProjection = hMeasureHistogram->ProjectionX("hMeasuredProjection",1,1);
    hMeasuredProjection->SetMarkerStyle(22);
    hMeasuredProjection->SetMarkerColor(kBlue);
    hMeasuredProjection->SetLineColor(kBlue);
    hMeasuredProjection->Draw("HIST SAME");
    gPad->BuildLegend(0.6,.75,.9,.9);
    cTest_->SaveAs("./result/TEST/Check1.pdf");
    hSecondConditional->SetMarkerStyle(25);
    hSecondConditional->SetMarkerColor(kRed);
    hSecondConditional->SetLineColor(kRed);
    
    hSecondConditional->SetMinimum(0);
    hSecondConditional->SetMaximum(0.5);
    hSecondConditional->Draw("HIST");
    hSecondMeasured->SetMarkerStyle(22);
    hSecondMeasured->SetMarkerColor(kBlue);
    hSecondMeasured->SetLineColor(kBlue);
    hSecondMeasured->Draw("SAME HIST");
    gPad->BuildLegend(0.6,.75,.9,.9);
    cTest_->SaveAs("./result/TEST/Check2.pdf");
    cout << "TEST INT: " << hMeasureHistogram->Integral()/2 << " TRUTH: " << fPairCount << endl;
    delete cTest_;
    
    cTest_ = new TCanvas();
    
    
    return;
    /*
    //---------------------//
    //  Setting up input   //
    //---------------------//
    
    //-// OPTIONS
    
    // Silencing warnings for smoother
    if ( true )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(true);
    }
    fChooseOption("");
    gROOT                               ->  ProcessLine(Form(".! mkdir -p ./result/yield/ExtrapolateCheck/1D/"));
    gROOT                               ->  ProcessLine(Form(".! mkdir -p ./result/yield/ExtrapolateCheck/2D/"));
    
    // Retrieving PreProcessed data histograms
    TFile*  insFile_DT_Yield            =   new TFile   (Form(kASigExtr_FitCheckRst,"Yield"));
    TFile*  insFile_EF_Yield            =   new TFile   (Form(kAnalysis_MCTruthHist,"Yield"));
    
    
    // Recovering the histograms-------------------------------------------------------------------------------

    // >-> GENERAL ANALYSIS //
    //
    TH1D       *hEvntEff;
    TH1D       *hEvntMlt;
    //
    hName       =   "fQC_Event_Enumerate";
    hEvntEff    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "fQC_Event_Enum_Mult";
    hEvntMlt    =   (TH1D*)(insFile_DT_Yield->Get(hName));
    //
    //
    // >-> YIELD ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F       *hRAW_1D;
    TH1F       *hREC_1D;
    TH1F       *hGEN_1D;
    TH1F       *hTRU_RECVTX_1D;
    TH1F       *hTRU_ALLVTX_1D;
    TH1F       *hREC_1D_in_2D_bin;
    TH1F       *hGEN_1D_in_2D_bin;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hRAW_1D";
    hRAW_1D     =   (TH1F*)(insFile_DT_Yield->Get(hName));
    //
    hName       =   "hREC_1D";
    hREC_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_1D";
    hGEN_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hTRU_RECVTX_1D";
    hTRU_RECVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hTRU_ALLVTX_1D";
    hTRU_ALLVTX_1D     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hREC_1D_in_2D_bin";
    hREC_1D_in_2D_bin     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    hName       =   "hGEN_1D_in_2D_bin";
    hGEN_1D_in_2D_bin     =   (TH1F*)(insFile_EF_Yield->Get(hName));
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH2F       *hRAW_2D;
    //
    //  Defining cumulative histogram over measurable pT
    //
    hName       =   "hRAW_2D";
    hRAW_2D     =   (TH2F*)(insFile_DT_Yield->Get(hName));
    //
    //---------------------//
    //  Setting up output  //
    //---------------------//
    //
    // Generating the binning array--------------------------------------------------------------------------
    //
    fSetAllBins();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept;
    //
    // Creating the histograms-------------------------------------------------------------------------------
    //
    //---------------------//
    // Preprocessing input //
    //---------------------//
    //
    gROOT                       ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_SigExtr_Dir,(TString("TEST/")+TString("Yield")).Data())));
    TFile*  outFile_DT_Yield            =   new TFile   (Form(kASigExtr_FitCheckRst,(TString("TEST/")+TString("Yield")).Data()),"recreate");
    //
    //                 N_raw            f_norm X f_vtx X f_SL
    // N_res = --------------------- X -----------------------
    //          EXA X DpT X Dy X BR             N_MB
    //
    // Scaling in pT [Done in PreProcessing]
    //
    // Scaling for efficiencies
    //
    auto        kN_Trg          =   (hEvntEff->GetBinContent(kEventCount::kTrigger));
    auto        kN_Vtx          =   (hEvntEff->GetBinContent(kEventCount::kVertex));
    auto        kN_MB           =   (hEvntEff->GetBinContent(kEventCount::kVertex10));
    Double_t    f1DCorrection   =   (1./kBR)        *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(kTriggerEff/kN_MB)*(kN_Vtx/kN_Trg);
    //
    // >-> YIELD ANALYSIS //
    //
    // >->-->-> 1-Dimension analysis //
    //
    //  Declaring all histograms
    //
    TH1F   *hRES_1D_Stat    =   fEfficiencycorrection   ( fEfficiencycorrection(hRAW_1D,hREC_1D,hGEN_1D,f1DCorrection),hTRU_RECVTX_1D,hTRU_ALLVTX_1D );
    TH1F   *hRES_1D_Syst    =   fSetSystErrors          ( hRES_1D_Stat );
    //
    // >->-->-> 2-Dimension analysis //
    //
    //  Declaring all histograms
    //
    std::vector<TH1F*>  hRES_2D_Cond1_Stat  =   fEfficiencycorrection   ( hRAW_2D,hREC_1D_in_2D_bin,hGEN_1D_in_2D_bin,f2DCorrection );
    std::vector<TH1F*>  hRES_2D_Cond1_Syst  =   fSetSystErrors          ( hRES_2D_Cond1_Stat );
    
    uIntegralError(hRES_2D_Cond1_Stat,hRES_2D_Cond1_Syst);

    outFile_DT_Yield->Close();
    return;
    //*/
    //(double) 1.6997588e-05
    //
    
    /*
    TFile* f1=TFile::Open("/Users/nikolajal/alidock/AliAnalysisPhiCount/result/yield/SCHistograms.root");//open the file that contains your measured histogram
    TGraphAsymmErrors * gStatistics     =   (TGraphAsymmErrors*) f1->Get("gRES_1D_Stat");
    TGraphAsymmErrors * gSystematics    =   (TGraphAsymmErrors*) f1->Get("gRES_1D_Syst");
    TGraphAsymmErrors * gTotal          =   new TGraphAsymmErrors(*(fSumGraphErrors(gStatistics,gSystematics)));
    TH1F* M = fMakeMeTH1F(gTotal);
    //Note that you may not be able to use a TF1 stored in a file for this purpose.  Please read note 6 above.
    
    TFile* f2=TFile::Open("/Users/nikolajal/alidock/AliAnalysisPhiCount/result/yield/PPReference.root");//open the file that contains your simulated histograms
    TH1F* G=(TH1F*) f2->Get("hGEN_Rw_1D");//get your generated histogram
    TH1F* R=(TH1F*) f2->Get("hREC_Rw_1D");//get your reconstructed histogram
    //get your measured histogram, the systematic uncertainties should be the sum in quadrature of the statistical and systematic uncertainties (but you may want to exclude sources of systematic uncertainty that are correlated between pT bins)

    
    fSetAllFunctions();
    fSetFunction(fLevyTsallis);
    M->Fit(fLevyTsallis);
    TCanvas * c1 = new TCanvas();
    c1->SetLogy();
    M->Draw();
    fLevyTsallis->Draw("same");
    c1->SaveAs("ee.pdf");
    c1->SetLogy(false);
    c1->SetLogx();
    c1->SaveAs("eeee.pdf");
    c1->SetLogx(false);
    c1->SaveAs("eeerre.pdf");
    delete c1;
    
    TFile* f3=new TFile("output_file.root","RECREATE","HistoFile");//open the new file that will contain your output
    ReweightEfficiency(M,fLevyTsallis,G,R,f3,0,2,1.e-3);//run the macro

    //Take the ratio R/G to get the reweighted efficiency; store it in f3.
    //The other objects stored in f3 can be useful, especially for plotting, but are not required.

    f1->Close();
    f2->Close();
    f3->Close();
     */
