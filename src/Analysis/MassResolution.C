#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: All Set!

void MassResolution ( TString fOption = "", bool fSilent = true )   {
    //---------------------//
    //  Setting up input   //
    //---------------------//
    //
    // Silencing warnings for smoother running
    if ( fSilent )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    fChooseOption(fOption);
    //
    // Retrieving Data
    //
    TFile *insFileDT;
    if ( kDoYield )         insFileDT =   new TFile   (Form(kMassResolution_Prod,"Yield"));
    if ( kDoMultiplicity )  insFileDT =   new TFile   (Form(kMassResolution_Prod,"Multiplicity"));
    //
    //---------------------//
    //  Setting up input   //
    //---------------------//
    //
    // Generating the binning array--------------------------------------------------------------------------
    //
    fSetAllBins();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept,  U_nAccep2;
    //
    // Creating the histograms-------------------------------------------------------------------------------

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
        hMRS_1D[iTer]   =   (TH1F*)(insFileDT->Get(hName));
        hName   =   Form("hMDS_1D_%i",iTer);
        hMDS_1D[iTer]   =   (TH1F*)(insFileDT->Get(hName));
        hName   =   Form("hTMD_1D_%i",iTer);
        hTMD_1D[iTer]   =   (TH1F*)(insFileDT->Get(hName));
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
        hMRS_1D_in_2D_bin[iTer] =   (TH1F*)(insFileDT->Get(hName));
        hName   =   Form("hMDS_1D_in_2D_bin_%i",iTer);
        hMDS_1D_in_2D_bin[iTer] =   (TH1F*)(insFileDT->Get(hName));
        hName   =   Form("hTMD_1D_in_2D_bin_%i",iTer);
        hTMD_1D_in_2D_bin[iTer] =   (TH1F*)(insFileDT->Get(hName));
        //
        hMRS_2D[iTer]           =   new TH2F   *[nBinPT2D];
        hMDS_2D[iTer]           =   new TH2F   *[nBinPT2D];
        hTMD_2D[iTer]           =   new TH2F   *[nBinPT2D];
        for ( Int_t jTer = 0; jTer < nBinPT2D; jTer++ )  {
            hName   =   Form("hMRS_2D_%i_%i",iTer,jTer);
            hMRS_2D[iTer][jTer] =   (TH2F*)(insFileDT->Get(hName));
            hName   =   Form("hMDS_2D_%i_%i",iTer,jTer);
            hMDS_2D[iTer][jTer] =   (TH2F*)(insFileDT->Get(hName));
            hName   =   Form("hTMD_2D_%i_%i",iTer,jTer);
            hTMD_2D[iTer][jTer] =   (TH2F*)(insFileDT->Get(hName));
        }
    }
    //
    //---------------------//
    //  Setting up output  //
    //---------------------//
    //
    // Creating the histograms-------------------------------------------------------------------------------
    //
    TH1F   *hSigmaLow_1D;
    TH1F   *hSigmaCnt_1D;
    TH1F   *hSigmaHig_1D;
    TH1F   *hSigmaLow_1D_in_2D_bin;
    TH1F   *hSigmaCnt_1D_in_2D_bin;
    TH1F   *hSigmaHig_1D_in_2D_bin;
    //
    //>>    Defining all histograms
    //
    hName           =   Form("hSigmaLow_1D");
    hTitle          =   Form("hSigmaLow_1D");
    hSigmaLow_1D    =   new TH1F(hName,hTitle,nBinPT1D,fArrPT1D);
    //
    hName           =   Form("hSigmaCnt_1D");
    hTitle          =   Form("hSigmaCnt_1D");
    hSigmaCnt_1D    =   new TH1F(hName,hTitle,nBinPT1D,fArrPT1D);
    //
    hName           =   Form("hSigmaHig_1D");
    hTitle          =   Form("hSigmaHig_1D");
    hSigmaHig_1D    =   new TH1F(hName,hTitle,nBinPT1D,fArrPT1D);
    //
    hName           =   Form("hSigmaLow_1D_in_2D_bin");
    hTitle          =   Form("hSigmaLow_1D_in_2D_bin");
    hSigmaLow_1D_in_2D_bin    =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName           =   Form("hSigmaCnt_1D_in_2D_bin");
    hTitle          =   Form("hSigmaCnt_1D_in_2D_bin");
    hSigmaCnt_1D_in_2D_bin    =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    hName           =   Form("hSigmaHig_1D_in_2D_bin");
    hTitle          =   Form("hSigmaHig_1D_in_2D_bin");
    hSigmaHig_1D_in_2D_bin    =   new TH1F(hName,hTitle,nBinPT2D,fArrPT2D);
    //
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    //
    RooRealVar      InvMass     =   RooRealVar        ("InvMass",   "m_{k^{+}k^{-}}^{REC}",      1.005, 1.0335                  );
    
    RooRealVar sMass, sWidt, sSlop;
                                sMass   =   RooRealVar      ("bMass",   "bMass",    kPhiMesonMass_, kPhiMesonMass_*0.5, kPhiMesonMass_*1.5);
                                sWidt   =   RooRealVar      ("bWidt",   "bWidt",    kPhiMesonWidth);
                                sSlop   =   RooRealVar      ("bSlop",   "bSlop",    0.0015,0.0000,0.0100);
    RooVoigtian                 fSig    =   RooVoigtian     ("fSig",    "fSig",     InvMass,    sMass,  sWidt,  sSlop);
    //
    gROOT->SetBatch(kTRUE);
    for ( Int_t iFit = 0; iFit < nBinPT1D; iFit++ ) {
        //
        //  Recovering Full Histogram Mean and STDV
        auto    kFull_Mean  =   hMRS_1D[iFit]->GetMean();
        auto    kFull_STDV  =   hMRS_1D[iFit]->GetRMS();
        //
        //  Restricting in +-3 Sigma Range
        hMRS_1D[iFit]       ->  GetXaxis()->SetRangeUser( kFull_Mean-(3)*kFull_STDV, kFull_Mean+(3)*kFull_STDV  );
        //
        //  Assigning Central Value
        hSigmaCnt_1D        ->  SetBinContent   (   iFit+1, hMRS_1D[iFit]->GetRMS()                             );
        hSigmaCnt_1D        ->  SetBinError     (   iFit+1, hMRS_1D[iFit]->GetRMSError()                        );
        //
        //  Fitting Gaussian in +-2 Sigma Range
        hMRS_1D[iFit]       ->  GetXaxis()->SetRangeUser( kFull_Mean-(2)*kFull_STDV, kFull_Mean+(2)*kFull_STDV  );
        hMRS_1D[iFit]       ->  Fit("gaus","IMREQ0S","");
        auto    fGaussFitRslt   =   hMRS_1D[iFit]->GetFunction("gaus");
        //
        //  Assigning Low Value
        hSigmaLow_1D        ->  SetBinContent   (   iFit+1, fGaussFitRslt->GetParameter(2)  );
        hSigmaLow_1D        ->  SetBinError     (   iFit+1, fGaussFitRslt->GetParError(2)   );
        //
        //  Fitting Voigtian
        sMass.setVal(kPhiMesonMass_-kFull_Mean);
        RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         InvMass,        Import(*hMDS_1D[iFit])                 );
        auto fFitResults = fSig.fitTo(*data,Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
        auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("bSlop"));
        //
        //  Assigning High Value
        hSigmaHig_1D        ->  SetBinContent   (   iFit+1, N_Raw->getVal()     );
        hSigmaHig_1D        ->  SetBinError     (   iFit+1, N_Raw->getError()   );
        //
        TCanvas            *cDrawPlot   =   new TCanvas();
        gStyle              ->  SetOptStat(0);
        hMRS_1D[iFit]       ->  GetXaxis()->SetRangeUser( -1, -1 );
        hMRS_1D[iFit]       ->  GetXaxis()->SetTitle( "m_{gen}-m_{rec} (GeV/c^{2})" );
        hMRS_1D[iFit]       ->  GetYaxis()->SetTitle( "Entries" );
        hMRS_1D[iFit]       ->  SetTitle( "" );
        hMRS_1D[iFit]       ->  Draw("SAME");
        fGaussFitRslt       ->  SetRange( kFull_Mean-(2)*kFull_STDV, kFull_Mean+(2)*kFull_STDV  );
        fGaussFitRslt       ->  Draw("SAME");
        if ( kDoYield )         cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Yield")+TString(Form("3SigGuassFit_1D_%i.pdf",iFit)));
        if ( kDoMultiplicity )  cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString(Form("3SigGuassFit_1D_%i.pdf",iFit)));
        delete              cDrawPlot;
        //
                            cDrawPlot   =   new TCanvas();
        auto fSaveToFrame   =   InvMass.frame(Name(""),Title(""));
        data                ->plotOn(fSaveToFrame,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
        fSig                .plotOn (fSaveToFrame,      LineColor(4),                   LineStyle(kSolid), Name("RooMod"));
        fSaveToFrame        ->  SetTitle("");
        fSaveToFrame        ->  Draw();
        if ( kDoYield )         cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Yield")+TString(Form("MDSVoigtFit_1D_%i.pdf",iFit)));
        if ( kDoMultiplicity )  cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString(Form("MDSVoigtFit_1D_%i.pdf",iFit)));
        delete              cDrawPlot;
    }
    for ( Int_t iFit = 0; iFit < nBinPT2D; iFit++ ) {
        //
        //  Recovering Full Histogram Mean and STDV
        auto    kFull_Mean  =   hMRS_1D_in_2D_bin[iFit]->GetMean();
        auto    kFull_STDV  =   hMRS_1D_in_2D_bin[iFit]->GetRMS();
        //
        //  Restricting in +-3 Sigma Range
        hMRS_1D_in_2D_bin[iFit] ->  GetXaxis()->SetRangeUser( kFull_Mean-(3)*kFull_STDV, kFull_Mean+(3)*kFull_STDV );
        //
        //  Assigning Central Value
        hSigmaCnt_1D_in_2D_bin               ->  SetBinContent   (   iFit+1, hMRS_1D_in_2D_bin[iFit]->GetRMS()      );
        hSigmaCnt_1D_in_2D_bin               ->  SetBinError     (   iFit+1, hMRS_1D_in_2D_bin[iFit]->GetRMSError() );
        //
        //  Fitting Gaussian in +-2 Sigma Range
        hMRS_1D_in_2D_bin[iFit] ->  GetXaxis()->SetRangeUser( kFull_Mean-(2)*kFull_STDV, kFull_Mean+(2)*kFull_STDV );
        hMRS_1D_in_2D_bin[iFit] ->  Fit("gaus","IMREQ0S","");
        auto    fGaussFitRslt   =   hMRS_1D_in_2D_bin[iFit]->GetFunction("gaus");
        //
        //  Assigning Low Value
        hSigmaLow_1D_in_2D_bin               ->  SetBinContent   (   iFit+1, fGaussFitRslt->GetParameter(2)  );
        hSigmaLow_1D_in_2D_bin               ->  SetBinError     (   iFit+1, fGaussFitRslt->GetParError(2)   );
        //
        //  Fitting Voigtian
        sMass.setVal(kPhiMesonMass_-kFull_Mean);
        RooDataHist*    data        =   new RooDataHist   ("Data",      "Data",         InvMass,        Import(*hMDS_1D_in_2D_bin[iFit])                 );
        auto fFitResults = fSig.fitTo(*data,Save(),NumCPU(kCPU_use,kCPUStrategy),Offset(kFitOffset),Strategy(kFitMinuitStrategy),InitialHesse(kFitInitHesse),Minos(kFitMinos));
        auto N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("bSlop"));
        //
        //  Assigning High Value
        hSigmaHig_1D_in_2D_bin        ->  SetBinContent   (   iFit+1, N_Raw->getVal()     );
        hSigmaHig_1D_in_2D_bin        ->  SetBinError     (   iFit+1, N_Raw->getError()   );
        //
        TCanvas            *cDrawPlot   =   new TCanvas();
        gStyle              ->  SetOptStat(0);
        hMRS_1D_in_2D_bin[iFit]       ->  GetXaxis()->SetRangeUser( -1, -1 );
        hMRS_1D_in_2D_bin[iFit]       ->  GetXaxis()->SetTitle( "m_{gen}-m_{rec} (GeV/c^{2})" );
        hMRS_1D_in_2D_bin[iFit]       ->  GetYaxis()->SetTitle( "Entries" );
        hMRS_1D_in_2D_bin[iFit]       ->  SetTitle( "" );
        hMRS_1D_in_2D_bin[iFit]       ->  Draw("SAME");
        fGaussFitRslt       ->  SetRange( kFull_Mean-(2)*kFull_STDV, kFull_Mean+(2)*kFull_STDV  );
        fGaussFitRslt       ->  Draw("SAME");
        if ( kDoYield )         cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Yield")+TString(Form("3SigGuassFit_2D_%i.pdf",iFit)));
        if ( kDoMultiplicity )  cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Mutliplicity")+TString(Form("3SigGuassFit_2D_%i.pdf",iFit)));
        delete              cDrawPlot;
        //
                            cDrawPlot   =   new TCanvas();
        auto fSaveToFrame   =   InvMass.frame(Name(""),Title(""));
        data                ->plotOn(fSaveToFrame,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
        fSig                .plotOn (fSaveToFrame,      LineColor(4),                   LineStyle(kSolid), Name("RooMod"));
        fSaveToFrame        ->  Draw();
        if ( kDoYield )         cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Yield")+TString(Form("MDSVoigtFit_2D_%i.pdf",iFit)));
        if ( kDoMultiplicity )  cDrawPlot           ->  SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString(Form("MDSVoigtFit_2D_%i.pdf",iFit)));
        delete              cDrawPlot;
    }
    //
    //--------------------------//
    // PostProcessin output obj //
    //--------------------------//
    //
    hSigmaLow_1D->Scale(1000);
    hSigmaCnt_1D->Scale(1000);
    hSigmaHig_1D->Scale(1000);
    hSigmaLow_1D_in_2D_bin->Scale(1000);
    hSigmaCnt_1D_in_2D_bin->Scale(1000);
    hSigmaHig_1D_in_2D_bin->Scale(1000);
    //
    if ( kDoYield )         gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Plot,"Yield")));
    if ( kDoMultiplicity )  gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Plot,"Multiplicity")));
    TCanvas        *cAllResolutions_1D  =   new TCanvas();
    gPad            ->  SetLogx();
    //
    TLegend        *lLegend     =   new TLegend(0.9,0.1,0.7,0.2);
    lLegend         ->  AddEntry(hSigmaCnt_1D,"3\sigma Cut RMS","EP");
    lLegend         ->  AddEntry(hSigmaLow_1D,"2\sigma Cut Gauss Fit","EP");
    lLegend         ->  AddEntry(hSigmaHig_1D,"Rec Mass Dist. Voigtian Fit","EP");
    //
    //  Central Resolution
    hSigmaCnt_1D->SetMaximum(2.2);
    hSigmaCnt_1D->SetMinimum(0.9);
    hSigmaCnt_1D->SetTitle("");
    hSigmaCnt_1D->GetXaxis()->SetTitle("Transverse Momentum (GeV/c^{2})");
    hSigmaCnt_1D->GetYaxis()->SetTitle("#sigma (MeV/c^{2})");
    hSigmaCnt_1D->SetMarkerStyle(22);
    hSigmaCnt_1D->SetMarkerColor(kRed);
    hSigmaCnt_1D->Draw();
    //
    //  Low Resolution
    hSigmaLow_1D->SetMarkerStyle(23);
    hSigmaLow_1D->SetMarkerColor(kBlue);
    hSigmaLow_1D->Draw("SAME");
    //
    //  High Resolution
    hSigmaHig_1D->SetMarkerStyle(21);
    hSigmaHig_1D->SetMarkerColor(kGreen-2);
    hSigmaHig_1D->Draw("SAME");
    //
    lLegend         ->Draw("SAME");
    if ( kDoYield )         cAllResolutions_1D->SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString("cAllResolutions_1D.pdf"));
    if ( kDoMultiplicity )  cAllResolutions_1D->SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString("cAllResolutions_1D.pdf"));
    delete      cAllResolutions_1D;
    //
                cAllResolutions_1D  =   new TCanvas();
    gPad            ->  SetLogx();
    //
    auto        hSigmaCnt_1D_N  =   (TH1F*)hSigmaCnt_1D->Clone();
    auto        hSigmaLow_1D_N  =   (TH1F*)hSigmaLow_1D->Clone();
    auto        hSigmaHig_1D_N  =   (TH1F*)hSigmaHig_1D->Clone();
    //
    hSigmaCnt_1D_N->Divide(hSigmaCnt_1D,hSigmaCnt_1D);
    hSigmaLow_1D_N->Divide(hSigmaLow_1D,hSigmaCnt_1D);
    hSigmaHig_1D_N->Divide(hSigmaHig_1D,hSigmaCnt_1D);
    //
    hSigmaCnt_1D_N->SetMaximum(1.25);
    hSigmaCnt_1D_N->SetMinimum(0.75);
    hSigmaCnt_1D_N->GetYaxis()->SetTitle("#sigma/#sigma_{cnt} (MeV/c^{2})");
    hSigmaCnt_1D_N->Draw();
    hSigmaLow_1D_N->Draw("SAME");
    hSigmaHig_1D_N->Draw("SAME");
    //
    lLegend         ->Draw("SAME");
    if ( kDoYield )         cAllResolutions_1D->SaveAs(Form(kMassResolution_Plot,"Yield")+TString("cAllResolutions_1D_N.pdf"));
    if ( kDoMultiplicity )  cAllResolutions_1D->SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString("cAllResolutions_1D_N.pdf"));
    delete      cAllResolutions_1D;
    //
    if ( kDoYield )         gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Plot,"Yield")));
    if ( kDoMultiplicity )  gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Plot,"Multiplicity")));
    TCanvas    *cAllResolutions_1D_in_2D_bin  =   new TCanvas();
    gPad            ->  SetLogx();
    //
    //  Central Resolution
    hSigmaCnt_1D_in_2D_bin->SetMaximum(2.2);
    hSigmaCnt_1D_in_2D_bin->SetMinimum(0.9);
    hSigmaCnt_1D_in_2D_bin->SetTitle("");
    hSigmaCnt_1D_in_2D_bin->GetXaxis()->SetTitle("Transverse Momentum (GeV/c^{2})");
    hSigmaCnt_1D_in_2D_bin->GetYaxis()->SetTitle("#sigma (MeV/c^{2})");
    hSigmaCnt_1D_in_2D_bin->SetMarkerStyle(22);
    hSigmaCnt_1D_in_2D_bin->SetMarkerColor(kRed);
    hSigmaCnt_1D_in_2D_bin->Draw();
    //
    //  Low Resolution
    hSigmaLow_1D_in_2D_bin->SetMarkerStyle(23);
    hSigmaLow_1D_in_2D_bin->SetMarkerColor(kBlue);
    hSigmaLow_1D_in_2D_bin->Draw("SAME");
    //
    //  High Resolution
    hSigmaHig_1D_in_2D_bin->SetMarkerStyle(21);
    hSigmaHig_1D_in_2D_bin->SetMarkerColor(kGreen-2);
    hSigmaHig_1D_in_2D_bin->Draw("SAME");
    //
    lLegend         ->Draw("SAME");
    if ( kDoYield )         cAllResolutions_1D_in_2D_bin->SaveAs(Form(kMassResolution_Plot,"Yield")+TString("cAllResolutions_1D_in_2D_bin.pdf"));
    if ( kDoMultiplicity )  cAllResolutions_1D_in_2D_bin->SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString("cAllResolutions_1D_in_2D_bin.pdf"));
    delete      cAllResolutions_1D_in_2D_bin;
    //
                cAllResolutions_1D_in_2D_bin  =   new TCanvas();
    gPad            ->  SetLogx();
    //
    auto        hSigmaCnt_1D_in_2D_bin_N  =   (TH1F*)hSigmaCnt_1D_in_2D_bin->Clone();
    auto        hSigmaLow_1D_in_2D_bin_N  =   (TH1F*)hSigmaLow_1D_in_2D_bin->Clone();
    auto        hSigmaHig_1D_in_2D_bin_N  =   (TH1F*)hSigmaHig_1D_in_2D_bin->Clone();
    //
    hSigmaCnt_1D_in_2D_bin_N->Divide(hSigmaCnt_1D_in_2D_bin,hSigmaCnt_1D_in_2D_bin);
    hSigmaLow_1D_in_2D_bin_N->Divide(hSigmaLow_1D_in_2D_bin,hSigmaCnt_1D_in_2D_bin);
    hSigmaHig_1D_in_2D_bin_N->Divide(hSigmaHig_1D_in_2D_bin,hSigmaCnt_1D_in_2D_bin);
    //
    hSigmaCnt_1D_in_2D_bin_N->SetMaximum(1.25);
    hSigmaCnt_1D_in_2D_bin_N->SetMinimum(0.75);
    hSigmaCnt_1D_in_2D_bin_N->GetYaxis()->SetTitle("#sigma/#sigma_{cnt} (MeV/c^{2})");
    hSigmaCnt_1D_in_2D_bin_N->Draw();
    hSigmaLow_1D_in_2D_bin_N->Draw("SAME");
    hSigmaHig_1D_in_2D_bin_N->Draw("SAME");
    //
    lLegend         ->Draw("SAME");
    if ( kDoYield )         cAllResolutions_1D_in_2D_bin->SaveAs(Form(kMassResolution_Plot,"Yield")+TString("cAllResolutions_1D_in_2D_bin_N.pdf"));
    if ( kDoMultiplicity )  cAllResolutions_1D_in_2D_bin->SaveAs(Form(kMassResolution_Plot,"Multiplicity")+TString("cAllResolutions_1D_in_2D_bin_N.pdf"));
    delete      cAllResolutions_1D_in_2D_bin;
    //
    gROOT->SetBatch(kFALSE);
    //
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    //
    // >> All Analysis Utility
    //
    if ( kDoYield ) {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Dir_,"Yield")));
        TFile *outFil0  =   new TFile   (Form(kMassResolution_Anal,"Yield"),"recreate");
        //
        hSigmaLow_1D->Write();
        hSigmaCnt_1D->Write();
        hSigmaHig_1D->Write();
        hSigmaLow_1D_in_2D_bin->Write();
        hSigmaCnt_1D_in_2D_bin->Write();
        hSigmaHig_1D_in_2D_bin->Write();
        //
        outFil0->Close();
    }
    if ( kDoMultiplicity ) {
        gROOT           ->  ProcessLine(Form(".! mkdir -p %s",Form(kMassResolution_Dir_,"Multiplicity")));
        TFile *outFil0  =   new TFile   (Form(kMassResolution_Anal,"Multiplicity"),"recreate");
        //
        hSigmaLow_1D->Write();
        hSigmaCnt_1D->Write();
        hSigmaHig_1D->Write();
        hSigmaLow_1D_in_2D_bin->Write();
        hSigmaCnt_1D_in_2D_bin->Write();
        hSigmaHig_1D_in_2D_bin->Write();
        //
        outFil0->Close();
    }
    
    //
    // >-> Close input File
    //
    insFileDT->Close();
    //    
}
