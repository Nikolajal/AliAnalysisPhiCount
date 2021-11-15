#include "../inc/AliAnalysisPhiPair.h"
#include "./Analysis/SignalExtraction.C"
#include "./PreProcessing.C"
// !TODO: Can there really be a TODO in a Test Macro?

void
TestMacro
()  {
    fSetAllBins();
    //
    //
    //
    //
    
    TFile*  kINFILE_HI   =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Reference_HI_MB.root");
    TProfile * hTest = (TProfile*)(kINFILE_HI->Get("hHISysErr"));
    TH1F* hTest2 = new TH1F("","",nBinPT1D,fArrPT1D);
    
    uRebin(hTest2,hTest);
    
    hTest2->Draw();
    hTest->Draw("SAME");
    return;
    
    std::vector<TString> fInputFiles;
    for ( auto kOption : kSyst_SEX_1D_Options ) {
        fInputFiles.push_back( TString( Form( "/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/SignalExtraction/Systematics/ExtractionCheck/%s/1D/FitResults_%s.root" , kOption.Data(), kOption.Data() ) ) );
    }
    uCompareResultsTH1F( fInputFiles, "hRAW_1D", kSyst_SEX_1D_Legend );
    
    return;
    
    TFile* kINFILE_1  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/PreProcessing/IM_MonteCarloTruth.root");
    TFile* kINFILE_2  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/Systematics/MaterialBudget/PreProcessing/IM_MonteCarloTruth_902.root");
    TFile* kINFILE_3  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/Systematics/MaterialBudget/PreProcessing/IM_MonteCarloTruth_900.root");
    TFile* kINFILE_4  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/Yield/Systematics/MaterialBudget/PreProcessing/IM_MonteCarloTruth_900.root");
    
    TH1F*   hEFF_FLL    =   (TH1F*)(kINFILE_1->Get("hEFF_1D"));
    TH1F*   hEFF_200    =   (TH1F*)(kINFILE_2->Get("hEFF_1D"));
    TH1F*   hEFF_201    =   (TH1F*)(kINFILE_3->Get("hEFF_1D"));
    
    TH1F*   hREC_FLL    =   (TH1F*)(kINFILE_1->Get("hREC_Rw_1D"));
    TH1F*   hREC_200    =   (TH1F*)(kINFILE_2->Get("hREC_Rw_1D"));
    TH1F*   hREC_201    =   (TH1F*)(kINFILE_3->Get("hREC_Rw_1D"));
    TH1F*   hREC_202    =   (TH1F*)(kINFILE_4->Get("hREC_Rw_1D"));
    TH1F*   hGEN_FLL    =   (TH1F*)(kINFILE_1->Get("hGEN_Rw_1D"));
    TH1F*   hGEN_200    =   (TH1F*)(kINFILE_2->Get("hGEN_Rw_1D"));
    TH1F*   hGEN_201    =   (TH1F*)(kINFILE_3->Get("hGEN_Rw_1D"));
    TH1F*   hGEN_202    =   (TH1F*)(kINFILE_4->Get("hGEN_Rw_1D"));
    
    auto hRatio1 = (TH1F*)hEFF_FLL->Clone();
    auto hRatio2 = (TH1F*)hEFF_FLL->Clone();
    auto hRatio3 = (TH1F*)hEFF_FLL->Clone();
    
    auto hRatio4 = (TH1F*)hREC_FLL->Clone();
    auto hRatio5 = (TH1F*)hREC_FLL->Clone();
    auto hRatio6 = (TH1F*)hREC_FLL->Clone();
    auto hRatio7 = (TH1F*)hREC_FLL->Clone();
    auto hRatio8 = (TH1F*)hREC_FLL->Clone();
    auto hRatio9 = (TH1F*)hREC_FLL->Clone();
    
    hRatio1->Divide(hEFF_FLL,hEFF_200,1.,1.,"b");
    hRatio2->Divide(hEFF_FLL,hEFF_201,1.,1.,"b");
    hRatio3->Divide(hEFF_200,hEFF_201,1.,1.,"b");
    
    hRatio4->Divide(hREC_FLL,hGEN_FLL,1.,1.,"b");
    hRatio5->Divide(hREC_200,hGEN_200,1.,1.,"b");
    hRatio6->Divide(hREC_201,hGEN_201,1.,1.,"b");
    
    hRatio7->Divide(hRatio4,hRatio5,1.,1.,"b");
    hRatio8->Divide(hRatio4,hRatio6,1.,1.,"b");
    hRatio9->Divide(hRatio5,hRatio6,1.,1.,"b");

    uOffset(hRatio1,-1,true);
    uOffset(hRatio2,-1,true);
    //uOffset(hRatio3,-1,true);
    
    TH1F* hReference1D = new TH1F("hReference1D","",nBinPT1D,fArrPT1D);
    TH1F* hReference2D = new TH1F("hReference2D","",nBinPT2D,fArrPT2D);
    
    
    TF1* fCorr = new TF1("fCorr","1 + [0] * TMath::Exp( [1] + x * [2] )",0,100);
    fCorr->SetParameter(0,+1.00);
    fCorr->SetParameter(1,-0.01);
    fCorr->SetParameter(2,+0.03);
    
    TCanvas* c1 = new TCanvas("","",1600,400);
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogx();
    hRatio7->Fit(fCorr,"IMEQS","R",0.,5.);
    hRatio7->DrawCopy("HIST");
    fCorr->DrawCopy("SAME");
    for ( Int_t i = 1; i <= nBinPT1D; i++ ) {
        auto kRefErr = fCorr->Integral( fArrPT1D[i-1], fArrPT1D[i] ) / ( fArrPT1D[i] - fArrPT1D[i-1] );
        hReference1D->SetBinContent( i, kRefErr-1 );
    }
    c1->cd(4);
    gPad->SetLogx();
    hReference1D->DrawCopy("HIST");
    fCorr->SetParameter(0,+1.00);
    fCorr->SetParameter(1,-0.01);
    fCorr->SetParameter(2,+0.03);
    c1->cd(2);
    gPad->SetLogx();
    hRatio8->Fit(fCorr,"IMEQS","R",0.,5.);
    hRatio8->DrawCopy("HIST");
    fCorr->DrawCopy("SAME");
    for ( Int_t i = 1; i <= nBinPT1D; i++ ) {
        auto kRefErr = fCorr->Integral( fArrPT1D[i-1], fArrPT1D[i] ) / ( fArrPT1D[i] - fArrPT1D[i-1] );
        hReference1D->SetBinContent( i, kRefErr-1 );
    }
    c1->cd(5);
    gPad->SetLogx();
    hReference1D->DrawCopy("HIST");
    fCorr->SetParameter(0,+1.00);
    fCorr->SetParameter(1,-0.01);
    fCorr->SetParameter(2,+0.03);
    c1->cd(3);
    gPad->SetLogx();
    hRatio9->Fit(fCorr,"IMEQS","R",0.,5.);
    hRatio9->DrawCopy("HIST");
    fCorr->DrawCopy("SAME");
    for ( Int_t i = 1; i <= nBinPT1D; i++ ) {
        auto kRefErr = fCorr->Integral( fArrPT1D[i-1], fArrPT1D[i] ) / ( fArrPT1D[i] - fArrPT1D[i-1] );
        hReference1D->SetBinContent( i, kRefErr-1 );
    }
    c1->cd(6);
    gPad->SetLogx();
    hReference1D->DrawCopy("HIST");
    
    return;
    
    gROOT->ProcessLine(".! $ROOT_SYS hadd -f ./result/MCout.root ./result/MC_Production/outGeneratorMC_00*.root");
    auto kREBIN = 1;
    TFile* kINFILE2  =   new TFile   ("/Users/nikolajal/alice/AliAnalysisPhiCount/result/MCout.root");
    auto hPhi1 = (TProfile*)(kINFILE2->Get("hYPhiPro1"));
    auto hPhi2 = (TProfile*)(kINFILE2->Get("hYPhiPro2"));
    hPhi1->Rebin(kREBIN);
    hPhi2->Rebin(kREBIN);
    cout << hPhi1->Integral("width") << endl;
    auto hPhiGamma  = new TH1F("hPhiGamma", "hPhiGamma",    (int)(100./kREBIN),0.,100);
    auto hPhiSigma  = new TH1F("hPhiSigma", "hPhiSigma",    (int)(100./kREBIN),0.,100);
    auto hPhiR2     = new TH1F("hPhiR2",    "hPhiR2",       (int)(100./kREBIN),0.,100);
    auto hPhiR1     = new TH1F("hPhiR1",    "hPhiR1",       (int)(100./kREBIN),0.,100);
    
    for ( Int_t iTer = 1; iTer <= (int)(100./kREBIN); iTer++ )  {
        auto nPhi1 = hPhi1->GetBinContent(iTer);
        auto nPhi2 = hPhi2->GetBinContent(iTer);
        auto nPhE1 = hPhi1->GetBinError(iTer);
        auto nPhE2 = hPhi2->GetBinError(iTer);
        if ( nPhi1 == 0 ) continue;
        if ( nPhi2 == 0 ) continue;
        hPhiGamma->SetBinContent( iTer, fGammaPhiValue( nPhi1, nPhi2 ) );
        hPhiGamma->SetBinError  ( iTer, fGammaPhiError( nPhi1, nPhi2, nPhE1, nPhE2 ) );
        hPhiSigma->SetBinContent( iTer, fSigmaPhiValue( nPhi1, nPhi2 ) );
        hPhiSigma->SetBinError  ( iTer, fSigmaPhiError( nPhi1, nPhi2, nPhE1, nPhE2 ) );
        hPhiR2->SetBinContent   ( iTer, nPhi2/(nPhi1*nPhi1) );
        hPhiR2->SetBinError     ( iTer, (nPhi2/(nPhi1*nPhi1))*sqrt( nPhE2*nPhE2/(nPhi2*nPhi2) + 4*nPhE1*nPhE1/(nPhi1*nPhi1) ) );
        hPhiR1->SetBinContent   ( iTer, nPhi2/(nPhi1) );
        hPhiR1->SetBinError     ( iTer, (nPhi2/(nPhi1))*sqrt( nPhE2*nPhE2/(nPhi2*nPhi2) + nPhE1*nPhE1/(nPhi1*nPhi1) ) );
    }
    
    hPhi1->SetLineColor(2);
    hPhi1->SetLineStyle(5);
    hPhi1->SetLineWidth(3);
    
    hPhi2->SetLineColor(2);
    hPhi2->SetLineStyle(5);
    hPhi2->SetLineWidth(3);
    
    hPhiR1->SetLineColor(2);
    hPhiR1->SetLineStyle(5);
    hPhiR1->SetLineWidth(3);
    
    hPhiR2->SetLineColor(2);
    hPhiR2->SetLineStyle(5);
    hPhiR2->SetLineWidth(3);
    
    hPhiGamma->SetLineColor(2);
    hPhiGamma->SetLineStyle(5);
    hPhiGamma->SetLineWidth(3);
    
    hPhiSigma->SetLineColor(2);
    hPhiSigma->SetLineStyle(5);
    hPhiSigma->SetLineWidth(3);
    
    TFile* kINFILE  =   new TFile   (Form(kASigExtp_FitCheckRst,"Multiplicity"));
    float    kNCH[]  =   {5.398,9.44,13.13,17.58,24};
    char*    kName[] =   {"1D","2D","R1","R2","P1","P2"};
    char*    kLabel[]=    {"#frac{dN_{#phi}}{dy}","#frac{dN_{#phi#phi}}{dy}","#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT}","#frac{#LT Y_{#phi#phi} #GT}{#LT Y_{#phi} #GT^{2}}","#sigma^{2}_{#phi}","#gamma_{#phi}"};
    
    
    std::vector<TGraphErrors*> fOutput;
    for ( int i = 0; i < 6; i++ )   {
        auto fCurrentTH1F = (TH1D*)(kINFILE->Get(Form("hShow%s",kName[i])));
        fOutput.push_back(new TGraphErrors());
        for ( int j = 0; j < 5; j++ )   {
            fOutput.at(i)->SetPoint(j,kNCH[j],fCurrentTH1F->GetBinContent(j+1));
            fOutput.at(i)->SetPointError(j,0.,fCurrentTH1F->GetBinError(j+1));
        }
        fOutput.at(i)->SetPoint(5,-1,fCurrentTH1F->GetBinContent(6));
        fOutput.at(i)->SetPointError(5,0.,fCurrentTH1F->GetBinError(6));
        fOutput.at(i)->SetMinimum(0);
        fOutput.at(i)->SetMarkerStyle(4);
        fOutput.at(i)->SetMarkerColor(4);
        fOutput.at(i)->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        if ( i == 5 ) fOutput.at(i)->SetMaximum(0.06);
        if ( i == 2 ) fOutput.at(i)->SetMaximum(0.08);
        if ( i == 3 ) fOutput.at(i)->SetMinimum(0.40);
        if ( i == 3 ) fOutput.at(i)->SetMaximum(1.50);
    }
    /*
    fOutput.push_back(new TGraphErrors());
    for ( int j = 0; j < 5; j++ )   {
        double fY1D, fY2D, fDump;
        fOutput.at(0)->GetPoint(j,fDump,fY1D);
        fOutput.at(4)->GetPoint(j,fDump,fY2D);
        
        fOutput.at(6)->SetPoint(j,kNCH[j],fY1D/fY2D);
        fOutput.at(6)->SetPointError(j,0.,0);
        fOutput.at(6)->SetMarkerStyle(markers[2+(6%4)]);
        fOutput.at(6)->SetMarkerColor(colors[1+(6%3)]);
    }
     */
    
    TCanvas* cTest = new TCanvas("","",1000,700);
    cTest->Divide(3,2);
    auto iTer = 1;
    for ( auto kGraph : fOutput )    {
        cTest->cd(iTer);
        kGraph->Draw("APE");
        uLatex->DrawLatexNDC(0.18,0.80,kLabel[iTer-1]);
        if ( iTer == 1 )    {
            hPhi1->DrawCopy("same HIST L ");
            hPhi1->SetFillColorAlpha(3,0.33);
            hPhi1->DrawCopy("same E3L ");
        }
        if ( iTer == 2 )    {
            hPhi2->DrawCopy("same HIST L ");
            hPhi2->SetFillColorAlpha(3,0.33);
            hPhi2->DrawCopy("same E3L ");
        }
        if ( iTer == 3 )    {
            hPhiR1->DrawCopy("same HIST L ");
            hPhiR1->SetFillColorAlpha(3,0.33);
            hPhiR1->DrawCopy("same E3L ");
        }
        if ( iTer == 4 )    {
            hPhiR2->DrawCopy("same HIST L ");
            hPhiR2->SetFillColorAlpha(3,0.33);
            hPhiR2->DrawCopy("same E3L ");
        }
        if ( iTer == 5 )    {
            hPhiSigma->DrawCopy("same HIST L ");
            hPhiSigma->SetFillColorAlpha(3,0.33);
            hPhiSigma->DrawCopy("same E3L ");
        }
        if ( iTer == 6 )    {
            hPhiGamma->DrawCopy("same HIST L ");
            hPhiGamma->SetFillColorAlpha(3,0.33);
            hPhiGamma->DrawCopy("same E3L ");
        }
        iTer++;
    }
    /*
    auto fResult = fOutput.at(0)->Fit("pol1","IMEQ0S");
    auto k1D_0  =   fResult->GetParams()[0];
    auto k1D_1  =   fResult->GetParams()[1];
    fResult = fOutput.at(1)->Fit("pol2","IMEQ0S");
    auto k2D_0  =   fResult->GetParams()[0];
    auto k2D_1  =   fResult->GetParams()[1];
    auto k2D_2  =   fResult->GetParams()[2];
    TH1F* cPrediction = new TH1F("","",30,0,30);
    for ( int i = 1; i <= 30; i++ )  {
        cPrediction->SetBinContent(i,fGammaPhiValue( k1D_0 + k1D_1*i, k2D_0 + k2D_1*i + k2D_2*i*i ) );
    }
    cPrediction->SetLineColor(colors[1]);
    cTest->cd(6);
    cPrediction->SetMinimum(0);
    cPrediction->Draw("HIST MIN0");
    fOutput.at(5)->Draw("EP SAME");
    cTest->cd(1);
    TMultiGraph* cmulti = new TMultiGraph();
    cmulti->Add(fOutput.at(0),"EP");
    cmulti->Add(fOutput.at(4),"EP");
    cmulti->Draw("APE");
     */
}
