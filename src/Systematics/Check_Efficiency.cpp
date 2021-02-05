#include "../../inc/AliAnalysisPhiPair.h"
// !TODO: Split

void Check_Efficiency ( TString fFileName1 = "", TString fFileName2 = "", TString fFileName3 = "", TString fFileName4 = "", TString fFileName5 = "", TString fFileName6 = "", TString fFileName7 = "" )
{
    //---------------------//
    //  Setting up input   //
    //---------------------//
    //
    // >-> OPTIONS
    //
    TFile     **fInputFiles =   new TFile*  [7];
    fInputFiles[0]  =   new TFile   (fFileName1.Data());
    fInputFiles[1]  =   new TFile   (fFileName2.Data());
    fInputFiles[2]  =   new TFile   (fFileName3.Data());
    fInputFiles[3]  =   new TFile   (fFileName4.Data());
    fInputFiles[4]  =   new TFile   (fFileName5.Data());
    fInputFiles[5]  =   new TFile   (fFileName6.Data());
    fInputFiles[6]  =   new TFile   (fFileName7.Data());
    
    Int_t   iFile       =   0;
    TH1F  **hTest       =   new TH1F*   [7+1];
    for ( Int_t iFil2 = 0; iFil2 < 4; iFil2++ ) {
        hTest[iFile] =   (TH1F*)(fInputFiles[iFil2]->Get("hEFF_1D"));
        if ( !hTest[iFile] ) continue;
        iFile++;
    }
    
    Int_t   fValidFiles =   iFile;
    
    if ( fValidFiles <= 1)  {
        cout << "[ERROR] Must Specify at least two valid input root file" << endl;
        cout << "[INFO] Usage PreProcessing(\"file1.root\",\"file2.root\",...)" << endl;
        return;
    }
    
    iFile   =   1;
    TH1D      **hREC_1D =   new TH1D*   [fValidFiles+1];
    TH1D      **hGEN_1D =   new TH1D*   [fValidFiles+1];
    for ( Int_t iFil2 = 1; iFil2 <= 7; iFil2++ ) {
        hREC_1D[iFile] =   (TH1D*)(fInputFiles[iFil2-1]->Get("hREC_1D"));
        hGEN_1D[iFile] =   (TH1D*)(fInputFiles[iFil2-1]->Get("hGEN_1D"));
        if ( !hREC_1D[iFile] ) continue;
        hREC_1D[0] =   (TH1D*)(fInputFiles[iFil2-1]->Get("hREC_1D"));
        hGEN_1D[0] =   (TH1D*)(fInputFiles[iFil2-1]->Get("hGEN_1D"));
        iFile++;
    }
    
    hREC_1D[0]->Clear();
    hGEN_1D[0]->Clear();
    
    for ( Int_t iFil2 = 1; iFil2 <= fValidFiles; iFil2++ ) {
        hREC_1D[0]->Add(hREC_1D[iFil2]);
        hGEN_1D[0]->Add(hGEN_1D[iFil2]);
    }
    
    TH1D  **hEFF_1D =   new TH1D*   [fValidFiles];
    for ( Int_t iFil2 = 0; iFil2 <= fValidFiles; iFil2++ ) {
        hEFF_1D[iFil2]  =   new TH1D(*hREC_1D[iFil2]);
        hEFF_1D[iFil2]->Divide(hREC_1D[iFil2],hGEN_1D[iFil2],1.,1.,"b");
    }
    
    TCanvas * cCheckEff = new TCanvas();
    
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1, 1.0);
    gStyle  ->SetOptStat(0);
    pad1    ->SetBottomMargin(0); // Upper and lower plot are joined
    pad1    ->Draw();             // Draw the upper pad: pad1
    pad1    ->cd();
    
    hEFF_1D[0]->SetMaximum(1.);
    hEFF_1D[0]->SetMinimum(0.);
    hEFF_1D[0]->SetLineWidth(1);
    hEFF_1D[0]->SetLineColor(1);
    hEFF_1D[0]->SetLineStyle(1);
    hEFF_1D[0]->SetFillColorAlpha(2,0.5);
    hEFF_1D[0]->GetXaxis()->SetTitle("");
    hEFF_1D[0]->SetOption("LE3");
    
    hEFF_1D[0]->Draw();
    for ( Int_t iFil2 = 1; iFil2 <= fValidFiles; iFil2++ ) {
        hEFF_1D[iFil2]  ->SetMarkerStyle(27);
        hEFF_1D[iFil2]  ->SetMarkerColor(41+iFil2);
        hEFF_1D[iFil2]  ->SetFillColorAlpha(0,0.);
        hEFF_1D[iFil2]  ->SetLineColorAlpha(41+iFil2,1.);
        hEFF_1D[iFil2]  ->SetLineWidth(1.);
        hEFF_1D[iFil2]  ->Draw("SAME L P E1 X0");
    }
    
    TH1D  **hEFF_1D_Norm =   new TH1D*   [fValidFiles+1];
    for ( Int_t iFil2 = 0; iFil2 <= fValidFiles; iFil2++ ) {
        hEFF_1D_Norm[iFil2]  =   new TH1D(*hEFF_1D[0]);
        hEFF_1D_Norm[iFil2]->Divide(hEFF_1D[iFil2],hEFF_1D[0]);
    }
    
    cCheckEff->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0., 0.05, 1, 0.3);
    gStyle  ->SetOptStat(0);
    pad2    ->SetTopMargin(0);
    pad2    ->SetBottomMargin(0.2);
    pad2    ->Draw();
    pad2    ->cd();
    
    // - // XAxis
    hEFF_1D_Norm[0]->GetXaxis()->SetTitle("P_{T}#phi_{1} (GeV/c)");
    hEFF_1D_Norm[0]->GetXaxis()->SetTitleSize(20);
    hEFF_1D_Norm[0]->GetXaxis()->SetTitleFont(43);
    hEFF_1D_Norm[0]->GetXaxis()->SetTitleOffset(3.);
    hEFF_1D_Norm[0]->GetXaxis()->SetLabelFont(43);
    hEFF_1D_Norm[0]->GetXaxis()->SetLabelSize(15);
    hEFF_1D_Norm[0]->SetMaximum(1.25);
    hEFF_1D_Norm[0]->SetMinimum(0.75);
    // - // YAxis
    hEFF_1D_Norm[0]->GetYaxis()->SetTitle("MODEL / DATA");
    hEFF_1D_Norm[0]->GetYaxis()->SetNdivisions(505);
    hEFF_1D_Norm[0]->GetYaxis()->SetTitleSize(20);
    hEFF_1D_Norm[0]->GetYaxis()->SetTitleFont(43);
    hEFF_1D_Norm[0]->GetYaxis()->SetTitleOffset(2.);
    hEFF_1D_Norm[0]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hEFF_1D_Norm[0]->GetYaxis()->SetLabelSize(15);
    hEFF_1D_Norm[0]->Draw();
    
    for ( Int_t iFil2 = 1; iFil2 <= fValidFiles; iFil2++ ) {
        hEFF_1D_Norm[iFil2]  ->Draw("SAME");
    }

    cCheckEff->SaveAs("cCheckEff.pdf");
    
    return;
}

