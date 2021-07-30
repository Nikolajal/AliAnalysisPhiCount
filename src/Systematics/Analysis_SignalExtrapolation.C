// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "RooMsgService.h"

void Analysis_SignalExtrapolation ( bool fSilent = false)
{
    TFile *insFileH1D  =   new TFile (fYldSigCorr);
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    // Creating the histograms-------------------------------------------------------------------------------
    //
    hName               =   Form("h1D_Syst_Bin_StSy");
    hTitle              =   Form("Percentage variation of raw yield for bin of PT [%.2f#;%.2f]",fArrPT1D[0],fArrPT1D[nBinPT1D]);
    auto fTest          =   (TGraphAsymmErrors*)insFileH1D->Get("gRES_1D_Stat");
    std::vector<TGraphAsymmErrors*> fPT2D_Graphs;
    for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ )  {
        fPT2D_Graphs.push_back((TGraphAsymmErrors*)insFileH1D->Get(Form("gRES_2D_Stat_%i",iPT2D)));
    }
    
    //------------//
    //  ANALYSIS  //
    //------------//
    
    fSetAllFunctions();
    auto iTer = 0;
    
    for ( auto fTestFit : fSystFitFunctions )   {
        fSetFunction(fTestFit);
        cout << endl;
        cout << endl;
        cout << fTestFit->GetName() << endl;
        cout << endl;
        cout << endl;
        fTest->Fit(fTestFit,"IMRE0S","",0.4,1.8);
        TCanvas *c1 = new TCanvas();
        
        fTest->Draw();
        fTestFit->Draw("same");
        TText * fText = new TText();
        fText -> DrawTextNDC(0.5,0.5,fTestFit->GetName());
        
        c1->SaveAs(Form("c1_%i.pdf",iTer));
        delete c1;
        
        auto jTer = 0;
        for ( auto fGraph2D : fPT2D_Graphs )    {
            
            fGraph2D->Fit(fTestFit,"IMRE0S","",0.4,1.8);
            TCanvas *c2 = new TCanvas();
            
            fGraph2D->Draw();
            fTestFit->Draw("same");
            TText * fText = new TText();
            fText -> DrawTextNDC(0.5,0.5,fTestFit->GetName());
            
            c2->SaveAs(Form("c1_%i_%i.pdf",iTer,jTer));
            delete c2;
            jTer++;
        }
        
        iTer++;
    }
    
    insFileH1D->Close();
}
