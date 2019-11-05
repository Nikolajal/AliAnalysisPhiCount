#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

//int main ()
int Analysis1D ()
{
    //Retrieving PreProcessing
    TFile *iFile_PP = new TFile(oFilePreP1D);
    TFile *iFile_Ef = new TFile(oFileEffici);
    
    //Target variables
    RooRealVar  xInvMass1D ("xInvMass1D","xInvMass1D",fMinIM1D,fMaxIM1D);
    
    //Recovering histograms in roofit
    RooDataHist ** hdM_dpT_Tot_Rec     = new RooDataHist * [nBinPT1D];
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        auto hName  = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        hdM_dpT_Tot_Rec[iHisto] = new RooDataHist (hName,hName,xInvMass1D,Import(*(TH1F*)(iFile_PP->Get(hName))));
    }

    //Final Histograms
    TH1F * hPhiEff_1D      = (TH1F*)(iFile_Ef->Get("hPhiEff_1D"));
    TH1F * hPhiTru_1D      = (TH1F*)(iFile_Ef->Get("hPhiTru_1D"));
    TH1F * hPhiRaw_1D      = new TH1F ("hPhiRaw_1D","Number of #varphi found in Fit per |y|<5 in K_{+}K_{-} Decay mode, 1D analysis",nBinPT1D,fMinPT1D,fMaxPT1D);
    hPhiRaw_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");
    TH1F * hPhiRes_1D      = new TH1F ("hPhiRes_1D","Number of #varphi reconstructed from Fit per |y|<5, 1D analysis",nBinPT1D,fMinPT1D,fMaxPT1D);
    hPhiRes_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");
    TH1F * hPhiChk_1D      = new TH1F ("hPhiChk_1D","Ratio of #varphi found in Fit and #varphi generated in MC per |y|<5, 1D analysis",nBinPT1D,fMinPT1D,fMaxPT1D);
    hPhiChk_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");

    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    // Fit Results and PlotOn object
    RooFitResult *  Results [nBinPT1D];
    
    // Fit  to  data, N_raw
    for (int iFit = 0; iFit < nBinPT1D; iFit++)
    {
        Results[iFit] = FitModel1D(hdM_dpT_Tot_Rec[iFit],xInvMass1D);
        
        /* Building N_Raw histogram */
        auto N_Raw      = static_cast<RooRealVar*>(Results[iFit]->floatParsFinal().at(inSig));
        hPhiRaw_1D->SetBinContent      (iFit+1,N_Raw->getVal());
        hPhiRaw_1D->SetBinError        (iFit+1,N_Raw->getError());
    }
    
    //N_res
    hPhiRes_1D->Divide(hPhiRaw_1D,hPhiEff_1D,1.,0.489,"");
    
    //Check on results
    hPhiChk_1D->Divide(hPhiRes_1D,hPhiTru_1D,1.,1.,"");
    
    TFile *oFile_A1 = new TFile(oFileAnal1D,"recreate");
    hPhiEff_1D->Write();
    hPhiRaw_1D->Write();
    hPhiRes_1D->Write();
    hPhiChk_1D->Write();
    hPhiTru_1D->Write();
    iFile_PP  -> Close();
    iFile_Ef  -> Close();
    oFile_A1  -> Close();
    TFile *oFile_Ht = new TFile(oFileHist1D,"recreate");
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        HistoModel(Results[iHisto],xInvMass1D,Form("Model1D_%f_%f",fBoundPT1D(iHisto),fBoundPT1D(iHisto+1)))->Write();
    }
    oFile_Ht  -> Close();
    return 0;
}
