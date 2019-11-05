#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

//int main ()
int Analysis2D ()
{
    //Retrieving PreProcessing
    TFile *iFile_PP = new TFile(oFilePreP2D);
    TFile *iFile_Ef = new TFile(oFileEffici);
    
    //Target variables
    RooRealVar  xInvMass2D ("xInvMass2D","xInvMass2D",fMinIM2D,fMaxIM2D);
    RooRealVar  yInvMass2D ("yInvMass2D","yInvMass2D",fMinIM2D,fMaxIM2D);
    
    //Recovering histograms in roofit
    RooDataHist **  hdM_dpT_Tot_Rec    = new RooDataHist *  [nBinPT2D];
    RooDataHist *** hdM_dpT_Tot_Rec2D  = new RooDataHist ** [nBinPT2D];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        auto hName = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        hdM_dpT_Tot_Rec[iHisto] = new RooDataHist (hName,hName,xInvMass2D,Import(*(TH1F*)(iFile_PP->Get(hName))));
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Tot_Rec2D[iHisto] = new RooDataHist * [nBinPT2D];
        for (int iHist2 = 0; iHist2 < nBinPT2D; iHist2++)
        {
            auto hName  = Form("hdM_dpT_Tot_Rec2D_%i_%i",iHisto,iHist2);
            hdM_dpT_Tot_Rec2D[iHisto][iHist2] = new RooDataHist(hName,hName,RooArgList(xInvMass2D,yInvMass2D),Import(*(TH2F*)(iFile_PP->Get(hName))));
        }
    }

    //Final Histograms
    TH1F * hPhiEff_2D      = (TH1F*)(iFile_Ef->Get("hPhiEff_2D"));
    TH1F * hPhiTru_2D      = (TH1F*)(iFile_Ef->Get("hPhiTru_2D"));
    TH2F * hPhiRaw_2D      = new TH2F ("hPhiRaw_2D","hPhiRaw_2D",        nBinPT2D,fMinPT2D,fMaxPT2D,nBinPT2D,fMinPT2D,fMaxPT2D);
    hPhiRaw_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiRaw_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    TH2F * hPhiRes_2D      = new TH2F ("hPhiRes_2D","hPhiRes_2D",        nBinPT2D,fMinPT2D,fMaxPT2D,nBinPT2D,fMinPT2D,fMaxPT2D);
    hPhiRes_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiRes_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    TH2F * hPhiChk_2D      = new TH2F ("hPhiChk_2D","hPhiChk_2D",        nBinPT2D,fMinPT2D,fMaxPT2D,nBinPT2D,fMinPT2D,fMaxPT2D);
    hPhiChk_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiChk_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    //Preprocessing
    RooFitResult *** Results = new RooFitResult **  [nBinPT2D];
    RooFitResult **  utility = new RooFitResult *   [nBinPT2D];
    
    // Fit  to  data, N_raw
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        Results[iFit]   = new RooFitResult * [nBinPT2D];
        utility[iFit]   = FitModel1D(hdM_dpT_Tot_Rec[iFit],xInvMass2D);
    }
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        for (int jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            Results[iFit][jFit] = FitModel2D(utility[iFit],utility[jFit],hdM_dpT_Tot_Rec2D[iFit][jFit],xInvMass2D,yInvMass2D);
            
            /* Building N_Raw histogram */
            auto N_Raw      = static_cast<RooRealVar*>(Results[iFit][jFit]->floatParsFinal().at(inSS));
            hPhiRaw_2D->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            hPhiRaw_2D->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
    
    //N_res
    hPhiRes_2D->Divide(hPhiRaw_2D,hPhiEff_2D,1.,0.489*0.489,"");
    
    //Check on results
    hPhiChk_2D->Divide(hPhiRes_2D,hPhiTru_2D,1.,1.,"");
    
    TFile *oFile_A1 = new TFile(oFileAnal2D,"recreate");
    hPhiEff_2D->Write();
    hPhiRaw_2D->Write();
    hPhiRes_2D->Write();
    hPhiChk_2D->Write();
    hPhiTru_2D->Write();
    iFile_PP  -> Close();
    iFile_Ef  -> Close();
    oFile_A1  -> Close();
    TFile *oFile_Ht = new TFile(oFileHist2D,"recreate");
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            HistoModel(utility[iHisto],utility[jHisto],Results[iHisto][jHisto],xInvMass2D,yInvMass2D,Form("Model2D_%f_%f_%f_%f",fBoundPT2D(iHisto),fBoundPT2D(iHisto+1),fBoundPT2D(jHisto),fBoundPT2D(jHisto+1)))->Write();
        }
    }
    oFile_Ht  -> Close();
    return 0;
}
