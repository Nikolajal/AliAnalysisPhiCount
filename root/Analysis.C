#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

//int main ()
int Analysis ()
{
    //Retrieving PreProcessing
    TFile *inFil1 = new TFile(oFilePreP1D);
    TFile *inFil2 = new TFile(oFilePreP2D);
    TFile *inFilE = new TFile(oFileEffici);
    TFile *outFile = new TFile(oFileHistog,"recreate");
    
    //Target variables
    RooRealVar  xInvMass1D ("xInvMass1D","xInvMass1D",minBound,maxBound);
    RooRealVar  xInvMass2D ("xInvMass2D","xInvMass2D",minBound,maxBound);
    RooRealVar  yInvMass2D ("yInvMass2D","yInvMass2D",minBound,maxBound);
    
    //Recovering histograms in roofit
    RooDataHist ** hdM_dpT_Tot_Rec     = new RooDataHist * [nBinPT1D];
    RooDataHist ** hdM_dpT_Tot_RecPP   = new RooDataHist * [nBinPT1D];
    RooDataHist *** hdM_dpT_Tot_Rec2D  = new RooDataHist ** [nBinPT2D];
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++)
    {
        auto hName  = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        hdM_dpT_Tot_Rec[iHisto] = new RooDataHist (hName,hName,xInvMass1D,Import(*(TH1F*)(inFil1->Get(hName))));
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        auto hName  = Form("hdM_dpT_Tot_RecPP_%i",iHisto);
        hdM_dpT_Tot_RecPP[iHisto] = new RooDataHist (hName,hName,xInvMass1D,Import(*(TH1F*)(inFil2->Get(Form("hdM_dpT_Tot_Rec_%i",iHisto)))));
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Tot_Rec2D[iHisto] = new RooDataHist * [nBinPT2D];
        for (int iHist2 = 0; iHist2 < nBinPT2D; iHist2++)
        {
            auto hName  = Form("hdM_dpT_Tot_Rec2D_%i_%i",iHisto,iHist2);
            hdM_dpT_Tot_Rec2D[iHisto][iHist2] = new RooDataHist(hName,hName,RooArgList(xInvMass2D,yInvMass2D),Import(*(TH2F*)(inFil2->Get(hName))));
        }
    }

    //Final Histograms
    TH1F * hPhiEff_1D      = (TH1F*)(inFilE->Get("hPhiEff_1D"));
    TH1F * hPhiEff_2D      = (TH1F*)(inFilE->Get("hPhiEff_2D"));
    TH1F * hPhiTru_1D      = (TH1F*)(inFilE->Get("hPhiTru_1D"));
    TH1F * hPhiTru_2D      = (TH1F*)(inFilE->Get("hPhiTru_2D"));
    TH1F * hPhiRaw_1D      = new TH1F ("hPhiRaw_1D","hPhiRaw_1D",        nBin_pT,nMin_pT,nMax_pT);
    TH2F * hPhiRaw_2D      = new TH2F ("hPhiRaw_2D","hPhiRaw_2D",        nBinPT2D,fMinPT2D,fMaxPT2D,nBinPT2D,fMinPT2D,fMaxPT2D);
    TH1F * hPhiRes_1D      = new TH1F ("hPhiRes_1D","hPhiRes_1D",        nBin_pT,nMin_pT,nMax_pT);
    TH2F * hPhiRes_2D      = new TH2F ("hPhiRes_2D","hPhiRes_2D",        nBinPT2D,fMinPT2D,fMaxPT2D,nBinPT2D,fMinPT2D,fMaxPT2D);
    TH1F * hPhiChk_1D      = new TH1F ("hPhiChk_1D","hPhiChk_1D",        nBin_pT,nMin_pT,nMax_pT);
    TH2F * hPhiChk_2D      = new TH2F ("hPhiChk_2D","hPhiChk_2D",        nBinPT2D,fMinPT2D,fMaxPT2D,nBinPT2D,fMinPT2D,fMaxPT2D);
    
    // Set-up model
    
    // Global Variables
    RooRealVar      KaonThreshold   ("KaonThreshold","KaonThreshold",0.9874,0.975,0.995);
    
    /*------*/
    /*  1D  */
    /*------*/
    
    // Fit Results and PlotOn object
    FIT1D_RESULT_ARRAY Result1D;
    
    // Fit  to  data, N_raw
    for (int iFit = 0; iFit < nBinPT1D; iFit++)
    {
        if (iFit == 0)  xInvMass1D.setRange(1.033,maxBound);
        if (iFit == 1)  xInvMass1D.setRange(1.018,maxBound);
        if (iFit == 2)  xInvMass1D.setRange(0.995,maxBound);
        if (iFit == 3)  xInvMass1D.setRange(minBound,maxBound);
        ModelFit1D(hdM_dpT_Tot_Rec[iFit],xInvMass1D,Result1D.Array[iFit],iFit);
        
        /* Building N_Raw histogram */
        auto N_Raw      = static_cast<RooRealVar*>(Result1D.Array[iFit].FitRes->floatParsFinal().at(inSig));
        hPhiRaw_1D->SetBinContent      (iFit+1,N_Raw->getVal());
        hPhiRaw_1D->SetBinError        (iFit+1,N_Raw->getError());
    }
    
    //N_res
    hPhiRes_1D->Divide(hPhiRaw_1D,hPhiEff_1D,1.,0.489,"");
    
    //Check on results
    hPhiChk_1D->Divide(hPhiRes_1D,hPhiTru_1D,1.,1.,"");
    
    /*------*/
    /*  2D  */
    /*------*/
    
    // Fit Results and PlotOn object
    FIT1D_RESULT_ARRAY * Result2D = new FIT1D_RESULT_ARRAY[nBinPT2D];
    
    //Preprocessing
    FIT1D_RESULT_ARRAY preprocess;
    ModelFit2D_Preprocess (hdM_dpT_Tot_RecPP,xInvMass1D,preprocess);
    
    // Fit  to  data, N_raw
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        for (int jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            ModelFit2D(hdM_dpT_Tot_Rec2D[iFit][jFit],xInvMass2D,yInvMass2D,Result2D[iFit].Array[jFit],preprocess,iFit,jFit);
            
            /* Building N_Raw histogram */
            auto N_Raw      = static_cast<RooRealVar*>(Result2D[iFit].Array[jFit].FitRes->floatParsFinal().at(inSS));
            hPhiRaw_2D->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            hPhiRaw_2D->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
    
    //N_res
    hPhiRes_2D->Divide(hPhiRaw_2D,hPhiEff_2D,1.,0.489,"");
    
    //Check on results
    hPhiChk_2D->Divide(hPhiRes_2D,hPhiTru_2D,1.,1.,"");
    
    hPhiEff_1D->Write();
    hPhiEff_2D->Write();
    hPhiRaw_1D->Write();
    hPhiRaw_2D->Write();
    hPhiRes_1D->Write();
    hPhiRes_2D->Write();
    hPhiChk_1D->Write();
    hPhiChk_2D->Write();
    hPhiTru_1D->Write();
    hPhiTru_2D->Write();
    inFil1    -> Close();
    inFil2    -> Close();
    inFilE    -> Close();
    outFile   -> Close();
    
    return 0;
}
