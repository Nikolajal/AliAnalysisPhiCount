#include "../inc/SetValues.h"

//int main ()
int Analysis ()
{
    //Retrieving PreProcessing
    TFile *inFil1 = new TFile(oFilePreP1D);
    TFile *inFil2 = new TFile(oFilePreP2D);
    TFile *inFilE = new TFile(oFileEffici);
    
    //Target variables
    RooRealVar  xInvMass1D ("xInvMass1D","xInvMass1D",minBound,maxBound);
    RooRealVar  xInvMass2D ("xInvMass2D","xInvMass2D",minBound,maxBound);
    RooRealVar  yInvMass2D ("yInvMass2D","yInvMass2D",minBound,maxBound);
    
    //Recovering histograms in roofit
    RooDataHist ** hdM_dpT_Tot_Rec     = new RooDataHist * [nBin_pT];
    RooDataHist *** hdM_dpT_Tot_Rec2D  = new RooDataHist ** [nBinPT2D];
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++)
    {
        auto hName  = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        hdM_dpT_Tot_Rec[iHisto] = new RooDataHist (hName,hName,xInvMass1D,Import(*(TH1F*)(inFil1->Get(hName))));
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
    
    // Parameters
    RooRealVar    ** ch0            = new RooRealVar * [nBin_pT];
    RooRealVar    ** ch1            = new RooRealVar * [nBin_pT];
    RooRealVar    ** ch2            = new RooRealVar * [nBin_pT];
    RooRealVar    ** ch3            = new RooRealVar * [nBin_pT];
    RooRealVar    ** ch4            = new RooRealVar * [nBin_pT];
    RooRealVar    ** m1Exp1D        = new RooRealVar * [nBin_pT];
    RooRealVar    ** m2Exp1D        = new RooRealVar * [nBin_pT];
    RooRealVar    ** m3Exp1D        = new RooRealVar * [nBin_pT];
    RooRealVar    ** m1Off1D        = new RooRealVar * [nBin_pT];
    RooRealVar    ** nSig1D         = new RooRealVar * [nBin_pT];
    RooRealVar    ** nBkg1D         = new RooRealVar * [nBin_pT];
    RooRealVar    ** phiMass        = new RooRealVar * [nBin_pT];
    RooRealVar    ** phiWidth       = new RooRealVar * [nBin_pT];
    Int_t         * nEntries1D      = new Int_t [nBin_pT];
    for (int iRVar = 0; iRVar < nBin_pT; iRVar++)
    {
        nEntries1D[iRVar] = hdM_dpT_Tot_Rec[iRVar]->sumEntries();
        
        auto hName  = Form("ch0_%i",iRVar);
        ch0[iRVar] = new RooRealVar (hName,hName,0.,-1.,1.);
        hName  = Form("ch1_%i",iRVar);
        ch1[iRVar] = new RooRealVar (hName,hName,0.,-1.,1.);
        hName  = Form("ch2_%i",iRVar);
        ch2[iRVar] = new RooRealVar (hName,hName,0.,-1.,1.);
        hName  = Form("ch3_%i",iRVar);
        ch3[iRVar] = new RooRealVar (hName,hName,0.,-1.,1.);
        hName  = Form("ch4_%i",iRVar);
        ch4[iRVar] = new RooRealVar (hName,hName,0.,-1.,1.);
        
        hName       = Form("nBkg1D_%i",iRVar);
        nBkg1D[iRVar] = new RooRealVar (hName,hName,0.8*nEntries1D[iRVar],0.,nEntries1D[iRVar]);
        hName       = Form("m1Exp1D_%i",iRVar);
        m1Exp1D[iRVar] = new RooRealVar (hName,hName,1.,-100.,100.);
        hName       = Form("m2Exp1D_%i",iRVar);
        m2Exp1D[iRVar] = new RooRealVar (hName,hName,0.,-100.,100.);
        hName       = Form("m3Exp1D_%i",iRVar);
        m3Exp1D[iRVar] = new RooRealVar (hName,hName,0.,-100.,100.);
        hName       = Form("m1Off1D_%i",iRVar);
        m1Off1D[iRVar] = new RooRealVar (hName,hName,1.027,-100.,100.);
        
        hName       = Form("nSig1D_%i",iRVar);
        nSig1D[iRVar] = new RooRealVar (hName,hName,0.2*nEntries1D[iRVar],0.,nEntries1D[iRVar]);
        
        hName       = Form("PhiMass_%i",iRVar);
        phiMass[iRVar] = new RooRealVar (hName,hName,1.020,1.010,1.030);
        hName       = Form("PhiWidth_%i",iRVar);
        phiWidth[iRVar] = new RooRealVar (hName,hName,0.005,0.0005,0.01);
    }
     
    // Signal, Background and Model
    RooGenericPdf   ** Bkg1D        = new RooGenericPdf * [nBin_pT];
    RooChebychev    ** Bkg1D_alt    = new RooChebychev * [nBin_pT];
    RooAddPdf       ** Model1D      = new RooAddPdf * [nBin_pT];
    RooBreitWigner  ** Sig1D        = new RooBreitWigner * [nBin_pT];
    for (int iFun = 0; iFun < nBin_pT; iFun++)
    {
        auto hName  = Form("Bkg1D_%i",iFun);
        //if (iFun == 2) Bkg1D[iFun] = new RooGenericPdf (hName,hName,"max(0,@1+@2*pow(@0,1)+@3*pow(@0,2)+@4*pow(@0,3))",RooArgSet(xInvMass1D,*m1Exp1D[iFun],*m1Off1D[iFun],*m2Exp1D[iFun],*m3Exp1D[iFun]));
        if (iFun <= 2)
        {
            Bkg1D[iFun] = new RooGenericPdf (hName,hName,"max(0,@2*(@0-@1)*abs(@0-@1)+@0*@3)",RooArgSet(xInvMass1D,*m1Exp1D[iFun],*m1Off1D[iFun],*m2Exp1D[iFun],*m3Exp1D[iFun]));
        }
        Bkg1D_alt[iFun] = new RooChebychev (hName,hName,xInvMass1D,RooArgSet(*ch0[iFun],*ch1[iFun],*ch2[iFun],*ch3[iFun],*ch4[iFun]));
        
        hName  = Form("Sig1D_%i",iFun);
        Sig1D[iFun] = new RooBreitWigner (hName,hName,xInvMass1D,*phiMass[iFun],*phiWidth[iFun]);
        
        hName       = Form("Model1D_%i",iFun);
        if (iFun <= 2) Model1D[iFun] = new
        RooAddPdf (hName,hName,RooArgList(*Sig1D[iFun],*Bkg1D[iFun]),RooArgList(*nSig1D[iFun],*nBkg1D[iFun]));
        if (iFun >= 3) Model1D[iFun] = new
        RooAddPdf (hName,hName,RooArgList(*Sig1D[iFun],*Bkg1D_alt[iFun]),RooArgList(*nSig1D[iFun],*nBkg1D[iFun]));
    }
    
    // Fit  to  data
    for (int iFit = 0; iFit < nBin_pT; iFit++)
    {
        if (iFit == 0)  xInvMass1D.setRange(1.033,maxBound);
        if (iFit == 1)  xInvMass1D.setRange(1.018,maxBound);
        if (iFit == 2)  xInvMass1D.setRange(0.995,maxBound);
        if (iFit == 3)  xInvMass1D.setRange(minBound,maxBound);
        (*Model1D[iFit]).fitTo((*hdM_dpT_Tot_Rec[iFit]));
    }
    
    //N_raw
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++ )
    {
        auto N_Raw      = ((*nSig1D[iHisto]).getVal());
        auto N_RawE     = ((*nSig1D[iHisto]).getError());

        hPhiRaw_1D->SetBinContent      (iHisto+1,N_Raw);
        hPhiRaw_1D->SetBinError        (iHisto+1,N_RawE);
    }
    
    //N_res
    hPhiRes_1D->Divide(hPhiRaw_1D,hPhiEff_1D,1.,0.489,"");
    
    //Check on results
    hPhiChk_1D->Divide(hPhiRes_1D,hPhiTru_1D,1.,1.,"");
    
    /*------*/
    /*  2D  */
    /*------*/

    // Parameters
    RooRealVar   *** ch02D          = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** ch12D          = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** ch22D          = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** ch32D          = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** ch42D          = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** nSig2D         = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** nBkg2D         = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** phiMass2D      = new RooRealVar ** [nBinPT2D];
    RooRealVar   *** phiWidth2D     = new RooRealVar ** [nBinPT2D];
    Int_t         ** nEntries2D     = new Int_t * [nBinPT2D];
    for (int iRVar = 0; iRVar < nBinPT2D; iRVar++)
    {
        nEntries2D[iRVar]   = new Int_t [nBinPT2D];
        ch02D[iRVar]        = new RooRealVar * [nBinPT2D];
        ch12D[iRVar]        = new RooRealVar * [nBinPT2D];
        ch22D[iRVar]        = new RooRealVar * [nBinPT2D];
        ch32D[iRVar]        = new RooRealVar * [nBinPT2D];
        ch42D[iRVar]        = new RooRealVar * [nBinPT2D];
        nSig2D[iRVar]       = new RooRealVar * [nBinPT2D];
        nBkg2D[iRVar]       = new RooRealVar * [nBinPT2D];
        phiMass2D[iRVar]    = new RooRealVar * [nBinPT2D];
        phiWidth2D[iRVar]   = new RooRealVar * [nBinPT2D];
        for (int iRVa2 = 0; iRVa2 < nBinPT2D; iRVa2++)
        {
            nEntries2D[iRVar][iRVa2] = hdM_dpT_Tot_Rec2D[iRVar][iRVa2]->sumEntries();
            auto hName  = Form("ch02D_%i_%i",iRVar,iRVa2);
            ch02D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.,-1.,1.);
            hName  = Form("ch12D_%i_%i",iRVar,iRVa2);
            ch12D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.,-1.,1.);
            hName  = Form("ch22D_%i_%i",iRVar,iRVa2);
            ch22D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.,-1.,1.);
            hName  = Form("ch32D_%i_%i",iRVar,iRVa2);
            ch32D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.,-1.,1.);
            hName  = Form("ch42D_%i_%i",iRVar,iRVa2);
            ch42D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.,-1.,1.);
            
            hName       = Form("nSig2D_%i_%i",iRVar,iRVa2);
            nSig2D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.01*nEntries2D[iRVar][iRVa2],0.,nEntries2D[iRVar][iRVa2]);
            
            hName       = Form("nBkg2D_%i_%i",iRVar,iRVa2);
            nBkg2D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.99*nEntries2D[iRVar][iRVa2],0.,nEntries2D[iRVar][iRVa2]);
            
            hName       = Form("PhiMass2D_%i_%i",iRVar,iRVa2);
            phiMass2D[iRVar][iRVa2] = new RooRealVar (hName,hName,1.020,1.010,1.030);
            
            hName       = Form("PhiWidth2D_%i_%i",iRVar,iRVa2);
            phiWidth2D[iRVar][iRVa2] = new RooRealVar (hName,hName,0.005,0.003,0.01);
        }
    }
     
    // Signal, Background and Model
    RooChebychev   *** xBkg2D_alt   = new RooChebychev      ** [nBinPT2D];
    RooChebychev   *** yBkg2D_alt   = new RooChebychev      ** [nBinPT2D];
    RooBreitWigner *** xSig2D       = new RooBreitWigner    ** [nBinPT2D];
    RooBreitWigner *** ySig2D       = new RooBreitWigner    ** [nBinPT2D];
    RooAddPdf      *** Model2D      = new RooAddPdf         ** [nBinPT2D];
    for (int iFun = 0; iFun < nBinPT2D; iFun++)
    {
        xBkg2D_alt[iFun]   = new RooChebychev      * [nBinPT2D];
        yBkg2D_alt[iFun]   = new RooChebychev      * [nBinPT2D];
        xSig2D[iFun]       = new RooBreitWigner    * [nBinPT2D];
        ySig2D[iFun]       = new RooBreitWigner    * [nBinPT2D];
        Model2D[iFun]      = new RooAddPdf         * [nBinPT2D];
        for (int iFu2 = 0; iFu2 < nBinPT2D; iFu2++)
        {
            auto hName  = Form("xBkg2D_%i_%i",iFun,iFu2);
            xBkg2D_alt[iFun][iFu2] = new RooChebychev (hName,hName,xInvMass2D,RooArgSet(*ch02D[iFun][iFu2],*ch12D[iFun][iFu2],*ch22D[iFun][iFu2],*ch32D[iFun][iFu2],*ch42D[iFun][iFu2]));
            
            hName  = Form("yBkg2D_%i_%i",iFun,iFu2);
            yBkg2D_alt[iFun][iFu2] = new RooChebychev (hName,hName,yInvMass2D,RooArgSet(*ch02D[iFun][iFu2],*ch12D[iFun][iFu2],*ch22D[iFun][iFu2],*ch32D[iFun][iFu2],*ch42D[iFun][iFu2]));
            
            hName  = Form("xSig2D_%i_%i",iFun,iFu2);
            xSig2D[iFun][iFu2] = new RooBreitWigner (hName,hName,xInvMass2D,*phiMass2D[iFun][iFu2],*phiWidth2D[iFun][iFu2]);
            
            hName  = Form("ySig2D_%i_%i",iFun,iFu2);
            ySig2D[iFun][iFu2] = new RooBreitWigner (hName,hName,yInvMass2D,*phiMass2D[iFun][iFu2],*phiWidth2D[iFun][iFu2]);
            
            hName  = Form("Model2D_%i_%i",iFun,iFu2);
            Model2D[iFun][iFu2] = new
            RooAddPdf (hName,hName,RooArgList(*xSig2D[iFun][iFu2],*ySig2D[iFun][iFu2],*xBkg2D_alt[iFun][iFu2],*yBkg2D_alt[iFun][iFu2]),RooArgList(*nSig2D[iFun][iFu2],*nSig2D[iFun][iFu2],*nBkg2D[iFun][iFu2],*nBkg2D[iFun][iFu2]));
        }
    }
    
    // Fit  to  data
    for (int iFit = 0; iFit < nBinPT2D; iFit++)
    {
        for (int iFi2 = 0; iFi2 < nBinPT2D; iFi2++)
        {
            (*Model2D[iFit][iFi2]).fitTo((*hdM_dpT_Tot_Rec2D[iFit][iFi2]));
        }
    }
    
    //N_raw
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        for (int iHist2 = 0; iHist2 < nBinPT2D; iHist2++ )
        {
            //x.setRange(“oneSigma”,-3,3) ;
            //RooAbsReal* fracInt = g.createIntegral(x,NormSet(x),Range(“cut”))
            
            auto N_Raw      = ((nSig2D[iHisto][iHist2])->getVal());
            auto N_RawE     = ((nSig2D[iHisto][iHist2])->getError());

            hPhiRaw_2D->SetBinContent      (iHisto+1,iHist2+1,N_Raw);
            hPhiRaw_2D->SetBinError        (iHisto+1,iHist2+1,N_RawE);
        }
    }
    
    //N_res
    hPhiRes_2D->Divide(hPhiRaw_2D,hPhiEff_2D,1.,0.489,"");
    
    //Check on results
    hPhiChk_2D->Divide(hPhiRes_2D,hPhiTru_2D,1.,1.,"");
    
    TFile *outFile = new TFile(oFileAnalys,"recreate");
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
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++ )
    {
        RooPlot * frame = xInvMass1D.frame();
        frame->SetNameTitle(Form("pT1D %f-%f",fBound_pT(iHisto),fBound_pT(iHisto+1)),Form("bb%i",iHisto));
        hdM_dpT_Tot_Rec[iHisto]->plotOn(frame);
        Model1D[iHisto]->plotOn(frame);
        frame->Write();
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        for (int iHist2 = 0; iHist2 < nBinPT2D; iHist2++ )
        {
            RooPlot * framex = xInvMass2D.frame();
            framex->SetNameTitle(Form("pT2D %f-%f %f-%f",fBound2D_pT(iHisto),fBound2D_pT(iHisto+1),fBound2D_pT(iHist2),fBound2D_pT(iHist2+1)),Form("bb2D%i%i",iHisto,iHist2));
            hdM_dpT_Tot_Rec2D[iHisto][iHist2]->plotOn(framex);
            Model2D[iHisto][iHist2]->plotOn(framex);
            framex->Write();
            RooPlot * framey = yInvMass2D.frame();
            framey->SetNameTitle(Form("pT2D %f-%f %f-%f",fBound2D_pT(iHisto),fBound2D_pT(iHisto+1),fBound2D_pT(iHist2),fBound2D_pT(iHist2+1)),Form("bb2D%i%i",iHisto,iHist2));
            hdM_dpT_Tot_Rec2D[iHisto][iHist2]->plotOn(framey);
            Model2D[iHisto][iHist2]->plotOn(framey);
            framey->Write();
        }
    }
    outFile     -> Close();
    
    return 0;
}
