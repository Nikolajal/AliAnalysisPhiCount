#include "../inc/SetValues.h"

int main ()
{
    //Retrieving Event/MC data
    TFile *inFile = new TFile(oFileMonCar);
        
    //Retrieving Event/MC TTree
    TTree *PtreeK2  = (TTree*)inFile->Get(PTreeNameK2);
    TTree *PtreePhi = (TTree*)inFile->Get(PTreeNamePhi);
        
    // Define some simple data structures
    EVKAONCOUPLE evKaonCouple;
    EVPHI evPhi;
        
    //Setting Branch Addresses
    // Filling TTree
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.nKaonCouple",&evKaonCouple.nKaonCouple);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.iKaon",&evKaonCouple.iKaon);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.jKaon",&evKaonCouple.jKaon);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.bPhi",&evKaonCouple.bPhi);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.bRec",&evKaonCouple.bRec);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.InvMass",&evKaonCouple.InvMass);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.pT",&evKaonCouple.pT);
    
    PtreePhi  ->SetBranchAddress    ("evPhi.nPhi",&evPhi.nPhi);
    PtreePhi  ->SetBranchAddress    ("evPhi.pT",&evPhi.pT);
    
    /*------*/
    /*  1D  */
    /*------*/
    TH1F * hPhiGen_1D   = new TH1F ("hPhiGen_1D","hPhiGen_1D",nBinPhi1D,fMinPhi1D,fMaxPhi1D);
    TH1F * hPhiRec_1D   = new TH1F ("hPhiRec_1D","hPhiRec_1D",nBinPhi1D,fMinPhi1D,fMaxPhi1D);
    TH1F * hPhiTru_1D   = new TH1F ("hPhiTru_1D","hPhiTru_1D",nBinPhi1D,fMinPhi1D,fMaxPhi1D);
    TH1F * hPhiEff_1D   = new TH1F ("hPhiEff_1D","hPhiEff_1D",nBinPhi1D,fMinPhi1D,fMaxPhi1D);
    
    //N_Gen and N_rec built
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++)
    {
        auto hCuts = Form("evKaonCouple.pT >= %f && evKaonCouple.pT < %f && evKaonCouple.bPhi",fBound_pT(iHisto),fBound_pT(iHisto+1));
        TH1F hdM_dpT ("hdM_dpT","hdM_dpT",nBins,minBound,maxBound);
        PtreeK2->Draw("evKaonCouple.InvMass>>hdM_dpT",hCuts,"goff");
        auto N_Val      = (hdM_dpT.GetEntries());
        auto N_ValE     = sqrt(N_Val);
        hPhiGen_1D->SetBinContent      (iHisto+1,N_Val);
        hPhiGen_1D->SetBinError        (iHisto+1,N_ValE);
        
        hCuts = Form("evKaonCouple.pT >= %f && evKaonCouple.pT < %f && evKaonCouple.bRec && evKaonCouple.bPhi",fBound_pT(iHisto),fBound_pT(iHisto+1));
        TH1F hdM_dpT1 ("hdM_dpT1","hdM_dpT1",nBins,minBound,maxBound);
        PtreeK2->Draw("evKaonCouple.InvMass>>hdM_dpT",hCuts,"goff");
        N_Val      = (hdM_dpT.GetEntries());
        N_ValE     = sqrt(N_Val);
        hPhiRec_1D->SetBinContent      (iHisto+1,N_Val);
        hPhiRec_1D->SetBinError        (iHisto+1,N_ValE);
    }
    
    //N_Tru
    PtreePhi->Draw("evPhi.pT >> hPhiTru_1D","","goff");
    
    //Efficiency
    hPhiEff_1D->Divide(hPhiRec_1D,hPhiGen_1D,1,1,"b");
    
    /*------*/
    /*  2D  */
    /*------*/
    
    TH2F * hPhiGen_2D   = new TH2F ("hPhiGen_2D","hPhiGen_2D",nBinPhi2D,fMinPhi2D,fMaxPhi2D,nBinPhi2D,fMinPhi2D,fMaxPhi2D);
    TH2F * hPhiRec_2D   = new TH2F ("hPhiRec_2D","hPhiRec_2D",nBinPhi2D,fMinPhi2D,fMaxPhi2D,nBinPhi2D,fMinPhi2D,fMaxPhi2D);
    TH2F * hPhiTru_2D   = new TH2F ("hPhiTru_2D","hPhiTru_2D",nBinPhi2D,fMinPhi2D,fMaxPhi2D,nBinPhi2D,fMinPhi2D,fMaxPhi2D);
    TH2F * hPhiEff_2D   = new TH2F ("hPhiEff_2D","hPhiEff_2D",nBinPhi2D,fMinPhi2D,fMaxPhi2D,nBinPhi2D,fMinPhi2D,fMaxPhi2D);
    
    for (Int_t iEvent = 0; iEvent < PtreeK2->GetEntries(); iEvent++)
    {
        PtreeK2->GetEntry(iEvent);
        if (evKaonCouple.nKaonCouple <= 4) continue;
        for (Int_t iKaon = 0; iKaon < evKaonCouple.nKaonCouple; iKaon++)
        {
            for (Int_t jKaon = iKaon+1; jKaon < evKaonCouple.nKaonCouple; jKaon++)
            {
                if ( !(evKaonCouple.bPhi[iKaon] && evKaonCouple.bPhi[jKaon]) ) continue;
                hPhiGen_2D->Fill(evKaonCouple.pT[iKaon],evKaonCouple.pT[jKaon]);
                if ( !(evKaonCouple.bRec[iKaon] && evKaonCouple.bRec[jKaon]) ) continue;
                hPhiRec_2D->Fill(evKaonCouple.pT[iKaon],evKaonCouple.pT[jKaon]);
            }
        }
    }
    
    //N_Tru
    for (Int_t iEvent = 0; iEvent < PtreePhi->GetEntries(); iEvent++)
    {
        PtreePhi->GetEntry(iEvent);
        for (Int_t iPhi = 0; iPhi < evPhi.nPhi; iPhi++)
        {
            for (Int_t jPhi = 1+iPhi; jPhi < evPhi.nPhi; jPhi++)
            {
                if (iPhi == jPhi) continue;
                hPhiTru_2D->Fill(evPhi.pT[iPhi],evPhi.pT[jPhi]);
            }
        }
    }
    
    //Efficiency
    hPhiEff_2D->Divide(hPhiRec_2D,hPhiGen_2D,1,1,"b");
    
    // Writing results out to file
    TFile * outFile = new TFile(oFileEffici,"recreate");
    hPhiTru_1D      ->Write();
    hPhiTru_2D      ->Write();
    hPhiGen_1D      ->Write();
    hPhiGen_2D      ->Write();
    hPhiRec_1D      ->Write();
    hPhiRec_2D      ->Write();
    hPhiEff_1D      ->Write();
    hPhiEff_2D      ->Write();
    outFile         ->Close();
    return 0;
}
