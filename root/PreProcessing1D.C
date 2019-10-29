#include "../inc/SetValues.h"

int main ()
{
    //Retrieving Event/MC data
    TFile *inFile = new TFile(oFileMonCar);
    //TFile *inFile = new TFile(outDF);
    
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
    
    PtreePhi  ->SetBranchAddress   ("evPhi.nPhi",&evPhi.nPhi);
    PtreePhi  ->SetBranchAddress   ("evPhi.pT",&evPhi.pT);
    
    TH1F ** hdM_dpT_Tot_Rec     = new TH1F * [nBin_pT];
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++)
    {
        auto hName = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        auto hFill = Form("evKaonCouple.InvMass>>hdM_dpT_Tot_Rec_%i",iHisto);
        auto hCuts = Form("evKaonCouple.pT >= %f && evKaonCouple.pT < %f && evKaonCouple.bRec == 1",fBound_pT(iHisto),fBound_pT(iHisto+1));
        hdM_dpT_Tot_Rec[iHisto] = new TH1F (hName,hName,nBins,minBound,maxBound);
        PtreeK2->Draw(hFill,hCuts,"goff");
    }
    
    TFile *outFile = new TFile(oFilePreP1D,"recreate");
    for (int iHisto = 0; iHisto < nBin_pT; iHisto++)
    {
        hdM_dpT_Tot_Rec[iHisto] -> Write();
    }
    outFile     -> Close();
    return 0;
}
