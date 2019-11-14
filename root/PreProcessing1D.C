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
    // Filling TTree of candidate Phis
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.nKaonCouple",&evKaonCouple.nKaonCouple);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.iKaon",&evKaonCouple.iKaon);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.jKaon",&evKaonCouple.jKaon);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.bPhi",&evKaonCouple.bPhi);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.bRec",&evKaonCouple.bRec);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.InvMass",&evKaonCouple.InvMass);
    PtreeK2  ->SetBranchAddress    ("evKaonCouple.pT",&evKaonCouple.pT);
    
    // Filling TTree of true Phis
    PtreePhi  ->SetBranchAddress   ("evPhi.nPhi",&evPhi.nPhi);
    PtreePhi  ->SetBranchAddress   ("evPhi.pT",&evPhi.pT);
    
    // Setting Bins for Invarian Mass and pT
    vSetBinsPT1D();
    vSetBinsIM1D();
    
    // Generating histograms
    auto hName  = "Name";
    auto hTitle = "Title";
    TH1F ** hdM_dpT_Tot_Rec     = new TH1F * [nBinPT1D];
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        // Setting up 1D Histogram
        hName = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        hTitle= Form("m_{K_{+}K_{-}} in p_{T} range %f to %f",fArrIM1D[iHisto],fArrIM1D[iHisto+1]);
        hdM_dpT_Tot_Rec[iHisto] = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
    }
    
    // Filling 1D histograms
    for (Int_t iEvent = 0; iEvent < PtreeK2->GetEntries(); iEvent++)
    {
        PtreeK2->GetEntry(iEvent);
        for (int iPhi = 0; iPhi < evKaonCouple.nKaonCouple; iPhi++ )
        {
            // Only Recordable Phi Candidates
            if (evKaonCouple.bRec[iPhi] == false) continue;
            
            // Information on pT based on defined bins
            Int_t ipT = 0;
            if (evKaonCouple.pT[iPhi] >  fMaxPT1D ) break;
            if (evKaonCouple.pT[iPhi] <  fMinPT1D ) break;
            for (Int_t ipT_ = 0; ipT_ <= nBinPT1D; ipT_++ )
            {
                ipT = ipT_;
                if (evKaonCouple.pT[iPhi] <= fArrPT1D[ipT_+1]) break;
            }
            hdM_dpT_Tot_Rec[ipT]->Fill(evKaonCouple.InvMass[iPhi]);
        }
    }
    TFile *outFile = new TFile(oFilePreP1D,"recreate");
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hdM_dpT_Tot_Rec[iHisto] -> Write();
    }
    inFile->Close();
    outFile->Close();
    
    return 0;
}
