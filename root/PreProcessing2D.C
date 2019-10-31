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
    
    TH1F ** hdM_dpT_Tot_Rec     = new TH1F * [nBin_pT];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        auto hName = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        auto hFill = Form("evKaonCouple.InvMass>>hdM_dpT_Tot_Rec_%i",iHisto);
        auto hCuts = Form("evKaonCouple.pT >= %f && evKaonCouple.pT < %f && evKaonCouple.bRec == 1",fBound2D_pT(iHisto),fBound2D_pT(iHisto+1));
        hdM_dpT_Tot_Rec[iHisto] = new TH1F (hName,hName,nBins,minBound,maxBound);
        PtreeK2->Draw(hFill,hCuts,"goff");
    }
    
    TH2F ***hdM_dpT_Tot_Rec2D   = new TH2F **[nBinPT2D];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Tot_Rec2D[iHisto] = new TH2F * [nBinPT2D];
        for (int iHist2 = 0; iHist2 < nBinPT2D; iHist2++)
        {
            auto hName = Form("hdM_dpT_Tot_Rec2D_%i_%i",iHisto,iHist2);
            hdM_dpT_Tot_Rec2D[iHisto][iHist2] = new TH2F (hName,hName,nBinIM2D,fMinIM2D,fMaxIM2D,nBinIM2D,fMinIM2D,fMaxIM2D);
        }
    }
    
    for (Int_t iEvent = 0; iEvent < PtreeK2->GetEntries(); iEvent++)
    {
        PtreeK2->GetEntry(iEvent);
        for (int iKaon = 0; iKaon < evKaonCouple.nKaonCouple; iKaon++ )
        {
            if (evKaonCouple.bRec[iKaon] == 0) continue;
            for (int jKaon = 0; jKaon < evKaonCouple.nKaonCouple; jKaon++ )
            {
                if (evKaonCouple.bRec[jKaon] == 0)  continue;
                if (iKaon == jKaon)  continue;
                auto ipT = 0;
                auto jpT = 0;
                for (int pT = 0; pT < nBinPT2D; pT++ )
                {
                    if (evKaonCouple.pT[iKaon] >= fBound2D_pT(pT) && evKaonCouple.pT[iKaon] < fBound2D_pT(pT+1)) ipT = pT;
                    if (evKaonCouple.pT[jKaon] >= fBound2D_pT(pT) && evKaonCouple.pT[jKaon] < fBound2D_pT(pT+1)) jpT = pT;
                }
                hdM_dpT_Tot_Rec2D[ipT][jpT]->Fill(evKaonCouple.InvMass[iKaon],evKaonCouple.InvMass[jKaon]);
            }
        }
    }
    
    TFile *outFile = new TFile(oFilePreP2D,"recreate");
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        for (int iHist2 = 0; iHist2 < nBinPT2D; iHist2++)
        {
            hdM_dpT_Tot_Rec2D[iHisto][iHist2]->Write();
        }
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Tot_Rec[iHisto] -> Write();
    }
    outFile     ->Close();
    
    return 0;
}
