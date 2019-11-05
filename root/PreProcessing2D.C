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
    
    TH1F ** hdM_dpT_Tot_Rec     = new TH1F * [nBinIM1D];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        auto hName = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        auto hFill = Form("evKaonCouple.InvMass>>hdM_dpT_Tot_Rec_%i",iHisto);
        auto hCuts = Form("evKaonCouple.pT >= %f && evKaonCouple.pT < %f && evKaonCouple.bRec == 1",fBoundPT2D(iHisto),fBoundPT2D(iHisto+1));
        hdM_dpT_Tot_Rec[iHisto] = new TH1F (hName,hName,nBinIM1D,fMinIM1D,fMaxIM1D);
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
    Int_t nPhiCouple = 0;
    Int_t nPhiCouple2 = 0;
    for (Int_t iEvent = 0; iEvent < PtreeK2->GetEntries(); iEvent++)
    {
        PtreeK2->GetEntry(iEvent);
        Int_t nPhi__ = 0;
        for (int iPhi = 0; iPhi < evKaonCouple.nKaonCouple; iPhi++ )
        {
            //if (evKaonCouple.bRec[iPhi] == false) continue;
            //if (evKaonCouple.bPhi[iPhi] == true) continue;
            nPhi__++;
            for (int jPhi = 0; jPhi < evKaonCouple.nKaonCouple; jPhi++ )
            {
                //if (evKaonCouple.bRec[jPhi] == false) continue;
                //if (evKaonCouple.bPhi[jPhi] == true) continue;
                if ( evKaonCouple.iKaon[iPhi] == evKaonCouple.iKaon[jPhi] ) continue;
                if ( evKaonCouple.iKaon[jPhi] == evKaonCouple.iKaon[iPhi] ) continue;
                if ( evKaonCouple.jKaon[iPhi] == evKaonCouple.jKaon[jPhi] ) continue;
                if ( evKaonCouple.jKaon[jPhi] == evKaonCouple.jKaon[iPhi] ) continue;
                nPhiCouple2++;
                auto ipT = 0;
                auto jpT = 0;
                for (int pT = 0; pT < nBinPT2D; pT++ )
                {
                    if (evKaonCouple.pT[iPhi] >= fBoundPT2D(pT) && evKaonCouple.pT[iPhi] < fBoundPT2D(pT+1)) ipT = pT;
                    if (evKaonCouple.pT[jPhi] >= fBoundPT2D(pT) && evKaonCouple.pT[jPhi] < fBoundPT2D(pT+1)) jpT = pT;
                }
                hdM_dpT_Tot_Rec2D[ipT][jpT]->Fill(evKaonCouple.InvMass[iPhi],evKaonCouple.InvMass[jPhi],0.5);
            }
        }
        if (nPhi__ < 2) continue;
        nPhiCouple += TMath::Binomial(nPhi__,2);
    }
    cout << nPhiCouple << endl;
    cout << nPhiCouple2 << endl;
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
