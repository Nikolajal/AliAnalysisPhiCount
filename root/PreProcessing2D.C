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
    
    PtreePhi ->SetBranchAddress    ("evPhi.nPhi",&evPhi.nPhi);
    PtreePhi ->SetBranchAddress    ("evPhi.pT",&evPhi.pT);
    
    vSetBinsPT1D();
    vSetBinsIM1D();
    vSetBinsPT2D();
    vSetBinsIM2D();
    
    auto hName  = "Name";
    auto hTitle = "Title";
    TH1F ** hdM_dpT_Tot_Rec     = new TH1F * [nBinIM1D];
    TH2F ***hdM_dpT_Tot_Rec2D   = new TH2F **[nBinPT2D];
    TH2F ***hdM_dpT_Tot_Bkg2D   = new TH2F **[nBinPT2D];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Tot_Rec2D[iHisto] = new TH2F * [nBinPT2D];
        hdM_dpT_Tot_Bkg2D[iHisto] = new TH2F * [nBinPT2D];
        
        // Setting up 1D Histogram
        hName = Form("hdM_dpT_Tot_Rec_%i",iHisto);
        hTitle= Form("m_{K_{+}K_{-}} in p_{T} range %f to %f",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        hdM_dpT_Tot_Rec[iHisto] = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            // Setting up 2D Histogram
            hName = Form("hdM_dpT_Tot_Rec2D_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K_{+}K_{-}} in p_{T} range %f to %f and %f to %f",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hdM_dpT_Tot_Rec2D[iHisto][jHisto] = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
        }
    }
    
    Int_t ipT, jpT;
    Int_t nBkgComb = 0;
    for (Int_t iEvent = 0; iEvent < PtreeK2->GetEntries(); iEvent++)
    {
        PtreeK2->GetEntry(iEvent);
        for (int iPhi = 0; iPhi < evKaonCouple.nKaonCouple; iPhi++ )
        {
            // Only Recordable Phi Candidates
            if (evKaonCouple.bRec[iPhi] == false) continue;
            
            // Only True Phi Candidates
            //if (evKaonCouple.bPhi[iPhi] == false) continue;
            
            // Information on pT based on defined bins
            if (evKaonCouple.pT[iPhi] >  fMaxPT2D ) break;
            if (evKaonCouple.pT[iPhi] <  fMinPT2D ) break;
            for (Int_t ipT_ = 0; ipT_ <= nBinPT2D; ipT_++ )
            {
                ipT = ipT_;
                if (evKaonCouple.pT[iPhi] <= fArrPT2D[ipT_+1]) break;
            }
            hdM_dpT_Tot_Rec[ipT]->Fill(evKaonCouple.InvMass[iPhi]);
            
            for (int jPhi = 0; jPhi < evKaonCouple.nKaonCouple; jPhi++ )
            {
                // Only Recordable Phi Candidates
                if (evKaonCouple.bRec[jPhi] == false) continue;
                
                // Only True Phi Candidates
                //if (evKaonCouple.bPhi[jPhi] == false) continue;
                
                // Information on pT based on defined bins
                if (evKaonCouple.pT[jPhi] >  fMaxPT2D ) break;
                if (evKaonCouple.pT[jPhi] <  fMinPT2D ) break;
                for (Int_t jpT_ = 0; jpT_ <= nBinPT2D; jpT_++ )
                {
                    jpT = jpT_;
                    if (evKaonCouple.pT[jPhi] <= fArrPT2D[jpT_+1]) break;
                }
                
                // Only non overlapping couples of Kaons
                if ( evKaonCouple.iKaon[iPhi] == evKaonCouple.iKaon[jPhi] ) continue;
                if ( evKaonCouple.iKaon[iPhi] == evKaonCouple.jKaon[jPhi] ) continue;
                if ( evKaonCouple.jKaon[iPhi] == evKaonCouple.iKaon[jPhi] ) continue;
                if ( evKaonCouple.jKaon[iPhi] == evKaonCouple.jKaon[jPhi] ) continue;
                
                hdM_dpT_Tot_Rec2D[ipT][jpT]->Fill(evKaonCouple.InvMass[iPhi],evKaonCouple.InvMass[jPhi],0.5);
                
            }
        }
    }

    TFile *outFile = new TFile(oFilePreP2D,"recreate");
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hdM_dpT_Tot_Rec2D[iHisto][jHisto]->Write();
            hdM_dpT_Tot_Bkg2D[iHisto][jHisto]->Write();
        }
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Tot_Rec[iHisto] -> Write();
    }
    outFile     ->Close();
    
    return 0;
}
