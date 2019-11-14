#include "../inc/SetValues.h"

int main ()
{
    //Retrieving Event/MC data
    TFile *inFile = new TFile(oFilePreBKG);
        
    //Retrieving Event/MC TTree
    TTree *PTreeKS  = (TTree*)inFile->Get(PTreeNameKS);
    TTree *PTreeKD  = (TTree*)inFile->Get(PTreeNameKD);
        
    // Define some simple data structures
    EVKAONCOUPLE evKaonS;
    EVKAONCOUPLE evKaonD;
        
    //Setting Branch Addresses
    // Filling TTree
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.nKaonCouple",&evKaonD.nKaonCouple);
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.iKaon",&evKaonD.iKaon);
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.jKaon",&evKaonD.jKaon);
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.bPhi",&evKaonD.bPhi);
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.bRec",&evKaonD.bRec);
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.InvMass",&evKaonD.InvMass);
    PTreeKD  ->SetBranchAddress    ("evKaonCouple.pT",&evKaonD.pT);
    
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.nKaonCouple",&evKaonS.nKaonCouple);
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.iKaon",&evKaonS.iKaon);
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.jKaon",&evKaonS.jKaon);
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.bPhi",&evKaonS.bPhi);
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.bRec",&evKaonS.bRec);
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.InvMass",&evKaonS.InvMass);
    PTreeKS  ->SetBranchAddress    ("evKaonCouple.pT",&evKaonS.pT);
    
    vSetBinsPT1D();
    vSetBinsIM1D();
    vSetBinsPT2D();
    vSetBinsIM2D();
    
    auto hName  = "Name";
    auto hTitle = "Title";
    TH1F ** hdM_dpT_Rec_BB1D        = new TH1F *    [nBinPT2D];
    TH2F ***hdM_dpT_Rec_BB2D        = new TH2F **   [nBinPT2D];
    TH2F ***hdM_dpT_Rec_SB2D        = new TH2F **   [nBinPT2D];
    Float_t** IMS                   = new Float_t*  [nBinPT2D];
    Float_t** IMD                   = new Float_t*  [nBinPT2D];
    Int_t * indexIMS                = new Int_t     [nBinPT2D];
    Int_t * indexIMD                = new Int_t     [nBinPT2D];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Rec_BB2D[iHisto]    = new TH2F * [nBinPT2D];
        hdM_dpT_Rec_SB2D[iHisto]    = new TH2F * [nBinPT2D];
        IMS[iHisto]                 = new Float_t [PEvents*25];
        IMD[iHisto]                 = new Float_t [PEvents*25];
        indexIMS[iHisto]            = 0;
        indexIMD[iHisto]            = 0;
        
        // Setting up 1D Histogram
        hName = Form("hdM_dpT_Rec_BB1D_%i",iHisto);
        hTitle= Form("m_{K_{#pm}K_{#pm}} in p_{T} range %f to %f",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        hdM_dpT_Rec_BB1D[iHisto] = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            // Setting up 2D Histogram
            hName = Form("hdM_dpT_Rec_BB2D_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K_{#pm}K_{#pm}} in p_{T} range %f to %f and %f to %f",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hdM_dpT_Rec_BB2D[iHisto][jHisto] = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            
            hName = Form("hdM_dpT_Rec_SB2D_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K_{+}K_{-}} in p_{T} range %f to %f and m_{K_{#pm}K_{#pm}} in p_{T} range %f to %f",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hdM_dpT_Rec_SB2D[iHisto][jHisto] = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
        }
    }
    
    Int_t ipT, jpT;
    for (Int_t iEvent = 0; iEvent < PTreeKS->GetEntries(); iEvent++)
    {
        PTreeKS->GetEntry(iEvent);
        for (Int_t iPhi = 0; iPhi < evKaonS.nKaonCouple; iPhi++ )
        {
            // Only Recordable Phi Candidates
            //if (evKaonS.bRec[iPhi] == false) continue;
            
            // Only True Phi Candidates
            //if (evKaonS.bPhi[iPhi] == false) continue;
            
            // Information on pT based on defined bins
            if (evKaonS.pT[iPhi] >  fMaxPT2D ) continue;
            if (evKaonS.pT[iPhi] <  fMinPT2D ) continue;
            for (Int_t ipT_ = 0; ipT_ <= nBinPT2D; ipT_++ )
            {
                if (evKaonS.pT[iPhi] < fArrPT2D[ipT_+1] && evKaonS.pT[iPhi] >= fArrPT2D[ipT_])
                {
                    ipT = ipT_;
                }
            }
            hdM_dpT_Rec_BB1D[ipT]->Fill(evKaonS.InvMass[iPhi]);
            IMS[ipT][indexIMS[ipT]] = evKaonS.InvMass[iPhi];
            indexIMS[ipT]++;
        }
        PTreeKD->GetEntry(iEvent);
        for (Int_t iPhi = 0; iPhi < evKaonD.nKaonCouple; iPhi++ )
        {
            // Only Recordable Phi Candidates
            //if (evKaonS.bRec[iPhi] == false) continue;
            
            // Only True Phi Candidates
            //if (evKaonS.bPhi[iPhi] == false) continue;
            
            // Information on pT based on defined bins
            if (evKaonD.pT[iPhi] >  fMaxPT2D ) continue;
            if (evKaonD.pT[iPhi] <  fMinPT2D ) continue;
            for (Int_t ipT_ = 0; ipT_ <= nBinPT2D; ipT_++ )
            {
                if (evKaonD.pT[iPhi] < fArrPT2D[ipT_+1] && evKaonD.pT[iPhi] > fArrPT2D[ipT_])
                {
                    ipT = ipT_;
                }
            }
            IMD[ipT][indexIMD[ipT]] = evKaonD.InvMass[iPhi];
            indexIMD[ipT]++;
        }
    }
    
    for (ipT = 0; ipT < nBinPT2D; ipT++ )
    {
        for (jpT = 0; jpT < nBinPT2D; jpT++ )
        {
            for (Int_t iPhi = 0; iPhi < indexIMS[ipT]; iPhi++ )
            {
                for (Int_t jPhi = 0; jPhi < indexIMS[jpT]; jPhi++ )
                {
                    hdM_dpT_Rec_BB2D[ipT][jpT]->Fill(IMS[ipT][iPhi],IMS[jpT][jPhi],0.5);
                }
                for (Int_t jPhi = 0; jPhi < indexIMD[jpT]; jPhi++ )
                {
                    hdM_dpT_Rec_SB2D[ipT][jpT]->Fill(IMS[ipT][iPhi],IMD[jpT][jPhi],0.5);
                    hdM_dpT_Rec_SB2D[ipT][jpT]->Fill(IMD[ipT][iPhi],IMS[jpT][jPhi],0.5);
                }
            }
        }
    }
    
    TFile *outFile = new TFile(oFilePrBKG2,"recreate");
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hdM_dpT_Rec_BB2D[iHisto][jHisto]->Write();
            hdM_dpT_Rec_SB2D[iHisto][jHisto]->Write();
        }
    }
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hdM_dpT_Rec_BB1D[iHisto] -> Write();
    }
    inFile->Close();
    outFile->Close();
    
    return 0;
}
