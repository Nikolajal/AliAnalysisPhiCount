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
    
    vSetBinsPT1D();
    vSetBinsIM1D();
    TH1F * hPhiGen_1D   = new TH1F ("hPhiGen_1D","Generated #varphi per |y|<5 in K_{+}K_{-} Decay mode, 1D analysis",nBinPT1D,fArrPT1D);
    hPhiGen_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");
    TH1F * hPhiRec_1D   = new TH1F ("hPhiRec_1D","Generated #varphi per |y|<5 in K_{+}K_{-} Decay mode, recordable in ALICE exp, 1D analysis",nBinPT1D,fArrPT1D);
    hPhiRec_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");
    TH1F * hPhiTru_1D   = new TH1F ("hPhiTru_1D","Generated #varphi per |y|<5, 1D analysis",nBinPT1D,fArrPT1D);
    hPhiTru_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");
    TH1F * hPhiEff_1D   = new TH1F ("hPhiEff_1D","#varphi reconstruction efficiency in ALICE exp (|y|<5), 1D analysis",nBinPT1D,fArrPT1D);
    hPhiEff_1D->GetXaxis()->SetTitle("pT #varphi (GeV)");
    
    /*------*/
    /*  2D  */
    /*------*/
    
    vSetBinsPT2D();
    vSetBinsIM2D();
    TH2F * hPhiEff_1D_2D    = new TH2F ("hPhiEff_1D_2D","hPhiEff_1D_2D",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2F * hPhiChk_xD   = new TH2F ("hPhiChk_xD","hPhiChk_xD",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH1F * hPhiEff_1D2D_= new TH1F ("hPhiEff_1D2D_","hPhiEff_1D2D_",nBinPT2D,fArrPT2D);
    
    
    TH2F * hPhiGen_2D   = new TH2F ("hPhiGen_2D","Generated #varphi per |y|<5 in K_{+}K_{-} Decay mode, 2D analysis",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hPhiGen_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiGen_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    TH2F * hPhiRec_2D   = new TH2F ("hPhiRec_2D","Generated #varphi per |y|<5 in K_{+}K_{-} Decay mode recordable in ALICE exp, 2D analysis",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hPhiRec_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiRec_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    TH2F * hPhiTru_2D   = new TH2F ("hPhiTru_2D","Generated #varphi per |y|<5, 2D analysis",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hPhiTru_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiTru_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    TH2F * hPhiEff_2D   = new TH2F ("hPhiEff_2D","#varphi reconstruction efficiency in ALICE exp (|y|<5), 2D analysis",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hPhiEff_2D->GetXaxis()->SetTitle("pT #varphi_{1} (GeV)");
    hPhiEff_2D->GetYaxis()->SetTitle("pT #varphi_{2} (GeV)");
    
    //N_Gen & N_Rec
    for (Int_t iEvent = 0; iEvent < PtreeK2->GetEntries(); iEvent++)
    {
        PtreeK2->GetEntry(iEvent);
        for (int iPhi = 0; iPhi < evKaonCouple.nKaonCouple; iPhi++ )
        {
            if (evKaonCouple.bPhi[iPhi] == false) continue;
            hPhiGen_1D->Fill(evKaonCouple.pT[iPhi]);
            if (evKaonCouple.bRec[iPhi] == false) continue;
            hPhiRec_1D->Fill(evKaonCouple.pT[iPhi]);
        }
        for (int iPhi = 0; iPhi < evKaonCouple.nKaonCouple; iPhi++ )
        {
            for (int jPhi = 0; jPhi < evKaonCouple.nKaonCouple; jPhi++ )
            {
                // Only non overlapping couples of Kaons
                if ( evKaonCouple.iKaon[iPhi] == evKaonCouple.iKaon[jPhi] ) continue;
                if ( evKaonCouple.iKaon[iPhi] == evKaonCouple.jKaon[jPhi] ) continue;
                if ( evKaonCouple.jKaon[iPhi] == evKaonCouple.iKaon[jPhi] ) continue;
                if ( evKaonCouple.jKaon[iPhi] == evKaonCouple.jKaon[jPhi] ) continue;
                
                // Only True Phi Candidates
                if (evKaonCouple.bPhi[iPhi] == false) continue;
                if (evKaonCouple.bPhi[jPhi] == false) continue;
                
                hPhiGen_2D->Fill(evKaonCouple.pT[iPhi],evKaonCouple.pT[jPhi],0.5);
                
                // Only Recordable Phi Candidates
                if (evKaonCouple.bRec[iPhi] == false) continue;
                if (evKaonCouple.bRec[jPhi] == false) continue;
                
                hPhiRec_2D->Fill(evKaonCouple.pT[iPhi],evKaonCouple.pT[jPhi],0.5);
            }
        }
    }
    
    //N_Tru
    for (Int_t iEvent = 0; iEvent < PtreePhi->GetEntries(); iEvent++)
    {
        PtreePhi->GetEntry(iEvent);
        for (Int_t iPhi = 0; iPhi < evPhi.nPhi; iPhi++)
        {
            hPhiTru_1D->Fill(evPhi.pT[iPhi]);
            for (Int_t jPhi = 0; jPhi < evPhi.nPhi; jPhi++)
            {
                if (iPhi == jPhi) continue;
                hPhiTru_2D->Fill(evPhi.pT[iPhi],evPhi.pT[jPhi],0.5);
            }
        }
    }
    
    //Efficiency
    hPhiEff_1D->Divide(hPhiRec_1D,hPhiGen_1D,1,1,"b");
    hPhiEff_2D->Divide(hPhiRec_2D,hPhiGen_2D,1,1,"b");
    
    //Check on Efficiency
    hPhiEff_1D2D_           ->Divide((hPhiRec_1D),(hPhiGen_1D),1,1,"b"); //->Rebin(1,"h2")
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hPhiEff_1D_2D->SetBinContent(iHisto+1,jHisto+1,(hPhiEff_1D2D_->GetBinContent(iHisto+1))*(hPhiEff_1D2D_->GetBinContent(jHisto+1)));
            hPhiEff_1D_2D->SetBinError(iHisto+1,jHisto+1,(hPhiEff_1D2D_->GetBinContent(iHisto+1))*(hPhiEff_1D2D_->GetBinContent(jHisto+1))*((hPhiEff_1D2D_->GetBinError(iHisto+1))/(hPhiEff_1D2D_->GetBinContent(iHisto+1))+(hPhiEff_1D2D_->GetBinError(jHisto+1))/(hPhiEff_1D2D_->GetBinContent(jHisto+1))));
        }
    }
    hPhiChk_xD              ->Divide(hPhiEff_2D,hPhiEff_1D_2D,1,1,"b");
    
    // Writing results out to file
    TFile * outFile = new TFile(oFileEffici,"recreate");
    hPhiTru_1D      ->Write();
    hPhiTru_2D      ->Write();
    hPhiGen_1D      ->Write();
    hPhiGen_2D      ->Write();
    hPhiRec_1D      ->Write();
    hPhiRec_2D      ->Write();
    /*
    hPhiEff_1D      ->SetMaximum(1.1);
    hPhiEff_1D_2D   ->SetMaximum(1.1);
    hPhiEff_1D2D_   ->SetMaximum(1.1);
    hPhiEff_2D      ->SetMaximum(1.1);
    hPhiChk_xD      ->SetMaximum(1.5);
    
    hPhiEff_1D      ->SetMinimum(0.);
    hPhiEff_1D_2D   ->SetMinimum(0.);
    hPhiEff_1D2D_   ->SetMinimum(0.);
    hPhiEff_2D      ->SetMinimum(0.);
    hPhiChk_xD      ->SetMinimum(0.5);
    */
    hPhiEff_1D      ->Write();
    hPhiEff_1D_2D   ->Write();
    hPhiEff_1D2D_   ->Write();
    hPhiEff_2D      ->Write();
    hPhiChk_xD      ->Write();
    
    outFile         ->Close();
    return 0;
}
