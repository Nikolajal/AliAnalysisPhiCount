void Process ( TString fInputName = "out.root",  TString fOutputName = "out2.root" ) {
    
    TFile      *fInputFile  =   new TFile (fInputName);
    TH2F       *hMPI        =   (TH2F*)(fInputFile->Get("hMPI"));
    TH2F       *hMult       =   (TH2F*)(fInputFile->Get("hMult"));
    
    auto nMPI_Entries = hMPI->GetEntries();
    hMPI->Scale(1./nMPI_Entries);
    auto nMultEntries = hMult->GetEntries();
    hMult->Scale(1./nMultEntries);
    
    TH1F * hMPI_1D          =   new TH1F ("hMPI_1D","hMPI_1D", 100, 0., 100.);
    TH1F * hMult1D          =   new TH1F ("hMult1D","hMult1D", 100, 0., 500.);
    TH2F * hMPI_12          =   new TH2F ("hMPI_12","hMPI_12", 100, 0., 100., 100, 0., 100.);
    TH2F * hMult12          =   new TH2F ("hMult12","hMult12", 100, 0., 500., 100, 0., 100.);
    TH1F * hMPI_13          =   new TH1F ("hMPI_13","hMPI_13", 100, 0., 100.);
    TH1F * hMult13          =   new TH1F ("hMult13","hMult13", 100, 0., 500.);
    TProfile *hMPI_PR       =   hMPI->ProfileX("hMPI_PR",-1,10000,"s");
    TProfile *hMultPR       =   hMult->ProfileX("hMultPR",-1,10000,"s");
    
    for ( Int_t i = 0; i < 100; i++ )
    {
        auto    MPI_Proj    =   hMPI->ProfileY(Form("projMPI__%i",i),i+1,i+1);
        
        auto    MPI_Mean    =   hMPI_PR->GetBinContent(i)-0.5;
        auto    MPI_MErr    =   0.;//(MPI_Proj    ->GetMeanError())/(MPI_Proj    ->GetMean());
        auto    MPI_StDv    =   hMPI_PR->GetBinError(i);
        auto    MPI_SErr    =   0.;//2*(MPI_Proj    ->GetStdDevError());
        
        if ( MPI_Mean != 0.5 && MPI_Mean != 0 )   {
            hMPI_1D->SetBinContent  (i+1,MPI_StDv*MPI_StDv/MPI_Mean-1);
            //hMPI_1D->SetBinError    (i+1,(MPI_MErr+MPI_SErr)*(MPI_StDv/MPI_Mean));
        }
        
        auto    MultProj    =   hMult->ProjectionY(Form("projMult_%i",i),i+1,i+1);
        
        auto    MultMean    =   hMultPR->GetBinContent(i)-0.5;
        auto    MultMErr    =   0.;//(MultProj    ->GetMeanError())/(MultProj    ->GetMean());
        auto    MultStDv    =   hMultPR->GetBinError(i);
        auto    MultSErr    =   0.;//2*(MultProj    ->GetStdDevError());
        
        if ( MultMean != -0.5 && MultMean != -0 )   {
            hMult1D->SetBinContent  (i+1,MultStDv*MultStDv/MultMean-1);
            //hMult1D->SetBinError    (i+1,(MultMErr+MultSErr)*(MultStDv/MultMean));
        }
        
        for ( Int_t iBin = 1; iBin <= 100; iBin++ )
        {
            for ( Int_t jBin = iBin; jBin < 100; jBin++ )
            {
                for ( Int_t kBin = 0; kBin < 100; kBin++ )
                {
                    auto fMultCenter = 5.*kBin;
                    hMPI_12->Fill(kBin,iBin,TMath::Binomial(jBin,iBin)*hMPI->GetBinContent(kBin,jBin+1));
                    hMult12->Fill(fMultCenter,iBin,TMath::Binomial(jBin,iBin)*hMult->GetBinContent(kBin,jBin+1));
                }
            }
        }
        
        auto    MPI_Yld1    =   hMPI_12->GetBinContent(i,2);
        auto    MPI_Err1    =   0.;//hMPI_12
        auto    MPI_Yld2    =   hMPI_12->GetBinContent(i,3);
        auto    MPI_Err2    =   0.;//hMPI_12
        
        if ( MPI_Yld1 != 0 )   {
            hMPI_13->SetBinContent  (i+1,2*MPI_Yld2/MPI_Yld1 - MPI_Yld1);
            //hMPI_1D->SetBinError    (i+1,(MPI_MErr+MPI_SErr)*(MPI_StDv/MPI_Mean));
        }
        
        auto    MultYld1    =   hMult12->GetBinContent(i,2);
        auto    MultErr1    =   0.;//hMult12
        auto    MultYld2    =   hMult12->GetBinContent(i,3);
        auto    MultErr2    =   0.;//hMult12
        
        if ( MultYld1 != 0 )   {
            hMult13->SetBinContent  (i+1,2*MultYld2/MultYld1 - MultYld1);
            //hMult1D->SetBinError    (i+1,(MultMErr+MultSErr)*(MultStDv/MultMean));
        }

    }
    
    TFile      *fOutputFile =   new TFile (fOutputName,"recreate");
    
    hMPI    ->Write();
    hMult   ->Write();
    hMPI_1D ->Write();
    hMult1D ->Write();
    hMPI_12 ->Write();
    hMult12 ->Write();
    hMPI_13 ->Write();
    hMult13 ->Write();
    hMPI_PR->Write();
    hMultPR->Write();
    
    fOutputFile->Close();
    fInputFile->Close();
}
