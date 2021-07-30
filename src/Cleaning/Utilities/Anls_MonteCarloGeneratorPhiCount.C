// Pythia
#include "Pythia8/Pythia.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"

using namespace Pythia8;

void Anls_MonteCarloGeneratorPhiCount(TString,Int_t,Int_t);
void fOutputAnalysis (TString,TString);

int main (int argc, char *argv[])
{
    // Check everything is good
    if (argc < 4)
    {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./GeneratorMC [filename] [nevents] [seed]" << endl;
        return -1;
    }
    
    if ( atoi(argv[2]) != -1 ) Anls_MonteCarloGeneratorPhiCount(argv[1],atoi(argv[2]),atoi(argv[3]));
    else fOutputAnalysis (argv[1],argv[3]);
    
    return 0;
}

void Anls_MonteCarloGeneratorPhiCount ( TString fOutputName = "out.root", Int_t nEvents = 0, Int_t nSeed = 0)
{
    // Utilities
    int   fnPhi, fMult, fnMPI;

    // Output File
    TFile * outFile     = new   TFile   (fOutputName,       "recreate","",101);
    
    TH2F * hMPI         =   new TH2F ("hMPI","hMPI",    20, 0., 50., 5, -0.5, 4.5);
    TH2F * hMult        =   new TH2F ("hMult","hMult",  20, 0., 500., 5, -0.5, 4.5);
    
    cout << "[INFO] Initialisation of Pythia8" << endl;
    
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    
    //Settings
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %i",nSeed));
    pythia.init();
    
    cout << "[INFO] Starting loops" << endl;
    
    // Cycling through events
    for ( int iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Next event
        pythia.next();
        
        // Resetting counters
        fnPhi             = 0;
        fMult             = 0;
        fnMPI             = 0;
        
        // Starting cycling through event particles
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            const auto particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( particle.id() == 333 )
            {
                if ( fabs(particle.p().rap()) < 0.5 ) fnPhi++;
            }
            
            //Skipping non-final particles
            if ( !particle.isFinal() )       continue;
            fMult++;
        }
        fnMPI   =   pythia.info.nMPI();
        
        // Fill
        hMPI    -> Fill(fnMPI,fnPhi);
        hMult   -> Fill(fMult,fnPhi);
    }
    
    hMPI    ->  Write();
    hMult   ->  Write();
    outFile     ->Close();
    return;
}

void fOutputAnalysis ( TString fInputName = "out.root",  TString fOutputName = "out2.root" ) {
    
    TFile      *fInputFile  =   new TFile (fInputName);
    TH2F       *hMPI        =   (TH2F*)(fInputFile->Get("hMPI"));
    TH2F       *hMult       =   (TH2F*)(fInputFile->Get("hMult"));
    
    auto nMPI_Entries = hMPI->GetEntries();
    hMPI->Scale(1./nMPI_Entries);
    auto nMultEntries = hMult->GetEntries();
    hMult->Scale(1./nMultEntries);
    
    TH1F       *hMPI_1D     =   new TH1F    ("hMPI_1D","hMPI_1D",20, 0., 50.);
    TH1F       *hMult1D     =   new TH1F    ("hMult1D","hMult1D",20, 0., 500.);
    
    for ( Int_t i = 0; i < 20; i++ )
    {
        auto    MPI_Proj    =   hMPI->ProjectionY(Form("projMPI__%i",i),i+1,i+1);
        
        auto    MPI_Mean    =   MPI_Proj    ->GetMean();
        auto    MPI_MErr    =   (MPI_Proj    ->GetMeanError())/(MPI_Proj    ->GetMean());
        auto    MPI_StDv    =   pow(MPI_Proj    ->GetStdDev(),2);
        auto    MPI_SErr    =   2*(MPI_Proj    ->GetStdDevError());
        
        if ( MPI_Mean != 0 )   {
            hMPI_1D->SetBinContent  (i,MPI_StDv/MPI_Mean-1);
            hMPI_1D->SetBinError    (i,(MPI_MErr+MPI_SErr)*(MPI_StDv/MPI_Mean));
        }
        
        auto    MultProj    =   hMult->ProjectionY(Form("projMult_%i",i),i+1,i+1);
        
        auto    MultMean    =   MultProj    ->GetMean();
        auto    MultMErr    =   (MultProj    ->GetMeanError())/(MultProj    ->GetMean());
        auto    MultStDv    =   pow(MultProj    ->GetStdDev(),2);
        auto    MultSErr    =   2*(MultProj    ->GetStdDevError());
        
        if ( MultMean != 0 )   {
            hMult1D->SetBinContent  (i,MultStDv/MultMean-1);
            hMult1D->SetBinError    (i,(MultMErr+MultSErr)*(MultStDv/MultMean));
        }

        
        /*
        auto fMPI_bi1   =   hMPI->GetBin(i,2);
        auto fMPI_bi2   =   hMPI->GetBin(i,3);
        
        auto fMultbi1   =   hMult->GetBin(i,2);
        auto fMultbi2   =   hMult->GetBin(i,3);
        
        auto fMPI_Yield1    =   hMPI->GetBinContent (fMPI_bi1);
        auto fMPI_Error1    =   hMPI->GetBinError   (fMPI_bi1);
        auto fMPI_Yield2    =   hMPI->GetBinContent (fMPI_bi2);
        auto fMPI_Error2    =   hMPI->GetBinError   (fMPI_bi2);
        
        auto fMultYield1    =   hMPI->GetBinContent (fMultbi1);
        auto fMultError1    =   hMPI->GetBinError   (fMultbi1);
        auto fMultYield2    =   hMPI->GetBinContent (fMultbi2);
        auto fMultError2    =   hMPI->GetBinError   (fMultbi2);
        
        auto fMPI_gammap    =   2*fMPI_Yield2/fMPI_Yield1 -fMPI_Yield1;
        auto fMPI_gammaE    =   (2/fMPI_Yield1)*sqrt(1+((fMPI_Error1*fMPI_Error1)/(fMPI_Yield1*fMPI_Yield1))+((fMPI_Yield1*fMPI_Yield1)/(4.))+(fMPI_Error2*fMPI_Error2));
        
        auto fMultgammap    =   2*fMultYield2/fMultYield1 -fMultYield1;
        auto fMultgammaE    =   (2/fMultYield1)*sqrt(1+((fMultError1*fMultError1)/(fMultYield1*fMultYield1))+((fMultYield1*fMultYield1)/(4.))+(fMultError2*fMultError2));
         */
    }
    
    TFile      *fOutputFile =   new TFile (fOutputName,"recreate");
    
    hMPI    ->Write();
    hMult   ->Write();
    hMPI_1D ->Write();
    hMult1D ->Write();
    
    fOutputFile->Close();
    fInputFile->Close();
}
