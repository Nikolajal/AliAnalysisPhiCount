typedef struct  {
    Int_t       nPrt,          nMult05,        nMult08,    nMult10,    fPDG[1024];
    Float_t     fPx[1024],      fPy[1024],      fPz[1024];
} SParticleData;

std::map<Int_t,Float_t> kPDGMass;
void
SetAllMasses
(){
    kPDGMass[333] = 1.019455;   //  GeV/c^2
}

void
QuickAnalysis
 ( TString fInputFileName, TString fInputTreeName, TString fOutputFolder = "", Long_t fParticlePDG = 333 ) {
    //
    //  --- Pre-Processing
    //
    //  --- --- Recover Tree
    TFile*  fInputFile  =   new TFile(fInputFileName);
    TTree*  fInputTree  =   (TTree*)(fInputFile->Get(fInputTreeName));
    //  --- --- Recover Histograms
    TH1D*   hEventCount =   (TH1D*)(fInputFile->Get("kEventCount"));
    //  --- --- Create Data Structure
    SParticleData   kCurrentParticle;
    fInputTree  ->  SetBranchAddress("nPrt",       &kCurrentParticle.nPrt);
    fInputTree  ->  SetBranchAddress("nMult05",     &kCurrentParticle.nMult05);
    fInputTree  ->  SetBranchAddress("nMult08",     &kCurrentParticle.nMult08);
    fInputTree  ->  SetBranchAddress("nMult10",     &kCurrentParticle.nMult10);
    fInputTree  ->  SetBranchAddress("fPDG",        &kCurrentParticle.fPDG);
    fInputTree  ->  SetBranchAddress("fPx",         &kCurrentParticle.fPx);
    fInputTree  ->  SetBranchAddress("fPy",         &kCurrentParticle.fPy);
    fInputTree  ->  SetBranchAddress("fPz",         &kCurrentParticle.fPz);
    //  --- --- Recover Entries
    Int_t   kTreeEntries    =   fInputTree  -> GetEntries();
    //  --- --- Set Used Masses
    SetAllMasses();
    //
    //
    //  --- Output Structure
    TH1F*   hPrtPTSpectrum  =   new TH1F( "hPrtPTSpectrum", "hPrtPTSpectrum", 200, 0, 20 );
    hPrtPTSpectrum  ->  GetXaxis()  ->SetTitle( "p_{T} ( GeV/c )" );
    hPrtPTSpectrum  ->  GetYaxis()  ->SetTitle( "dN/dy" );
    TH1F*   hYield  =   new TH1F( "hYield", "hYield", 6, 0, 6 );
    hYield          ->  GetXaxis()  ->SetTitle( "" );
    hYield          ->  GetYaxis()  ->SetTitle( "" );
    //
    //
    //  --- Analysis
    //
    //  --- --- Looping Events
    for ( Int_t iEv = 0; iEv < kTreeEntries; iEv++ ) {
        fInputTree  -> GetEvent( iEv );
        //  --- --- --- Looping Particles
        Int_t   nParticles  =   0;
        for ( Int_t iPrt = 0; iPrt < kCurrentParticle.nPrt; iPrt++ ) {
            //  --- --- --- --- Select chosen particle
            if ( kCurrentParticle.fPDG[iPrt] != fParticlePDG ) continue;
            //  --- --- --- --- Build TLorentz for kinematics
            TLorentzVector  kCurrentParticleTLV;
            kCurrentParticleTLV.SetXYZM( kCurrentParticle.fPx[iPrt], kCurrentParticle.fPy[iPrt], kCurrentParticle.fPz[iPrt], kPDGMass[kCurrentParticle.fPDG[iPrt]] );
            //  --- --- --- --- Physics Cuts
            if ( fabs( kCurrentParticleTLV.Rapidity() ) > 0.5 ) continue;
            hPrtPTSpectrum  -> Fill( kCurrentParticleTLV.Pt() );
            nParticles++;
        }
        hYield  -> Fill( nParticles );
    }
    //  --- --- Post-Processing Histograms
    Int_t   kEventNorm  =   hEventCount -> GetEntries();
    hPrtPTSpectrum  ->  Scale( 1., "width" );
    hPrtPTSpectrum  ->  Scale( 1./kEventNorm );
    //
    //  --- Output
    //
    //  --- --- Save To File
    TFile*  fOutputFile =   new TFile( fOutputFolder + TString("./QuickAnalysisResults.root"), "recreate" );
    hPrtPTSpectrum  ->  Write();
    //  --- --- Save To PDF
    gROOT->SetBatch( kTRUE );
    TCanvas*    cDrawQuickResults   = new TCanvas("cDrawQuickResults","cDrawQuickResults",1000,1000);
    gStyle              -> SetOptStat( 0 );
    gPad                -> SetLogy();
    hPrtPTSpectrum      -> Draw();
    cDrawQuickResults   -> SaveAs( fOutputFolder + TString("./hPrtPTSpectrum.pdf") );
    delete cDrawQuickResults;
    gROOT->SetBatch( kFALSE );
}
