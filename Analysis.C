#include "SetValues.h"

//int main ()
void Analysis ()
{
    //Retrieving PreProcessing
    TFile *inFile = new TFile(outPP);
    //TFile *inFile = new TFile(outDF);
    
    //Target variables
    RooRealVar  xInvMass1D ("xInvMass1D","xInvMass1D",0.9874,maxBound);
    RooRealVar  xInvMass2D ("xInvMass2D","xInvMass2D",0.9874,maxBound);
    RooRealVar  yInvMass2D ("yInvMass2D","yInvMass2D",0.9874,maxBound);
    
    RooDataHist hK1D ("hK1D","hK1D",xInvMass1D,Import(*(TH1F*)inFile->Get("Kaon_InvMass1D")));
    RooDataHist hK2D ("hK2D","hK2D",RooArgList(xInvMass2D,yInvMass2D),Import(*(TH2F*)inFile->Get("Kaon_InvMass2D")));
    //RooDataHist hK3D ("hK3D","hK3D",mass,Import(*(TH1F*)inFile->Get("Kaon_InvM3D")));
    
    // Set-up model
    
    //Constants
    RooRealVar      KaonThreshold   ("KaonThreshold","KaonThreshold",0.9874);
    RooRealVar      phiMass         ("phiMass","phiMass",1.020);
    RooRealVar      phiWidth1D      ("phiWidth1D","phiWidth1D",0.0045);
    RooRealVar      phiWidth2D      ("phiWidth2D","phiWidth2D",0.0045);
    
    /* 1D */
    
    // Parameters
    RooRealVar      mExp1D          ("mExp1D","mExp1D",12.);
    RooRealVar      fracSig1D       ("fracSig1D","fracSig1D",0.084267,0.,0.1);
    
    // Functions
    RooGenericPdf   Bkg1D           ("Bkg1D","Bkg1D","sqrt(1-exp(-@2*(@0-@1)))",RooArgSet(xInvMass1D,KaonThreshold,mExp1D));
    
    RooBreitWigner  Sig1D           ("xSig1D","xSig1D",xInvMass1D,phiMass,phiWidth1D);
    
    RooAddPdf       Model1D         ("model1D","model1D",RooArgList(Sig1D,Bkg1D),fracSig1D);
    
    /* 2D */
    // Parameters
    RooRealVar      xmExp2D         ("mModExp","mModExp",12.);
    RooRealVar      ymExp2D         ("mModExp","mModExp",12.);
    RooRealVar      fracBkg2D_      ("fracBkg2D_","fracBkg2D_",0.084267,0.,1.0);
    RooRealVar      fracSig2D_      ("fracSig2D_","fracSig2D_",0.084267,0.,1.0);
    RooRealVar      fracSig2D       ("fracSig2D","fracSig2D",0.682861,0.,0.1);
    RooRealVar      offSetx         ("offSetx","offSetx",0.);
    RooRealVar      offSety         ("offSety","offSety",0.);
    
    // Functions
    RooGenericPdf   xBkg2D          ("xBkg2D","xBkg2D","sqrt(1-exp(-@2*(@0-@1)))-@3",RooArgSet(xInvMass2D,KaonThreshold,xmExp2D,offSetx));
    RooGenericPdf   yBkg2D          ("yBkg2D","yBkg2D","sqrt(1-exp(-@2*(@0-@1)))-@3",RooArgSet(yInvMass2D,KaonThreshold,ymExp2D,offSety));
    
    RooBreitWigner  xSig2D          ("xSig2D","xSig2D",xInvMass2D,phiMass,phiWidth2D);
    RooBreitWigner  ySig2D          ("ySig2D","ySig2D",yInvMass2D,phiMass,phiWidth2D);
    
    RooAddPdf       Bkg2D           ("Bkg2D","Bkg2D",RooArgList(xSig2D,xBkg2D),fracBkg2D_);
    RooAddPdf       Sig2D           ("Bkg2D","Bkg2D",RooArgList(ySig2D,yBkg2D),fracSig2D_);
    
    RooAddPdf       Model2D         ("model2D","model2D",RooArgList(Sig2D,Bkg2D),fracSig2D);
    
    
    
    // Fit  to  data
    Model1D.fitTo(hK1D);
    Model2D.fitTo(hK2D);
    
    Float_t fracSIG     = hK1D.sumEntries()*(fracSig1D.getVal());
    Float_t fracSIGhigh = hK1D.sumEntries()*(fracSig1D.getAsymErrorHi());
    Float_t fracSIGlow  = hK1D.sumEntries()*(fracSig1D.getAsymErrorLo());
    
    Float_t fracSIG2    = hK2D.sumEntries()*(fracSig2D.getVal());
    Float_t fracSIGhigh2= hK2D.sumEntries()*(fracSig2D.getAsymErrorHi());
    Float_t fracSIGlow2 = hK2D.sumEntries()*(fracSig2D.getAsymErrorLo());
    
     
    xInvMass2D.setRange("sigx",1.034,1.006);
    RooAbsReal* DPhiIntx= Sig2D.createIntegral(xInvMass2D,NormSet(xInvMass2D),Range("sigx"));
    yInvMass2D.setRange("sigy",1.034,1.006);
    RooAbsReal* DPhiInty= Sig2D.createIntegral(yInvMass2D,NormSet(yInvMass2D),Range("sigy"));
    Float_t DPhiCorrection  = (DPhiIntx->getVal()+DPhiInty->getVal());
    
    
    // Visual
    /*
    RooPlot* frame = xInvMass1D.frame();
    hK1D.plotOn(frame);
    Model1D.plotOn(frame);
    Model1D.plotOn(frame,Components("Bkg1D"),LineStyle(kDashed));
    frame->Draw();
    */
    
    RooPlot* frame2D = yInvMass2D.frame();
    hK2D.plotOn(frame2D,Range("sigy"));
    Model2D.plotOn(frame2D,Range("sigy"));
    frame2D->Draw();
    TCanvas *c2 = new TCanvas();
    //((TH2F*)inFile->Get("Kaon_InvMass2D"))->Draw("SURF2");
    (Model2D.createHistogram("Model2D",xInvMass2D,Binning(100), YVar(yInvMass2D,Binning(100))))->Draw("SURF2");
    
    cout << "--------------------------------" << endl;
    cout << "|  Time for some results!      |" << endl;
    cout << "________________________________" << endl;
    cout << "| 1 D  C A S E" << endl;
    cout << "| " << endl;
    cout << "| Total entries    : " << hK1D.sumEntries() << endl;
    cout << "| " << endl;
    cout << "| Expected signal  : 20329" << endl;
    cout << "| Measured signal  : " << fracSIG          << " + (" << fracSIGhigh << ") - (" << fracSIGlow << ")"<< endl;
    cout << "| Difference       : " << fracSIG-20329    << " + (" << fracSIGhigh << ") - (" << fracSIGlow << ")"<< endl;
    cout << "| " << endl;
    cout << "| Expected frac    : " << 20329/hK1D.sumEntries()   << endl;
    cout << "| fracSig1D        : " << fracSig1D.getVal()   << endl;
    cout << "| " << endl;
    cout << "| mExp1D           : " << mExp1D.getVal()      << endl;
    cout << "| phiWidth1D       : " << phiWidth1D.getVal()      << endl;
    cout << "________________________________" << endl;
    cout << "| 2 D  C A S E" << endl;
    cout << "| " << endl;
    cout << "| S I N G L E  P H I " << endl;
    cout << "| " << endl;
    cout << "| Total entries    : " << hK2D.sumEntries() << endl;
    cout << "| " << endl;
    cout << "| Expected signal  : " << 894608 << endl;
    cout << "| Measured signal  : " << fracSIG2         << " + (" << fracSIGhigh2 << ") - (" << fracSIGlow2 << ")"<< endl;
    cout << "| Difference       : " << fracSIG2-894608   << " + (" << fracSIGhigh2 << ") - (" << fracSIGlow2 << ")"<< endl;
    cout << "| " << endl;
    cout << "| Expected frac    : " << 894608/hK2D.sumEntries()   << endl;
    cout << "| fracSig2D        : " << fracSig2D.getVal()   << endl;
    cout << "| " << endl;
    cout << "| mxExp2D          : " << xmExp2D.getVal()         << endl;
    cout << "| myExp2D          : " << ymExp2D.getVal()         << endl;
    cout << "| fracBkg2D_ (x)   : " << fracBkg2D_.getVal()      << endl;
    cout << "| fracSig2D_ (x)   : " << fracSig2D_.getVal()      << endl;
    cout << "| offSetx          : " << offSetx.getVal()         << endl;
    cout << "| offSety          : " << offSety.getVal()         << endl;
    cout << "| phiWidth2D       : " << phiWidth2D.getVal()      << endl;
    cout << "| " << endl;
    cout << "| D O U B L E  P H I " << endl;
    cout << "| " << endl;
    cout << "| Expected signal  : " << 6802 << endl;
    cout << "| Measured signal  : " << fracSIG2*DPhiCorrection         << " + (" << fracSIGhigh2*DPhiCorrection << ") - (" << fracSIGlow2*DPhiCorrection << ")"<< endl;
    cout << "| Difference       : " << fracSIG2*DPhiCorrection-6802    << " + (" << fracSIGhigh2*DPhiCorrection << ") - (" << fracSIGlow2*DPhiCorrection << ")"<< endl;
    cout << "| " << endl;
    cout << "| Expected frac    : " << 6802/hK2D.sumEntries()   << endl;
    cout << "| fracSig2D        : " << fracSIG2*DPhiCorrection/hK2D.sumEntries()   << endl;
    cout << "| " << endl;
    cout << "| DPhiIntx         : " << DPhiIntx->getVal()   << endl;
    cout << "| DPhiInty         : " << DPhiInty->getVal()   << endl;
    cout << "________________________________" << endl;
    //return 0;
}
