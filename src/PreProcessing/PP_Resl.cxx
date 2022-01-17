#include "../../inc/AliAnalysisPhiPair.h"

void PP_Resl ( TString fOption = "all", TString kFolder = "", Bool_t fSilent = kTRUE )    {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    // --- Silencing warnings for smoother running
    if ( fSilent )  {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
        gErrorIgnoreLevel   =   kWarning;
    }
    //
    // --- INFO on Set-up variables
    fChooseOption(fOption);
    //
    fSetAllBins();
    auto    kPlotFolder1D   =   TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/1D/");
    auto    kPlotFolder2D   =   TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/2D/");
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s/Resolution/",kPlotFolder1D.Data()));
    gROOT                   ->  ProcessLine(Form(".! mkdir -p %s/Resolution/",kPlotFolder2D.Data()));
    TFile*  insFile_Data        =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()));
    //
    // --- Retrieving Histograms
    std::vector<TH1F*>  h1D_TruInvMass      =   uLoadHistograms<1,TH1F>( insFile_Data, "h1D_TruInvMass_PT_%i" );
    std::vector<TH1F*>  h1D_RecInvMass      =   uLoadHistograms<1,TH1F>( insFile_Data, "h1D_RecInvMass_PT_%i" );
    std::vector<TH1F*>  h1D_InvMassRes      =   uLoadHistograms<1,TH1F>( insFile_Data, "h1D_InvMassRes_PT_%i" );
    std::vector<TH1F*>  h1D_TruInvMass_2Db  =   uLoadHistograms<1,TH1F>( insFile_Data, "h1D_TruInvMass_2Db_PT_%i" );
    std::vector<TH1F*>  h1D_RecInvMass_2Db  =   uLoadHistograms<1,TH1F>( insFile_Data, "h1D_RecInvMass_2Db_PT_%i" );
    std::vector<TH1F*>  h1D_InvMassRes_2Db  =   uLoadHistograms<1,TH1F>( insFile_Data, "h1D_InvMassRes_2Db_PT_%i" );
    //
    //  TODO: create folders before here so that the functions can plot control histos
    gROOT->SetBatch( kTRUE );
    TH1F*   h1D_Template        =   new TH1F( "", "", nBinPT1D, fArrPT1D );
    TH1F*   h1D_Template_2Db    =   new TH1F( "", "", nBinPT2D, fArrPT2D );
    //
    
    auto h1D_FullResolutions    =   uCalculateResolution ( h1D_InvMassRes,      h1D_RecInvMass,     h1D_TruInvMass,     h1D_Template,       { kPhiMesonMass_, 1.005,  1.0335 }, { kPhiMesonWidth, kPhiMesonWidth*0.9, kPhiMesonWidth*1.1 }, 0.001, kPlotFolder1D+TString("/Resolution/") );
    cout << "2D" << endl;
    auto h2Db_FullResolutions   =   uCalculateResolution ( h1D_InvMassRes_2Db,  h1D_RecInvMass_2Db, h1D_TruInvMass_2Db, h1D_Template_2Db,   { kPhiMesonMass_, 1.005,  1.0335 }, { kPhiMesonWidth, kPhiMesonWidth*0.9, kPhiMesonWidth*1.1 }, 0.001, kPlotFolder2D+TString("/Resolution/") );
    gROOT->SetBatch( kFALSE );
    //
    for ( auto kSave : h1D_FullResolutions  )   kSave->SetName(Form("%s_1D",    kSave->GetName()));
    for ( auto kSave : h2Db_FullResolutions )   kSave->SetName(Form("%s_2Db",   kSave->GetName()));
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- --- YIELD ANALYSIS
    if ( kDoYield ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kMassResolution_Anal,(TString("Yield")+kFolder).Data()),"recreate");
        //
        // --- Saving to File
        for ( auto kSave : h1D_FullResolutions  )   kSave->Write();
        for ( auto kSave : h2Db_FullResolutions )   kSave->Write();
        //
        // --- Printing to Plots
        gROOT->SetBatch( kTRUE );
        auto    kDrawCanvases   =   uPlotResolution( h1D_FullResolutions );
        for ( auto kDraw : kDrawCanvases )  kDraw   ->  SaveAs( kPlotFolder1D + TString(Form("/%s_1D.pdf",kDraw->GetName())) );
                kDrawCanvases   =   uPlotResolution( h2Db_FullResolutions );
        for ( auto kDraw : kDrawCanvases )  kDraw   ->  SaveAs( kPlotFolder2D + TString(Form("/%s_2Db.pdf",kDraw->GetName())) );
        gROOT->SetBatch( kFALSE );
        //
        outFile_Yield   ->  Close();
    }
    //  --- --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kMassResolution_Anal,(TString("Multiplicity")+kFolder).Data()),"recreate");
        auto    kPlotFolder1D   =   TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/1D/");
        auto    kPlotFolder2D   =   TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/2D/");
        //
        // --- Saving to File
        for ( auto kSave : h1D_FullResolutions  )   kSave->Write();
        for ( auto kSave : h2Db_FullResolutions )   kSave->Write();
        //
        // --- Printing to Plots
        gROOT->SetBatch( kTRUE );
        auto    kDrawCanvases   =   uPlotResolution( h1D_FullResolutions );
        for ( auto kDraw : kDrawCanvases )  kDraw   ->  SaveAs( kPlotFolder1D + TString(Form("/%s_1D.pdf",kDraw->GetName())) );
                kDrawCanvases   =   uPlotResolution( h2Db_FullResolutions );
        for ( auto kDraw : kDrawCanvases )  kDraw   ->  SaveAs( kPlotFolder2D + TString(Form("/%s_2Db.pdf",kDraw->GetName())) );
        gROOT->SetBatch( kFALSE );
        //
        outFile_Yield   ->  Close();
    }
    //  --- --- CORRELATION ANALYSIS
    if ( kDoCorrelation ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile *outFile_Yield    =   new TFile   (Form(kMassResolution_Anal,(TString("Correlation")+kFolder).Data()),"recreate");
        auto    kPlotFolder1D   =   TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/1D/");
        auto    kPlotFolder2D   =   TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/2D/");
        //
        // --- Saving to File
        for ( auto kSave : h1D_FullResolutions  )   kSave->Write();
        for ( auto kSave : h2Db_FullResolutions )   kSave->Write();
        //
        outFile_Yield->Close();
        //
        // --- Printing to Plots
        gROOT->SetBatch( kTRUE );
        auto    kDrawCanvases   =   uPlotResolution( h1D_FullResolutions );
        for ( auto kDraw : kDrawCanvases )  kDraw   ->  SaveAs( kPlotFolder1D + TString(Form("/%s_1D.pdf",kDraw->GetName())) );
                kDrawCanvases   =   uPlotResolution( h2Db_FullResolutions );
        for ( auto kDraw : kDrawCanvases )  kDraw   ->  SaveAs( kPlotFolder2D + TString(Form("/%s_2Db.pdf",kDraw->GetName())) );
        gROOT->SetBatch( kFALSE );
        //
    }
}
