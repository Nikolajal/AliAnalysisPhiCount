#include "../../inc/AliAnalysisPhiPair.h"
#include "PP_MC_Partial.cxx"
#include "PP_Resl.cxx"

void PP_MC ( TString fFolderDT = "/Volumes/\[HD\]\[Nikolajal\]_Toshiba\ 2/Dataset/_Data/20210805_LHC10X/Partial/", TString fFolderMC = "/Volumes/\[HD\]\[Nikolajal\]_Toshiba\ 2/Dataset/_Sim/20210805_LHC14j4X/Partial/", std::vector<pair<TString,TString>> fDataset = kpp7TeVDataset, TString fOption = "yield", Int_t nEventsCut = -1., TString kFolder = "", Bool_t kReEvaluate = true )   {
    // --- --- --- --- --- --- --- SET-UP --- --- --- --- --- --- --- --- --- --- ---
    //
    //  Check Correct Syntax
    if ( fFolderMC == "" )  {
        cout << "[ERROR] Must Specify an input root file" << endl;
        return;
    }
    //
    // --- INFO on Set-up variables
    if ( nEventsCut != -1 ) cout << "[INFO] Choosing to limit the datasample to " << nEventsCut << " events" <<endl;
    fChooseOption(fOption);
    //
    // --- Separating MC and Data datasets
    std::vector<TString>    fDTDataset;
    for ( auto kDataset : fDataset ) fDTDataset.push_back( kDataset.first );
    std::vector<TString>    fMCDataset;
    for ( auto kDataset : fDataset ) fMCDataset.push_back( kDataset.second );
    //
    // --- Starting MC Partial Analysis
    auto    iFil = 0;
    if ( kReEvaluate )  {
        cout << "[INFO] Re-preprocessing of single datasets requested" << endl;
        for ( auto kDataSet : fMCDataset )    {
            cout << "[INFO] Starting preprocessing of " << kDataSet.Data() << endl;
            PP_MC_Partial( fFolderMC + TString("/") + fMCDataset.at(iFil) + TString("/") + fMCDataset.at(iFil) + TString("_STD.root"), fOption, nEventsCut, kFolder + TString("/PreProcessing/") + fMCDataset.at(iFil) );
            iFil++;
        }
    }
    //
    std::vector<TFile*> insFile_Data;
    for ( auto kDataSet : fDTDataset )    insFile_Data.push_back( new TFile ( fFolderDT+kDataSet+TString("/")+kDataSet+TString("_STD.root") ) );
    std::vector<TFile*> insFile_MntC;
    if ( kDoMultiplicity )  for ( auto kDataSet : fMCDataset )  insFile_MntC.push_back( new TFile ( Form( kAnalysis_MCTruthHist, (TString("Multiplicity")+kFolder+TString("/PreProcessing/")+kDataSet+TString("/")).Data()  ) ) );
    if ( kDoYield )         for ( auto kDataSet : fMCDataset )  insFile_MntC.push_back( new TFile ( Form( kAnalysis_MCTruthHist, (TString("Yield")+kFolder+TString("/PreProcessing/")+kDataSet+TString("/")).Data() ) ) );
    //
    // --- Setting the input datastructure
    fSetAllBins();
    //
    std::vector<TH1D*>  hEvCountHistoUtility;
    std::vector<TH1D*>  hEvCountHistoUtilityMult;
    std::vector<TH1F*>  h1D_Nrec;
    std::vector<TH1F*>  h1D_Ngen;
    std::vector<TH1F*>  h1D_Nrec_2Db;
    std::vector<TH1F*>  h1D_Ngen_2Db;
    std::vector<TH1F*>  h1D_Nrec_Fin;
    std::vector<TH1F*>  h1D_Ngen_Fin;
    std::vector<std::vector<TH1F*>> h1D_TruInvMass;
    std::vector<std::vector<TH1F*>> h1D_RecInvMass;
    std::vector<std::vector<TH1F*>> h1D_InvMassRes;
    std::vector<std::vector<TH1F*>> h1D_TruInvMass_2Db;
    std::vector<std::vector<TH1F*>> h1D_RecInvMass_2Db;
    std::vector<std::vector<TH1F*>> h1D_InvMassRes_2Db;
    //
    iFil = 0;
    for ( auto kFile : insFile_MntC )   {
        // --- Single Histograms
        hEvCountHistoUtility        .push_back( new TH1D ( *(TH1D*)(( insFile_Data.at(iFil) )->Get("fQCOutputList")->FindObject("fQC_Event_Enum_FLL")) ) );
        hEvCountHistoUtilityMult    .push_back( new TH1D ( *(TH1D*)(( insFile_Data.at(iFil) )->Get("fQCOutputList")->FindObject("fQC_Event_Enum_V0M")) ) );
        h1D_Nrec                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC.at(iFil) )->Get("h1D_Nrec")) ) );
        h1D_Ngen                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC.at(iFil) )->Get("h1D_Ngen")) ) );
        h1D_Nrec_2Db                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC.at(iFil) )->Get("h1D_Nrec_2Db")) ) );
        h1D_Ngen_2Db                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC.at(iFil) )->Get("h1D_Ngen_2Db")) ) );
        h1D_Nrec_Fin                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC.at(iFil) )->Get("h1D_Nrec_Fin")) ) );
        h1D_Ngen_Fin                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC.at(iFil) )->Get("h1D_Ngen_Fin")) ) );
        hEvCountHistoUtility        .at(iFil)   ->  SetName( TString("fQC_Event_Enum_FLL_") +   fMCDataset.at(iFil) );
        hEvCountHistoUtilityMult    .at(iFil)   ->  SetName( TString("fQC_Event_Enum_V0M_") +   fMCDataset.at(iFil) );
        h1D_Nrec                    .at(iFil)   ->  SetName( TString("h1D_Nrec_")           +   fMCDataset.at(iFil) );
        h1D_Ngen                    .at(iFil)   ->  SetName( TString("h1D_Ngen_")           +   fMCDataset.at(iFil) );
        h1D_Nrec_2Db                .at(iFil)   ->  SetName( TString("h1D_Nrec_2Db_")       +   fMCDataset.at(iFil) );
        h1D_Ngen_2Db                .at(iFil)   ->  SetName( TString("h1D_Ngen_2Db_")       +   fMCDataset.at(iFil) );
        h1D_Nrec_Fin                .at(iFil)   ->  SetName( TString("h1D_Nrec_Fin_")       +   fMCDataset.at(iFil) );
        h1D_Ngen_Fin                .at(iFil)   ->  SetName( TString("h1D_Ngen_Fin_")       +   fMCDataset.at(iFil) );
        // --- std::vector Histograms
        h1D_TruInvMass      .push_back( uLoadHistograms<1,TH1F>( insFile_MntC.at(iFil), "h1D_TruInvMass_PT_%i",     TString("h1D_TruInvMass_PT_%i_")    + fMCDataset.at(iFil) ) );
        h1D_RecInvMass      .push_back( uLoadHistograms<1,TH1F>( insFile_MntC.at(iFil), "h1D_RecInvMass_PT_%i",     TString("h1D_RecInvMass_PT_%i_")    + fMCDataset.at(iFil) ) );
        h1D_InvMassRes      .push_back( uLoadHistograms<1,TH1F>( insFile_MntC.at(iFil), "h1D_InvMassRes_PT_%i",     TString("h1D_InvMassRes_PT_%i_")    + fMCDataset.at(iFil) ) );
        h1D_TruInvMass_2Db  .push_back( uLoadHistograms<1,TH1F>( insFile_MntC.at(iFil), "h1D_TruInvMass_2Db_PT_%i", TString("h1D_TruInvMass_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
        h1D_RecInvMass_2Db  .push_back( uLoadHistograms<1,TH1F>( insFile_MntC.at(iFil), "h1D_RecInvMass_2Db_PT_%i", TString("h1D_RecInvMass_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
        h1D_InvMassRes_2Db  .push_back( uLoadHistograms<1,TH1F>( insFile_MntC.at(iFil), "h1D_InvMassRes_2Db_PT_%i", TString("h1D_InvMassRes_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
        //
        iFil++;
    }
    //
    // --- Setting the output datastructure
    //
    std::vector<TH1F*>  h1D_Eff;
    std::vector<TH1F*>  h1D_Eff_2Db;
    std::vector<TH2F*>  h2D_Eff_2Db;
    std::vector<Float_t>    kPeriodWeights;
    //
    // --- Evaluating Periods Weighting
    for ( auto kEvCountHist : hEvCountHistoUtility )    kPeriodWeights.push_back( kEvCountHist->GetBinContent(kEventCount::kVertex10) );
    uAddSumHistogram( hEvCountHistoUtility, "fQC_Event_Enum_FLL" );
    for ( auto&& kWeight : kPeriodWeights )             kWeight     /=  hEvCountHistoUtility.at(0)->GetBinContent(kEventCount::kVertex10);
    //
    // --- --- --- --- --- --- --- ANALYSIS --- --- --- --- --- --- --- --- --- --- -
    //
    iFil = 0;
    for ( auto kFile : insFile_MntC )   {
        // --- Single Dataset Efficiency
        auto    kEffUtil        =   (TH1F*)(h1D_Nrec.at(iFil)->Clone());
        kEffUtil                ->  SetName( TString("h1D_Eff_")  +   fMCDataset.at(iFil) );
        kEffUtil                ->  Divide( h1D_Nrec.at(iFil), h1D_Ngen.at(iFil), 1., 1., "b" );
        h1D_Eff                 .push_back( kEffUtil );
        //
        iFil++;
    }
    //
    uAddSumHistogram( h1D_Eff,              "h1D_Eff",                  kPeriodWeights );
    uAddSumHistogram( h1D_Nrec,             "h1D_Nrec",                 kPeriodWeights );
    uAddSumHistogram( h1D_Ngen,             "h1D_Ngen",                 kPeriodWeights );
    uAddSumHistogram( h1D_Nrec_2Db,         "h1D_Nrec_2Db",             kPeriodWeights );
    uAddSumHistogram( h1D_Ngen_2Db,         "h1D_Ngen_2Db",             kPeriodWeights );
    uAddSumHistogram( h1D_Nrec_Fin,         "h1D_Nrec_Fin",             kPeriodWeights );
    uAddSumHistogram( h1D_Ngen_Fin,         "h1D_Ngen_Fin",             kPeriodWeights );
    uAddSumHistogram( h1D_TruInvMass,       "h1D_TruInvMass_PT_%i",     kPeriodWeights );
    uAddSumHistogram( h1D_RecInvMass,       "h1D_RecInvMass_PT_%i",     kPeriodWeights );
    uAddSumHistogram( h1D_InvMassRes,       "h1D_InvMassRes_PT_%i",     kPeriodWeights );
    uAddSumHistogram( h1D_TruInvMass_2Db,   "h1D_TruInvMass_2Db_PT_%i", kPeriodWeights );
    uAddSumHistogram( h1D_RecInvMass_2Db,   "h1D_RecInvMass_2Db_PT_%i", kPeriodWeights );
    uAddSumHistogram( h1D_InvMassRes_2Db,   "h1D_InvMassRes_2Db_PT_%i", kPeriodWeights );
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- --- YIELD ANALYSIS
    if ( kDoYield ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/Full/")).Data()));
        TFile*  outFile_Yield   =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Yield")+kFolder).Data()),"recreate");
        auto    kPlotDirectory  =   TString(Form(kAnalysis_PreProc_Dir,(TString("Yield")+kFolder).Data()))+TString("/Plots/");
        //
        // --- Saving to File
        for ( auto kSave    : h1D_Eff )                                             kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_2Db )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_2Db )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_Fin )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_Fin )                                        kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        // --- Printing to Plots
        SetStyle();
        auto    fLegend =   fMCDataset;
        fLegend.insert( fLegend.begin(), "Inclusive" );
        gROOT           ->  SetBatch(kTRUE);
        //
        auto    cDrawEff    =   uPlotEfficiencies( h1D_Eff, fLegend );
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        // TODO: Ratio plot under the efficiencies to see how they compare
        cDrawEff        ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficiencies1D.pdf") );
        delete cDrawEff;
        //
        gROOT           ->  SetBatch(kFALSE);
        //
        outFile_Yield   ->  Close();
    }
    //  --- --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity/")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity/")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity/")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile*  outFile_Yield   =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity/")+kFolder).Data()),"recreate");
        auto    kPlotDirectory  =   TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/");
        //
        // --- Saving to File
        for ( auto kSave    : h1D_Eff )                                             kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_2Db )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_2Db )                                        kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        // --- Printing to Plots
        SetStyle();
        auto    fLegend =   fMCDataset;
        fLegend.insert( fLegend.begin(), "Inclusive" );
        gROOT           ->  SetBatch(kTRUE);
        //
        auto    cDrawEff    =   uPlotEfficiencies( h1D_Eff, fLegend );
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        cDrawEff        ->  SaveAs( kPlotDirectory + TString("hFullEfficiencies.pdf") );
        delete cDrawEff;
        //
        gROOT           ->  SetBatch(kFALSE);
        //
        outFile_Yield   ->  Close();
    }
    //  --- --- CORRELATION ANALYSIS
    if ( kDoCorrelation ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Correlation/")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation/")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation/")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile*  outFile_Yield   =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Correlation/")+kFolder).Data()),"recreate");
        auto    kPlotDirectory  =   TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/");
        //
        // --- Saving to File
        for ( auto kSave    : h1D_Eff )                                             kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_2Db )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_2Db )                                        kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        // --- Printing to Plots
        SetStyle();
        auto    fLegend =   fMCDataset;
        fLegend.insert( fLegend.begin(), "Inclusive" );
        gROOT           ->  SetBatch(kTRUE);
        //
        auto    cDrawEff    =   uPlotEfficiencies( h1D_Eff, fLegend );
        //
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        cDrawEff        ->  SaveAs( kPlotDirectory + TString("hFullEfficiencies.pdf") );
        delete cDrawEff;
        //
        gROOT           ->  SetBatch(kFALSE);
        //
        outFile_Yield   ->  Close();
    }
    //
    for ( auto kFile : insFile_MntC ) kFile->Close();
}
