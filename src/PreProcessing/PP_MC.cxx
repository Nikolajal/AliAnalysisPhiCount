#include "../../inc/AliAnalysisPhiPair.h"
#include "PP_MC_Partial.cxx"
#include "PP_Resl.cxx"

void PP_MC ( TString fFolderDT, TString fFolderMC, std::vector<pair<TString,TString>> fDataset, TString fOption, Int_t nEventsCut, TString kFolder, Bool_t kReEvaluate, TString kType )   {
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
            PP_MC_Partial( fFolderMC + TString("/") + fMCDataset.at(iFil) + TString("/") + fMCDataset.at(iFil) + TString("_") + kType + TString(".root"), fOption, nEventsCut, kFolder + TString("/PreProcessing/") + fMCDataset.at(iFil) );
            iFil++;
        }
    }
    //
    std::vector<TFile*> insFile_Data;
    for ( auto kDataSet : fDTDataset )    insFile_Data.push_back( new TFile ( fFolderDT+kDataSet+TString("/")+kDataSet + TString("_") + kType + TString(".root") ) );
    std::vector<TFile*> insFile_MntC_YLD;
    std::vector<TFile*> insFile_MntC_MLT;
    if ( kDoMultiplicity )  for ( auto kDataSet : fMCDataset )  insFile_MntC_MLT.push_back( new TFile ( Form( kAnalysis_MCTruthHist, (TString("Multiplicity")+kFolder+TString("/PreProcessing/")+kDataSet+TString("/")).Data()  ) ) );
    if ( kDoYield )         for ( auto kDataSet : fMCDataset )  insFile_MntC_YLD.push_back( new TFile ( Form( kAnalysis_MCTruthHist, (TString("Yield")+kFolder+TString("/PreProcessing/")+kDataSet+TString("/")).Data() ) ) );
    //
    // --- Setting the input datastructure
    fSetAllBins();
    //
    // --- Histograms for Yield Analysis
    std::vector<TH1D*>  hEvCountHistoUtility;
    std::vector<TH1D*>  hEvCountHistoUtilityMult;
    std::vector<TH1F*>  h1D_Nrec;
    std::vector<TH1F*>  h1D_Ngen;
    std::vector<TH1F*>  h1D_Ntru;
    std::vector<TH2F*>  h2D_Nrec;
    std::vector<TH2F*>  h2D_Ngen;
    std::vector<TH2F*>  h2D_Ntru;
    std::vector<TH1F*>  h1D_Nrec_2Db;
    std::vector<TH1F*>  h1D_Ngen_2Db;
    std::vector<TH1F*>  h1D_Nrec_Fin;
    std::vector<TH1F*>  h1D_Ngen_Fin;
    std::vector<TH1F*>  h1D_Ngen_VTX10;
    std::vector<TH1F*>  h1D_Ngen_VTX10_2Db;
    std::vector<TH2F*>  h2D_Ngen_VTX10;
    std::vector<std::vector<TH1F*>> h1D_TruInvMass;
    std::vector<std::vector<TH1F*>> h1D_RecInvMass;
    std::vector<std::vector<TH1F*>> h1D_InvMassRes;
    std::vector<std::vector<TH1F*>> h1D_TruInvMass_2Db;
    std::vector<std::vector<TH1F*>> h1D_RecInvMass_2Db;
    std::vector<std::vector<TH1F*>> h1D_InvMassRes_2Db;
    //
    // --- Histograms for Multiplicity Analysis
    std::vector<std::vector<TH1F*>>  h1D_Mult_Nrec;
    std::vector<std::vector<TH1F*>>  h1D_Mult_Ngen;
    std::vector<std::vector<TH1F*>>  h1D_Mult_Nrec_2Db;
    std::vector<std::vector<TH1F*>>  h1D_Mult_Ngen_2Db;
    std::vector<std::vector<TH1F*>>  h1D_Mult_Nrec_Fin;
    std::vector<std::vector<TH1F*>>  h1D_Mult_Ngen_Fin;
    std::vector<std::vector<TH1F*>>  h1D_Mult_Ngen_VTX10;
    //
    iFil = 0;
    if ( kDoYield ) {
        for ( auto kFile : insFile_MntC_YLD )   {
            // --- Single Histograms
            hEvCountHistoUtility        .push_back( new TH1D ( *(TH1D*)(( insFile_Data.at(iFil) )->Get("fQCOutputList")->FindObject("fQC_Event_Enum_FLL")) ) );
            hEvCountHistoUtilityMult    .push_back( new TH1D ( *(TH1D*)(( insFile_Data.at(iFil) )->Get("fQCOutputList")->FindObject("fQC_Event_Enum_V0M")) ) );
            h1D_Nrec                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Nrec")) ) );
            h1D_Ngen                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Ngen")) ) );
            h1D_Ntru                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Ntru")) ) );
            h2D_Nrec                    .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_YLD.at(iFil) )->Get("h2D_Nrec")) ) );
            h2D_Ngen                    .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_YLD.at(iFil) )->Get("h2D_Ngen")) ) );
            h2D_Ntru                    .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_YLD.at(iFil) )->Get("h2D_Ntru")) ) );
            h1D_Nrec_2Db                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Nrec_2Db")) ) );
            h1D_Ngen_2Db                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Ngen_2Db")) ) );
            h1D_Nrec_Fin                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Nrec_Fin")) ) );
            h1D_Ngen_Fin                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Ngen_Fin")) ) );
            h1D_Ngen_VTX10              .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Ngen_VTX10")) ) );
            h1D_Ngen_VTX10_2Db          .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_YLD.at(iFil) )->Get("h1D_Ngen_VTX10_2Db")) ) );
            h2D_Ngen_VTX10              .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_YLD.at(iFil) )->Get("h2D_Ngen_VTX10")) ) );
            hEvCountHistoUtility        .at(iFil)   ->  SetName( TString("fQC_Event_Enum_FLL_") +   fMCDataset.at(iFil) );
            hEvCountHistoUtilityMult    .at(iFil)   ->  SetName( TString("fQC_Event_Enum_V0M_") +   fMCDataset.at(iFil) );
            h1D_Nrec                    .at(iFil)   ->  SetName( TString("h1D_Nrec_")           +   fMCDataset.at(iFil) );
            h1D_Ngen                    .at(iFil)   ->  SetName( TString("h1D_Ngen_")           +   fMCDataset.at(iFil) );
            h1D_Ntru                    .at(iFil)   ->  SetName( TString("h1D_Ntru_")           +   fMCDataset.at(iFil) );
            h2D_Nrec                    .at(iFil)   ->  SetName( TString("h2D_Nrec_")           +   fMCDataset.at(iFil) );
            h2D_Ngen                    .at(iFil)   ->  SetName( TString("h2D_Ngen_")           +   fMCDataset.at(iFil) );
            h2D_Ntru                    .at(iFil)   ->  SetName( TString("h2D_Ntru_")           +   fMCDataset.at(iFil) );
            h1D_Nrec_2Db                .at(iFil)   ->  SetName( TString("h1D_Nrec_2Db_")       +   fMCDataset.at(iFil) );
            h1D_Ngen_2Db                .at(iFil)   ->  SetName( TString("h1D_Ngen_2Db_")       +   fMCDataset.at(iFil) );
            h1D_Nrec_Fin                .at(iFil)   ->  SetName( TString("h1D_Nrec_Fin_")       +   fMCDataset.at(iFil) );
            h1D_Ngen_Fin                .at(iFil)   ->  SetName( TString("h1D_Ngen_Fin_")       +   fMCDataset.at(iFil) );
            h1D_Ngen_VTX10              .at(iFil)   ->  SetName( TString("h1D_Ngen_VTX10_")     +   fMCDataset.at(iFil) );
            h1D_Ngen_VTX10_2Db          .at(iFil)   ->  SetName( TString("h1D_Ngen_VTX10_2Db_") +   fMCDataset.at(iFil) );
            h2D_Ngen_VTX10              .at(iFil)   ->  SetName( TString("h2D_Ngen_VTX10_")     +   fMCDataset.at(iFil) );
            // --- std::vector Histograms
            h1D_TruInvMass              .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_YLD.at(iFil), "h1D_TruInvMass_PT_%i",     TString("h1D_TruInvMass_PT_%i_")    + fMCDataset.at(iFil) ) );
            h1D_RecInvMass              .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_YLD.at(iFil), "h1D_RecInvMass_PT_%i",     TString("h1D_RecInvMass_PT_%i_")    + fMCDataset.at(iFil) ) );
            h1D_InvMassRes              .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_YLD.at(iFil), "h1D_InvMassRes_PT_%i",     TString("h1D_InvMassRes_PT_%i_")    + fMCDataset.at(iFil) ) );
            h1D_TruInvMass_2Db          .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_YLD.at(iFil), "h1D_TruInvMass_2Db_PT_%i", TString("h1D_TruInvMass_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
            h1D_RecInvMass_2Db          .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_YLD.at(iFil), "h1D_RecInvMass_2Db_PT_%i", TString("h1D_RecInvMass_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
            h1D_InvMassRes_2Db          .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_YLD.at(iFil), "h1D_InvMassRes_2Db_PT_%i", TString("h1D_InvMassRes_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
            //
            iFil++;
        }
    }
    if ( kDoMultiplicity ) {
        iFil = 0;
        for ( auto kFile : insFile_MntC_MLT )   {
            // --- Single Histograms
            if ( !kDoYield )    {
                hEvCountHistoUtility        .push_back( new TH1D ( *(TH1D*)(( insFile_Data.at(iFil) )->Get("fQCOutputList")->FindObject("fQC_Event_Enum_FLL")) ) );
                hEvCountHistoUtilityMult    .push_back( new TH1D ( *(TH1D*)(( insFile_Data.at(iFil) )->Get("fQCOutputList")->FindObject("fQC_Event_Enum_V0M")) ) );
                h1D_Nrec                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Nrec")) ) );
                h1D_Ngen                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Ngen")) ) );
                h1D_Ntru                    .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Ntru")) ) );
                h2D_Nrec                    .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_MLT.at(iFil) )->Get("h2D_Nrec")) ) );
                h2D_Ngen                    .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_MLT.at(iFil) )->Get("h2D_Ngen")) ) );
                h2D_Ntru                    .push_back( new TH2F ( *(TH2F*)(( insFile_MntC_MLT.at(iFil) )->Get("h2D_Ntru")) ) );
                h1D_Nrec_2Db                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Nrec_2Db")) ) );
                h1D_Ngen_2Db                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Ngen_2Db")) ) );
                h1D_Nrec_Fin                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Nrec_Fin")) ) );
                h1D_Ngen_Fin                .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Ngen_Fin")) ) );
                h1D_Ngen_VTX10              .push_back( new TH1F ( *(TH1F*)(( insFile_MntC_MLT.at(iFil) )->Get("h1D_Ngen_VTX10")) ) );
                hEvCountHistoUtility        .at(iFil)   ->  SetName( TString("fQC_Event_Enum_FLL_") +   fMCDataset.at(iFil) );
                hEvCountHistoUtilityMult    .at(iFil)   ->  SetName( TString("fQC_Event_Enum_V0M_") +   fMCDataset.at(iFil) );
                h1D_Nrec                    .at(iFil)   ->  SetName( TString("h1D_Nrec_")           +   fMCDataset.at(iFil) );
                h1D_Ngen                    .at(iFil)   ->  SetName( TString("h1D_Ngen_")           +   fMCDataset.at(iFil) );
                h1D_Ntru                    .at(iFil)   ->  SetName( TString("h1D_Ntru_")           +   fMCDataset.at(iFil) );
                h2D_Nrec                    .at(iFil)   ->  SetName( TString("h2D_Nrec_")           +   fMCDataset.at(iFil) );
                h2D_Ngen                    .at(iFil)   ->  SetName( TString("h2D_Ngen_")           +   fMCDataset.at(iFil) );
                h2D_Ntru                    .at(iFil)   ->  SetName( TString("h2D_Ntru_")           +   fMCDataset.at(iFil) );
                h1D_Nrec_2Db                .at(iFil)   ->  SetName( TString("h1D_Nrec_2Db_")       +   fMCDataset.at(iFil) );
                h1D_Ngen_2Db                .at(iFil)   ->  SetName( TString("h1D_Ngen_2Db_")       +   fMCDataset.at(iFil) );
                h1D_Nrec_Fin                .at(iFil)   ->  SetName( TString("h1D_Nrec_Fin_")       +   fMCDataset.at(iFil) );
                h1D_Ngen_Fin                .at(iFil)   ->  SetName( TString("h1D_Ngen_Fin_")       +   fMCDataset.at(iFil) );
                h1D_Ngen_VTX10              .at(iFil)   ->  SetName( TString("h1D_Ngen_VTX10_")     +   fMCDataset.at(iFil) );
                // --- std::vector Histograms
                h1D_TruInvMass              .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_TruInvMass_PT_%i",     TString("h1D_TruInvMass_PT_%i_")    + fMCDataset.at(iFil) ) );
                h1D_RecInvMass              .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_RecInvMass_PT_%i",     TString("h1D_RecInvMass_PT_%i_")    + fMCDataset.at(iFil) ) );
                h1D_InvMassRes              .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_InvMassRes_PT_%i",     TString("h1D_InvMassRes_PT_%i_")    + fMCDataset.at(iFil) ) );
                h1D_TruInvMass_2Db          .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_TruInvMass_2Db_PT_%i", TString("h1D_TruInvMass_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
                h1D_RecInvMass_2Db          .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_RecInvMass_2Db_PT_%i", TString("h1D_RecInvMass_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
                h1D_InvMassRes_2Db          .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_InvMassRes_2Db_PT_%i", TString("h1D_InvMassRes_2Db_PT_%i_")+ fMCDataset.at(iFil) ) );
            }
            h1D_Mult_Nrec                   .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Nrec_%i",         TString("h1D_Mult_Nrec_%i_")        + fMCDataset.at(iFil) ) );
            h1D_Mult_Ngen                   .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Ngen_%i",         TString("h1D_Mult_Ngen_%i_")        + fMCDataset.at(iFil) ) );
            h1D_Mult_Nrec_2Db               .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Nrec_2Db_%i",     TString("h1D_Mult_Nrec_2Db_%i_")    + fMCDataset.at(iFil) ) );
            h1D_Mult_Ngen_2Db               .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Ngen_2Db_%i",     TString("h1D_Mult_Ngen_2Db_%i_")    + fMCDataset.at(iFil) ) );
            h1D_Mult_Nrec_Fin               .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Nrec_Fin_%i",     TString("h1D_Mult_Nrec_Fin_%i_")    + fMCDataset.at(iFil) ) );
            h1D_Mult_Ngen_Fin               .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Ngen_Fin_%i",     TString("h1D_Mult_Ngen_Fin_%i_")    + fMCDataset.at(iFil) ) );
            h1D_Mult_Ngen_VTX10             .push_back( uLoadHistograms<1,TH1F>( insFile_MntC_MLT.at(iFil), "h1D_Mult_Ngen_VTX10_%i",   TString("h1D_Mult_Ngen_VTX10_%i_")  + fMCDataset.at(iFil) ) );
            //
            iFil++;
        }
    }
    //
    // --- Setting the output datastructure
    //
    std::vector<TH1F*>              h1D_Eff;
    std::vector<TH1F*>              h1D_SLC;
    std::vector<std::vector<TH1F*>> h1D_Mult_Eff;
    std::vector<std::vector<TH1F*>> h1D_Mult_SLC;
    std::vector<Float_t>            kPeriodWeights;
    std::vector<Float_t>            kUniformWeight;
    //
    // --- Evaluating Periods Weighting
    for ( auto kEvCountHist : hEvCountHistoUtility )    kPeriodWeights.push_back( kEvCountHist->GetBinContent(kEventCount::kVertex10) );
    uAddSumHistogram( hEvCountHistoUtility, "fQC_Event_Enum_FLL" );
    for ( auto&& kWeight : kPeriodWeights )             kWeight     /=  hEvCountHistoUtility.at(0)->GetBinContent(kEventCount::kVertex10);
    //
    // --- --- --- --- --- --- --- ANALYSIS --- --- --- --- --- --- --- --- --- --- -
    //
    if ( kDoYield ) {
        iFil = 0;
        for ( auto kFile : insFile_MntC_YLD )   {
            // --- Single Dataset Efficiency
            auto    kCurrentPeriod_Efficiency   =   (TH1F*)( h1D_Nrec.at(iFil)->Clone() );
                    kCurrentPeriod_Efficiency   ->  SetName( TString("h1D_Eff_")  +   fMCDataset.at(iFil) );
                    kCurrentPeriod_Efficiency   ->  Divide( h1D_Nrec.at(iFil), h1D_Ngen.at(iFil), 1., 1., "b" );
            h1D_Eff .push_back( kCurrentPeriod_Efficiency );
            //
            // --- Single Dataset Signal Loss Efficiency
            auto    kCurrentTarget_SigLossCor   =   (TH1F*) ( h1D_Ngen.at(iFil)->Clone() );
                    kCurrentTarget_SigLossCor   ->  SetName ( TString("h1D_SLC_")  +   fMCDataset.at(iFil) );
                    kCurrentTarget_SigLossCor   ->  Divide  ( h1D_Ngen_VTX10.at(iFil), h1D_Ngen.at(iFil), 1., 1., "b" );
            h1D_SLC .push_back( kCurrentTarget_SigLossCor );
            //
            iFil++;
        }
    }
    if ( kDoMultiplicity ) {
        iFil = 0;
        for ( auto kFile : insFile_MntC_MLT )   {
            if ( !kDoYield ) {
                // --- Single Dataset Efficiency
                auto    kCurrentPeriod_Efficiency   =   (TH1F*)( h1D_Nrec.at(iFil)->Clone() );
                        kCurrentPeriod_Efficiency   ->  SetName( TString("h1D_Eff_") + fMCDataset.at(iFil) );
                        kCurrentPeriod_Efficiency   ->  Divide( h1D_Nrec.at(iFil), h1D_Ngen.at(iFil), 1., 1., "b" );
                h1D_Eff .push_back( kCurrentPeriod_Efficiency );
                //
                // --- Single Dataset Signal Loss Efficiency
                auto    kCurrentTarget_SigLossCor   =   (TH1F*) ( h1D_Ngen.at(iFil)->Clone() );
                        kCurrentTarget_SigLossCor   ->  SetName ( TString("h1D_SLC_")  +   fMCDataset.at(iFil) );
                        kCurrentTarget_SigLossCor   ->  Divide  ( h1D_Ngen_VTX10.at(iFil), h1D_Ngen.at(iFil), 1., 1., "b" );
                h1D_SLC .push_back( kCurrentTarget_SigLossCor );
            }
            auto iMult = 0;
            std::vector<TH1F*>  kMultiplicityEfficiency;
            std::vector<TH1F*>  kMultiplicitySignalLoss;
            for ( auto kCurrentRec : h1D_Mult_Nrec.at(iFil) ) {
                // --- Single Dataset-Multiplicity Efficiency
                auto    kCurrentTarget_Efficiency   =   (TH1F*)( kCurrentRec->Clone() );
                        kCurrentTarget_Efficiency   ->  SetName( TString( Form("h1D_Mutl_Eff_%i_",iMult) ) + fMCDataset.at(iFil) );
                        kCurrentTarget_Efficiency   ->  Divide( kCurrentRec, h1D_Mult_Ngen.at(iFil).at(iMult), 1., 1., "b" );
                //
                kMultiplicityEfficiency             .push_back( kCurrentTarget_Efficiency );
                //
                // --- Single Dataset-Multiplicity Signal Loss Efficiency
                auto    kCurrentTarget_SigLossCor   =   (TH1F*)( kCurrentRec->Clone() );
                        kCurrentTarget_SigLossCor   ->  SetName( TString( Form("h1D_Mutl_SLC_%i_",iMult) ) + fMCDataset.at(iFil) );
                        kCurrentTarget_SigLossCor   ->  Divide( h1D_Mult_Ngen_VTX10.at(iFil).at(iMult), h1D_Mult_Ngen.at(iFil).at(iMult), 1., 1., "b" );
                //
                kMultiplicitySignalLoss             .push_back( kCurrentTarget_SigLossCor );
                iMult++;
            }
            h1D_Mult_Eff                  .push_back( kMultiplicityEfficiency );
            h1D_Mult_SLC                  .push_back( kMultiplicitySignalLoss );
            iFil++;
        }
    }
    //
    uAddSumHistogram( h1D_Eff,              "h1D_Eff",                  kPeriodWeights );
    uAddSumHistogram( h1D_SLC,              "h1D_SLC",                  kPeriodWeights );
    uAddSumHistogram( h1D_Nrec,             "h1D_Nrec",                 kPeriodWeights );
    uAddSumHistogram( h1D_Ngen,             "h1D_Ngen",                 kPeriodWeights );
    uAddSumHistogram( h1D_Ntru,             "h1D_Ntru",                 kPeriodWeights );
    uAddSumHistogram( h2D_Nrec,             "h2D_Nrec",                 kPeriodWeights );
    uAddSumHistogram( h2D_Ngen,             "h2D_Ngen",                 kPeriodWeights );
    uAddSumHistogram( h2D_Ntru,             "h2D_Ntru",                 kPeriodWeights );
    uAddSumHistogram( h1D_Nrec_2Db,         "h1D_Nrec_2Db",             kPeriodWeights );
    uAddSumHistogram( h1D_Ngen_2Db,         "h1D_Ngen_2Db",             kPeriodWeights );
    uAddSumHistogram( h1D_Nrec_Fin,         "h1D_Nrec_Fin",             kPeriodWeights );
    uAddSumHistogram( h1D_Ngen_Fin,         "h1D_Ngen_Fin",             kPeriodWeights );
    uAddSumHistogram( h1D_Ngen_VTX10,       "h1D_Ngen_VTX10",           kPeriodWeights );
    uAddSumHistogram( h1D_Ngen_VTX10_2Db,   "h1D_Ngen_VTX10_2Db",       kPeriodWeights );
    uAddSumHistogram( h2D_Ngen_VTX10,       "h2D_Ngen_VTX10",           kPeriodWeights );
    uAddSumHistogram( h1D_TruInvMass,       "h1D_TruInvMass_PT_%i",     kPeriodWeights );
    uAddSumHistogram( h1D_RecInvMass,       "h1D_RecInvMass_PT_%i",     kPeriodWeights );
    uAddSumHistogram( h1D_InvMassRes,       "h1D_InvMassRes_PT_%i",     kPeriodWeights );
    uAddSumHistogram( h1D_TruInvMass_2Db,   "h1D_TruInvMass_2Db_PT_%i", kPeriodWeights );
    uAddSumHistogram( h1D_RecInvMass_2Db,   "h1D_RecInvMass_2Db_PT_%i", kPeriodWeights );
    uAddSumHistogram( h1D_InvMassRes_2Db,   "h1D_InvMassRes_2Db_PT_%i", kPeriodWeights );
    //
    // --- --- --- --- --- --- --- OUTPUT --- --- --- --- --- --- --- --- --- --- ---
    //
    //  --- --- GENERAL
    //
    SetStyle();
    gROOT   ->  SetBatch(kTRUE);
    auto    kPeriodLegend   =   fMCDataset;
    push_to_front( kPeriodLegend, TString("Inclusive") );
    //
    //  --- --- --- 1D Period Efficiency
    auto    cDrawPeriodEfficiency   =   uPlotEfficiencies( h1D_Eff, kPeriodLegend, "cDrawPeriodEfficiency" );
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
    if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
    uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    //
    //  --- --- --- 1D Period Signal Loss
    auto    cDrawPeriodSignalLoss   =   uPlotEfficiencies( h1D_SLC, kPeriodLegend, "cDrawPeriodSignalLoss", true );
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
    if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
    uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    //
    //  --- --- --- 2D vs 1D Full Efficiency
    auto    cDraw1Dvs2DEfficiency   =   uEfficiencyCompare_1D_2D( h1D_Nrec_2Db.at(0), h1D_Ngen_2Db.at(0), h2D_Nrec.at(0), h2D_Ngen.at(0)  );
    cDraw1Dvs2DEfficiency->cd();
    cDraw1Dvs2DEfficiency->cd(1);
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.08);
    uLatex->DrawLatexNDC(0.14, 0.88,"ALICE Performance");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.06);
    if ( is_pp_anl ) uLatex->DrawLatexNDC(0.14, 0.80, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
    if ( is_pb_anl ) uLatex->DrawLatexNDC(0.14, 0.80, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
    uLatex->DrawLatexNDC(0.14, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    //
    //  --- --- --- 2D vs 1D Full Signal Loss
    auto    cDraw1Dvs2DSignalLoss   =   uEfficiencyCompare_1D_2D( h1D_Ngen_VTX10_2Db.at(0), h1D_Ngen_2Db.at(0), h2D_Ngen_VTX10.at(0), h2D_Ngen.at(0), true );
    cDraw1Dvs2DSignalLoss->cd();
    cDraw1Dvs2DSignalLoss->cd(1);
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.08);
    uLatex->DrawLatexNDC(0.14, 0.88,"ALICE Performance");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.06);
    if ( is_pp_anl ) uLatex->DrawLatexNDC(0.14, 0.80, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
    if ( is_pb_anl ) uLatex->DrawLatexNDC(0.14, 0.80, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
    uLatex->DrawLatexNDC(0.14, 0.73,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
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
        for ( auto kSave    : h1D_SLC )                                             kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Ntru )                                            kSave   ->  Write();
        for ( auto kSave    : h2D_Nrec )                                            kSave   ->  Write();
        for ( auto kSave    : h2D_Ngen )                                            kSave   ->  Write();
        for ( auto kSave    : h2D_Ntru )                                            kSave   ->  Write();
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
        //  --- Printing to Plots
        //
        //  --- --- General QA Plots
        cDrawPeriodEfficiency   ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficiencies1D.pdf") );
        cDrawPeriodSignalLoss   ->  SaveAs( kPlotDirectory + TString("/Full/hFullSignalLosses1D.pdf") );
        cDraw1Dvs2DEfficiency   ->  SaveAs( kPlotDirectory + TString("/Full/hCompareEfficiency_1D_2D.pdf") );
        cDraw1Dvs2DSignalLoss   ->  SaveAs( kPlotDirectory + TString("/Full/hCompareSignalLoss_1D_2D.pdf") );
        //
        //  --- --- Specific QA Plots
        //
        outFile_Yield   ->  Close();
    }
    //  --- --- MULTIPLICITY ANALYSIS
    if ( kDoMultiplicity ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/Full/")).Data()));
        TFile*  outFile_Yield   =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Multiplicity")+kFolder).Data()),"recreate");
        auto    kPlotDirectory  =   TString(Form(kAnalysis_PreProc_Dir,(TString("Multiplicity")+kFolder).Data()))+TString("/Plots/");
        //
        // --- Saving to File
        for ( auto kSave    : h1D_Nrec )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_Fin )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Nrec_2Db )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen )                                            kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_2Db )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_Fin )                                        kSave   ->  Write();
        for ( auto kSave    : h1D_Ngen_VTX10 )                                      kSave   ->  Write();
        for ( auto kSave    : h1D_Eff )                                             kSave   ->  Write();
        for ( auto kSave    : h1D_SLC )                                             kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes )      for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_TruInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_RecInvMass_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_InvMassRes_2Db )  for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        for ( auto kVecSave : h1D_Mult_Nrec )       for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_Nrec_2Db )   for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_Nrec_Fin )   for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_Ngen )       for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_Ngen_2Db )   for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_Ngen_Fin )   for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_Eff )        for ( auto kSave : kVecSave )   kSave   ->  Write();
        for ( auto kVecSave : h1D_Mult_SLC )        for ( auto kSave : kVecSave )   kSave   ->  Write();
        //
        //  --- Printing to Plots
        //
        //  --- --- General QA Plots
        cDrawPeriodEfficiency   ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficiencies1D.pdf") );
        cDrawPeriodSignalLoss   ->  SaveAs( kPlotDirectory + TString("/Full/hFullSignalLosses1D.pdf") );
        cDraw1Dvs2DEfficiency   ->  SaveAs( kPlotDirectory + TString("/Full/hCompareEfficiency_1D_2D.pdf") );
        cDraw1Dvs2DSignalLoss   ->  SaveAs( kPlotDirectory + TString("/Full/hCompareSignalLoss_1D_2D.pdf") );
        //
        //  --- --- Specific QA Plots
        //  --- --- --- Extract Full Period Efficiencies
        std::vector<TH1F*>  kMultEfficiencies;
        auto iTer = 0;
        for ( auto kCurrent_Multiplicity : h1D_Mult_Eff.at(0) ) {
            if ( iTer == 0 )    kMultEfficiencies.push_back( (TH1F*)( kCurrent_Multiplicity->Clone(Form("V0M [0.00-100.]")) ) );
            else                kMultEfficiencies.push_back( (TH1F*)( kCurrent_Multiplicity->Clone(Form("V0M [%3.2f-%3.2f]",   fArrMult[iTer-1], fArrMult[iTer] )) ) );
            iTer++;
        }
        //uAddSumHistogram( kMultEfficiencies,   "V0M [00.00-100.00]", kPeriodWeights );
        //
        std::vector<TH1F*>  kMultSignalLosses;
        iTer = 0;
        for ( auto kCurrent_SignalLoss : h1D_Mult_SLC.at(0) )   { kMultSignalLosses.push_back( (TH1F*)( kCurrent_SignalLoss->Clone(Form("V0M [%3.2f-%3.2f]",   fArrMult[iTer-1], fArrMult[iTer] )) ) ); iTer; }
        uAddSumHistogram( kMultSignalLosses,   "V0M [00.00-100.00]", kPeriodWeights );
        //
        auto    cDrawMultipEfficiency    =   uPlotEfficiencies( kMultEfficiencies, {}, "cDrawMultipEfficiency" );
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        cDrawMultipEfficiency        ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficienciesMult.pdf") );
        delete cDrawMultipEfficiency;
        //
        /*
        auto    cDrawMultipSignalLoss    =   uPlotEfficiencies( kMultEfficiencies, {}, "cDrawMultipEfficiency" );
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        cDrawMultipEfficiency        ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficienciesMult.pdf") );
        delete cDrawMultipSignalLoss;
         */
        //
        /*  !TODO: Add in the Analysis Task information about multiplicity when the event is discarded
        auto    cDrawMultipSignalLoss    =   uPlotEfficiencies( kMultSignalLosses, {}, "cDrawMultipSignalLoss", true );
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        cDrawMultipSignalLoss        ->  SaveAs( kPlotDirectory + TString("/Full/hFullSignalLossesMult.pdf") );
        delete cDrawMultipSignalLoss;
         */
        //
        /*
        cDrawEff        =   uEfficiencyCompare_1D_2D( h1D_Nrec_2Db.at(0), h1D_Ngen_2Db.at(0), h2D_Nrec.at(0), h2D_Ngen.at(0)  );
        cDrawEff        ->  SaveAs( kPlotDirectory + TString("/Full/hCompareEfficiency_1D_2D.pdf") );
        delete cDrawEff;
        //
        cDrawEff    =   uPlotEfficiencies( kMultEfficiencies );
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.19, 0.83,"ALICE Performance");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        if ( is_pp_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pp #sqrt{#it{s}}= %.2f TeV",  kEnergy ) );
        if ( is_pb_anl ) uLatex->DrawLatexNDC(0.19, 0.77, Form( "pPb #sqrt{#it{s}}= %.2f TeV", kEnergy ) );
        uLatex->DrawLatexNDC(0.19, 0.71,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
        //
        cDrawEff        ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficienciesMult.pdf") );
        delete cDrawEff;
        //
        gROOT           ->  SetBatch(kFALSE);
        //
         */
        outFile_Yield   ->  Close();
    }
    //  --- --- CORRELATION ANALYSIS
    if ( kDoCorrelation ) {
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data())));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/1D/")).Data()));
        gROOT                   ->  ProcessLine(Form(".! mkdir -p %s",(TString(Form(kAnalysis_PreProc_Dir,(TString("Correlation")+kFolder).Data()))+TString("/Plots/2D/")).Data()));
        TFile*  outFile_Yield   =   new TFile   (Form(kAnalysis_MCTruthHist,(TString("Correlation")+kFolder).Data()),"recreate");
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
        //  --- Printing to Plots
        //
        //  --- --- General QA Plots
        cDrawPeriodEfficiency   ->  SaveAs( kPlotDirectory + TString("/Full/hFullEfficiencies1D.pdf") );
        cDrawPeriodSignalLoss   ->  SaveAs( kPlotDirectory + TString("/Full/hFullSignalLosses1D.pdf") );
        //
        outFile_Yield   ->  Close();
    }
    //
    gROOT   ->  SetBatch(kFALSE);
    delete cDrawPeriodEfficiency;
    delete cDrawPeriodSignalLoss;
    for ( auto kFile : insFile_MntC_YLD ) kFile->Close();
    for ( auto kFile : insFile_MntC_MLT ) kFile->Close();
}
