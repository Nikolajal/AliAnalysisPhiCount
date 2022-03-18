#include "../inc/AliAnalysisPhiPair.h"
#include "./PreProcessing/PP_Data.cxx"
#include "./PreProcessing/PP_MC.cxx"
#include "./Analysis/AN_SigExtraction.cxx"
#include "./Analysis/AN_SigCorrections.cxx"
#include "./Analysis/AN_ShowPlots.cxx"
#include "./Analysis/AN_dQuantities_FinalPlots.cxx"
// !TODO: All Set!

void
runAnalysis
 ( TString fFolderMC, TString fFolderDT, TString kDataFile, std::vector<pair<TString,TString>> kDataSet, TString fOption, Int_t nEventsCut, TString kFolder, TString kType ) {
    //
    //  --- Pre-Processing
    //PP_Data                     ( fFolderDT + kDataFile, fOption, nEventsCut, kFolder );
    //PP_MC                       ( fFolderDT, fFolderMC, kDataSet, fOption, nEventsCut, kFolder, true, kType );
    //PP_Resl                     ( fOption, kFolder, true );
    //
    //  --- Signal Processing
    //AN_SigExtraction            ( fOption, kFolder, true );
    AN_SigCorrections           ( fOption, kFolder, true );
    //
    //  --- Final Plots for Show
    //AN_ShowPlots                ( fOption, kFolder, true );
    //  --- Results Re-Fining
    //AN_dQuantities_FinalPlots   ( );
}
//
//  --- FULL COMMANDS LIST
//
//  --- --- pp 7TeV
//.x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__7TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__7TeV/","/LHC10_STD.root",{{"LHC10b","LHC14j4b"},{"LHC10c","LHC14j4c"},{"LHC10d","LHC14j4d"},{"LHC10e","LHC14j4e"}},"yield",-1,"_p_p__7TeV","STD")
//.x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__7TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__7TeV/","/HighPileup/HighPileup_STD.root",{{"HighPileup","HighPileup"}},"yield",-1,"_p_p__7TeV/PileUpTest/HighRate/")
//.x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__7TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__7TeV/","/Low_Pileup/Low_Pileup_STD.root",{{"Low_Pileup","Low_Pileup"}},"yield",-1,"_p_p__7TeV/PileUpTest/Low_Rate/")
//.x src/Systematics/SY_Production.cxx("PID","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__7TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__7TeV/","/LHC10",{{"LHC10b","LHC14j4b"},{"LHC10c","LHC14j4c"},{"LHC10d","LHC14j4d"},{"LHC10e","LHC14j4e"}},"yield",-1.,"_p_p__7TeV/Systematics/PID/")
//.x src/Systematics/SY_Production.cxx("TRK","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__7TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__7TeV/","/LHC10",{{"LHC10b","LHC14j4b"},{"LHC10c","LHC14j4c"},{"LHC10d","LHC14j4d"},{"LHC10e","LHC14j4e"}},"yield",-1.,"_p_p__7TeV/Systematics/TRK/")
//
//  --- --- pp 5TeV
//.x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__5TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__5TeV/","LHC15n17p_STD.root",{ {"LHC15n","LHC18j3"}, {"LHC17p_FAST","LHC18j2_FAST"}, {"LHC17p_CENT_woSDD","LHC18j2_CENT_woSDD"} }, "all", -1, "_p_p__5TeV", "STD")
//.x src/Systematics/SY_Production.cxx("PID","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__5TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__5TeV/","LHC15n17p",{ {"LHC15n","LHC17e2"}, {"LHC17p_FAST","LHC18j2_FAST"}, {"LHC17p_CENT_woSDD","LHC18j2_CENT_woSDD"} },"all",-1.,"_p_p__5TeV/Systematics/PID/")
//.x src/Systematics/SY_Production.cxx("TRK","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__5TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__5TeV/","LHC15n17p",{ {"LHC15n","LHC17e2"}, {"LHC17p_FAST","LHC18j2_FAST"}, {"LHC17p_CENT_woSDD","LHC18j2_CENT_woSDD"} },"all",-1.,"_p_p__5TeV/Systematics/TRK/")
//
//  --- --- pp 13TeV
// .x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_2/Dataset/_Sim/p_p__13TeV/","/Volumes/NRUBINI_DATASTASH_2/Dataset/_Data/p_p__13TeV/","LHC13TeV_STD.root", { {"LHC16d","LHC16d_LHC17f6"}, {"LHC16e","LHC16e_LHC17f9"}, {"LHC16g","LHC16g_LHC17d17"}, {"LHC16h","LHC16h_LHC17f5"}, {"LHC16i","LHC16i_LHC17d3"}, {"LHC16j","LHC16j_LHC17e5"}, {"LHC16k","LHC16k_LHC18f1"}, {"LHC16l","LHC16l_LHC18d8"}, {"LHC16o","LHC16o_LHC17d16"}, {"LHC16p","LHC16p_LHC17d18"}, {"LHC17c","LHC17c_LHC18d3"}, {"LHC17e","LHC17e_LHC17h1"}, {"LHC17f","LHC17i_LHC17k4"}, {"LHC17g","LHC17i_LHC17k4"}, {"LHC17i","LHC17i_LHC17k4"}, {"LHC17j","LHC17j_LHC17h11"}, {"LHC17k","LHC17k_LHC18c13"}, {"LHC17l","LHC17l_LHC18a8"}, {"LHC17m","LHC17m_LHC17l5"}, {"LHC17o","LHC17m_LHC17l5"}, {"LHC17r","LHC17r_LHC18a1"}, {"LHC18b","LHC18b_LHC18g4"}, {"LHC18d","LHC18d_LHC18g5"}, {"LHC18e","LHC18e_LHC18g6"}, {"LHC18f","LHC18f_LHC18h2"}, {"LHC18g","LHC18ghijk_LHC18h4"}, {"LHC18h","LHC18ghijk_LHC18h4"}, {"LHC18i","LHC18ghijk_LHC18h4"}, {"LHC18k","LHC18ghijk_LHC18h4"}, {"LHC18l","LHC18l_LHC18j1"}, {"LHC18m","LHC18m_LHC18j4"}, {"LHC18n","LHC18n_LHC18k1"}, {"LHC18o","LHC18o_LHC18k2"}, {"LHC18p","LHC18p_LHC18k3" } }, "correlation", -1, "_p_p__13TeV", "STD")
//
//  --- --- pPb 5TeV
//.x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_Pb_5TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_Pb_5TeV/","LHC16q_STD.root",{ {"LHC16q_CENT_woSDD","LHC17f2b_CENT_woSDD"}, {"LHC16q_FAST","LHC17f2b_FAST"} }, "all", -1, "_p_Pb_5TeV", "STD")
//
//  --- --- TEST
//
// .x src/runAnalysis.C("/Volumes/NRUBINI_DATASTASH_1/Dataset/_Sim/p_p__5TeV/","/Volumes/NRUBINI_DATASTASH_1/Dataset/_Data/p_p__5TeV/","LHC15n/LHC15n_STD.root",{ {"LHC15n","LHC18j3"} }, "all", -1, "_p_p__5TeV_test", "STD")
