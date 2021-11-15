// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../../inc/AliAnalysisPhiPair.h"
#include "../GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SystAnalysis_SEX
( TString fOption = "Yield" )    {
    //  Setting binnings
    fSetAllBins();
    //  Option selection
    fChooseOption(fOption);
    if ( kDoYield ) {
        //
        //  YIELD Efficiency
        TFile      *insFile_EF_Yield    =   new TFile   ( Form(kAnalysis_MCTruthHist,"yield") ); // /Systematics/Standard/") );
        TFile      *insFile_DT_Yield    =   new TFile   ( Form(kASigExtr_FitCheckRst,"yield") ); // /Systematics/Standard/") );
        //
        TH1F       *hEFF_1D;
        hName                           =   Form("hEFF_1D");
        hEFF_1D                         =   (TH1F*)((insFile_EF_Yield)->Get(hName));
        //
        TH1F       *hEFF_2D_fr_1D;
        hName                           =   Form("hEFF_2D_fr_1D");
        hEFF_2D_fr_1D                   =   (TH1F*)((insFile_EF_Yield)->Get(hName));
        //
        TH1F       *h1DStandard;
        hName                           =   Form("hRAW_1D");
        h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
        //
        hName                           =   Form("fQC_Event_Enum_FLL");
        auto        kN_Trg          =   ((TH1F*)((insFile_DT_Yield)->Get(hName)))->GetBinContent(kEventCount::kTrigger);
        auto        kN_Vtx          =   ((TH1F*)((insFile_DT_Yield)->Get(hName)))->GetBinContent(kEventCount::kVertex);
        auto        kN_MB           =   ((TH1F*)((insFile_DT_Yield)->Get(hName)))->GetBinContent(kEventCount::kVertex10);
        Double_t    f1DCorrection   =   (1./kBR)        *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        Double_t    f2DCorrection   =   (1./(kBR*kBR))  *(1./kN_MB) *(kTriggerEff/1.)   *(kN_Vtx/kN_Trg);
        //
        h1DStandard->Scale(f1DCorrection);
        h1DStandard->Divide(h1DStandard,hEFF_1D,1.,1.,"b");
        h1DStandard->SetName("hSystematicUncertainty_SEX");
        //
        std::vector<TH1F*>  f1DVars;
        std::vector<TFile*> f1DInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
            hName           =   Form("hRAW_1D");
            f1DInFiles.insert(f1DInFiles.begin() , new TFile   (Form("%s/ExtractionCheck/%s/1D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard/*/),kCurrent_Syst.Data(),kCurrent_Syst.Data())));
            auto    fTarget =   (TH1F*)((f1DInFiles.at(0))->Get(hName));
            fTarget->Scale(f1DCorrection);
            fTarget->Divide(fTarget,hEFF_1D,1.,1.,"b");
            fTarget->SetName(kCurrent_Syst.Data());
            f1DVars.push_back(fTarget);
        }
        //
        GeneralAnalysis(h1DStandard,f1DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard/*/)));
        //
        TH2F       *h2DStandard;
        hName                           =   Form("hRAW_2D");
        h2DStandard                     =   (TH2F*)((insFile_DT_Yield)->Get(hName));
        //
        h2DStandard->Scale(f2DCorrection);
        h2DStandard->Divide(h2DStandard,hEFF_2D_fr_1D,1.,1.,"b");
        h2DStandard->SetName("hSystematicUncertainty_SEX");
        //
        std::vector<TH2F*>  f2DVars;
        std::vector<TFile*> f2DInFiles;
        for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
            hName           =   Form("hRAW_2D");
            f2DInFiles.insert(f2DInFiles.begin() , new TFile   (Form("%s/ExtractionCheck/%s/2D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard/*/),kCurrent_Syst.Data(),kCurrent_Syst.Data())));
            auto    fTarget =   (TH2F*)((f2DInFiles.at(0))->Get(hName));
            fTarget->Scale(f2DCorrection);
            fTarget->Divide(fTarget,hEFF_2D_fr_1D,1.,1.,"b");
            fTarget->SetName(kCurrent_Syst.Data());
            f2DVars.push_back(fTarget);
        }
        //
        GeneralAnalysis(h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard/*/)));
        //
        //GeneralAnalysis(h1DStandard,f1DVars,h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard/*/)));
        //
        for ( auto CurrentFile : f1DInFiles ) CurrentFile->Close();
        for ( auto CurrentFile : f2DInFiles ) CurrentFile->Close();
        //
        insFile_DT_Yield->Close();
        insFile_EF_Yield->Close();
        //
        //GeneralAnalysis();
        //
    }
    if ( kDoMultiplicity ) {
        //
        //  YIELD Efficiency
        TFile*      insFile_EF_Yield        =   new TFile   ( Form(kAnalysis_MCTruthHist,"multiplicity") ); // /Systematics/Standard/") );
        TFile*      insFile_DT_Yield        =   new TFile   ( Form(kASigExtr_FitCheckRst,"multiplicity") ); // /Systematics/Standard/") );
        //
        //  Event Count for Normalisation
        TH1D*       hUtilEventMultiplicity  =   (TH1D*)( insFile_DT_Yield -> Get ( "fQC_Event_Enum_V0M" ) );
        //
        TH1F       *hEFF_1D;
        hName                           =   Form("hEFF_1D");
        hEFF_1D                         =   (TH1F*)((insFile_EF_Yield)->Get(hName));
        //
        TH1F       *hEFF_2D_fr_1D;
        hName                           =   Form("hEFF_2D_fr_1D");
        hEFF_2D_fr_1D                   =   (TH1F*)((insFile_EF_Yield)->Get(hName));
        //
        for ( Int_t iMult = 0; iMult <= nBinMult; iMult++ )    {
            TH1F       *h1DStandard;
            hName                           =   Form("hRAW_1D_in_Mlt_%i",iMult);
            h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
            //
            auto    kNormalisation1D    =   (1.)/(fEvaluateINELgt0(iMult-1,hUtilEventMultiplicity) * kBR );
            auto    kNormalisation2D    =   (1.)/(fEvaluateINELgt0(iMult-1,hUtilEventMultiplicity) * kBR * kBR );
            //
            h1DStandard->Scale(kNormalisation1D);
            h1DStandard->Divide(h1DStandard,hEFF_1D);
            h1DStandard->SetName("hSystematicUncertainty_SEX");
            //
            std::vector<TH1F*>  f1DVars;
            std::vector<TFile*> f1DInFiles;
            for ( auto kCurrent_Syst : kSyst_SEX_1D_Options ) {
                hName           =   Form("hRAW_1D_MLT_%i",iMult);
                f1DInFiles.insert(f1DInFiles.begin() , new TFile   (Form("%s/ExtractionCheck/%s/1D/MLT_%i/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"multiplicity"/*yield/Systematics/Standard/*/),kCurrent_Syst.Data(),iMult,kCurrent_Syst.Data())));
                auto    fTarget =   (TH1F*)((f1DInFiles.at(0))->Get(hName));
                fTarget->Scale(kNormalisation1D);
                fTarget->Divide(fTarget,hEFF_1D);
                fTarget->SetName(kCurrent_Syst.Data());
                f1DVars.push_back(fTarget);
            }
            //
            GeneralAnalysis(h1DStandard,f1DVars,TString(Form(kAnalysis_SgExSys_Dir,"multiplicity"/*yield/Systematics/Standard/*/)),Form("MLT_%i",iMult));
            //
            TH2F       *h2DStandard;
            hName                           =   Form("hRAW_2D_in_Mlt_%i",iMult);
            h2DStandard                     =   (TH2F*)((insFile_DT_Yield)->Get(hName));
            //
            h2DStandard->Scale(kNormalisation2D);
            h2DStandard->Divide(h2DStandard,hEFF_2D_fr_1D);
            h2DStandard->SetName("hSystematicUncertainty_SEX");
            //
            std::vector<TH2F*>  f2DVars;
            std::vector<TFile*> f2DInFiles;
            for ( auto kCurrent_Syst : kSyst_SEX_2D_Options ) {
                hName           =   Form("hRAW_2D_MLT_%i",iMult);
                f2DInFiles.insert(f2DInFiles.begin() , new TFile   (Form("%s/ExtractionCheck/%s/2D/MLT_%i/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"multiplicity"),kCurrent_Syst.Data(),iMult,kCurrent_Syst.Data())));
                auto    fTarget =   (TH2F*)((f2DInFiles.at(0))->Get(hName));
                fTarget->Scale(kNormalisation2D);
                fTarget->Divide(fTarget,hEFF_2D_fr_1D);
                fTarget->SetName(kCurrent_Syst.Data());
                f2DVars.push_back(fTarget);
            }
            //
            GeneralAnalysis(h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"multiplicity"/*yield/Systematics/Standard/*/)),Form("MLT_%i",iMult));
            //
            //GeneralAnalysis(h1DStandard,f1DVars,h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield"/*yield/Systematics/Standard/*/)));
            //
            for ( auto CurrentFile : f1DInFiles ) CurrentFile->Close();
            for ( auto CurrentFile : f2DInFiles ) CurrentFile->Close();
        }
        //
        insFile_DT_Yield->Close();
        insFile_EF_Yield->Close();
        //
        //GeneralAnalysis();
        //
    }
}
