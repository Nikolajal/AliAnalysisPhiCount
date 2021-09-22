// File for 1-Dimensional Analysis:
// !TODO: All Set!
#include "../../inc/AliAnalysisPhiPair.h"
#include "GeneralAnalysis.cxx"
#include "RooMsgService.h"

void
SystAnalysis
( TString fType = "PID", Int_t iTot = 0, Int_t iInit = 0 )    {
    if ( fType.Contains("SEX") ) {
        //
        //  YIELD Efficiency
        TFile      *insFile_EF_Yield    =   new TFile   ( Form(kAnalysis_MCTruthHist,"yield/Systematics/Standard/") );
        TFile      *insFile_DT_Yield    =   new TFile   ( Form(kASigExtr_FitCheckRst,"yield/Systematics/Standard/") );
        //
        TH1F       *h1DStandard;
        hName                           =   Form("hRAW_1D");
        h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
        //
        std::vector<TH1F*>  f1DVars;
        std::vector<TFile*> f1DInFiles;
        for ( Int_t iTer = 0; iTer < nOptions; iTer++ )    {
            hName           =   Form("hRAW_1D");
            f1DInFiles.push_back(new TFile   (Form("%s/ExtractionCheck/%s/1D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),sOptions.at(iTer).Data(),sOptions.at(iTer).Data())));
            auto    fTarget =   (TH1F*)((f1DInFiles.at(iTer))->Get(hName));
            fTarget->SetName(sOptions.at(iTer).Data());
            f1DVars.push_back(fTarget);
        }
        //
        GeneralAnalysis(h1DStandard,f1DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield")));
        //
        TH2F       *h2DStandard;
        hName                           =   Form("hRAW_2D");
        h2DStandard                     =   (TH2F*)((insFile_DT_Yield)->Get(hName));
        //
        std::vector<TH2F*>  f2DVars;
        std::vector<TFile*> f2DInFiles;
        for ( Int_t iTer = 0; iTer < nOption2; iTer++ )    {
            hName           =   Form("hRAW_2D");
            f2DInFiles.push_back(new TFile   (Form("%s/ExtractionCheck/%s/2D/FitResults_%s.root",Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/"),sOption2.at(iTer).Data(),sOption2.at(iTer).Data())));
            auto    fTarget =   (TH2F*)((f2DInFiles.at(iTer))->Get(hName));
            fTarget->SetName(sOption2.at(iTer).Data());
            f2DVars.push_back(fTarget);
        }
        //
        GeneralAnalysis(h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/")));
        //
        uEvaluateRatioError(h1DStandard,f1DVars,h2DStandard,f2DVars,TString(Form(kAnalysis_SgExSys_Dir,"yield/Systematics/Standard/")),(TH1F*)(insFile_EF_Yield->Get("hEFF_1D")),(TH2F*)(insFile_EF_Yield->Get("hEFF_2D_fr_1D")));
        //
        for ( auto CurrentFile : f1DInFiles ) CurrentFile->Close();
        for ( auto CurrentFile : f2DInFiles ) CurrentFile->Close();
        //
        insFile_DT_Yield->Close();
        insFile_EF_Yield->Close();
        //
        GeneralAnalysis();
        //
    } else {
        //
        //  YIELD Efficiency
        TFile      *insFile_DT_Yield    =   new TFile   ( Form(kASigExtp_FitCheckRst,"Yield/Systematics/Standard/") );
        //
        TH1F       *h1DStandard;
        hName                           =   Form("hRES_1D_Stat");
        h1DStandard                     =   (TH1F*)((insFile_DT_Yield)->Get(hName));
        //
        fSetAllBins();
        //
        auto iVec = 0;
        std::vector<TH1F*>  f1DVars;
        std::vector<TFile*> f1DInFiles;
        for ( Int_t iTer = iInit; iTer <  iTot; iTer++ )    {
            hName           =   Form("hRES_1D_Stat");
            f1DInFiles.push_back(new TFile  (Form("%s/%s/%s/SignalExtrapolation/FitResults.root",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data(),Form("%s_%i",fType.Data(),iTer+1))));
            auto    fTarget =   (TH1F*)((f1DInFiles.at(iVec))->Get(hName));
            fTarget->SetName(Form("%s%i",fType.Data(),iTer+1));
            f1DVars.push_back(fTarget);
            iVec++;
        }
        //
        GeneralAnalysis(h1DStandard,f1DVars,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        //
        TH2F       *h2DStandard;
        hName                           =   Form("hRES_2D_Stat");
        h2DStandard                     =   new TH2F    ("h2DStandard","h2DStandard",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
        for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
            auto fInput = (TH1F*)(insFile_DT_Yield->Get(Form("hRES_2D_Cond1_Stat_%i",iPT2D)));
            for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                h2DStandard->SetBinContent  (iPT2D+1,jPT2D+1,fInput->GetBinContent(jPT2D+1));
                h2DStandard->SetBinError    (iPT2D+1,jPT2D+1,fInput->GetBinError(jPT2D+1));
            }
        }
        //
        iVec = 0;
        std::vector<TH2F*>  f2DVars;
        std::vector<TFile*> f2DInFiles;
        for ( Int_t iTer = iInit; iTer < iTot; iTer++ )    {
            hName           =   Form("hRES_2D_Stat");
            f2DInFiles.push_back(new TFile   (Form("%s/%s/%s/SignalExtrapolation/FitResults.root",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data(),Form("%s_%i",fType.Data(),iTer+1))));
            auto    fTarget =   new TH2F    (Form("%s%i",fType.Data(),iTer+1),Form("%s%i",fType.Data(),iTer),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
            for ( Int_t iPT2D = 0; iPT2D < nBinPT2D; iPT2D++ ) {
                auto fInput = (TH1F*)(f2DInFiles.at(iVec)->Get(Form("hRES_2D_Cond1_Stat_%i",iPT2D)));
                for ( Int_t jPT2D = 0; jPT2D < nBinPT2D; jPT2D++ ) {
                    fTarget->SetBinContent  (iPT2D+1,jPT2D+1,fInput->GetBinContent(jPT2D+1));
                    fTarget->SetBinError    (iPT2D+1,jPT2D+1,fInput->GetBinError(jPT2D+1));
                }
            }
            f2DVars.push_back(fTarget);
            iVec++;
        }
        //
        GeneralAnalysis(h2DStandard,f2DVars,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        //
        uEvaluateRatioError(h1DStandard,f1DVars,h2DStandard,f2DVars,TString(Form("%s/%s/",Form(kAnalysis_Systemt_Dir,"yield"),fType.Data())));
        //
        for ( auto CurrentFile : f1DInFiles ) CurrentFile->Close();
        for ( auto CurrentFile : f2DInFiles ) CurrentFile->Close();
        //
        insFile_DT_Yield->Close();
        //
        GeneralAnalysis();
        //
    }
}
