#include "../../inc/AliAnalysisPhiPair.h"

void convert_EXT_sys(TString fOption = "Yield mult", TString kFolder = "_p_p__5TeV")
{
    //!
    //! Set-Up
    //! --- Option chosing
    if (!fChooseOption(fOption))
        return;
    //! --- Style
    SetStyle();
    //! --- Binning
    fSetAllBins();
    //!
    //! Yield analysis
    if (kDoYield)
    {
        //!
        //! Building output directory
        gROOT->ProcessLine(Form(".! mkdir -p %s", Form(kAnalysis_Systemt_Dir + TString("EXT"), (TString("Yield") + kFolder).Data())));
        //!
        //! Input file
        TFile *input_file = new TFile(Form(kASigExtp_FitCheckRst, (TString("Yield") + kFolder).Data()));
        //  Recover sys
        auto yl_ext_syst = (TH1F *)(input_file->Get("hxD_YL_Extrapol_Syst"));
        auto yl_full_val = (TH1F *)(input_file->Get("hXD_Nfqs_stat"));
        auto yl_full_pt_ = (TH1F *)(input_file->Get("h2D_MeanPT_stat"));
        auto yl_full_1d_ = (TH1F *)(input_file->Get("h1D_Nres_stat"));
        auto yl_full_2d_ = (TH1F *)(input_file->Get("h2D_Nres_stat"));
        uScale(yl_full_pt_,0);
        uScale(yl_full_1d_,0);
        uScale(yl_full_2d_,0);
        yl_full_pt_->SetName("kPTSystematics");
        yl_full_1d_->SetName("k1DSystematics");
        yl_full_2d_->SetName("k2DSystematics");
        //  Build final plotter
        TH1F *kRatioSystematics = new TH1F("kRatioSystematics", "kRatioSystematics", 6, 0.5, 6.5);
        uSetHisto(kRatioSystematics, "FNL");
        //
        auto d1_val = yl_full_val->GetBinContent(1);
        auto d2_val = yl_full_val->GetBinContent(2);
        //
        auto d1_err = yl_ext_syst->GetBinError(1) / d1_val;
        auto d2_err_tmp = 0.;
        auto d2_val_temp = 0.;
        for (int i = 3; i < 12; i++)
        {
            d2_err_tmp += yl_ext_syst->GetBinError(i) * yl_ext_syst->GetBinError(i);
            d2_val += yl_ext_syst->GetBinContent(i);
        }
        auto d2_err = sqrt(d2_err_tmp) / d2_val;
        //
        kRatioSystematics->SetBinContent(1, d1_err);
        kRatioSystematics->SetBinContent(2, d2_err);
        kRatioSystematics->SetBinContent(3, sqrt(d1_err*d1_err+d2_err*d2_err));
        kRatioSystematics->SetBinContent(4, sqrt(4*d1_err*d1_err+d2_err*d2_err));
        kRatioSystematics->SetBinContent(5, fSigmaPhiError(d1_val,d2_val,d1_val*d1_err,d2_val*d2_err) / fSigmaPhiValue(d1_val,d2_val) );
        kRatioSystematics->SetBinContent(6, fGammaPhiError(d1_val,d2_val,d1_val*d1_err,d2_val*d2_err) / fGammaPhiValue(d1_val,d2_val) );
        //
        //! Output file
        TFile *output_file = new TFile(Form(kAnalysis_Systemt_Dir + TString("EXT/FullSystematics.root"), (TString("Yield") + kFolder).Data()), "RECREATE");
        yl_ext_syst->Write();
        kRatioSystematics->Write();
        yl_full_pt_->Write();
        yl_full_1d_->Write();
        yl_full_2d_->Write();
    }
}

void delete_me()
{

    convert_EXT_sys();
}