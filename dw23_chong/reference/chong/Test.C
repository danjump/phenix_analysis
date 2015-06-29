#include "HbgExtpolTools.C"

void Test(int ARM = 1, int CHG = 1, bool Show = true, bool Print = false)
{
	TFile *F = TFile::Open("./stobg_high.root");
	TTree *T = (TTree*)F->Get("stobg_data510");

	TF1 *WnessPol4 = F1WnessPol4(T, ARM, CHG, 10, Show);

	TH1F *H1dw23Slice = H1GetDw23Slice(T, ARM, CHG, 30, 0.1, 0.1, Show);

	TH1F *H1dw23DAT = H1GetSignalSample(0, "dw23", ARM, CHG, 32, -0.1, 0.1);
	TH1F *H1dw23SIG = H1GetSignalSample(1, "dw23", ARM, CHG, 32, -0.1, 0.1);
	TH1F *H1dw23MBG = H1GetSignalSample(2, "dw23", ARM, CHG, 32, -0.1, 0.1);

	TH2F *H2Had[2][2];
        //fill with dw23 vs wness from hadronic background simulation
	GetHadMCSummed(H2Had, 50, 30, Show);

	TF1 *F1Lo = new TF1(Form("F1Lo_%s", H1dw23MBG->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
	TF1 *F1Hi = new TF1(Form("F1Hi_%s", H1dw23MBG->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
	DecomposeDw23(H1dw23MBG, F1Lo, F1Hi, true,  Show);
	//DecomposeDw23(H1dw23MBG, F1Lo, F1Hi, false, Show);
	//DecomposeDw23Manual(H1dw23MBG, F1Lo, F1Hi, -5.e-2, 5.e-2);

	float HbgDw23SignalLo = GetHbgDw23SignalLoWidth(ARM, CHG, Show);
	cout <<HbgDw23SignalLo <<endl;

	TCanvas *cNL = new TCanvas("cNL", "", 1);
	TLegend *lNL = new TLegend(0.15, 0.15, 0.35, 0.5);
	for (int a=1; a<9; a++)
	{
		TF1 *F1 = new TF1(Form("F1Pol1NL_%i", a), Pol1NLFcn, 0.1, 1.0, 3);
		F1->SetParameter(0,  0.135);
		F1->SetParameter(1, -0.1);
		F1->SetParameter(2,  0.5*a);
		F1->SetLineColor(a+1);
		lNL->AddEntry(F1, Form("Drop p = %2.1f", F1->GetParameter(2)), "l");
		cNL->cd();
		F1->Draw(a==1?"":"same");
		if (a==8) lNL->Draw();
	}

	float FitSeed[4] = {0};
	GetHbgDw23Pattern(T, ARM, CHG, FitSeed, Show, Print);

	return;
}
