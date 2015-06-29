#include "HbgExtpolTools.C"

#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

RooRealVar dw23("dw23", "dw23", -0.3, 0.3);

//2D Fit function
//---------------------------------
double Gaus2D(double *x, double *p)
{
	double X = x[0];	//Wness
	double Y = x[1];	//dw23

	double INT = p[0] + p[1]*X + p[2]*pow(X, 2) + p[3]*pow(X, 3) + p[4]*pow(X, 4);
	double AMP = p[5];

	double Center1 = p[6];
	double Center2 = p[7];
	double Width1 = p[8]  + p[9] *X + NLFcn(X, p[10]);
	double Width2 = p[11] + p[12]*(X - 0.96);
	double Frac = p[13] + p[14]*X;

	double R = ((INT*AMP) / (Width1*(1-Frac) + Width2*Frac)) *
				(TMath::Gaus(Y, Center1, Width1) * (1-Frac) + TMath::Gaus(Y, Center1-Center2, Width2) * Frac);

	return R;
};

//2D fit + extrapol, then return RooGenericPdf
//--------------------------------------------------------------------------------------------------------------------
RooGenericPdf *HbgExtpol(TTree *T, int ARM, int CHG, int BinWness, int BinDw23, bool Show = false, bool Print = false)
{
	//Get pol4 function for Wness axis
	//-----------------------------------------------------
	TF1 *F1Wness = F1WnessPol4(T, ARM, CHG, BinWness, Show);

	//Get Width1 (high + narrow) parameters (pol1 from data, NL drop from mcHad) and Width2 (low + wide) slope
	//--------------------------------------------------------------------------------------------------------
	float FitSeed[4] = {0};
	GetHbgDw23Pattern(T, ARM, CHG, FitSeed, Show, Print);

	//Get Width2's intercept in signal region (Wness > 0.92)
	//----------------------------------------------------------------
	const float Width2IntVal = GetHbgDw23SignalLoWidth(ARM, CHG, Show);

	//Perform 2D fit: Width1's NL and Width2's Int is fixed
	//----------------------------------------------------------------
	TCut CutWnessBg = "Wness > 0.1 && Wness < 0.9";
	TCut CutFin		= CutAC[ARM][CHG] + CutEta + CutFvtxCone + CutWnessBg;

	TH2F *H2 = new TH2F(
			Form("H2_%i%i", ARM, CHG), Form("%s;Wness;dw23", AC[ARM][CHG]), BinWness, 0, 1, BinDw23, -0.3, 0.3
			);
	T->Project(H2->GetName(), "dw23:Wness", CutFin);

	enum
	{
		INT0, INT1, INT2, INT3, INT4, AMP,
		Center1, Center2, Width1Int, Width1Slo, Width1NL, Width2Int, Width2Slo, FracInt, FracSlo
	};
	TF2 *F2 = new TF2(Form("F2_%i%i", ARM, CHG), Gaus2D, 0.1, 0.9, -0.3, 0.3, 15);

	F2->SetParNames("pol4_0", "pol4_1", "pol4_2", "pol4_3", "pol4_4", "AMP", "Center1", "Center2");
	F2->SetParName(Width1Int, "Width1_Int");
	F2->SetParName(Width1Slo, "Width1_Slo");
	F2->SetParName(Width1NL, "Width1_NL");
	F2->SetParName(Width2Int, "Width2_Int");
	F2->SetParName(Width2Slo, "Width2_Slo");
	F2->SetParName(FracInt, "Frac_Int");
	F2->SetParName(FracSlo, "Frac_Slo");

	F2->FixParameter(INT0, F1Wness->GetParameter(0));
	F2->FixParameter(INT1, F1Wness->GetParameter(1));
	F2->FixParameter(INT2, F1Wness->GetParameter(2));
	F2->FixParameter(INT3, F1Wness->GetParameter(3));
	F2->FixParameter(INT4, F1Wness->GetParameter(4));

	F2->SetParameter(AMP,		5.e-3);
	F2->SetParameter(Center1,	-0.02 + 0.04*CHG);
	F2->SetParameter(Center2,	0);
	F2->SetParameter(Width1Int, FitSeed[0]);
	F2->SetParameter(Width1Slo, FitSeed[1]);
	F2->FixParameter(Width1NL,  FitSeed[2]);		//Fixed
	F2->FixParameter(Width2Int, Width2IntVal);		//Fixed
	F2->SetParameter(Width2Slo, FitSeed[3]);
	F2->SetParameter(FracInt,	0.3);
	F2->SetParameter(FracSlo,	0);

	H2->Fit(F2->GetName(), "ER0");

	//Extrapolate into signal region. Set fine binning to increase statistics, then rebin to the similar binning
	//----------------------------------------------------------------------------------------------------------
	TH2F *H2Ep = new TH2F(Form("%s_Ep", H2->GetName()), "Extrapolated;Wness;dw23", 10*5, 0, 1, 90*10, -0.3, 0.3);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(H2Ep);
	H2Ep->Rebin2D(5, 10);

	//Get dw23 in target signal region
	//-----------------------------------------
	TH1F *H1EpSlice = (TH1F*)H2Ep->ProjectionY(
			Form("H1EpSlice_%i%i", ARM, CHG),
			H2Ep->GetXaxis()->FindBin(0.96 + 1.e-3),
			H2Ep->GetXaxis()->FindBin(0.97 - 1.e-3)
			);

	//Decompose extrapolated dw23, then compose RooGenericPdf
	//-----------------------------------------------------------------------------------------------------
	TF1 *F1Lo = new TF1(Form("F1LoEp_%s", H1EpSlice->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
	TF1 *F1Hi = new TF1(Form("F1HiEp_%s", H1EpSlice->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);

	if (ARM==0 && CHG==1) DecomposeDw23Manual(H1EpSlice, F1Lo, F1Hi, 1.6e-2, 5.1e-2);
	else DecomposeDw23(H1EpSlice, F1Lo, F1Hi, true);

	RooGenericPdf *GpDw23 = new RooGenericPdf(
			Form("GpDw23_%i%i", ARM, CHG),
			Form("%f*TMath::Gaus(dw23, %f, %f) + %f*TMath::Gaus(dw23, %f, %f)",
				F1Hi->GetParameter(0)/(F1Hi->GetParameter(0) + F1Lo->GetParameter(0)),
				F1Hi->GetParameter(1),
				F1Hi->GetParameter(2),
				F1Lo->GetParameter(0)/(F1Hi->GetParameter(0) + F1Lo->GetParameter(0)),
				F1Lo->GetParameter(1),
				F1Lo->GetParameter(2)
				),
			RooArgSet(dw23)
			);

	//if (Show)
	{
		TCanvas *c1 = new TCanvas(Form("c1HbgExtpol_%i%i", ARM, CHG), "", 1280, 800);
		c1->Divide(2, 2);
		c1->cd(1);
		H2->DrawCopy();
		F2->Draw("cont1 same");
		c1->cd(2);
		H2Ep->DrawCopy("cont1");
		c1->cd(3);
		TH1F *H1Sig = H1GetSignalSample(0, "dw23", ARM, CHG, 30, -0.1, 0.1);
		H1Sig->Scale(1/H1Sig->Integral());
		H1Sig->SetLineColor(1);
		H1Sig->Draw("hist");
		H1EpSlice->Scale(1/H1EpSlice->Integral());
		H1EpSlice->SetLineColor(4);
		H1EpSlice->DrawCopy("hist same");
		c1->cd(4);
		H1EpSlice->GetXaxis()->SetRangeUser(-0.1, 0.1);
		H1EpSlice->DrawCopy("hist e");
		RooPlot *p1 = dw23.frame();
		GpDw23->plotOn(p1);
		p1->Draw("same");

		if (Print) c1->Print(Form("Fit2DEp_%i%i.png", ARM, CHG));
	}

	return GpDw23;
}
