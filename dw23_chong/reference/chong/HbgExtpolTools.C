#include <Math/ProbFuncMathCore.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <TBranch.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualFitter.h>

//=================================================================================================

//Common cut
//-------------------------------------------------
const char *AC[2][2] = {{"SM", "SP"}, {"NM", "NP"}};
const char *CutAC[2][2] =
{
	{"eta < 0 && charge < 0", "eta < 0 && charge > 0"}, {"eta > 0 && charge < 0", "eta > 0 && charge > 0"}
};

TCut CutEta		 = "TMath::Abs(eta) > 1.0 && TMath::Abs(eta) < 2.6";	//!!!!!
TCut CutFvtxCone = "fvtx_cone > 0.";
TCut CutWnessFin = "Wness > 0.92";

//=================================================================================================

//Get pol4 function for Wness axis
//----------------------------------------------------------------------------
TF1 *F1WnessPol4(TTree *T, int ARM, int CHG, int Bins = 10, bool Show = false)
{
	const float FitRange[2] = {0.1, 0.95};

	TH1F *H1 = new TH1F(Form("H1WnessPol4_%i%i", ARM, CHG), Form("%s;Wness;Entries", AC[ARM][CHG]),	Bins, 0, 1);
	T->Project(H1->GetName(), "Wness", CutAC[ARM][CHG] + CutEta + CutFvtxCone);

	TF1 *F1 = new TF1(Form("F1WnessPol4_%i%i", ARM, CHG), "pol4", 0, 1);
	H1->Fit(F1->GetName(), "EQR0", "kRed", FitRange[0], FitRange[1]);

	if (Show)
	{
		gStyle->SetOptStat(0);
		TCanvas *c1 = new TCanvas(Form("c1WnessPol4_%i%i", ARM, CHG), AC[ARM][CHG], 1);
		c1->cd()->SetLogy();

		H1->SetLineColor(1);
		H1->SetMarkerColor(1);
		H1->SetMarkerSize(0.8);
		H1->SetMarkerStyle(20);
		H1->DrawCopy("pe");
		F1->SetLineStyle(2);
		F1->Draw("same");

		TLine *L1 = new TLine(FitRange[0], 0, FitRange[0], H1->GetMaximum());
		TLine *L2 = new TLine(FitRange[1], 0, FitRange[1], H1->GetMaximum());
		L1->SetLineStyle(2);
		L2->SetLineStyle(2);
		L1->Draw();
		L2->Draw();

		return F1;
	}
	else
	{
		H1->Delete();
		return F1;
	}
}//F1WnessPol4

//Get HBG dw23 slice in target region
//-------------------------------------------------------------------------------------------------------------
TH1F *H1GetDw23Slice(TTree *T, int ARM, int CHG, int BinDw23, float CutVal, float StepWidth, bool Show = false)
{
	TCut CutWnessSlice = Form("Wness > %f && Wness < %f", CutVal, CutVal + StepWidth);

	TH1F *H1 = new TH1F(
			Form("H1dw23Slice_%i%i_%i", ARM, CHG, (int)(CutVal*100)),
			Form("%s, %3.2f < Wness < %3.2f;dw23;Entries", AC[ARM][CHG], CutVal, CutVal + StepWidth),
			BinDw23, -0.3, 0.3
			);
	T->Project(H1->GetName(), "dw23", CutAC[ARM][CHG] + CutEta + CutFvtxCone + CutWnessSlice);

	if (Show)
	{
		TCanvas *c1 = new TCanvas(Form("c1dw23Slice_%i%i_%i", ARM, CHG, (int)(CutVal*100)), AC[ARM][CHG], 1);
		c1->cd()->SetGrid();
		H1->SetLineColor(1);
		H1->SetMarkerColor(1);
		H1->SetMarkerSize(0.8);
		H1->SetMarkerStyle(20);
		H1->DrawCopy("pe");
	}

	return H1;
}//H1dw23Slice

//Get sample of signal region
//--------------------------------------------------------------------------------------------------
TH1F *H1GetSignalSample(int Type, const char *Var, int ARM, int CHG, int Bins, float AxS, float AxE)
{
	enum {data510, w, onlyz, dy, onium, openbottom, opencharm, whad, wtau, z};

	const char *TitlePrs[10] =
	{
		"data510", "w", "onlyz", "dy", "onium", "openbottom", "opencharm", "whad", "wtau", "z"
	};

	TFile *F = TFile::Open("stobg_high.root");		//!!!!!
	TTree *T[10];
	for (int a=0; a<10; a++) T[a] = (TTree*)F->Get(Form("stobg_%s", TitlePrs[a]));

	const float LumiPerPb[10] =		//high condition luminosity		//!!!!!
	{
		53.10,		//data, Run12 pp510, by Sanghwa (2014 May)
		104457.83,	//w
		800751.88,	//onlyz
		120.30,		//dy
		410.89,		//onium
		548.36,		//openbottom
		235.06,		//opencharm
		48795.18,	//whad
		49397.59,	//wtau
		15421.38	//z
	};

	//Get temporary histogram for each process
	//----------------------------------------
	TH1F *H1Temp[10];
	for (int a=0; a<10; a++)
	{
		H1Temp[a] = new TH1F(Form("H1Temp%s_%i_%i%i_%s", Var, Type, ARM, CHG, TitlePrs[a]), "", Bins, AxS, AxE);
		T[a]->Project(H1Temp[a]->GetName(), Var, CutAC[ARM][CHG] + CutEta + CutWnessFin);

		H1Temp[a]->Sumw2();
		H1Temp[a]->Scale(LumiPerPb[data510]/LumiPerPb[a]);

		//2nd scale by dimuon study results (Sanghwa)
		//-----------------------------------------------
		if		(a == onlyz)	  H1Temp[a]->Scale(1.445);		//+- 0.353,
		else if (a == onium)	  H1Temp[a]->Scale(0.634);		//+- 0.079, Upsilon
		else if (a == openbottom) H1Temp[a]->Scale(3.186);		//+- 1.064
		else if (a == opencharm)  H1Temp[a]->Scale(3.539);		//+- 0.823
		else if (a == z)		  H1Temp[a]->Scale(1.445);		//+- 0.353
	}

	//Get summed histograms
	//------------------------------------------------------------------------------------------------------
	TH1F *H1 = new TH1F(Form("H1%s_%i_%i%i", Var, Type, ARM, CHG), Form(";%s;Entries", Var), Bins, AxS, AxE);
	H1->Sumw2();
	H1->SetLineColor(Type+1);

	enum {DAT, SIG, MBG};
	if		(Type == DAT) H1->Add(H1Temp[data510]);
	else if (Type == SIG)
	{
		H1->Add(H1Temp[w]);
		H1->Add(H1Temp[onlyz]);
	}
	else if (Type == MBG)
	{
		H1->Add(H1Temp[onlyz], -1);
		H1->Add(H1Temp[dy]);
		H1->Add(H1Temp[onium]);
		H1->Add(H1Temp[openbottom]);
		H1->Add(H1Temp[opencharm]);
		H1->Add(H1Temp[whad]);
		H1->Add(H1Temp[wtau]);
		H1->Add(H1Temp[z]);

		//Scale by Muon efficiency (overall trigger + hit)	//!!!!!
		//-----------------------------------------------------------------------
		if (ARM==0) H1->Scale(0.47 * 0.9 * 0.9 * (0.77 - 0.432*0.32) * 0.9*0.798);
		else		H1->Scale(0.57 * 0.9 * 0.9 * (0.77 - 0.432*0.32) * 0.9*0.774);		
	}

	return H1;
}//H1GetSignalSample

//=================================================================================================

//Get Hadron MC luminosity by UA1 parameters
//---------------------------------------------------------------
double GetHadMCLumi(int ARM, int pT, int Type, bool Show = false)
{
	const long GenEventsPerK[2][9][4] =
	{
		{//South
			{		0,		  0,		0,		  0},	//0, pT bin, pi+/pi-/K+/K-
			{14737000, 14642000, 14351000, 14378000},	//1, from revents
			{24980000, 24765000, 24734000, 24664000},	//2, from revents
			{ 5898550,  5881250,  5886650,  5881100},	//3
			{ 3854850,  3856450,  1310250,  1314900},	//4
			{  902900,   902900,   311200,   309900},	//5
			{  259400,   258400,	91100,	  88200},	//6
			{	85800,	  86300,	30600,	  30400},	//7
			{	32800,	  32600,	30700,	  30300}	//8
		},
		{//North
			{		0,		  0,		0,		  0},
			{14750000, 14641000, 14512000, 14450000},	//1, from revents
			{25033000, 24952000, 24830000, 24879000},	//2, from revents
			{ 5905000,  5893450,  5881300,  5889700},	//3
			{ 3855700,  3851200,  1314150,  1316500},	//4
			{  903600,   903400,   308900,   309400},	//5
			{  259600,   260700,    88000,    87500},	//6
			{	87200,	  85600,	30600,	  30000},	//7
			{	32900,	  33100,	30300,	  30200}	//8
		}
	};

	double pT_offset[9] = {0, 0.2909, 0.3279, 0.3560, 0.3759, 0.3924, 0.4048, 0.4159, 0.4225};
	double pT_real	    = pT + pT_offset[pT];

	double XsecPerPb = 0;
	if (Type < 2) XsecPerPb = 1.800 * TMath::Power(1 + pT_real/1.61, -10.64) * pT_real * 1.e12;		//pi
	else		  XsecPerPb = 0.616 * TMath::Power(1 + pT_real/1.61, -10.64) * pT_real * 1.e12;		//K

	double LumiPerPb = (GenEventsPerK[ARM][pT][Type] * 1.e3) / XsecPerPb;
	if (Show) cout <<Form("%i %i %i | %7.2f M | %10.5f", ARM, pT, Type, XsecPerPb/1.e6, LumiPerPb) <<endl;

	return LumiPerPb;
}//GetHadMCLumi

//Get summed (pi, K) Hadron MC distribution
//------------------------------------------------------------------------------------------
void GetHadMCSummed(TH2F *(*H2)[2], int BinWness = 100, int BinDw23 = 30, bool Show = false)
{
	const int TypeNum[4] = {8, 9, 11, 12};	//pi+, pi-, K+, and K-
	const float LumiData = 53.1;	//!!!!!

	//Link trees
	//-----------------------------------------------
	TFile *F = TFile::Open("./stobg_high_HadMC.root");
	TTree *T[2][9][4];
	for (int i=0; i<2; i++) {	//Arm
	for (int a=1; a<9; a++) {	//pT bins
	for (int b=0; b<4; b++) {	//Type

		T[i][a][b] = (TTree*)F->Get(Form("HadMC_%i_%i_s%i", i, a, TypeNum[b]));

	}}}

	//Get temporary H2 (dw23:Wness) for given tree, then scale it by Luminosity
	//-------------------------------------------------------------------------
	TCut CutWnessAll = "Wness > 0.1 && Wness < 1.0";

	TH2F *H2Temp[2][2][9][4];
	for (int i=0; i<2; i++) {	//Arm
	for (int j=0; j<2; j++) {	//Charge
	for (int a=1; a<9; a++) {	//pT bins
	for (int b=0; b<4; b++) {	//Type

		H2Temp[i][j][a][b] = new TH2F(
				Form("H2Temp_%s_%s", T[i][a][b]->GetName(), j==0?"M":"P"), "", BinWness, 0, 1, BinDw23, -0.3, 0.3
				);
		T[i][a][b]->Project(
				H2Temp[i][j][a][b]->GetName(), "dw23:Wness", CutAC[1-i][j] + CutEta + CutFvtxCone + CutWnessAll
				);

		double LumiHadMC = GetHadMCLumi(i, a, b, (i==0&&j==0)?Show:0);
		H2Temp[i][j][a][b]->Scale(LumiData/LumiHadMC);

	}}}}

	//Fill up - Caution: Arm is swapped, Exclude pT bin = 1 as it stirs distribution
	//------------------------------------------------------------------------------
	for (int i=0; i<2; i++) {
	for (int j=0; j<2; j++) {

		H2[i][j] = new TH2F(Form("H2_HadMC_%i%i", i, j), "", BinWness, 0, 1, BinDw23, -0.3, 0.3);

		for (int a=2; a<9; a++) {	//! Exclude pT bin = 1
		for (int b=0; b<4; b++) {
			
			H2[i][j]->Add(H2Temp[1-i][j][a][b]);

		}}
	}}

	return;
}//GetHadMCSummed

//=================================================================================================

//Decompose dw23 into 'Low + Wide' part and 'High + Narrow' part
//------------------------------------------------------------------------------------------------------------------
void DecomposeDw23(TH1F *H1, TF1 *F1Lo, TF1 *F1Hi, bool SignalRegion = false, bool Show = false, bool Print = false)
{
	//Determine sidebands (outer sides of dw23)
	//-----------------------------------------
	float LimL = 0, LimR = 0;
	if (SignalRegion)
	{
		int CHG = -1;
		if		(H1->GetMean() < 0) CHG = 0;
		else if (H1->GetMean() > 0) CHG = 1;

		LimL = -0.04 + 0.03*CHG - 1.e-3;
		LimR = +0.01 + 0.03*CHG + 1.e-3;
	}
	else
	{
		LimL = H1->GetMean() - H1->GetRMS()*1.35 - 1.e-3;
		LimR = H1->GetMean() + H1->GetRMS()*1.35 + 1.e-3;
	}
	const int BinLimL = H1->FindBin(LimL), BinLimR = H1->FindBin(LimR);

	//Get sidebands only histogram
	//----------------------------
	TH1F *H1a = new TH1F(
			Form("%s_SBonly", H1->GetName()), "",
			H1->GetNbinsX(), H1->GetXaxis()->GetXmin(), H1->GetXaxis()->GetXmax()
			);
	for (int x=1;		x<=BinLimL;			x++) H1a->SetBinContent(x, H1->GetBinContent(x));
	for (int y=BinLimR; y<=H1->GetNbinsX(); y++) H1a->SetBinContent(y, H1->GetBinContent(y));

	//Fit on sidebands only histogram, then extrapolate it
	//----------------------------------------------------
	F1Lo->SetParameter(0, H1a->GetMaximum());
	F1Lo->SetParameter(1, H1a->GetMean());
	F1Lo->SetParameter(2, H1a->GetRMS());
	F1Lo->SetParLimits(0, 0, H1a->GetMaximum());	//Height cannot have negative value
	H1a->Fit(F1Lo->GetName(), "EQR0");

	TH1F *H1aEp = new TH1F(
			Form("%sEp", H1a->GetName()), "",
			H1a->GetNbinsX(), H1a->GetXaxis()->GetXmin(), H1a->GetXaxis()->GetXmax()
			);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(H1aEp);

	//Get remaining histograms after sidebands extrapolation subtracted
	//-----------------------------------------------------------------
	TH1F *H1b = (TH1F*)H1->Clone(Form("%sSub", H1aEp->GetName()));
	H1b->Add(H1aEp, -1);
	F1Hi->SetParameter(0, H1b->GetMaximum());
	F1Hi->SetParameter(1, H1b->GetMean());
	F1Hi->SetParameter(2, H1b->GetRMS());
	F1Hi->SetParLimits(0, 0, H1b->GetMaximum());
	H1b->Fit(F1Hi->GetName(), "EQR0");

	if (Show)
	{
		gStyle->SetOptStat();
		gStyle->SetOptFit();

		const int YMin = H1b->GetMinimum()*1.5, YMax = H1->GetMaximum()*1.25;
		H1   ->GetYaxis()->SetRangeUser(YMin, YMax);
		H1a  ->GetYaxis()->SetRangeUser(YMin, YMax);
		H1aEp->GetYaxis()->SetRangeUser(YMin, YMax);
		H1b  ->GetYaxis()->SetRangeUser(YMin, YMax);

		TLine *LineL = new TLine(LimL, 0, LimL, H1->GetMaximum());
		TLine *LineR = new TLine(LimR, 0, LimR, H1->GetMaximum());
		LineL->SetLineStyle(2);
		LineR->SetLineStyle(2);

		TCanvas *c1 = new TCanvas(Form("c1_%s", H1->GetName()), "", 1280, 800);
		c1->Divide(2, 2);
		c1->cd(1);
		H1->DrawCopy("hist e");
		LineL->Draw("same");
		LineR->Draw("same");
		c1->cd(2);
		H1->DrawCopy("hist e");
		TF1 *F1Sum = new TF1(
				Form("F1Sum_%s", H1->GetName()), Form("%s + %s", F1Lo->GetName(), F1Hi->GetName()), -0.3, 0.3
				);
		F1Sum->Draw("same");
		c1->cd(3);
		H1a->DrawCopy("hist e");
		F1Lo->Draw("same");
		H1aEp->SetLineColor(2);
		H1aEp->DrawCopy("le same");
		LineL->Draw("same");
		LineR->Draw("same");
		c1->cd(4);
		H1b->DrawCopy("hist 2");
		F1Hi->DrawCopy("same");
		LineL->Draw("same");
		LineR->Draw("same");

		if (Print) c1->Print(Form("%s.png", c1->GetName()));
	}

	return;
}//DecomposeDw23

//Decompose dw23, but set sidebands range manually
//------------------------------------------------------------------------------
void DecomposeDw23Manual(TH1F *H1, TF1 *F1Lo, TF1 *F1Hi, float LimL, float LimR)
{
	const int BinLimL = H1->FindBin(LimL);
	const int BinLimR = H1->FindBin(LimR);

	TH1F *H1a = new TH1F(
			Form("%s_SBonly", H1->GetName()), "",
			H1->GetNbinsX(), H1->GetXaxis()->GetXmin(), H1->GetXaxis()->GetXmax()
			);
	for (int x=1;		x<=BinLimL;			x++) H1a->SetBinContent(x, H1->GetBinContent(x));
	for (int y=BinLimR; y<=H1->GetNbinsX(); y++) H1a->SetBinContent(y, H1->GetBinContent(y));

	F1Lo->SetParameter(0, H1a->GetMaximum());
	F1Lo->SetParameter(1, H1a->GetMean());
	F1Lo->SetParameter(2, H1a->GetRMS());
	F1Lo->SetParLimits(0, 0, H1a->GetMaximum());
	H1a->Fit(F1Lo->GetName(), "EQR0");

	TH1F *H1aEp = new TH1F(
			Form("%sEp", H1a->GetName()), "",
			H1a->GetNbinsX(), H1a->GetXaxis()->GetXmin(), H1a->GetXaxis()->GetXmax()
			);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(H1aEp);

	TH1F *H1b = (TH1F*)H1->Clone(Form("%sSub", H1aEp->GetName()));
	H1b->Add(H1aEp, -1);
	F1Hi->SetParameter(0, H1b->GetMaximum());
	F1Hi->SetParameter(1, H1b->GetMean());
	F1Hi->SetParameter(2, H1b->GetRMS());
	F1Hi->SetParLimits(0, 0, H1b->GetMaximum());
	H1b->Fit(F1Hi->GetName(), "EQR0");

	return;
}

//=================================================================================================

//Get Low + Wide Gaussian's width in Signal sample
//----------------------------------------------------------------
float GetHbgDw23SignalLoWidth(int ARM, int CHG, bool Show = false)
{
	enum {DAT, SIG, MBG};

	TH1F *H1DAT = H1GetSignalSample(DAT, "dw23", ARM, CHG, 32, -0.1, 0.1);
	TH1F *H1MBG = H1GetSignalSample(MBG, "dw23", ARM, CHG, 32, -0.1, 0.1);

	H1DAT->Add(H1MBG, -1);
	H1DAT->SetTitle(Form("%s, MBG subtracted data (Wness > 0.92)", AC[ARM][CHG]));

	TF1 *F1Lo = new TF1(Form("F1Lo_%s", H1DAT->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
	TF1 *F1Hi = new TF1(Form("F1Hi_%s", H1DAT->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
	DecomposeDw23(H1DAT, F1Lo, F1Hi, true, Show);

	return F1Lo->GetParameter(2);
}//GetHbgDw23SignalLo

//Nonlinear function
//------------------------------
double NLFcn(double x, double p)
{
	return 0.03 * (TMath::Erf(10 * p * (1.01 - x)) -1);
};

//pol1 + Nonlinear
//------------------------------------
double Pol1NLFcn(double *x, double *p)
{
	return p[0] + p[1] * x[0] + 0.03 * (TMath::Erf(10 * p[2] * (1.01 - x[0])) -1);
}

//Get dw23 patterns for the 2D fit
//---------------------------------------------------------------------------------------------------------
void GetHbgDw23Pattern(TTree *T, int ARM, int CHG, float FitPars[4], bool Show = false, bool Print = false)
{
	enum {d510, mcHad};
	enum {Lo, Hi};

	const float CutSet[] = {0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.94, 0.98, 1.0};
	const int nCutSet = sizeof(CutSet)/sizeof(CutSet[0]) - 1;

	TF1 *F1[2][2][nCutSet];		//Sample, Low/High, Cut

	//Get dw23 parameters from data 510 
	//---------------------------------
	for (int a=0; a<nCutSet; a++)
	{
		TH1F *H1Temp = H1GetDw23Slice(T, ARM, CHG, 30, CutSet[a], CutSet[a+1] - CutSet[a]);

		F1[d510][Lo][a] = new TF1(Form("F1Lo_%s", H1Temp->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
		F1[d510][Hi][a] = new TF1(Form("F1Hi_%s", H1Temp->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
		DecomposeDw23(H1Temp, F1[d510][Lo][a], F1[d510][Hi][a], CutSet[a] > 0.92?true:false, Show, Print);
	}

	//Get parameters from Hadron MC
	//-----------------------------
	TH2F *H2Had[2][2];
	GetHadMCSummed(H2Had);

	for (int a=0; a<nCutSet; a++)
	{
		TH1F *H1HadTemp = (TH1F*)H2Had[ARM][CHG]->ProjectionY(
				Form("%s_%i", H2Had[ARM][CHG]->GetName(), a),
				H2Had[ARM][CHG]->GetXaxis()->FindBin(CutSet[a]	+ 1.e-3),
				H2Had[ARM][CHG]->GetXaxis()->FindBin(CutSet[a+1] - 1.e-3)
				);
		H1HadTemp->SetTitle(Form("MChad, %s, %3.2f < Wness < %3.2f", AC[ARM][CHG], CutSet[a], CutSet[a+1]));

		F1[mcHad][Lo][a] = new TF1(Form("F1Lo_%s", H1HadTemp->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
		F1[mcHad][Hi][a] = new TF1(Form("F1Hi_%s", H1HadTemp->GetName()), "[0]*TMath::Gaus(x, [1], [2])", -0.3, 0.3);
		DecomposeDw23(H1HadTemp, F1[mcHad][Lo][a], F1[mcHad][Hi][a], CutSet[a] > 0.92?true:false, Show, Print);
	}

	//Get Parameter Vs. Wness graphs
	//------------------------------
	TGraphErrors *G1[2][2][4];
	for (int a=0; a<2; a++) {		//data, Hadron MC
	for (int b=0; b<2; b++) {		//Lo, Hi
	for (int c=0; c<4; c++) {		//Height, Center, Width, chi2/NDF

		G1[a][b][c] = new TGraphErrors();
		G1[a][b][c]->SetName(Form("G1_%i%i%i", a, b, c));

		for (int d=0; d<nCutSet; d++)
		{
			if (c<3)	//Height, Center, Width
			{
				G1[a][b][c]->SetPoint(
						G1[a][b][c]->GetN(), (CutSet[d] + CutSet[d+1])/2, F1[a][b][d]->GetParameter(c)
						);
				G1[a][b][c]->SetPointError(
						G1[a][b][c]->GetN() - 1, (CutSet[d+1] - CutSet[d])/2, F1[a][b][d]->GetParError(c)
						);
			}
			else	//chi2/NDF
			{
				G1[a][b][c]->SetPoint(
						G1[a][b][c]->GetN(), (CutSet[d] + CutSet[d+1])/2,
						F1[a][b][d]->GetChisquare()/F1[a][b][d]->GetNDF()
						);			
			}
		}

	}}}

	//Perform pol1 fit on Width2 (Lo) Vs. Wness in data sample, then get its slope
	//----------------------------------------------------------------------------
	TF1 *F1LoWidth = new TF1(Form("F1LoWidth_%i%i", ARM, CHG), "pol1", 0.1, 1);
	F1LoWidth->SetParameter(0, ARM==0?0.1:0.15);
	F1LoWidth->SetParameter(1, -0.1);
	G1[d510][Lo][2]->Fit(F1LoWidth->GetName(), "EQR0", "", 0.1, 1);

	//Perform pol1 + NL fit on Width1 (Hi) Vs. Wness - Get Pol1 part from data, Get NL drop starting point in mcHad
	//-------------------------------------------------------------------------------------------------------------
	TF1 *F1HiWidth[2];
	for (int a=0; a<2; a++)
	{
		F1HiWidth[a] = new TF1(Form("F1HiWidth_%i%i", ARM, CHG), Pol1NLFcn, 0.1, 1, 3);
		F1HiWidth[a]->SetParameter(0,  0.1);
		F1HiWidth[a]->SetParameter(1, -0.1);
		F1HiWidth[a]->SetParameter(2,  1);
		G1[a][Hi][2]->Fit(F1HiWidth[a]->GetName(), "EQR0", "", 0.1, 1);
	}

	//if (Show)
	{
		const char *ParTitle[] = {"Height", "Center", "Width", "Chi2/NDF"};
		TCanvas *c1 = new TCanvas(Form("c1_dw23FitParms_%i%i", ARM, CHG), AC[ARM][CHG], 1280, 800);
		c1->Divide(4, 2);

		for (int a=0; a<2; a++) {	//Sample
		for (int b=0; b<2; b++) {	//Lo/Hi
		for (int c=0; c<4; c++) {

			c1->cd(b*4 + c + 1)->SetGrid();

			if (b==Hi && c==2)
			{
				gStyle->SetOptStat();
				gStyle->SetOptFit();
			}

			G1[a][b][c]->SetLineColor(a+1);
			G1[a][b][c]->SetMarkerColor(a+1);
			G1[a][b][c]->SetMarkerStyle(a+24);
			G1[a][b][c]->SetTitle(Form("%s, %s, %s;Wness", AC[ARM][CHG], b==0?"Lo":"Hi", ParTitle[c]));

			if		(c==0) G1[a][b][c]->GetYaxis()->SetRangeUser(0, G1[a][b][c]->GetHistogram()->GetMaximum()*1.1);
			else if (c==1) G1[a][b][c]->GetYaxis()->SetRangeUser(-0.05 + 0.05*CHG, 0 + 0.05*CHG);
			else if (c==2) G1[a][b][c]->GetYaxis()->SetRangeUser(0, b==0?0.15:0.075);
			else if (c==3) G1[a][b][c]->GetYaxis()->SetRangeUser(0, 3);

			G1[a][b][c]->Draw(a==0?"ap":"p sames");

			if (a==d510 && b==Lo && c==2)
			{
				F1LoWidth->SetLineColor(1);
				F1LoWidth->SetLineStyle(2);
				F1LoWidth->Draw("same");
			}
			if (b==Hi && c==2)
			{
				F1HiWidth[a]->SetLineColor(a+1);
				F1HiWidth[a]->SetLineStyle(2);
				F1HiWidth[a]->Draw("same");
			}

		}}}

		if (Print) c1->Print(Form("%s.png", c1->GetName()));		
	}//Show

	//Return parameters: Width1's pol1 part from data, Width1's NL part from mcHad, Width2's slope
	//--------------------------------------------------------------------------------------------
	FitPars[0] = F1HiWidth[d510]->GetParameter(0);
	FitPars[1] = F1HiWidth[d510]->GetParameter(1);
	FitPars[2] = F1HiWidth[mcHad]->GetParameter(2);
	FitPars[3] = F1LoWidth->GetParameter(1);

	return;
}//GetHbgDw23Pattern
