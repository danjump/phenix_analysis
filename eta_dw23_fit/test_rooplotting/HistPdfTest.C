#include <TBranch.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>

#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif//__CINT__
#include <RooAbsPdf.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>

void HistPdfTest(void)
{
	using namespace RooFit;

	TFile *fin = TFile::Open("./stobg_high.root");
	TTree *tin = (TTree*)fin->Get("stobg_data");

	TCut Cut_SC   = "eta < 0 && charge < 0";
	TCut Cut_DCA  = "Rpc1dca < 100 && Rpc3dca < 100";
	TCut Cut_FVTX = "((fvtx_dr * fvtx_dtheta) < 150) && (TMath::Abs(fvtx_dphi) < 1.5)";
	TCut Cut_Comb = (Cut_SC + Cut_DCA + Cut_FVTX);

	TH2F *TestH = new TH2F("DCA_r_chi2", "", 100, 0, 30, 100, 0, 20);
	tin->Project(TestH->GetName(), "chi2:DCA_r", Cut_Comb);

	RooRealVar DCA_r("DCA_r", "DCA_r", 0, 30.0);
	RooRealVar chi2 ("chi2",  "chi2",  0, 20.0);

	RooArgSet *ASX  = new RooArgSet(DCA_r);
	RooArgSet *ASY  = new RooArgSet(chi2);
	RooArgSet *ASXY = new RooArgSet(DCA_r, chi2);

	//Projected into each direction
	RooDataHist *D1x = new RooDataHist("D1x", "", *ASX, TestH->ProjectionX());
	RooDataHist *D1y = new RooDataHist("D1y", "", *ASY, TestH->ProjectionY());
	RooHistPdf  *P1x = new RooHistPdf ("P1x", "", *ASX, *D1x);
	RooHistPdf  *P1y = new RooHistPdf ("P1y", "", *ASY, *D1y);
	RooAbsPdf	*AP1 = new RooProdPdf ("AP1", "", RooArgSet(*P1x, *P1y));

	//Whole 2D histogram
	RooDataHist *D2  = new RooDataHist("D2",  "", *ASXY, TestH);
	RooHistPdf  *P2  = new RooHistPdf ("P2",  "", *ASXY, *D2);
	RooAbsPdf	*AP2 = new RooProdPdf ("AP2", "", RooArgSet(*P2));


	//! ------------------------------------------------------------------

	//Draw
	RooPlot *f1_x = DCA_r.frame();
	RooPlot *f1_y =  chi2.frame();
	AP1->plotOn(f1_x);
	AP1->plotOn(f1_y);

	RooPlot *f2_x = DCA_r.frame();
	RooPlot *f2_y =  chi2.frame();
	//D2 ->plotOn(f2_x, LineColor(2), LineWidth(2));
	//D2 ->plotOn(f2_y, LineColor(2), LineWidth(2));
	AP2->plotOn(f2_x, LineColor(4), LineWidth(2));
	AP2->plotOn(f2_y, LineColor(4), LineWidth(2));

	//! ------------------------------------------------------------------

	TCanvas *c1 = new TCanvas("c1", "", 1280, 800);
	c1->Divide(2, 2);
	c1->cd(1);
	f1_x->Draw();
	c1->cd(2);
	f1_y->Draw();
	c1->cd(3);
	f2_x->Draw();
	c1->cd(4);
	f2_y->Draw();

	DCA_r = 4.10000;

	double Val1 = AP1->getVal(DCA_r);
	double Val2 = AP2->getVal(DCA_r);

	cout <<endl <<Val1 <<" " <<Val2 <<endl;

	return;
}

