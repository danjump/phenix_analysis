#include "HbgExtpol.C"

void TestExtpol(int ARM = 1, int CHG = 1, bool Show = false, bool Print = false)
{
	TFile *F = TFile::Open("./stobg_high.root");
	TTree *T = (TTree*)F->Get("stobg_data510");

	RooGenericPdf *HbgDw23 = HbgExtpol(T, ARM, CHG, 10, 30, Show, Print);

	return;
}
