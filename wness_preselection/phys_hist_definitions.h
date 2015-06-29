#ifndef __phys_hist_definitions_h
#define __phys_hist_definitions_h
#include <TH1F.h>

// physics related variables and Names
const int n_arms     = 2;		//arm
const int n_stations = 2;		//sta
const int n_octants  = 8;		//oct
const int n_halfocts = 2;		//hal
const int n_radsegments = 3;	//rad
const int n_paddles = 9;		//pad
const int n_strips   = 64;		//str
const int n_projections = 4;	//pro
const int n_dists = 14;    //distributions
const int n_trig = 5;    //3 triggers? rpc3, rpc1 and only sg1
const char distchar[14][64] = {"pT","pz","phi","eta",
			       "dg0","ddg0","dg4","chi2",
			       "dcaz","dcar","dphi12","dphi23",
			       "rpc1dca","rpc3dca"};
const char disttitles[14][64] = {"pT [GeV]","pz [GeV]","#phi [rad]","#eta ",
			       "DG0 [cm]","DDG0  [cm]","DG4  [cm]","#chi^{2}",
			       "DCA_z [cm]","DCA_r  [cm]","#delta #phi_{12} [rad]","#delta #phi_{23} [rad]",
			       "Rpc1dca [cm]","Rpc3dca [cm]"};
const float distmin[14] = {0,-250,-3.1415,1.2,
			 0.,0.,0.,0.,
			 0.,0.,-0.1,-0.1,
			 0.,0.};
const float distmax[14] = {60.,250.,3.1415,2.6,
			 20.,9.,28.,20.,
			 55.,55.,0.1,0.1,
			 30.,30.};


TH1F * physdists[n_arms][n_trig][n_dists];


#endif
