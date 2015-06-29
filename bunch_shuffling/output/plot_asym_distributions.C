#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TFormula.h>
#include <TLine.h>
#include <TH1F.h>
#include <TCanvas.h>

void plot_asym_distributions(char * infilename = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/bunch_shuffling/output/shuffled_100000.root",char * plots_subfolder = "test") {
  TFile * input = new TFile(infilename);

  
  TH1F * h_single_asym_b[2][2][4]; // [arm][charge][eta bin]
  TH1F * h_single_asym_y[2][2][4]; // [arm][charge][eta bin]
  TH1F * h_double_asym[2][2][4];   // [arm][charge][eta bin]

  char * arm_label[2] = {"South","North"};
  char charge_label[2] = {'-','+'};

  char name[300];

  for(int arm=0; arm<2; arm++) {
    for(int charge=0; charge<2; charge++) {
      for(int eta_bin=0; eta_bin<4; eta_bin++) {
        sprintf(name,"h_single_asym_b_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        h_single_asym_b[arm][charge][eta_bin] = (TH1F*)input->Get(name);

        sprintf(name,"h_single_asym_y_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        h_single_asym_y[arm][charge][eta_bin] = (TH1F*)input->Get(name);

        sprintf(name,"h_double_asym_arm%d_charge%d_eta%d",arm,charge,eta_bin);
        h_double_asym[arm][charge][eta_bin] = (TH1F*)input->Get(name);
      }
    }
  }
  
  /*float raw_single_asymmetry_fit[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ y=0 b=1 ]
  float raw_double_asymmetry_fit[2][2][4];    // [ arm ] [ charge ] [eta bins]

  float raw_single_asymmetry_err_fit[2][2][4][2]; // [ arm ] [ charge ] [eta bins] [ y=0 b=1 ]
  float raw_double_asymmetry_err_fit[2][2][4];    // [ arm ] [ charge ] [eta bins]*/
  
  TTree * t_real_asym = (TTree*)input->Get("t_real_asym");
  int tree_arm, tree_charge, tree_eta, tree_which;
  float tree_asym, tree_asym_err;

  t_real_asym->SetBranchAddress("arm",&tree_arm);
  t_real_asym->SetBranchAddress("charge",&tree_charge);
  t_real_asym->SetBranchAddress("eta_bin",&tree_eta);
  t_real_asym->SetBranchAddress("which_asym",&tree_which);
  t_real_asym->SetBranchAddress("asym",&tree_asym);
  t_real_asym->SetBranchAddress("asym_err",&tree_asym_err);

  int entries = t_real_asym->GetEntries();
  float asym_b[2][2][4];
  float asym_err_b[2][2][4];
  float asym_y[2][2][4];
  float asym_err_y[2][2][4];
  float asym_ll[2][2][4];
  float asym_err_ll[2][2][4];
  for(int i=0; i<entries; i++) {
    t_real_asym->GetEntry(i);

    if(tree_which==0) {
      asym_y[tree_arm][tree_charge][tree_eta] = tree_asym;
      asym_err_y[tree_arm][tree_charge][tree_eta] = tree_asym_err;

    } else if(tree_which==1) {
      asym_b[tree_arm][tree_charge][tree_eta] = tree_asym;
      asym_err_b[tree_arm][tree_charge][tree_eta] = tree_asym_err;
      
    } else if(tree_which==2) {
      asym_ll[tree_arm][tree_charge][tree_eta] = tree_asym;
      asym_err_ll[tree_arm][tree_charge][tree_eta] = tree_asym_err;
      
    } else {
      printf("\nproblem!\n\n");
    }
  }

  TTree * t_rand_asym = (TTree*)input->Get("t_rand_asym");
  int tree_rand_index;
  t_rand_asym->SetBranchAddress("rand_index",&tree_rand_index);
  t_rand_asym->SetBranchAddress("arm",&tree_arm);
  t_rand_asym->SetBranchAddress("charge",&tree_charge);
  t_rand_asym->SetBranchAddress("eta_bin",&tree_eta);
  t_rand_asym->SetBranchAddress("which_asym",&tree_which);
  t_rand_asym->SetBranchAddress("asym",&tree_asym);
  t_rand_asym->SetBranchAddress("asym_err",&tree_asym_err);
 
  printf("rand entries:%d\n",(int)t_rand_asym->GetEntries());
  
  TCanvas *c_asym_y[4];
  TCanvas *c_asym_b[4];
  TCanvas *c_asym_ll[4];
  TCanvas *c_paper_plots[4];

  TH1F * h_plots[4][2][2][3];
  TH1F * h_asym_y[2][2][4];
  TH1F * h_asym_b[2][2][4];
  TH1F * h_asym_ll[2][2][4];
  
  //gStyle->SetOptStat(0);


  for(int eta_bin=0; eta_bin<1; eta_bin++) {
    sprintf(name,"c_asym_y_eta%d",eta_bin);
    c_asym_y[eta_bin] = new TCanvas(name,name,2000,1400);
    c_asym_y[eta_bin]->Divide(2,2);

    sprintf(name,"c_asym_b_eta%d",eta_bin);
    c_asym_b[eta_bin] = new TCanvas(name,name,2000,1400);
    c_asym_b[eta_bin]->Divide(2,2);

    sprintf(name,"c_asym_ll_eta%d",eta_bin);
    c_asym_ll[eta_bin] = new TCanvas(name,name,2000,1400);
    c_asym_ll[eta_bin]->Divide(2,2);
    
    sprintf(name,"c_paper_plots_eta%d",eta_bin);
    c_paper_plots[eta_bin] = new TCanvas(name,name,2200,2200);
    c_paper_plots[eta_bin]->Divide(3,4);
    
    char eta_cond[300];
    sprintf(eta_cond,"eta_bin == %d",eta_bin);

    for(int arm=0; arm<2; arm++) {
      
      char arm_cond[300];
      sprintf(arm_cond," && arm == %d",arm);

      for(int charge=0; charge<2; charge++) {
        
        for(int which=0; which<3; which++) {
          sprintf(name,"h_plot_eta%d_arm%d_charge%d_which%d",eta_bin,arm,charge,which);
          h_plots[eta_bin][arm][charge][which] = new TH1F(name,name,100,-.15,.15);
        }
        
        char charge_cond[300];
        sprintf(charge_cond," && charge == %d",charge);

        char which_asym_cond[300];
        char conditions[500];

        double pars[3];
        double x_min,x_max,x_width,y_max;
        //TText * t_asym;
        TLine * line_asym;
        TLine * line_real_1;
        TLine * line_real_2;
        TLine * line_fit_1;
        TLine * line_fit_2;
        TF1 * f_gaus;

        sprintf(which_asym_cond," && which_asym == %d",0);
        sprintf(conditions,"%s%s%s%s",eta_cond,arm_cond,charge_cond,which_asym_cond);
        c_asym_y[eta_bin]->cd(2*arm+charge+1);
        //h_plots[eta_bin][arm][charge][0]->Draw();
        t_rand_asym->Draw("scaled_asym>>htemp(100,-5,5)",conditions);//,"same");
        sprintf(name,"h_asym_y_arm%d_charge%d_eta%d",arm, charge, eta_bin);
        h_asym_y[arm][charge][eta_bin] = (TH1F*)gPad->GetPrimitive("htemp")->Clone(name);
        sprintf(name,"Scaled A_{Ly} Distribution %s W^{%c}",arm_label[arm],charge_label[charge]);
        h_asym_y[arm][charge][eta_bin]->SetTitle(name);
        h_asym_y[arm][charge][eta_bin]->GetXaxis()->SetTitle("Scaled Asymmetry");
        h_asym_y[arm][charge][eta_bin]->Draw();
        h_asym_y[arm][charge][eta_bin]->Fit("gaus");
        f_gaus = (TF1*)h_asym_y[arm][charge][eta_bin]->GetFunction("gaus");
        f_gaus->GetParameters((Double_t*)pars);
        line_fit_1 = new TLine(-pars[2]+pars[1],0,-pars[2]+pars[1],pars[0]);
        line_fit_2 = new TLine(pars[2]+pars[1],0,pars[2]+pars[1],pars[0]);
        line_fit_1->SetLineColor(2);
        line_fit_2->SetLineColor(2);
        line_fit_1->Draw();
        line_fit_2->Draw();
        line_real_1 = new TLine(-1+pars[1],0,-1+pars[1],pars[0]);
        line_real_2 = new TLine(1+pars[1],0,1+pars[1],pars[0]);
        line_real_1->SetLineStyle(2);
        line_real_2->SetLineStyle(2);
        line_real_1->SetLineColor(9);
        line_real_2->SetLineColor(9);
        line_real_1->Draw();
        line_real_2->Draw();
        x_min = h_asym_y[arm][charge][eta_bin]->GetXaxis()->GetXmin();
        x_max = h_asym_y[arm][charge][eta_bin]->GetXaxis()->GetXmax();
        x_width = x_max - x_min;
        y_max = h_asym_y[arm][charge][eta_bin]->GetMaximum();
        sprintf(name,"Fit Mean: %f +- %f",pars[1],f_gaus->GetParError(1));
        t_asym = new TText(x_min+.04*x_width,y_max,name);
        t_asym->SetTextSize(.045);
        t_asym->SetTextColor(9);
        t_asym->Draw();
        /*line_asym = new TLine(asym_y[arm][charge][eta_bin],0,asym_y[arm][charge][eta_bin],pars[0]);
        line_asym->SetLineColor(9);
        line_asym->SetLineStyle(3);
        line_asym->Draw();*/
        sprintf(name,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/bunch_shuffling/%s/asym_y_eta%d.png",plots_subfolder,eta_bin);
        c_asym_y[eta_bin]->SaveAs(name);
        c_paper_plots[eta_bin]->cd(1+arm*2*3+charge*3);
        sprintf(name,"A_{Ly} %s W^{%c}",arm_label[arm],charge_label[charge]);
        h_asym_y[arm][charge][eta_bin]->SetTitle(name);
        h_asym_y[arm][charge][eta_bin]->Draw();
        line_fit_1->Draw();
        line_fit_2->Draw();
        line_real_1->Draw();
        line_real_2->Draw();
        t_asym->Draw();

        sprintf(which_asym_cond," && which_asym == %d",1);
        sprintf(conditions,"%s%s%s%s",eta_cond,arm_cond,charge_cond,which_asym_cond);
        c_asym_b[eta_bin]->cd(2*arm+charge+1);
       // h_plots[eta_bin][arm][charge][1]->Draw();
        t_rand_asym->Draw("scaled_asym>>htemp(100,-5,5)",conditions);//,"same");
        sprintf(name,"h_asym_b_arm%d_charge%d_eta%d",arm, charge, eta_bin);
        h_asym_b[arm][charge][eta_bin] = (TH1F*)gPad->GetPrimitive("htemp")->Clone(name);
        sprintf(name,"Scaled A_{Lb} Distribution %s W^{%c}",arm_label[arm],charge_label[charge]);
        h_asym_b[arm][charge][eta_bin]->SetTitle(name);
        h_asym_b[arm][charge][eta_bin]->GetXaxis()->SetTitle("Scaled Asymmetry");
        //h_asym_b[arm][charge][eta_bin]->GetXaxis()->SetTitleSize(.01);
        h_asym_b[arm][charge][eta_bin]->Draw();
        h_asym_b[arm][charge][eta_bin]->Fit("gaus");
        f_gaus = (TF1*)h_asym_b[arm][charge][eta_bin]->GetFunction("gaus");
        f_gaus->GetParameters((Double_t*)pars);
        line_fit_1 = new TLine(-pars[2]+pars[1],0,-pars[2]+pars[1],pars[0]);
        line_fit_2 = new TLine(pars[2]+pars[1],0,pars[2]+pars[1],pars[0]);
        line_fit_1->SetLineColor(2);
        line_fit_2->SetLineColor(2);
        line_fit_1->Draw();
        line_fit_2->Draw();
        line_real_1 = new TLine(-1+pars[1],0,-1+pars[1],pars[0]);
        line_real_2 = new TLine(1+pars[1],0,1+pars[1],pars[0]);
        line_real_1->SetLineStyle(2);
        line_real_2->SetLineStyle(2);
        line_real_1->SetLineColor(9);
        line_real_2->SetLineColor(9);
        line_real_1->Draw();
        line_real_2->Draw();
        x_min = h_asym_b[arm][charge][eta_bin]->GetXaxis()->GetXmin();
        x_max = h_asym_b[arm][charge][eta_bin]->GetXaxis()->GetXmax();
        x_width = x_max - x_min;
        y_max = h_asym_b[arm][charge][eta_bin]->GetMaximum();
        sprintf(name,"Fit Mean: %f +- %f",pars[1],f_gaus->GetParError(1));
        t_asym = new TText(x_min+.04*x_width,y_max,name);
        t_asym->SetTextSize(.045);
        t_asym->SetTextColor(9);
        t_asym->Draw();
        /*line_asym = new TLine(asym_b[arm][charge][eta_bin],0,asym_b[arm][charge][eta_bin],pars[0]);
        line_asym->SetLineColor(9);
        line_asym->SetLineStyle(3);
        line_asym->Draw();*/
        sprintf(name,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/bunch_shuffling/%s/asym_b_eta%d.png",plots_subfolder,eta_bin);
        c_asym_b[eta_bin]->SaveAs(name);
        c_paper_plots[eta_bin]->cd(2+arm*2*3+charge*3);
        sprintf(name,"A_{Lb} %s W^{%c}",arm_label[arm],charge_label[charge]);
        h_asym_b[arm][charge][eta_bin]->SetTitle(name);
        gStyle->SetTitleSize(0.07,"t");
        gStyle->SetTitleOffset(0.6,"t");
        //h_asym_b[arm][charge][eta_bin]->SetTitleSize(.08,"title");
        h_asym_b[arm][charge][eta_bin]->SetTitleSize(.065,"x");
        h_asym_b[arm][charge][eta_bin]->SetTitleOffset(.6,"x");
        h_asym_b[arm][charge][eta_bin]->Draw();
        line_fit_1->Draw();
        line_fit_2->Draw();
        line_real_1->Draw();
        line_real_2->Draw();
        t_asym->Draw();

        sprintf(which_asym_cond," && which_asym == %d",2);
        sprintf(conditions,"%s%s%s%s",eta_cond,arm_cond,charge_cond,which_asym_cond);
        c_asym_ll[eta_bin]->cd(2*arm+charge+1);
        //h_plots[eta_bin][arm][charge][2]->Draw();
        t_rand_asym->Draw("scaled_asym>>htemp(100,-5,5)",conditions);//,"same");
        sprintf(name,"h_asym_ll_arm%d_charge%d_eta%d",arm, charge, eta_bin);
        h_asym_ll[arm][charge][eta_bin] = (TH1F*)gPad->GetPrimitive("htemp")->Clone(name);
        sprintf(name,"Scaled A_{LL} Distribution %s W^{%c}",arm_label[arm],charge_label[charge]);
        h_asym_ll[arm][charge][eta_bin]->SetTitle(name);
        h_asym_ll[arm][charge][eta_bin]->GetXaxis()->SetTitle("Scaled Asymmetry");
        h_asym_ll[arm][charge][eta_bin]->Draw();
        h_asym_ll[arm][charge][eta_bin]->Fit("gaus");
        f_gaus = (TF1*)h_asym_ll[arm][charge][eta_bin]->GetFunction("gaus");
        f_gaus->GetParameters((Double_t*)pars);
        f_gaus->GetParError(1);
        line_fit_1 = new TLine(-pars[2]+pars[1],0,-pars[2]+pars[1],pars[0]);
        line_fit_2 = new TLine(pars[2]+pars[1],0,pars[2]+pars[1],pars[0]);
        line_fit_1->SetLineColor(2);
        line_fit_2->SetLineColor(2);
        line_fit_1->Draw();
        line_fit_2->Draw();
        line_real_1 = new TLine(-1+pars[1],0,-1+pars[1],pars[0]);
        line_real_2 = new TLine(1+pars[1],0,1+pars[1],pars[0]);
        line_real_1->SetLineStyle(2);
        line_real_2->SetLineStyle(2);
        line_real_1->SetLineColor(9);
        line_real_2->SetLineColor(9);
        line_real_1->Draw();
        line_real_2->Draw();
        x_min = h_asym_ll[arm][charge][eta_bin]->GetXaxis()->GetXmin();
        x_max = h_asym_ll[arm][charge][eta_bin]->GetXaxis()->GetXmax();
        x_width = x_max - x_min;
        y_max = h_asym_ll[arm][charge][eta_bin]->GetMaximum();
        sprintf(name,"Fit Mean: %f +- %f",pars[1],f_gaus->GetParError(1));
        t_asym = new TText(x_min+.04*x_width,y_max,name);
        t_asym->SetTextSize(.045);
        t_asym->SetTextColor(9);
        t_asym->Draw();
        /*line_asym = new TLine(asym_ll[arm][charge][eta_bin],0,asym_ll[arm][charge][eta_bin],pars[0]);
        line_asym->SetLineColor(9);
        line_asym->SetLineStyle(3);
        line_asym->Draw();*/
        sprintf(name,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/bunch_shuffling/%s/asym_ll_eta%d.png",plots_subfolder,eta_bin);
        c_asym_ll[eta_bin]->SaveAs(name);
        c_paper_plots[eta_bin]->cd(3+arm*2*3+charge*3);
        sprintf(name,"A_{LL} %s W^{%c}",arm_label[arm],charge_label[charge]);
        h_asym_ll[arm][charge][eta_bin]->SetTitle(name);
        h_asym_ll[arm][charge][eta_bin]->Draw();
        line_fit_1->Draw();
        line_fit_2->Draw();
        line_real_1->Draw();
        line_real_2->Draw();
        t_asym->Draw();

        sprintf(name,"/direct/phenix+WWW/p/draft/danielj/analysis_plots/run13_w_analysis/bunch_shuffling/%s/paper_plots%d.png",plots_subfolder,eta_bin);
        c_paper_plots[eta_bin]->SaveAs(name);

      }
    }
  }

  
}
