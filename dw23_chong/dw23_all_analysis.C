#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <TF2.h>
#include <TF1.h>
#include <TH2F.h>
#include <stdio.h>
#include "dw23_common.h"
#include "dw23_decomposition.C"
#include "dw23_2d_fits.C"

void write_dw23_text_file(TH2F * h2_dw23_vs_wness[][2]);
void read_dw23_text_file(TH2F * h2_gpr_dw23_vs_wness[][2]);

void dw23_all_analysis() {
  TF1 * f_double_gaus_1d[10][2][2]; 
  TH2F * h2_dw23_vs_wness[2][2];
  TH1F * h_dw23_slices[2][2][10], *h_wness_slices[2][2][10];
  //TH2F * h2_dw23_vs_wness_w[2][2];

  //open files and fill 2d dw23 vs wness histograms
  get_dw23_from_file(0, // int fileindex
      0, // bool do_rpc_cluster_cut
      3, // int cluster_cut
      h2_dw23_vs_wness,h_wness_slices,h_dw23_slices); //get data dw23_vs_wness

  write_dw23_text_file(h2_dw23_vs_wness);
  
  
  float target_wness[2][2];
  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      target_wness[a][c] = h_wness_slices[a][c][9]->GetMean();
    }
  }

  // read in the distribution resulting form the GPR fit
  TH2F * h2_gpr_dw23_vs_wness[2][2];
  read_dw23_text_file(h2_gpr_dw23_vs_wness);
  
  //1d dw23 wness slice fits for extrapolation
  TF1 * f_double_gaus_decomp[2][2][10]; 
  TF1 * f_double_gaus_decomp_target[2][2][3];
  TF1 * f_decomp_params[20][6][2][2];
  dw23_decomposition(0,h2_dw23_vs_wness,h_wness_slices,f_double_gaus_decomp,f_decomp_params,target_wness,f_double_gaus_decomp_target);
  
  
  
  //2d dw23 vs wness fits for linear & non-linear  extraplation
  TF1 * f_double_gaus_coax_linear[2][2][10]; 
  TF1 * f_double_gaus_offset_linear[2][2][10]; 
  TF1 * f_tripple_gaus_linear[2][2][10];
  /*TF1 * f_double_gaus_coax_nonl[2][2][10]; 
  TF1 * f_double_gaus_offset_nonl[2][2][10]; 
  TF1 * f_tripple_gaus_nonl[2][2][10];*/

  TF2 * f2_double_gaus_coax_linear[2][2]; 
  TF2 * f2_double_gaus_offset_linear[2][2]; 
  TF2 * f2_tripple_gaus_linear[2][2];
  /*TF2 * f2_double_gaus_coax_nonl[2][2]; 
  TF2 * f2_double_gaus_offset_nonl[2][2]; 
  TF2 * f2_tripple_gaus_nonl[2][2];*/
  
  TF1 * f_double_gaus_coax_target[2][2]; 
  TF1 * f_double_gaus_offset_target[2][2]; 
  TF1 * f_tripple_gaus_target[2][2];
  
  //assumed linear dw23 parameter vs wness relation
  dw23_2d_fits(0,h2_dw23_vs_wness,h_wness_slices,h_dw23_slices,
      f_double_gaus_coax_linear,
      f_double_gaus_offset_linear,
      f_tripple_gaus_linear,
      f2_double_gaus_coax_linear,
      f2_double_gaus_offset_linear,
      f2_tripple_gaus_linear,
      target_wness,
      f_double_gaus_coax_target,
      f_double_gaus_offset_target,
      f_tripple_gaus_target);

  //----
  //TH1F * h_sbg_fit_dw23[2][2];
  TFile * f_sbg_had_dw23 = new TFile("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/dw23_chong/output.root","RECREATE");

  TH1D * h_gpr_target[2][2];

  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      //sprintf(name,"h_sbg_fit_dw23_source0_a%d_c%d",arm,charge);
      //h_sbg_fit_dw23[a][c] = new TH1F(
     printf("test1\n");
      f_double_gaus_coax_target[a][c]->Write(); 
      f_double_gaus_offset_target[a][c]->Write(); 
      //f_tripple_gaus_target[a][c]->Write();
      f_double_gaus_decomp_target[a][c][0]->Write();
      f_double_gaus_decomp_target[a][c][1]->Write();
      f_double_gaus_decomp_target[a][c][2]->Write();
      
     printf("test2\n");
      for(int i=0; i<10; i++) {
        //f_double_gaus_coax_linear[a][c][i]->Write(); 
        //f_double_gaus_offset_linear[a][c][i]->Write();
        h_wness_slices[a][c][i]->Write();
        h_dw23_slices[a][c][i]->Write();
      }
      
     printf("test3\n");
      int nwness = h2_gpr_dw23_vs_wness[a][c]->GetYaxis()->GetNbins();
      h_gpr_target[a][c] = h2_gpr_dw23_vs_wness[a][c]->ProjectionX(Form("h_gpr_target_a%d_c%d",a,c),nwness*.92,nwness);
      h_gpr_target[a][c]->Scale(1.0/h_gpr_target[a][c]->Integral("width"));
      h_gpr_target[a][c]->Write();
      h2_gpr_dw23_vs_wness[a][c]->Write();
      h2_dw23_vs_wness[a][c]->Write();
     printf("test4\n");

    }
  }
  f_sbg_had_dw23->Close();
}



void write_dw23_text_file(TH2F * h2_dw23_vs_wness[][2]) {
  FILE *f1 = fopen("dw23_vs_wness_list.txt","w");
  FILE *f2[2][2];

  float content;
  double xcenter,ycenter;

  fprintf(f1,"index arm charge wness_bin_center dw23_bin_center entries\n");
  
  int index = 0;
  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      int wmax = h2_dw23_vs_wness[a][c]->GetNbinsX();
      int dmax = h2_dw23_vs_wness[a][c]->GetNbinsY();
  
      f2[a][c] = fopen(Form("dw23_vs_wness_grid_a%d_c%d.txt",a,c),"w");
      
      for(int wbin=1; wbin<=wmax; wbin++) {
        for(int dbin=1; dbin<=dmax; dbin++) {
          content = h2_dw23_vs_wness[a][c]->GetBinContent(wbin,dbin);
          xcenter = h2_dw23_vs_wness[a][c]->GetXaxis()->GetBinCenter(wbin);
          ycenter = h2_dw23_vs_wness[a][c]->GetYaxis()->GetBinCenter(dbin);

          fprintf(f1,"%d %d %d %f %f %.0f\n",index,a,c,xcenter,ycenter,content);
          index++;
          if(content==0 || wbin<=wmax/10 || wbin >=wmax*.9) {
            if(dbin<dmax)
              fprintf(f2[a][c],"nan ");
            else
              fprintf(f2[a][c],"nan\n");
          } else {
            if(dbin<dmax)
              fprintf(f2[a][c],"%.0f ",content);
            else
              fprintf(f2[a][c],"%.0f\n",content);
          }
        }
      }

      fclose(f2[a][c]);
    }
  }
      
  fclose(f1);
}

struct gpr_fit_values {
  int index, arm, charge;
  float dw23_bin_center, wness_bin_center, value7, uncertainty7, value8, uncertainty8, value9, uncertainty9;
  std::string GetDataString() {
    std::stringstream ss;
    ss << index << " " << arm << " " << charge  
      << " " << dw23_bin_center << " " << wness_bin_center 
      << " " << value7
      << " " << uncertainty7
      << " " << value8
      << " " << uncertainty8
      << " " << value9
      << " " << uncertainty9;
    return ss.str();
  }
};

void read_dw23_text_file(TH2F * h2_gpr_dw23_vs_wness[][2]) {
  std::vector<gpr_fit_values> data[2][2];

  for(int a=0; a<2; a++) {
    for(int c=0; c<2; c++) {
      std::string input_file = 
        Form("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/dw23_chong/predicted_dw23_wness_points_a%d_c%d.txt",
            a,c);
      std::cout << std::endl << "Loading file: " << input_file << std::endl;
      std::ifstream in_file(input_file.c_str());

      std::string line = "";

      bool first_time=true;

      if(in_file) { //check that infile opened
        while(getline(in_file,line)) {
          if(first_time) {first_time=false; continue;}

          gpr_fit_values data_temp;

          std::stringstream ss;
          ss.str(line.c_str());
          ss >> data_temp.index >> data_temp.arm >> data_temp.charge
            >> data_temp.dw23_bin_center >> data_temp.wness_bin_center
            >> data_temp.value7 >> data_temp.uncertainty7
            >> data_temp.value8 >> data_temp.uncertainty8
            >> data_temp.value9 >> data_temp.uncertainty9;

          data[a][c].push_back(data_temp);
        }
        std::cout << " Lines successfully extracted: " << data[a][c].size() << std::endl;
      } else
        std::cout << "ERROR! GPR RESULT TEXT FILE FAILED TO OPEN" << std::endl;

      float dw23=data[a][c][1].dw23_bin_center;
      float dw23_previous=data[a][c][0].dw23_bin_center;
      int nwness=1;

      while(dw23==dw23_previous) {
        dw23_previous=dw23;
        dw23=data[a][c][nwness+1].dw23_bin_center;
        nwness++;
      }
      
      int ndw23 = data[a][c].size()/nwness;

      std::cout << "dw23:" << ndw23 << " wness:" << nwness << std::endl;
      
      char *name = Form("h2_gpr_dw23_vs_wness_a%d_c%d",a,c);
      h2_gpr_dw23_vs_wness[a][c] = new TH2F(name,name,ndw23,-.3,.3,nwness,0,1);
      for(int i=0; i<data[a][c].size(); i++) {
        h2_gpr_dw23_vs_wness[a][c]->Fill(data[a][c][i].dw23_bin_center,data[a][c][i].wness_bin_center,data[a][c][i].value9);
        h2_gpr_dw23_vs_wness[a][c]->SetBinError(data[a][c][i].dw23_bin_center,data[a][c][i].wness_bin_center,data[a][c][i].uncertainty9);
      }

    }
  }
  
}

