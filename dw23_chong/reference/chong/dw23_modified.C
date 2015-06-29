
int get_hists_from_file(TFile * infile[11],TH1F * dw23[2][][]);

void dw23_modified( const std::string bkg_had_infilename,
    const std::string sig_w_infilename,
    const std::string bkg_dy_infilename,
    const std::string bkg_light_infilename,
    const std::string bkg_onium_infilename,
    const std::string bkg_onlyz_infilename,
    const std::string bkg_openbottom_infilename,
    const std::string bkg_opencharm_infilename,
    const std::string bkg_whad_infilename,
    const std::string bkg_wtau_infilename,
    const std::string bkg_zsum_infilename) {

  char name[300]; sprintf(name,"prevent unused variable error");
  
  Int_t
  Run_Number,  Evt_Number,  triggerbit,
  clockcross,  fvtx_cone;
  
  Float_t 
  Evt_bbcZ,    Wness,
  tree_charge,      pT,          pz,
  phi,         eta,         DG0,
  DDG0,        DG4,         chi2,
  DCA_z,       DCA_r,       dphi23,
  dw23,        Rpc1dca,     Rpc3dca,
  fvtx_dphi,   fvtx_dr,     fvtx_dtheta, fvtx_dr_dtheta;
  
  TFile *infile[11];
  std::string temp_str[11];

  temp_str[0] = bkg_had_infilename;
  temp_str[1] = sig_w_infilename;
  temp_str[2] = bkg_dy_infilename;
  temp_str[3] = bkg_light_infilename;
  temp_str[4] = bkg_onium_infilename;
  temp_str[5] = bkg_onlyz_infilename;
  temp_str[6] = bkg_openbottom_infilename;
  temp_str[7] = bkg_opencharm_infilename;
  temp_str[8] = bkg_whad_infilename;
  temp_str[9] = bkg_wtau_infilename;
  temp_str[10] = bkg_zsum_infilename;
  
  infile[0]= new TFile(bkg_had_infilename.c_str());
  infile[1]= new TFile(sig_w_infilename.c_str());
  infile[2] = new TFile(bkg_dy_infilename.c_str());
  infile[3] = new TFile(bkg_light_infilename.c_str());
  infile[4] = new TFile(bkg_onium_infilename.c_str());
  infile[5] = new TFile(bkg_onlyz_infilename.c_str());
  infile[6] = new TFile(bkg_openbottom_infilename.c_str());
  infile[7] = new TFile(bkg_opencharm_infilename.c_str());
  infile[8] = new TFile(bkg_whad_infilename.c_str());
  infile[9] = new TFile(bkg_wtau_infilename.c_str());
  infile[10] = new TFile(bkg_zsum_infilename.c_str());



  //load trigger efficiencies 
  float trig_eff_factor[2][2][15]; // 15=eta range bins
  for(int ebin=0; ebin < 15; ebin++) {
    for(int arm=0; arm < 2; arm++) {
      for(int charge=0; charge < 2; charge++) {
        trig_eff_factor[arm][charge][ebin] = get_eta_trigeff(ebin,arm,charge); 
        printf("e%d a%d c%d ef%f\n",ebin,arm,charge,trig_eff_factor[arm][charge][ebin]);
      }
    }
  }
  
  define_dw23_vs_eta_hists();
  
  TH2F * dw23_vs_eta = new TH2F("dw23_vs_eta","dw23_vs_eta",100,distmin[5],distmax[5],100,distmin[14],distmax[14]);
  
  int in_count[11][2][2],over_count[11][2][2],under_count[11][2][2],raw_target_wness_count[11][2][2];
  float trigeff_target_wness_count[11][2][2];


  //loop over the different input files and all events within each file to get distributions
  for(int file_index=0; file_index<1; file_index++) {//read loop over files
    for(int arm=0; arm<2; arm++) {
      for(int charge_index=0; charge_index<2; charge_index++) {
        in_count[file_index][arm][charge_index] = 0;
        over_count[file_index][arm][charge_index] = 0;
        under_count[file_index][arm][charge_index] = 0;
        raw_target_wness_count[file_index][arm][charge_index] = 0;
        trigeff_target_wness_count[file_index][arm][charge_index] = 0;
      }
    }

    printf("test2\n");
    
    printf("Reading file:\n%s\n",temp_str[file_index].c_str());
    define_wness_tree(infile[file_index]);

    wness_tree->SetBranchAddress("Run_Number",&Run_Number);
    wness_tree->SetBranchAddress("Evt_Number",&Evt_Number);
    wness_tree->SetBranchAddress("triggerbit",&triggerbit);
    wness_tree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
    wness_tree->SetBranchAddress("clockcross",&clockcross);
    wness_tree->SetBranchAddress("Wness",&Wness);
    wness_tree->SetBranchAddress("charge",&tree_charge);
    wness_tree->SetBranchAddress("pT",&pT);
    wness_tree->SetBranchAddress("pz",&pz);
    wness_tree->SetBranchAddress("phi",&phi);
    wness_tree->SetBranchAddress("eta",&eta);
    wness_tree->SetBranchAddress("DG0",&DG0);
    wness_tree->SetBranchAddress("DDG0",&DDG0);
    wness_tree->SetBranchAddress("DG4",&DG4);
    wness_tree->SetBranchAddress("chi2",&chi2);
    wness_tree->SetBranchAddress("DCA_z",&DCA_z);
    wness_tree->SetBranchAddress("DCA_r",&DCA_r);
    wness_tree->SetBranchAddress("dphi23",&dphi23);
    wness_tree->SetBranchAddress("dw23",&dw23);
    wness_tree->SetBranchAddress("Rpc1dca",&Rpc1dca);
    wness_tree->SetBranchAddress("Rpc3dca",&Rpc3dca);
    wness_tree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
    wness_tree->SetBranchAddress("fvtx_dr",&fvtx_dr);
    wness_tree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
    wness_tree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
    wness_tree->SetBranchAddress("fvtx_cone",&fvtx_cone);


    int entries = wness_tree->GetEntries();

    int percent_done=0;
    int percent_incriment=20;
    int percent_done_previous;

    time_t rawtime;

    cout << "\nNumber of events:  " << entries << endl;
    cout << "Starting Main Event Loop...\n\n";
    time(&rawtime);
    printf("Start time:  %s",ctime(&rawtime));

    printf("test3\n");

    //Get events from file, fill necessary histograms
    for(int i=0; i<entries; i++) {//Events loop
      //loop progress command line output
      percent_done_previous=percent_done;
      percent_done=(int)floor((float)(i+1)/(float)entries*(float)100);
      if(percent_done%percent_incriment==0 && percent_done != percent_done_previous) {
        printf("%3i%% done",percent_done);
        if(percent_done==100) {
          cout << "^_^";
          time( &rawtime );
          printf(" %s",ctime(&rawtime));
        } else {
          cout << "...";
          time( &rawtime );
          printf(" %s",ctime(&rawtime));
        }
      }

      wness_tree->GetEntry(i);

      // Set arm (0 = South, 1 = North)
      int arm = 0;
      if		(pz < 0) arm = 0;
      else if (pz > 0) arm = 1;


      int charge_index = 0;
      if		(tree_charge < 0) charge_index = 0;
      else if (tree_charge > 0) charge_index = 1;

      int rpcdca_condition = -1;
      if((Rpc1dca < 100) && (Rpc3dca < 100)) {
        rpcdca_condition = 2;
      } else if(Rpc1dca < 100) {
        rpcdca_condition = 0;
        //continue;
      } else if(Rpc3dca < 100) {
        rpcdca_condition = 1;
        //continue;
      }

      int goodfvtx = 0;
      if(-1.5 < fvtx_dphi && fvtx_dphi < 1.5 &&
          0 < fvtx_dtheta && fvtx_dtheta < 1.5  &&
          0 < fvtx_dr && fvtx_dr < 100) {
        bool fvtx_outofbounds = false;
        if(-1.5 > fvtx_dphi || fvtx_dphi > distmax[19]) {
          printf("fvtx_dphi: %f  ",fvtx_dphi);
          fvtx_outofbounds = true;
        }
        if(distmin[20] > fvtx_dr || fvtx_dr > distmax[20]) {
          printf("fvtx_dr: %f  ",fvtx_dr);
          fvtx_outofbounds = true;
        }
        if(distmin[21] > fvtx_dtheta || fvtx_dtheta > distmax[21]) {
          printf("fvtx_dtheta: %f  ",fvtx_dtheta);
          fvtx_outofbounds = true;
        }
        //if(distmin[17] > fvtx_dr_dtheta || fvtx_dr_dtheta > distmax[17]) {
          //printf("fvtx_dr_dtheta: %f",fvtx_dr_dtheta);
          //fvtx_outofbounds = true;
       // }

        if(fvtx_outofbounds) {
          printf("\n");
          continue;
        }

        goodfvtx = 1;

      } else {
        goodfvtx = 0;
      }

      int wness_section;
      if(Wness==1.0) wness_section = 9;
      else wness_section = floor(Wness*10);

      float trig_scaling_factor;
      if(file_index==0)
        trig_scaling_factor = 1;
      else
        trig_scaling_factor = get_eta_trigeff(eta,trig_eff_factor[arm][charge_index]);
      

      
      if(wness_section==9) dw23_vs_eta->Fill(eta,dw23,trig_scaling_factor);
      h_wness[file_index][arm][charge_index]->Fill(Wness,trig_scaling_factor);
      h2_dw23_vs_wness[file_index][arm][charge_index]->Fill(Wness,dw23,trig_scaling_factor);
      
      h_dw23_sections[file_index][arm][charge_index][wness_section]->Fill(dw23,trig_scaling_factor);
      h_wness_sections[file_index][arm][charge_index][wness_section]->Fill(Wness,trig_scaling_factor);
      h_eta[file_index][arm][charge_index][wness_section]->Fill(eta,trig_scaling_factor);
      h2_dw23_vs_eta[file_index][arm][charge_index][wness_section]->Fill(eta,dw23,trig_scaling_factor);
      
      if(Wness>0 && Wness<1)
        in_count[file_index][arm][charge_index]++;
      else if(Wness>=1) {
        over_count[file_index][arm][charge_index]++;
      } else if(Wness <= 0) {
        under_count[file_index][arm][charge_index]++;
      }


    } //end main event loop
    printf("test1\n");
    printf("Num bad etas file %d: %d out of %d\n", file_index,bad_eta_count,entries); bad_eta_count = 0;
  }

}
