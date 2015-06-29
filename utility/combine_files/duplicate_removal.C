#include "DuplicateSearch.h"

int DuplicateSearch(std::string searchFileName="" , std::string uniqFileName="")
{
  std::cout << "Filtering out duplicate events from " << searchFileName << " and sending to " << uniqFileName << std::endl;
  TFile* searchFile = new TFile(searchFileName.c_str(),"READ");
  TTree* searchTree = (TTree*)searchFile->Get(wness_tree.c_str());

  Int_t
    Run_Number,  Evt_Number,  triggerbit,
    clockcross,  fvtx_cone,   fvtx_tracklcone;

  Float_t
    Evt_bbcZ,    Wness,
    charge,      pT,          px,           py,         pz,
    phi,         eta,         DG0,
    DDG0,        DG4,         chi2,
    DCA_z,       DCA_r,       dphi12,       dphi23,
    dw23,        Rpc1dca,     Rpc1time,     Rpc3dca,    Rpc3time,
    fvtx_dphi,   fvtx_dr,     fvtx_dtheta,
    fvtx_dr_dtheta, Rpc1x, Rpc1y, Rpc3x, Rpc3y;

  searchTree->SetBranchAddress("Run_Number",&Run_Number);
  searchTree->SetBranchAddress("Evt_Number",&Evt_Number);
  searchTree->SetBranchAddress("triggerbit",&triggerbit);
  searchTree->SetBranchAddress("Evt_bbcZ",&Evt_bbcZ);
  searchTree->SetBranchAddress("clockcross",&clockcross);
  searchTree->SetBranchAddress("Wness",&Wness);
  searchTree->SetBranchAddress("charge",&charge);
  searchTree->SetBranchAddress("pT",&pT);
  searchTree->SetBranchAddress("px",&px);
  searchTree->SetBranchAddress("pz",&py);
  searchTree->SetBranchAddress("pz",&pz);
  searchTree->SetBranchAddress("phi",&phi);
  searchTree->SetBranchAddress("eta",&eta);
  searchTree->SetBranchAddress("DG0",&DG0);
  searchTree->SetBranchAddress("DDG0",&DDG0);
  searchTree->SetBranchAddress("DG4",&DG4);
  searchTree->SetBranchAddress("chi2",&chi2);
  searchTree->SetBranchAddress("DCA_z",&DCA_z);
  searchTree->SetBranchAddress("DCA_r",&DCA_r);
  searchTree->SetBranchAddress("dphi12",&dphi12);
  searchTree->SetBranchAddress("dphi23",&dphi23);
  searchTree->SetBranchAddress("dw23",&dw23);
  searchTree->SetBranchAddress("Rpc1dca",&Rpc1dca);
  searchTree->SetBranchAddress("Rpc1time",&Rpc1time);
  searchTree->SetBranchAddress("Rpc1x",&Rpc1x);
  searchTree->SetBranchAddress("Rpc1y",&Rpc1y);
  searchTree->SetBranchAddress("Rpc3dca",&Rpc3dca);
  searchTree->SetBranchAddress("Rpc3time",&Rpc3time);
  searchTree->SetBranchAddress("Rpc3x",&Rpc3x);
  searchTree->SetBranchAddress("Rpc3y",&Rpc3y);
  searchTree->SetBranchAddress("fvtx_dphi",&fvtx_dphi);
  searchTree->SetBranchAddress("fvtx_dr",&fvtx_dr);
  searchTree->SetBranchAddress("fvtx_dtheta",&fvtx_dtheta);
  searchTree->SetBranchAddress("fvtx_dr_dtheta",&fvtx_dr_dtheta);
  searchTree->SetBranchAddress("fvtx_cone",&fvtx_cone);
  searchTree->SetBranchAddress("fvtx_tracklcone",&fvtx_tracklcone);
  
  TFile* uniqFile = new TFile(uniqFileName.c_str(),"RECREATE");
  TTree* uniqTree = new TTree(wness_tree.c_str(), "Wness Tree - Data");
  uniqTree->Branch("Run_Number",&Run_Number,"Run_Number/I");
  uniqTree->Branch("Evt_Number",&Evt_Number,"Evt_Number/I");
  uniqTree->Branch("triggerbit",&triggerbit,"triggerbit/I");
  uniqTree->Branch("Evt_bbcZ",&Evt_bbcZ,"Evt_bbcZ/F");
  uniqTree->Branch("clockcross",&clockcross,"clockcross/I");
  uniqTree->Branch("Wness",&Wness,"Wness/F");
  uniqTree->Branch("charge",&charge,"charge/F");
  uniqTree->Branch("pT",&pT,"pT/F");
  uniqTree->Branch("px",&px,"px/F");
  uniqTree->Branch("pz",&py,"pz/F");
  uniqTree->Branch("pz",&pz,"pz/F");
  uniqTree->Branch("phi",&phi,"phi/F");
  uniqTree->Branch("eta",&eta,"eta/F");
  uniqTree->Branch("DG0",&DG0,"DG0/F");
  uniqTree->Branch("DDG0",&DDG0,"DDG0/F");
  uniqTree->Branch("DG4",&DG4,"DG4/F");
  uniqTree->Branch("chi2",&chi2,"chi2/F");
  uniqTree->Branch("DCA_z",&DCA_z,"DCA_z/F");
  uniqTree->Branch("DCA_r",&DCA_r,"DCA_r/F");
  uniqTree->Branch("dphi12",&dphi12,"dphi12/F");
  uniqTree->Branch("dphi23",&dphi23,"dphi23/F");
  uniqTree->Branch("dw23",&dw23,"dw23/F");
  uniqTree->Branch("Rpc1dca",&Rpc1dca,"Rpc1dca/F");
  uniqTree->Branch("Rpc1time",&Rpc1time,"Rpc1time/F");
  uniqTree->Branch("Rpc1x",&Rpc1x,"Rpc1x/F");
  uniqTree->Branch("Rpc1y",&Rpc1y,"Rpc1y/F");
  uniqTree->Branch("Rpc3dca",&Rpc3dca,"Rpc3dca/F");
  uniqTree->Branch("Rpc3time",&Rpc3time,"Rpc3time/F");
  uniqTree->Branch("Rpc3x",&Rpc3x,"Rpc3x/F");
  uniqTree->Branch("Rpc3y",&Rpc3y,"Rpc3y/F");
  uniqTree->Branch("fvtx_dphi",&fvtx_dphi,"fvtx_dphi/F");
  uniqTree->Branch("fvtx_dr",&fvtx_dr,"fvtx_dr/F");
  uniqTree->Branch("fvtx_dtheta",&fvtx_dtheta,"fvtx_dtheta/F");
  uniqTree->Branch("fvtx_dr_dtheta",&fvtx_dr_dtheta,"fvtx_dr_dtheta/F");
  uniqTree->Branch("fvtx_cone",&fvtx_cone,"fvtx_cone/I");
  uniqTree->Branch("fvtx_tracklcone",&fvtx_tracklcone,"fvtx_tracklcone/I");
  
  std::map< int, std::set< int >  > evtMap;
  std::pair<std::set<int>::iterator,bool> ret; // insert returns pair - first being the iterator to the entry, second being whether or not entry was contained.
  std::map< int, std::map< int, std::vector<Long64_t> > > dupe_map;

  std::cout << "--Searching " << searchTree->GetName() << " for Duplicated Events--" << std::endl;
  long int entry_counter = 0;
  for(int i = 0; i < searchTree->GetEntries(); i++)
  {
    searchTree->GetEntry(i);
    ret = evtMap[Run_Number].insert(Evt_Number);
    dupe_map[Run_Number][Evt_Number].push_back(i);

    if( ret.second == false )
    {
      entry_counter++;
    }
  }

  std::map< int, std::map< int, std::vector<Long64_t> > >::iterator dupe_itr; 
  /*
  std::cout << "checking that events which are flagged as duplicates are REALLY duplicates" << std::endl;
  for(dupe_itr = dupe_map.begin(); dupe_itr != dupe_map.end(); ++dupe_itr)
  {
    std::map< int, std::vector<Long64_t> > evt_map = dupe_itr->second;
    std::map< int, std::vector<Long64_t> >::iterator evt_itr;
    for(evt_itr = evt_map.begin(); evt_itr != evt_map.end(); ++evt_itr)
    {
      if( (evt_itr->second).size() > 1 )
      {
        std::vector< Long64_t > evt_i = evt_itr->second;
        std::vector< Long64_t >::iterator evt_i_itr;
        searchTree->GetEntry(*(evt_i.begin()));
        double DG0_first = var.DG0;
        double DCA_first = var.DCA_r;
        for(evt_i_itr = evt_i.begin(); evt_i_itr != evt_i.end(); ++evt_i_itr)
        {
          searchTree->GetEntry(*evt_i_itr);
          double DG0_second = var.DG0;
          double DCA_second = var.DCA_r;
          if( fabs( abs(DG0_first - DG0_second) ) > 0.001 ) 
          {
            std::cout << "DGO IS NOT DUPLICATED FOR RUN[EVENT]: " << Run_Number << "[" << Evt_Number << "]" << std::endl;
          }
          if( fabs( abs(DCA_first - DCA_second) ) > 0.001 ) 
          {
            std::cout << "DCA_r IS NOT DUPLICATED FOR RUN[EVENT]: " << Run_Number << "[" << Evt_Number << "]" << std::endl;
          }
        }
      }
    }
  }
  */

  int runs_processed = 0;
  std::cout << "Sorting Tree For Faster Performance" << std::endl;
  // Sort the entries from the tree so it runs faster.
  std::set<Long64_t> entries; // fill with entry index we want to copy
  for(dupe_itr = dupe_map.begin(); dupe_itr != dupe_map.end(); ++dupe_itr)
  {
    runs_processed++;
    std::cout << "Progress: " << std::setprecision(4) << runs_processed << " runs               \r" << std::flush;
    std::map< int, std::vector<Long64_t> > evt_map = dupe_itr->second;
    std::map< int, std::vector<Long64_t> >::iterator evt_itr;
    for(evt_itr = evt_map.begin(); evt_itr != evt_map.end(); ++evt_itr)
    {
      std::vector< Long64_t > evt_i = evt_itr->second;
      entries.insert(evt_i[0]);
    }
  }
  std::cout << "We will now copy over " << entries.size() << " entries." << std::endl;
  
  std::cout << "Dumping uniqe tree from " << searchFileName << " to " << uniqFileName << std::endl;
  runs_processed = 0;
  for(std::set<Long64_t>::iterator entry_i = entries.begin(); entry_i != entries.end(); ++entry_i)
  {
    runs_processed++;
    std::cout << "Progress: " << *entry_i << " entries               \r" << std::flush;
    searchTree->GetEntry(*entry_i);
    uniqTree->Fill();
  }
  std::cout << entry_counter << " duplicates removed" << std::endl;
  searchFile->Close();
  delete searchFile;

  uniqFile->cd();
  uniqTree->Write();
  delete uniqTree;
  delete uniqFile;
	return 0;
}

