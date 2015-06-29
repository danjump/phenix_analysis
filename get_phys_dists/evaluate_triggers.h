bool didrpc3bctrigger ( int bit){
  
  bool ret = ( bit & 0x04000000 )
  || ( bit & 0x08000000 ) ; 
  //MUON_S_SG1_RPC3_1_B||C 
  //MUON_N_SG1_RPC3_1_B||C
  
  return ret;
}                       



bool didrpc3bcStrigger ( int bit){                 
  bool ret = ( bit & 0x04000000 ) ; 
  //MUON_S_SG1_RPC3_1_B||C 
  return ret;
}                       

bool didrpc3bcNtrigger ( int bit){
  bool ret = ( bit & 0x08000000 ) ; 
  //MUON_N_SG1_RPC3_1_B||C                  
  return ret;
}                       


 


bool didrpc3atrigger ( int bit){                 
  bool ret = ( bit &  0x00400000 ) 
  || ( bit & 0x00800000 );  
  //MUON_S_SG1_RPC3A&MUID_S1D
  //MUON_N_SG1_RPC3A&MUID_N1D
  return ret;
}
bool didrpc3aStrigger ( int bit){                 
  bool ret = ( bit &  0x00400000 ) ;
  //MUON_S_SG1_RPC3A&MUID_S1D
  return ret;
}
bool didrpc3aNtrigger ( int bit){                 
  bool ret =  ( bit & 0x00800000 );  
  //MUON_N_SG1_RPC3A&MUID_N1D
  return ret;
}

bool didrpc1trigger ( int bit){
  
  bool ret = bit & 0x00200000;  //SG1+RPC1(C)&MUIDLL1_N||S
  return ret;
}                       

bool didsg1trigger ( int bit){                 
  bool ret = bit &  0x01000000  //MUON_S_SG1&BBCLL1(noVtx) 
  || bit &  0x02000000; //MUON_N_SG1&BBCLL1(noVtx)                  
  return ret;
}                       
bool didsg1Strigger ( int bit){                 
  bool ret = bit &  0x01000000;  //MUON_S_SG1&BBCLL1(noVtx) 
  return ret;
}                       
bool didsg1Ntrigger ( int bit){                 
  bool ret = bit &  0x02000000; //MUON_N_SG1&BBCLL1(noVtx)  
  return ret;
}                       


bool didotherNtrigger ( int bit){                 
  bool ret = (  ! ( bit &  0x02000000 ) &&
  ! ( bit &  0x00200000 ) && 
  ! ( bit &  0x00800000 ) &&
  ! ( bit &  0x08000000 )) ; 
  //no rpc or sg1 trigger 
  return ret;
}                       


bool didotherStrigger ( int bit){                 
  bool ret = (  ! ( bit &  0x01000000 ) &&
  ! ( bit &  0x00200000 ) && 
  ! ( bit &  0x00400000 ) &&
  ! ( bit &  0x04000000 ) ); 
  //no rpc or sg1 trigger 
  return ret;
}       
