#include "W2eGetPol.h"
#include <odbc++/connection.h>
#include <odbc++/setup.h>
#include <odbc++/types.h>
#include <odbc++/errorhandler.h>
#include <sql.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <odbc++/preparedstatement.h>
#include <odbc++/databasemetadata.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>

using namespace odbc;
using namespace std;

W2eGetPol::W2eGetPol(){
 // spin_out.SetUserName("phnxrc");  // it seems this does not make any difference, Ralf is not using it
}

int W2eGetPol::InitPat(int run)
{
  
  runnumber=run;

  //spin_out.SetTableName("spin_oncal");
  //spin_out.StoreDBContent(run,run,65535);         
  spin_out.StoreDBContent(run,run);

  if(spin_out.CheckRunRowStore(run)!=1){return 5;}
  spin_out.GetDBContentStore(spin_cont,run);

  //xingshift=0;//spin_cont.GetCrossingShift(); // the correction for the shift is already applied in the funciotn:
						//GetSpinPatternBluePHENIX(i) and GetSpinPatternYellowPHENIX(i)
  //spin_cont.SetCrossingShift(0); questo toglie il xing shift

  xingshift=spin_cont.GetCrossingShift();

  qalevel = spin_cont.GetQALevel(); 

  //cout << "  Spin-database Runnumber " << run << " qalevel " << qalevel<<endl; 


  Pattern();
  Scaler(); 
  return 0;
}


int W2eGetPol::InitSc(int run)
{
  spin_out.SetTableName("spin_daq");
  spin_out.StoreDBContent(run,run,65535);         

  if(spin_out.CheckRunRowStore(run)!=1){return 5;}
  spin_out.GetDBContentStore(spin_cont,run);
  spin_cont.SetCrossingShift(0);

  fout.open("runs_bad_trig.lst",fstream::app);//this is a really bad way of doing things
  Scaler();
  fout.close();
  return 0;
}

void W2eGetPol::Pattern()
{
   double bpol,bpolerr,ypol,ypolerr;
   spin_cont.GetPolarizationBlue(0,bpol,bpolerr);
   spin_cont.GetPolarizationYellow(0,ypol,ypolerr);
   //spin_cont.GetCrossingShift();// dont use if PHENIX method 


  pol[0]=bpol;
  polerr[0]=bpolerr;
  pol[1]=ypol;
  polerr[1]=ypolerr;


  for(int i=0;i<120;i++)
    {
      //int pb=spin_cont.GetSpinPatternBluePHENIX(i);  //this function already apply the crossing shift 
      //int py=spin_cont.GetSpinPatternYellowPHENIX(i);//this function already apply the crossing shift
    
      int pb=spin_cont.GetSpinPatternBlue(i);  // this function need the correct xid
      int py=spin_cont.GetSpinPatternYellow(i);//
 
 
  /*    if     (pb== -1 && py== -1) pattern[i]=0;
      else if(pb==-1 && py== 1) pattern[i]=1;
      else if(pb== 1 && py==-1) pattern[i]=2;
      else if(pb== 1 && py==1) pattern[i]=3;
      else                      pattern[i]=4;
*/
      if     (pb== 1 && py== 1) pattern[i]=0;
      else if(pb==-1 && py== 1) pattern[i]=1;
      else if(pb== 1 && py==-1) pattern[i]=2;
      else if(pb== -1 && py==-1) pattern[i]=3;
      else                      pattern[i]=4;


      //if(pattern[i]==4)cout<<i<<":"<<pattern[i]<< " " << pb << " " << py<< endl;
    }
}


void W2eGetPol::Scaler()
{
//  long long GetScalerBbcVertexCut(int bunch);
//   long long GetScalerBbcNoCut(int bunch);
//   long long GetScalerZdcWide(int bunch);
//   long long GetScalerZdcNarrow(int bunch);
//   long long GetScaler(int channel,int bunch);

  for(int i=0;i<120;i++)
    {
      //scaler[pattern[(i+xingshift)%120]]+=spin_cont.GetScalerBbcVertexCut((i)%120)/10000;
      scaler[pattern[(i+xingshift)%120]]+=spin_cont.GetScalerBbcNoCut((i)%120)/10000;

      if(pattern[(i+xingshift)%120]==4 && spin_cont.GetScalerBbcNoCut((i)%120)/10000 > 2000 )
	{
	  fout<<runnumber<<" "<<(i+xingshift)%120<<" "<<spin_cont.GetScalerBbcVertexCut((i)%120)/10000<<endl;
	  //cout<<runnumber<<" "<<(i+xingshift)%120<<" "<<spin_cont.GetScalerBbcVertexCut((i)%120)/10000<<endl;
	}
    }
 }
