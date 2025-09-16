/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void corr(dmatrix& MM,adstring& s,dvar_fish_stock_history& fsh)
{
  int nsmp=MM.indexmax();
  dmatrix TM=trans(MM);
  int rmin=TM.indexmin();
  int rmax=TM.indexmax();
  dvector mn(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    TM(i)-=mean(TM(i));
  }
  
  dmatrix N(rmin,rmax,rmin,rmax);
  dmatrix C(rmin,rmax,rmin,rmax);
  C.initialize();
  N.initialize();
  for (int k=rmin;k<=rmax;k++)
  {
    for (int l=rmin;l<=k;l++)
    {
//      N(k,l)=TM(k)*TM(l)/(rmax-rmin+1);
      N(k,l)=TM(k)*TM(l)/(nsmp);  //NMD_jan18-19
      N(l,k)=N(k,l);
    }
  }
  for (int k=rmin;k<=rmax;k++)
  {
    for (int l=rmin;l<=k;l++)
    {
      C(k,l)=N(k,l)/(1.e-20+sqrt(N(k,k)*N(l,l)));
      C(l,k)=C(k,l);
    }
  }
  MY_DOUBLE_TYPE bound=0.995;
  if (fsh.parest_flags(194)/1000.>0)
  {
    bound=fsh.parest_flags(194)/1000.;
  }
  adstring ss= s + ".par_corr";
  ofstream ofs1(ss);
  for (int k=rmin;k<=rmax;k++)
  {
    for (int l=rmin;l<k;l++)
    {
      if (fabs(C(k,l))>bound)
      {
        ofs1 << k << " " << l << " " << C(k,l) << endl;
      }
    }
  }
}

//int main(int argc,char * argv[])
void correlation_report(adstring& s,dvar_fish_stock_history& fsh)
{
  adstring ss=s+".xbsamples";
  ifstream ifs(ss);
  int ncol = 0;
  int nrow = 0;
  string line;
  while (getline(ifs, line)) nrow++;
  cout << "Numbers of lines in the file : " << nrow << endl;
  ifs.clear();
  ifs.seekg(0);
  getline(ifs, line);
  istringstream iss(line);
  do
  {
      std::string sub;
      iss >> sub;
      if (sub.length())
          ++ncol;
  }
  while(iss);
  cout << "Numbers of columns in the file : " << ncol << endl;
  cout<< "nrows = " << nrow << " ncol = " << ncol << endl;
  dmatrix M(1,nrow,1,ncol);
  ifs.clear();
  ifs.seekg(0);
  ifs >> M;
  if (!ifs)
  {
   cerr << " ifs has problems" << endl;
   exit(1);
  }
  corr(M,s,fsh);
}
