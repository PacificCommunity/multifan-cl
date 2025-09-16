/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include <admodel.h>

dmatrix mod_gram_schmidt(dmatrix & M)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int n=M(mmin).indexmax();

  for (int i=mmin;i<=mmax;i++)
  {
    if (norm(M(i))<1.e-3)
    {
      M(i)=0;
    }
    else
    {
      M(i)/=norm(M(i));
    }
    for (int j=i+1;j<=mmax;j++)
    {
      M(j)-= M(j)*M(i)*M(i);
    }
  }
  int ii=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (norm(M(i))>1.e-6) ii++;
  }
//  cout << "There are "  << ii << " nonzero rows "
//       << " and " << n << " columns " << endl;
  if (ii != n)
  {
    cerr << "Error in movement gram schmidt. Not enough linearly independent"
         << endl << "vectors" << endl;
    ad_exit(1);
  }

  dmatrix MM(1,n,1,n);
  ii=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (norm(M(i))>1.e-6) 
    {
      ii++;
      MM(ii)=M(i);
    }
  }
  return MM;
}

dmatrix dvar_fish_stock_history::get_n_season_diffusion_matrix(void)
{
  int nslots=sum(Dflags);  //NMD_Jan14-19
  if (nslots==0)
  {
    dmatrix MM;
    return MM;
  }
  else
  {
    //int nslots=max(1,sum(Dflags));
    int nmp=xdiff_coffs.indexmax();
    int tnslots=nslots*nmp;
    int ns=1;
    int nrr=num_regions;
    if (pmsd)
    {
      ns=pmsd->num_species;
      nrr=pmsd->num_real_regions;
    }
    dmatrix M(1,3*tnslots,1,tnslots);
    M.initialize();  //NMD 11Feb2019

    ivector slots_by_species(1,ns);
    if (pmsd)
    {
      for (int is=1;is<=ns;is++)
      {
        int offset=pmsd->rowoffset(is);
        slots_by_species(is)=sum(Dflags.sub(1+offset,nrr+offset));   
      }
    }
    else
    {
      slots_by_species(1)=nslots;
    }
//    cout << slots_by_species << endl;
  
    int rowoffset=1;
    int coloffset=0;
    M(rowoffset)=1.0;   // overall mean
    rowoffset++;
    for (int is=1;is<ns;is++)
    {
      M(rowoffset).sub(1+coloffset,nmp*slots_by_species(is))  // species effect
         =1.0;
      rowoffset++;
      coloffset+=nmp*slots_by_species(is);
    }

    coloffset=0;
    for (int is=1;is<=ns;is++)
    {
      for (int ip=1;ip<=nmp;ip++)
      {
        if (ip<nmp)
        {
          M(rowoffset)(1+coloffset,slots_by_species(is)+coloffset) 
            =1.0;                               // movement period effect
          rowoffset++;
        }
        coloffset+=slots_by_species(is);
      }
    }
    for (int i=1;i<=nslots*nmp;i++)
    {
      M(rowoffset+i-1,i)=1.0;
    }
    dmatrix MM=mod_gram_schmidt(M);
    return MM;
  }
}
  


