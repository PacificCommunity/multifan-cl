/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"

void make_stddev_report_noeff(dvar_len_fish_stock_history& fsh)
{
  int nvar;
  int nvar1;
  int nnvar1;
  int num_deps;
  int nnum_deps;
  int i;

  adstring current_path;
  adstring root;

  //getpaths(current_path,root);

  adstring fname= ad_root+".hes";

  uistream ifs(fname);
  if (!ifs)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  ifs >> nvar;
  if (!ifs)
  {
    cerr << "Error reading nvar from " << fname << endl;
  }
  dmatrix hessinv(1,nvar,1,nvar);
  dvector row(1,nvar);

  ifs >> hessinv;
  if (!ifs)
  {
    cerr << "Error reading hessian from " << fname << endl;
  }
  // check for zero row
  for (i=1;i<=nvar;i++)
  {
    for (int j=1;j<=nvar;j++)
    {
      if (hessinv(i,j)!=0) row(i)=1;
    }
  }

  for (i=1;i<=nvar;i++)
  {
    if (row(i)==0)
    {
      cout << " Error there was a zero row " << endl;
      cout << " in row " << i << endl;
      ad_exit(1);
    }
  }
  {
    hessinv=inv(hessinv);
  }

  fname=ad_root+".dep";
  uistream ifs1(fname);
  if (!ifs1)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  ifs1 >> nvar1 >> num_deps;
  if (!ifs1) {
    cerr << "Error reading nvar1 from depgrad.rpt" << endl;
  }
  if (nvar!=nvar1) {
    cerr << "number of independent variables is not consistent"
	    " between files" << "hessinv.rpt and depgrad.rpt" << endl;
    ad_exit(1);
  }
  dmatrix depgrad(1,num_deps,1,nvar);
  fname=ad_root+".var";
  ofstream ofs(fname);
  ifstream ifs2("deplabel.tmp");
  line_adstring deplabel;
  dvector depvalue(1,num_deps);

  adstring_array dstrings(1,num_deps);
  print_identifier_stuff(ofs);
  for (i=1;i<=num_deps;i++)
  {
    ifs1 >> depgrad(i);
    if (!ifs1) {
      cerr << "Error reading gradient for dependent variable " << i
	   << "  from depgrad.rpt" << endl;
      ad_exit(1);
    }
    ifs2 >> depvalue(i);
    ifs2 >> deplabel;
    dstrings(i)=deplabel;
    MY_DOUBLE_TYPE variance=depgrad(i)*hessinv*depgrad(i);
    ofs << setw(3) << i; 
    if (variance>=0)
    {
      ofs << "  "  << setw(8) << setprecision(3) << sqrt(variance);
    }
    else
    {
      cout << "variance = " << variance << endl;
      ofs << "  ********";
    }
    ofs << "  "  << setw(9) << setprecision(4) 
         << depvalue(i) << "  "  << deplabel << endl;
  }
  adstring nfname=ad_root+".dp2";
  uistream nifs1(nfname);
  if (!nifs1)
  {
    cerr << "Error trying to open file " << nfname << endl;
  }
  nifs1 >> nnvar1 >> nnum_deps;
  dvector ndepvalue(1,nnum_deps);
  if (!nifs1) {
    cerr << "Error reading nnvar from depgrad.rpt" << endl;
  }
  if (nvar!=nnvar1) {
    cerr << "number of independent variables is not consistent"
	    " between files" << "hessinv.rpt and depgrad.rpt" << endl;
    ad_exit(1);
  }
  ifstream nifs2("deplabel_noeff.tmp");
  dmatrix ndepgrad(1,nnum_deps,1,nvar);
  adstring_array ndstrings(1,nnum_deps);
  for (i=1;i<=nnum_deps;i++)
  {
    nifs1 >> ndepgrad(i);
    if (!nifs1) {
      cerr << "Error reading gradient for dependent variable " << i
	   << "  from depgrad.rpt" << endl;
      ad_exit(1);
    }
    nifs2 >> ndepvalue(i);
    nifs2 >> deplabel;
    ndstrings(i)=deplabel;
  }
  for (i=1;i<=fsh.nyears*2+fsh.num_regions*fsh.nyears*2;i++)
  {
    if (i>nnum_deps) break;
    if (i>num_deps) break;
    dvector gdiff=depgrad(i)-ndepgrad(i);
    MY_DOUBLE_TYPE variance=gdiff*hessinv*gdiff;
    ofs << setw(3) << i; 
    if (variance>=0)
    {
      ofs << "  "  << setw(8) << setprecision(3) << sqrt(variance);
    }
    else
    {
      cout << "variance = " << variance << endl;
      ofs << "  ********";
    }
    ofs << "  "  << setw(9) << setprecision(4) 
         << depvalue(i) - ndepvalue(i) << "  "  
         << dstrings(i) << " - " << ndstrings(i) << endl;
  }
}
