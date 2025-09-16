/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
void dvar_fish_stock_history::do_fragment_stuff(void)
{
  get_movement_pointer();

  int nregions=num_regions;
  dvar4_array A1(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A2(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A3(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A4(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A5(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A6(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A7(1,4,1,nage,1,nregions,1,nregions);
  dvar4_array A8(1,4,1,nage,1,nregions,1,nregions);

  
  A1=Dad;   // original Dad ?
  movement_coffs_option_0_to_1();
  assign_movement_coffs_option_1();
  A7=Dad;
  movement_coffs_option_1_to_0();
  assign_movement_coffs_option_0();
  A8=Dad;

  old_assign_movement_coffs_option_0();
  A5=Dad;
  assign_movement_coffs_option_0();
  A2=Dad;
  movement_coffs_option_0_to_2();
  assign_movement_coffs_option_0();
  A3=Dad;
  movement_coffs_option_2_to_0();
  assign_movement_coffs_option_0();
  A4=Dad;
  assign_movement_coffs_option_2();
  A6=Dad;
  dvar_matrix MM=A5(1,1)-A2(1,1);

  ofstream ofs("MM");
  {
    ofs << setw(7) << setfixed() << setprecision(3) << MM << endl;
    ofs << endl << "Original Dad(1,1)  " << endl;
    ofs << endl << setw(7) << setfixed() << setprecision(3) << inv(A1(1,1)) << endl;
    ofs << endl << setw(7) << setfixed() << setprecision(3) << A2(1,1) << endl;
    ofs << endl << setw(7) << setfixed() << setprecision(3) << A5(1,1) << endl;
    ofs << endl << "A3" << endl
        << setw(7) << setfixed() << setprecision(3) << A3(1,1) << endl;
    ofs << endl << "A4" << endl 
        << setw(7) << setfixed() << setprecision(3) << A4(1,1) << endl;
    ofs << endl << "A6" << endl 
        << setw(7) << setfixed() << setprecision(3) << A6(1,1) << endl;
    ofs << endl << "A7" << endl 
        << setw(7) << setfixed() << setprecision(3) << A7(1,1) << endl;
    ofs << endl << "A8" << endl 
        << setw(7) << setfixed() << setprecision(3) << A8(1,1) << endl;
  }

  for (int mp=1;mp<=4;mp++)
  {
    ofs << "Movement period " << mp << endl;
    for (int i=1;i<=num_regions;i++)
    { 
      ofs << i << " 1 " << norm2(inv(A1(mp,1))-A2(mp,1)) << endl;
      ofs << i << " 2 " << norm2(A3(mp,1)-A2(mp,1)) << endl;
      ofs << i << " 3 " << norm2(A4(mp,1)-A2(mp,1)) << endl;
      ofs << i << " 0 " << norm2(A5(mp,1)-A2(mp,1)) << endl;
    }
  }
 
  ad_exit(1);
}


void dvar_fish_stock_history::old_assign_movement_coffs_option_0(void)
{
  int mmin=1;
  int mmax=mo.num_periods();
  //int mmax=xdiff_coffs.indexmax();
  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  Dad.initialize();
  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    dvar_matrix td=Dad(k,1).sub(rmin,rmax).shift(1);
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }
      ii=1;
      for (int i=1;i<=nrr;i++)
      {
        td(i,i)=1.0;
        for (int j=1;j<=nrr;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            dvariable tmp=diff_coffs(k,ii);
            td(i,i)+=tmp;
            td(j,i)-=tmp;
            ii++;
          }
        }
      }
    }
  }
}

void dvar_fish_stock_history::get_movement_pointer(void)
{
  int mmin=1;
  int mmax=mo.num_periods();
  //int mmax=xdiff_coffs.indexmax();
  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  if (!allocated(movement_ptr))
  {
    movement_ptr.allocate(mmin,mmax,1,nrr,1,nrr);
  }
  movement_ptr.initialize();
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  //Dad.initialize();
  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }
      ii=1;
      for (int i=1;i<=nrr;i++)
      {
        for (int j=1;j<=nrr;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            movement_ptr(k,j,i)=ii;
            ii++;
           /*
            dvariable tmp=0.0;
            tmp=diff_coffs(k,ii);
            td(i,i)+=tmp;
            td(j,i)-=tmp;
            ii++;
           */
          }
        }
      }
    }
  }
}

void dvar_fish_stock_history::assign_movement_coffs_option_1(void)
{
  movement_coffs_option_1_to_0();
  assign_movement_coffs_option_0();
}
void dvar_fish_stock_history::assign_movement_coffs_option_2(void)
{
  movement_coffs_option_2_to_0();
  assign_movement_coffs_option_0();
}

void dvar_fish_stock_history::assign_movement_coffs_option_0(void)
{
  int mmin=1;
  int mmax=mo.num_periods();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  Dad.initialize();

  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }
      dvar_matrix td=Dad(k,1).sub(rmin,rmax).shift(1);

      for (int i=1;i<=nrr;i++)
      {
        td(i,i)=1;
        for (int j=1;j<i;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            dvariable tmp=diff_coffs(k,movement_ptr(k,j,i));
            dvariable tmp1=diff_coffs(k,movement_ptr(k,i,j));
            td(i,i)+=tmp;
            td(j,i)=-tmp;
            td(j,j)+=tmp1;
            td(i,j)=-tmp1;
          }
        }
      }
    }
  }
}



void dvar_fish_stock_history::movement_coffs_option_2_to_0(void)
{
  int mmin=1;
  int mmax=mo.num_periods();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  //Dad.initialize();

  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }

      for (int i=1;i<=nrr;i++)
      {
        for (int j=1;j<i;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            dvariable z1=zdiff_coffs(k,movement_ptr(k,j,i));
            dvariable z2=zdiff_coffs(k,movement_ptr(k,i,j));
            diff_coffs(k,movement_ptr(k,j,i))=exp(0.5*(z1+z2));
            diff_coffs(k,movement_ptr(k,i,j))=exp(0.5*(z1-z2));
          }
        }
      }
    }
  }
}
dvariable dvar_fish_stock_history::movement_coffs_option_2_penalty(void)
{
  dvariable fpen=0.0;
  int mmin=1;
  int mmax=mo.num_periods();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  //Dad.initialize();

  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }

      for (int i=1;i<=nrr;i++)
      {
        for (int j=1;j<i;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            dvariable z1=zdiff_coffs(k,movement_ptr(k,j,i));
            dvariable z2=zdiff_coffs(k,movement_ptr(k,i,j));
            fpen+=10.0*square(z2);
            fpen+=10.0*square(z1-1.0);
            //diff_coffs(k,movement_ptr(k,j,i))=exp(0.5*(z1+z2));
            //diff_coffs(k,movement_ptr(k,i,j))=exp(0.5*(z1-z2));
          }
        }
      }
    }
  }
  return fpen;
}

void dvar_fish_stock_history::movement_coffs_option_1_to_0(void)
{
  if (sum(Dflags)>0)   //NMD_jan14-19
  {
    int mmin=xdiff_coffs.indexmin();
    int mmax=xdiff_coffs.indexmax();
    int nrows=mmax-mmin+1;
    dvar_vector v1=lineup(xdiff_coffs);  // lineup the rows of diff_coffs into
    dvar_vector tmp1=exp(v1*new_orthogonal_diffusion_matrix);
    dvar_matrix M=stack(tmp1,nrows);
    diff_coffs=M;
  }
  else
  {
    if (allocated(diff_coffs))
    {
      diff_coffs=0.0;
    }   //NMD_jan14-19
  }
}
void dvar_fish_stock_history::rationalize_all_coffs(void)
{
  rationalize_movement_coffs(0);
}

void dvar_fish_stock_history::rationalize_movement_coffs(int af184)
{
  if (!allocated(movement_ptr))
  {
    get_movement_pointer();
  }
  switch(af184)
  {
  case(0):
    movement_coffs_option_0_to_1();
    movement_coffs_option_0_to_2();
    break;
  case(1):
    movement_coffs_option_1_to_0();
    movement_coffs_option_0_to_2();
    break;
  case(2):
    movement_coffs_option_2_to_0();
    movement_coffs_option_0_to_1();
    break;
  default:
    cerr << "error in historical_age_flags(184)" << endl;
    ad_exit(1);
  }
}

void dvar_fish_stock_history::movement_coffs_option_0_to_1(void)
{
  if (sum(Dflags)>0)   //NMD_jan14-19
  {

    int mmin=xdiff_coffs.indexmin();
    int mmax=xdiff_coffs.indexmax();
    int nrows=mmax-mmin+1;
    int cmin=xdiff_coffs(1).indexmin();
    int cmax=xdiff_coffs(1).indexmax();
    dvar_vector ldc=lineup(diff_coffs);  // lineup the rows of diff_coffs into
    dvar_vector tmp=new_orthogonal_diffusion_matrix*log(ldc);
    xdiff_coffs=stack(tmp,nrows);  // break up vector into nc pieces
  }
  else
  {
    if (allocated(xdiff_coffs))
    {
      xdiff_coffs=0.0;
    }   //NMD_jan14-19
  }
}

void dvar_fish_stock_history::movement_coffs_option_0_to_2(void)
{
  int mmin=1;
  int mmax=mo.num_periods();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  //Dad.initialize();

  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }

      for (int i=1;i<=nrr;i++)
      {
        for (int j=1;j<i;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            //movement_ptr(k,j,i)=ii;
            dvariable z1=diff_coffs(k,movement_ptr(k,j,i));
            dvariable z2=diff_coffs(k,movement_ptr(k,i,j));
            zdiff_coffs(k,movement_ptr(k,j,i))=log(z1)+log(z2);
            zdiff_coffs(k,movement_ptr(k,i,j))=log(z1)-log(z2);
            //z1=log(u1)+log(u2);
            //z2=log(u1)-log(u2);
           /*
            dvariable tmp=0.0;
            tmp=diff_coffs(k,ii);
            td(i,i)+=tmp;
            td(j,i)=-tmp;
            ii++;
           */
          }
        }
      }
    }
  }
}

void dvar_fish_stock_history::test_movement_pointer_2(void)
{
  int mmin=1;
  int mmax=mo.num_periods();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  //Dad.initialize();

  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }
      dvar_matrix td=Dad(k,1).sub(rmin,rmax).shift(1);

      for (int i=1;i<=nrr;i++)
      {
        td(i,i)=1;
        for (int j=1;j<i;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            //movement_ptr(k,j,i)=ii;
            dvariable z1=zdiff_coffs(k,movement_ptr(k,i,j));
            dvariable z2=zdiff_coffs(k,movement_ptr(k,j,i));
            td(j,i)=exp(0.5*(z1+z2));
            td(i,j)=exp(0.5*(z1-z2));
           /*
            dvariable tmp=0.0;
            tmp=diff_coffs(k,ii);
            td(i,i)+=tmp;
            td(j,i)-=tmp;
            ii++;
           */
          }
        }
      }
    }
  }
}
void dvar_fish_stock_history::test_movement_pointer_2_inv(void)
{
  int mmin=1;
  int mmax=mo.num_periods();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  //Dad(1,num_move_periods,1,nage,1,nregions,1,nregions),
  //Dad.initialize();

  int offset=0;
  int rmin=1;
  int rmax=num_regions;
  for (int k=mmin;k<=mmax;k++)   // loop over movement periods
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      if (pmsd)
      {
        offset=pmsd->rowoffset(is);
        rmin=pmsd->region_bounds(is,1);
        rmax=pmsd->region_bounds(is,2);
      }

      for (int i=1;i<=nrr;i++)
      {
        for (int j=1;j<i;j++)
        {
          if (Dflags(i+offset,j)) 
          {
            //movement_ptr(k,j,i)=ii;
            dvariable x1=xdiff_coffs(k,movement_ptr(k,i,j));
            dvariable x2=xdiff_coffs(k,movement_ptr(k,j,i));
            //td(j,i)=exp(0.5*(x1+x2));
            //td(i,j)=exp(0.5*(x1-x2));

            xdiff_coffs(movement_ptr(k,j,i))=
              log(diff_coffs(movement_ptr(k,i,j)))
             +log(diff_coffs(movement_ptr(k,j,i)));
            xdiff_coffs(movement_ptr(k,i,j))=
              log(diff_coffs(movement_ptr(k,j,i)))
             -log(diff_coffs(movement_ptr(k,i,j)));

           /*
            dvariable tmp=0.0;
            tmp=diff_coffs(k,ii);
            td(i,i)+=tmp;
            td(j,i)-=tmp;
            ii++;
           */
          }
        }
      }
    }
  }
}

