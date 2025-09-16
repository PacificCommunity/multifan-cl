/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include "scbd.hpp"

int get_true_date(int i,int j,dvar_fish_stock_history& fsh);

  void read_selectivity_blocks(dvar_len_fish_stock_history& fsh,
    ivector& ff71)
  {
    ivector ff57=column(fsh.fish_flags,57);
    ivector ff61=column(fsh.fish_flags,66);
    ivector ff74=column(fsh.fish_flags,74);
    ivector ff3=column(fsh.fish_flags,3);
    ivector ff16=column(fsh.fish_flags,16);
    if(!allocated(fsh.block_ptr))
      fsh.block_ptr.allocate(1,fsh.num_fisheries);
    
    fsh.block_ptr.initialize();
    cifstream cif("selblocks.dat");
    if (!cif)
    {
      cerr << "Error trying to open selblocks.dat " << endl;
    }
    int ii=0;
    for (int i=1;i<=fsh.num_fisheries;i++)
    {
      if (ff71(i)>0) 
      {
        ii++;
        fsh.block_ptr(i)=ii;
      }
    }

    int mmin=1;
    int mmax=fsh.num_fisheries;
    ivector group=column(fsh.fish_flags,24);
    int sg=sum(group);
    if (sg)
    {
      ivector key(mmin,mmax);
      sort(group,key);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (ff71(key(i)) !=ff71(key(i-1)))
          {
            cerr << "unequal block flags ff71 for grouped fisheries " 
                 << key(i-1) << " and " << key(i) << endl;
            ad_exit(1);
          }
        }
      }
    }

    int num_blocked_fisheries=ii;
    fsh.sel_block_index.allocate(1,num_blocked_fisheries);
    fsh.sel_block_fisheries.allocate(1,num_blocked_fisheries);
    ii=0;
    ivector itmp(1,num_blocked_fisheries);
    
    for (int i=1;i<=fsh.num_fisheries;i++)
    {
      if (ff71(i)>0) 
      {
        ii++;
        // sel_block_index(i) is the number of breaks in blocked fisheries
        // sel_block_fisheries(i) is the i'th blocked fishery
        fsh.sel_block_index(ii)=ff71(i);
        fsh.sel_block_fisheries(ii)=i;
        switch (ff57(i))
        {
        case 0:
          if (ff3(i)==0)
          {
            fsh.pmsd_error();
            itmp(ii)=fsh.nage-1;
          }
          else
          {
	    itmp(ii)=ff3(i);
          }
          break;
        case 1:
	  itmp(ii)=2;
          break;
        case 2:
	  itmp(ii)=3;
          break;
        case 3:
          int num_nodes=ff61(i);
	  itmp(ii)=num_nodes;
          break;
        }
      }
    }
    fsh.sel_block_breaks.allocate(1,num_blocked_fisheries,
     1,fsh.sel_block_index);

    //fsh.blocked_sel_dev_coffs.allocate(1,num_blocked_fisheries,
      //1,fsh.sel_block_index,1,itmp);
    fsh.num_blocks.allocate(1,fsh.num_fisheries);
    fsh.num_breaks.allocate(1,fsh.num_fisheries);
    fsh.num_blocks=1;
    fsh.num_breaks=0;
    fsh.num_blocked_seasons.allocate(1,num_blocked_fisheries);
    for (int i=1;i<=num_blocked_fisheries;i++)
    {
      fsh.num_blocks(fsh.sel_block_fisheries(i))=fsh.sel_block_index(i)+1;
      fsh.num_breaks(fsh.sel_block_fisheries(i))=fsh.sel_block_index(i);
      fsh.num_blocked_seasons(i)=ff74(fsh.sel_block_fisheries(i));
    }

    fsh.better_sbb.allocate(1,fsh.num_fisheries,1,ff71+1);
    if (sum(ff71))
    {
      cif >> fsh.sel_block_breaks;
      if (!cif)
      {
        cerr << "Error reading sel_block_breaks from file selblocks.dat " 
             << endl;
        ad_exit(1);
      }

      int ic=0;
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        if (ff71(i)>0)
        {
          fsh.better_sbb(i)(1,ff71(i))=fsh.sel_block_breaks(++ic);
        }
        fsh.better_sbb(i)(ff71(i)+1)=fsh.nyears+fsh.year1+1;
      }
    }
    else
    {
      fsh.better_sbb=fsh.nyears+fsh.year1+1;
    }

    if(!allocated(fsh.selblks_ff16))
      fsh.selblks_ff16.allocate(1,fsh.num_fisheries);
    if(!allocated(fsh.selseas_ff16))
      fsh.selseas_ff16.allocate(1,fsh.num_fisheries);

    int chkff16=0;
    for (int i=1;i<=fsh.num_fisheries;i++)
    {
      if (ff16(i)==1 && ((ff71(i)!=0) || (ff74(i)>1)))
      {
        chkff16=1;
        break;
      }
    }
    if (chkff16)
    {
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        fsh.selblks_ff16(i).allocate(1,fsh.num_blocks(i));
        fsh.selseas_ff16(i).allocate(1,ff74(i));
      }
      fsh.selblks_ff16.initialize();
      fsh.selseas_ff16.initialize();
      cifstream cif1("selblocks_ff16.dat");
      if (!cif1)
      {
        cerr << "Error trying to open selblocks_ff16.dat " << endl;
        ad_exit(1);
      }
      cifstream cif2("selseas_ff16.dat");
      if (!cif2)
      {
        cerr << "Error trying to open selseas_ff16.dat " << endl;
        ad_exit(1);
      }
      cif1 >> fsh.selblks_ff16;
      if (!cif1)
      {
        cerr << "Error input of selblocks_ff16.dat" << endl;
        ad_exit(1);
      }
      cif2 >> fsh.selseas_ff16;
      if (!cif2)
      {
        cerr << "Error input of selseas_ff16.dat" << endl;
        ad_exit(1);
      }
// check sanity of fish_flags(16) inputs
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        int tmp_blks=sum(fsh.selblks_ff16(i));
        int tmp_seas=sum(fsh.selseas_ff16(i));
        if (ff16(i)==1) 
        {
          if ((tmp_blks<= 0) || (tmp_seas<= 0))
          {
            cerr << "Input error time-variant selectivity"
                 << "  and fish_flags(16) settings, fishery: " << i  << endl;
            ad_exit(1);
          }
        }
        else
        {
          if ((tmp_blks!= 0) || (tmp_seas!= 0))
          {
            cerr << "Input error time-variant selectivity"
                 << "  and fish_flags(16) settings, fishery: " << i  << endl;
            ad_exit(1);
          }
        }
      }
    }
    else
    {
      fsh.selblks_ff16.deallocate();
      fsh.selseas_ff16.deallocate();
    }
  }


void dvar_fish_stock_history::add_blocked_sel_devs(void)
{
  int nb=sel_block_index.indexmax();
  for (int i=1;i<=nb;i++)
  {
    int fi=sel_block_fisheries(i);
    int num_ages=min(fish_flags(fi,72),nage);
    for (int ii=1;ii<=num_fish_times(fi);ii++)
    {
      int bj=break_block(i,ii);
      int rr=realization_region(fi,ii);
      int rp=realization_period(fi,ii);
      int ri=realization_incident(fi,ii);
      //if (break_year<year(rr,rp)) break;
      //if (len_sample_size(rr,rp,ri)>0 && bj>0)
      if (bj>0)
      {
        //if (num_ages==nage)
          fish_mort(rr,rp,ri)+=sel_dev_coffs(fi,bj);
        //else
          //fish_mort(rr,rp,ri)(1,num_ages)+=sel_dev_coffs(fi,bj)(1,num_ages);
      }
    }
  }
}

void dvar_fish_stock_history::set_selectivity_breaks(void)
{
  cerr << "Too bad" << endl;
  ad_exit(1);
  int nb=sel_block_index.indexmax();
  break_block.allocate(1,nb);
  blocked_fishery_index.allocate(1,num_fisheries);
  blocked_fishery_index.initialize();
  for (int i=1;i<=nb;i++)
  {
    int fi=sel_block_fisheries(i);
    blocked_fishery_index(i)=fi;
    break_block(i).allocate(1,num_fish_times(fi));
    break_block(i).initialize();
    int rr=realization_region(fi,1);
    int jj=1;
    int br=0;
    int ub=sel_block_breaks(i).indexmax();
    int ii;
    for (ii=1;ii<=num_fish_times(fi);ii++)
    {
      int ip=realization_period(fi,ii);
      int yr=really_true_year(rr,ip)+year1-1;
      if (sel_block_breaks(i,jj)<=yr )
      {
        if (jj==ub)
          break;
        jj++;
      }
      break_block(i,ii)=jj-1;
    }    
    for (int i2=ii;i2<=num_fish_times(fi);i2++)
    {
      break_block(i,i2)=jj;
    }
  }
}

void dvar_len_fish_stock_history::new_set_selectivity_breaks(void)
{
  ivector ff71=column(fish_flags,71);
  ivector ff74=column(fish_flags,74);
  ivector ff75=column(fish_flags,75);
  ivector ff26=column(fish_flags,26);
  ivector ff57=column(fish_flags,57);
  ff263flag=check(column(fish_flags,26),3);
  ivector size_stuff=calculate_blocked_sel_sizes();
  ivector nd=column(fish_flags,3);
  for (int i=1;i<=num_fisheries;i++)
  {
    if (nd(i)==0) nd(i)=nage;
  }
  //dvar_matrix tlength(1,num_fisheries,1,nd);
  ivector tmplength(1,num_fisheries);
  ivector tmpweight(1,num_fisheries);
  for (int i=1;i<=num_fisheries;i++)
  {
    if (ff26(i)==3)
    {
      //cout << "YYY Need to fix this" << endl;
      //tmplength(i)=min(nlint,nage-1);
      tmplength(i)=nlint;
      tmpweight(i)=nwint;
    }
    else
    {  
      tmplength(i)=nd(i);
      tmpweight(i)=0;
    }  
  }
  for (int i=1;i<=num_fisheries;i++)
  {
    // the number of distict slots being parameterized
    int nslots=tmplength(i)-ff75(i);
    // the number of parameters for the parameterization
    int nparams=bs_selcoff(i,1,1).indexmax();
    switch(ff57(i))
    {
    case 0:
      if (nparams != nslots)
      {
        cerr << "Maybe this shouldn;t happen! "
             << " PARAMETER NUMBER MISMATCH nlsot mismatch for fishery "
             << i << endl << "Contact your supplier " << endl;
        ad_exit(1);
      }
      break;
    case 1:
    case 2:
    case 3:
      if (nparams >= nslots)
      {
        cerr << "Too many parameters for the selecitivity parameterization"
             << endl << "for fishery " << i << " selectivity (ff57) = " 
             << ff57(i) << endl;
        ad_exit(1);
      }
      break;
     default:
      break;
    }
  }  
  if (!ff263flag)
  {
    bstempsel.allocate(1,num_fisheries,1,ff74,1,ff71+1,ff75+1,tmplength);
  }
  else
  {
    bstempsel.allocate(1,num_fisheries,1,ff74,1,ff71+1,1,tmplength);
    bswtempsel.allocate(1,num_fisheries,1,ff74,1,ff71+1,1,tmpweight);
    bstempsel_afl.allocate(1,num_fisheries,1,ff74,1,ff71+1,1,nage);
    bswtempsel.initialize(); //NMD_25apr2025
    bstempsel_afl.initialize(); //NMD_25apr2025
  }
  bstempsel.initialize(); //NMD_25apr2025
  yearblock.allocate(1,num_fisheries,1,nyears);
  yearblock.initialize();
  for (int i=1;i<=num_fisheries;i++)
  {
    for (int j=1;j<=num_fish_times(i);j++)
    {
      int rr=realization_region(i,j);
      int rp=realization_period(i,j);
      int ri=realization_incident(i,j);
      int fi=parent(rr,rp,ri);
      int rtm=(really_true_month(rr,rp));    
      switch(months_used(fi,0))
      {
      case 1:
        {
          int fish_per_year=months_used(fi).indexmax();
          int num_per_group=fish_per_year/ff74(fi);
          int div=fish_per_year%ff74(fi);
          if (div)
          {
            cerr << "bad ff74 value for fishery " << fi << endl;
            cerr << "value is " << ff74 << "  fisheries per year is "
                 << fish_per_year << endl;
            ad_exit(1);
          }
          //int mod=months_used(fi,num_per_group+1)-months_used(fi,1);
          //sseason(rr,rp,ri)=(rtm-months_used(fi,1))/mod+1;
          if (sel_pointer(fi,rtm)==0)
          {
            cerr << "This cannot happen" << endl;
            ad_exit(1);
          }
          sseason(rr,rp,ri)=sel_pointer(fi,rtm);
        }
        break;
      case -1:
      case 0:
        if (ff74(fi)>1)
        {
           cerr << "bad ff74 value for fishery " << fi << "  "
                << " it must be 1 but is " << ff74(fi) << endl;
           ad_exit(1);
        }
        sseason(rr,rp,ri)=1;
        break;
      default:  
        cerr << "bad months used indicator for fishery " << fi << "  "
             << " it must be a -1 or 0 1 but is " << months_used(fi,0)
              << endl;
        ad_exit(1);
      }
      int yr=really_true_year(rr,rp)+year1-1;
      int rtyr=really_true_year(rr,rp);
      int myr=year(rr,rp);  //NMD_7Nov2017
      if (rtyr>nyears)
         rtyr=nyears;
      {
        for (int ii=1;ii<=ff71(i)+1;ii++)
        {
          if (yr<better_sbb(i,ii))
          {
//            yearblock(i,rtyr)=ii;
            yearblock(i,myr)=ii;    //NMD_7Nov2017
            bblock(rr,rp,ri)=ii;
            break;
          }
        }
      }
    }
  }
//  cout << bblock << endl;
}
          

void dvar_fish_stock_history::setup_sel_ptr(void)
{
  ivector ff74=column(fish_flags,74);
  if (allocated(sel_pointer))
  {
    sel_pointer.deallocate();
  }
  sel_pointer.allocate(1,num_fisheries);
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    sel_pointer(fi)=setup_group_selectivity_pointer(months_used(fi),month_1,
      ff74(fi));
  }
}
  

ivector setup_group_selectivity_pointer(ivector& months_used,int month_1,
  int ff74)
{
  int gmax=months_used.indexmax();
  int first_month=0;
  int rem=gmax%ff74;
  if (rem) 
  {
    cerr << "ff74 must divide the number of groups in fishery" << endl;
    ad_exit(1);
  }
  int num_in_group=gmax/ff74;
//  cout << "rem = " << rem << " num_in_group = " << num_in_group << endl;
  for (int i=1;i<=gmax;i++)
  {
    if( months_used(i)>=month_1)
    {
      first_month=i;
      break;
    }
  }
  if (!first_month) first_month=1; //NMD_26jan2023 - no month>month_1
  ivector tmp(0,gmax);
  tmp=months_used;
  for (int i=1;i<=gmax;i++)
  {
    months_used(i) = tmp((i+first_month-2)%gmax+1);
  }
//  cout << setw(5) << months_used << endl;

  ivector pointer(1,12);
  pointer.initialize();
  
  rem=12%ff74;
  int spread=12/ff74;
  if (rem) spread++;
  int minf=months_used(1);
  // 8    10    12     2     5     6

  int gp=1;
  for (int i=1;i<=gmax;i++)
  {
    pointer(months_used(i))=gp;
    if (i%num_in_group==0) gp++;
  }
    
    
  return pointer;
  
}

