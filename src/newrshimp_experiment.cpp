/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

#if !defined(ZFEV)
extern dvar_len_fish_stock_history * pcfsh;
#endif
extern dvar_vector * psv;


void check_sanity(ivector &t);

void check_sanity(ivector &t, dvar3_array& tg, ivector& rip,int it,
  dvar_fish_stock_history& fsh);
extern int sip;
extern adstring full_output_parfile_path;

static void set_xxx(int & ipold,int & ip)
{ 
  ipold=ip;
}

static void newxxx(void)
{
}
dvariable calculate_term(dvar_vector& y,dvector& x,dvector& I)
{
  dvariable alpha=((y-x)*I)/(I*I);
  return alpha;
}

dvariable dvar_len_fish_stock_history::implicit_catch_effort_relationship(void)
{
  dvariable tmp=old_build_implicit_catch_effort_design_matrix();
  return tmp;
}

group_flag_manager::group_flag_manager(ivector& group)
{
  int mmin=group.indexmin();
  int mmax=group.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (group(i)<0)
    {
      cerr << " group(" << i << ") <= 0 is a grouping flags error - exiting" 
           << endl;
      ad_exit(1);
    }
  }
  if (sum (group)==0)
  {
    numgroups=0;
  }
  else
  {
    ivector key(mmin,mmax);
    sort(group,key);
    check_group_for_holes(group);
    numgroups=max(group);
    group_members.allocate(1,numgroups);
    num_in_group.allocate(1,numgroups);
    int ir=0;
    int ic=1;
    for (int i=1;i<mmax;i++)
    {
      if (group(key(i+1))>group(key(i)))
      {
        ir++;
        num_in_group(ir)=ic;
        group_members(ir).allocate(1,num_in_group(ir));
        ic=1;
      }
      else
      {
        ic++;
      }
    }
    num_in_group(numgroups)=ic;
    group_members(numgroups).
      allocate(1,num_in_group(numgroups));

    int ii=0;
    for (int i=1;i<=numgroups;i++)
    {
      for (int j=1;j<=num_in_group(i);j++)
      {
        group_members(i,j)=key(++ii);
      }
    }
  }
}

void dvar_len_fish_stock_history::put_terminal_catchability(int ir,int ip)
{
  int nfi;
  nfi=num_fish_incidents(ir,ip);
  if (nfi>0)
  {
    for (int fi=1;fi<=nfi;fi++)
    {
      int i;
      i=parent(ir,ip,fi);
      if (eff_proj_fshry(i)>0 || catch_proj_fshry(i)>0)
      {
        int rmnth=really_true_month(ir,ip);
        for (int fnyri=1;fnyri<=fshtms_finyr(i);fnyri++)
        {
          int pmnth=ptr_fml_Mbs_mnth_finyr(i,fnyri);
          if (rmnth==pmnth)
          {
            MY_DOUBLE_TYPE leff=log_effort_by_fishery(i,fish_times(ir,ip,fi));
            q_level_finyr(i,fnyri)=value(log(fm_level(ir,ip,fi)))-leff;
            break;
          }
        }
      }
    }
  }
}


void dvar_len_fish_stock_history::do_row_bounds_count
  (int fi,int& irow)
{
  ivector ff92=column(fish_flags,92);
  ivector ff27=column(fish_flags,27);
  ivector ff73=column(fish_flags,73);
  //int icol;
  for (int nt=1;nt<=num_fish_times(fi);nt++)
  {
    int ir=realization_region(fi,nt);
    int ip=realization_period(fi,nt);
    int it=realization_incident(fi,nt);  
    //Need logic for number of columns
    
    if (ff92(fi)==0 && missing_catch_for_incident_flag(ir,ip,it)==0
     &&   missing_effort_by_region_flag(ir,ip,it)==0)
    {
      irow++;
    }
    else if (missing_catch_for_incident_flag(ir,ip,it))
    {
      //cout << "Missing obs number " << fi << "  " << nt << endl;
    }
    else if (missing_effort_by_region_flag(ir,ip,it))
    {
      //cout << "Missing effort " << fi << "  " << nt << endl;
    }
  }
}

void dvar_len_fish_stock_history::do_col_bounds_count(int fi,int& icol)
{
  ivector ff27=column(fish_flags,27);
  ivector ff73=column(fish_flags,73);
  //int icol;
  icol=1;
  if (ff27(fi)>0)
  {
    icol+=2;
  }  
  if (ff73(fi)>0)
  {
    icol+=ff73(fi);
  }  
}

void dvar_len_fish_stock_history::new_do_build_part_for_projections1(void)
{
  // Assignments to fml_designvec as needed for effort-conditioned
  // fishing incidents in the projection periods

  ivector nrft=num_real_fish_times;
  ivector nft=num_fish_times;
  fshtms_finyr.initialize();
  for (int fi=1;fi<=num_fisheries;fi++)
  {
/*
Steps:
- get fishery-specific pointer to final estimation period calendar year
- identify the quarters for which data exists and make pointers
- allocate ptr_fml_Mbs_mnth_finyr
- exclude CPUE index fisheries
 */
// - loop for obtaining pointers in final calendar year
    if (fish_flags(fi,92)==0)
    {
      int cnt_frst=0;
      int cnt_inc=0;
      int fin_cal_year;
      int ip_bfr;
      for (int nt=nrft(fi);nt>=nrft(fi)-100;nt--) //loop backwards over fish times
      {
        int ir=realization_region(fi,nt);
        int ip=realization_period(fi,nt);
        int it=realization_incident(fi,nt);
        if (cnt_frst==0)
        {
          if (age_flags(57) == 1)
          {
            fin_cal_year=year(ir,ip);
            if (fin_cal_year != last_real_year) //NMD_23Apr2024
            {
              cerr << "Fishery ends before last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              break;
            }
          }
          else
          {
            fin_cal_year=really_true_year(ir,ip);
            if (fin_cal_year != really_true_year(ir,last_real_year))
            {
              cerr << "Fishery ends before last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              break;
            }
          }
          cnt_frst=1;
        }
        if (missing_catch_for_incident_flag(ir,ip,it)==0
            &&   missing_effort_by_region_flag(ir,ip,it)==0)
        {
          cnt_inc++;
        }
        else if (missing_effort_by_region_flag(ir,ip,it)==0)
        {
          cnt_inc++;
        }

        ip_bfr=realization_period(fi,nt-1);
        int yr_cal;
        int yr_cal_bfr;
        if (age_flags(57) == 1)
        {
          yr_cal=year(ir,ip);
          yr_cal_bfr=year(ir,ip_bfr);
        }
        else
        {
          yr_cal=really_true_year(ir,ip);
          yr_cal_bfr=really_true_year(ir,ip_bfr);
        }
        if (yr_cal_bfr != yr_cal)
        {
          cnt_frst=0;
          break;
        }
      }
      fshtms_finyr(fi)=cnt_inc;
    }
  }

// - loop for checking projection fisheries data conditioned on effort or catch
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (!fish_flags(fi,92))
    {
      if (fshtms_finyr(fi) > 0)
      {
        for (int nt=nrft(fi)+1;nt<=nft(fi);nt++) //loop over projn fish times
        {
          int ir=realization_region(fi,nt);
          int ip=realization_period(fi,nt);
          int it=realization_incident(fi,nt);  
    
          if (missing_catch_for_incident_flag(ir,ip,it)==1
           &&   missing_effort_by_region_flag(ir,ip,it)==0)
          {
            eff_proj_fshry(fi)=1;
            if (fish_flags(fi,92) != 0)
            {
              cerr << "Error - effort-conditioned projection fishery: "
                   << fi << "  is a survey fishery. Exiting" << endl;
              ad_exit(1);
            }
            break;
          }
          else if (missing_catch_for_incident_flag(ir,ip,it)==0
           &&   missing_effort_by_region_flag(ir,ip,it)==1)
          {
            catch_proj_fshry(fi)=1;
          }
        }
      }
    }
  }
  if (!allocated(ptr_fml_Mbs_mnth_finyr))
    ptr_fml_Mbs_mnth_finyr.allocate(1,num_fisheries);

  for (int fi=1;fi<=num_fisheries;fi++)
  {    
    ptr_fml_Mbs_mnth_finyr(fi).deallocate();
    ptr_fml_Mbs_mnth_finyr(fi).allocate(1,fshtms_finyr(fi));
  }
  ptr_fml_Mbs_mnth_finyr.initialize();

  // Assign months to the pointer
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    int cnt_frst=0;
    int cnt_inc=0;
    int fin_cal_year;
    int ip_bfr;
    if (eff_proj_fshry(fi) || catch_proj_fshry(fi))
    {
      for (int nt=1;nt<=nrft(fi);nt++) //loop forwards over fish times
      {
        int ir=realization_region(fi,nt);
        int ip=realization_period(fi,nt);
        int it=realization_incident(fi,nt);

        if (age_flags(57) == 1)
        {
          fin_cal_year=year(ir,ip);
          if (fin_cal_year == last_real_year)
          {
            if (missing_catch_for_incident_flag(ir,ip,it)==0
              &&   missing_effort_by_region_flag(ir,ip,it)==0)
            {
              cnt_inc++;
              if (cnt_inc <= fshtms_finyr(fi))
              {
                ptr_fml_Mbs_mnth_finyr(fi,cnt_inc)=really_true_month(ir,ip);
              }
            }
          }
        }
        else
        {
          fin_cal_year=really_true_year(ir,ip);
          if (fin_cal_year == really_true_year(ir,last_real_year))
          {
            if (missing_catch_for_incident_flag(ir,ip,it)==0
              &&   missing_effort_by_region_flag(ir,ip,it)==0)
            {
              cnt_inc++;
              if (cnt_inc <= fshtms_finyr(fi))
              {
                ptr_fml_Mbs_mnth_finyr(fi,cnt_inc)=really_true_month(ir,ip);
              }
            }
          }
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::new_do_build_part_for_projections2(void)
{
  // Assignments to fml_designvec as needed for effort-conditioned
  // fishing incidents in the projection periods

  ivector nrft=num_real_fish_times;
  ivector nft=num_fish_times;

  /*
Steps:
- check corresponding months in projection years vs ptr_fml_Mbs_month_finyr
- assign fml_Mbs_finyr to fml_designvec
  */

  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (eff_proj_fshry(fi) == 1 || catch_proj_fshry(fi)==1)
    {
// - loop for assigning to fml_designvec     
      for (int nt=nrft(fi)+1;nt<=nft(fi);nt++) //loop over projection fish times
      {
        int ir=realization_region(fi,nt);
        int ip=realization_period(fi,nt);
        int it=realization_incident(fi,nt);  
        int rmnth=really_true_month(ir,ip);
        for (int fnyri=1;fnyri<=fshtms_finyr(fi);fnyri++)
        {
          int pmnth=ptr_fml_Mbs_mnth_finyr(fi,fnyri);
          if (rmnth==pmnth)
          {
            fml_designvec(ir,ip,it)=fml_Mbs_finyr(fi,fnyri);
            break;
          }
        }
      }
    }
  }
}



void dvar_len_fish_stock_history::remove_extra_columns_from_fml_designvec(int fi)
{
  int nft;     //NMD_27jun2022
  if (sum(data_fish_flags(2)))
  {
    nft=num_real_fish_times(fi);
  }
  else
  {
    nft=num_fish_times(fi);
  }     //NMD_27jun2022
//  for (int nt=1;nt<=num_fish_times(fi);nt++)
  for (int nt=1;nt<=nft;nt++)
  {
    int ir=realization_region(fi,nt);
    int ip=realization_period(fi,nt);
    int it=realization_incident(fi,nt);  
    
    if (missing_catch_for_incident_flag(ir,ip,it)==1
     &&   missing_effort_by_region_flag(ir,ip,it)==0)
    {
      dvector& v=fml_designvec(ir,ip,it);
      if (!allocated(v)) 
      {
        cerr  << "This can not hapen" << endl;
        ad_exit(1);
      }
     
      int nc=fml_columns(fi).indexmax();
      if (nc<v.indexmax())
      {
        dvector tmp(1,nc);
        for (int i=1;i<=nc;i++)
        {
          tmp(i)=v(fml_columns(fi,i));
        }
        fml_designvec(ir,ip,it).deallocate();
        fml_designvec(ir,ip,it)=tmp;
      }
    }
  }
}

void dvar_len_fish_stock_history::do_old_build_part(int fi,dvar_vector& O,
  int& irow)
{
  dvector leff=log_effort_by_fishery(fi);
  dvar_vector lfm=log(1.e-5+new_fm_level(fi));
  
  for (int nt=1;nt<=num_fish_times(fi);nt++)
  {
    int ir=realization_region(fi,nt);
    int ip=realization_period(fi,nt);
    int it=realization_incident(fi,nt);  
    if (missing_catch_for_incident_flag(ir,ip,it)==0
     &&   missing_effort_by_region_flag(ir,ip,it)==0)
    {
      irow++;
      O(irow)=lfm(nt)-leff(nt);
    }
    else if (missing_catch_for_incident_flag(ir,ip,it))
    {
      //cout << "Missing obs number " << nt << endl;
    }
    else if (missing_effort_by_region_flag(ir,ip,it))
    {
      //cout << "Missing effort number " << nt << endl;
    }
  }
}

void dvar_len_fish_stock_history::
  new_build_implicit_catch_effort_design_matrix(void)
{
  ivector ff29=column(fish_flags,29);
  ivector ff60=column(fish_flags,60);
  ivector ff27=column(fish_flags,27);
  ivector group=ff29; 
  dvariable fpen=0.0;
  if (sum(data_fish_flags(2)))
  {
    new_do_build_part_for_projections1();
  }

  implicit_fml_bounds=get_bounds_for_implicit_catch_effort_design_matrix();
  if (sum(group))
  {
    if (pfml_group)
    {
      cerr << "pfml_group already allocated. fix this." << endl;
      //ad_exit(1);
    }
    else
    {
      pfml_group=new group_flag_manager(group);
    }
    imatrix group_members=pfml_group->get_group_members();
    ivector num_in_group=pfml_group->get_num_in_group();
    int numgroups=pfml_group->get_numgroups();

    if (!allocated(fml_M))
      fml_M.allocate(1,numgroups);

    if (!allocated(fml_Q))
      fml_Q.allocate(1,numgroups);

    if (!allocated(fml_R))
      fml_R.allocate(1,numgroups);

    if (!allocated(fml_columns))
      fml_columns.allocate(1,num_fisheries);

    ivector nrows=implicit_fml_bounds(1);
    ivector ncols=implicit_fml_bounds(2);
    for (int ig=1;ig<=numgroups;ig++)
    {
      dmatrix tmp_M(1,nrows(ig),1,ncols(ig));
      tmp_M.initialize();
      if (allocated(tmp_M))
      {
        int irow=0;   //NMD_30aug2021
        int icol=0;
        for (int in=1;in<=num_in_group(ig);in++)
        {
          int fi=group_members(ig,in);
          new_do_build_part(fi,tmp_M,irow,icol);
        }            //NMD_30aug2021
        ivector columns=gram_schmidt_remove_extra_columns(tmp_M);
        for (int in=1;in<=num_in_group(ig);in++)
        {
          int fi=group_members(ig,in);
          fml_columns(fi)=columns;
        }
        int tmp_ncols=columns.indexmax();
        // Are there extra columns?
        if (tmp_ncols<ncols(ig))
        {
          ncols(ig)=tmp_ncols;
          for (int in=1;in<=num_in_group(ig);in++)
          {
            int fi=group_members(ig,in);
            implicit_fml_bounds(3,fi)=ncols(ig);
          }
          // remove redundant columns from  tmp_M
          dmatrix tmp(1,ncols(ig),1,nrows(ig));
          for (int i=1;i<=ncols(ig);i++)
          {
            tmp(i)=column(tmp_M,columns(i));
          }
          fml_M(ig)=trans(tmp);
        }
        else
        {
          fml_M(ig)=tmp_M;
        }
        if (!allocated(fml_Q(ig)))
        {
          fml_Q(ig).allocate(1,nrows(ig),1,tmp_ncols);
          fml_Q(ig).initialize();
        }
        if (!allocated(fml_R(ig)))
        {
          fml_R(ig).allocate(1,tmp_ncols,1,tmp_ncols);
          fml_R(ig).initialize();
        }
        // Do the QR dcomposition of the design matrix via
        // modified gram schmidt
        gram_schmidt_qr(fml_M(ig),fml_Q(ig),fml_R(ig));
      }
    }
  }
  else
  {
    if (!allocated(fml_M))
      fml_M.allocate(1,num_fisheries);
    fml_M.initialize();

    if (!allocated(fml_Q))
      fml_Q.allocate(1,num_fisheries);
    fml_Q.initialize();

    if (!allocated(fml_R))
      fml_R.allocate(1,num_fisheries);
    fml_R.initialize();

    if (!allocated(fml_columns))
      fml_columns.allocate(1,num_fisheries);

    ivector nrows=implicit_fml_bounds(1);
    ivector ncols=implicit_fml_bounds(2);

    for (int fi=1;fi<=num_fisheries;fi++)
    {
      //int nrows=implicit_fml_bounds(1,fi);
      //int ncols=implicit_fml_bounds(2,fi);
      dmatrix tmp_M(1,nrows(fi),1,ncols(fi));
      tmp_M.initialize();
      if (allocated(tmp_M))
      {
        int irow=0;
        int icol=0;
        new_do_build_part(fi,tmp_M,irow,icol);
        ivector columns=gram_schmidt_remove_extra_columns(tmp_M);
        fml_columns(fi)=columns;
        int tmp_ncols=columns.indexmax();
        if (tmp_ncols<ncols(fi))
        {
          ncols(fi)=tmp_ncols;
          // remove redundant columns from  tmp_M
//          dmatrix tmp(1,ncols(fi),1,nrows); //NMD_28mar2022
          dmatrix tmp(1,ncols(fi),1,nrows(fi));
          for (int i=1;i<=ncols(fi);i++)
          {
            tmp(i)=column(tmp_M,columns(i));
          }
          fml_M(fi)=trans(tmp);
        }
        else
        {
          fml_M(fi)=tmp_M;
        }
        if (!allocated(fml_Q(fi)))
        {
          fml_Q(fi).allocate(1,nrows(fi),1,ncols(fi));
          fml_Q(fi).initialize();
        }
        if (!allocated(fml_R(fi)))
        {
          fml_R(fi).allocate(1,ncols(fi),1,ncols(fi));
          fml_R(fi).initialize();
        }
        // Do the QR dcomposition of the design matrix via
        // modified gram schmidt
        gram_schmidt_qr(fml_M(fi),fml_Q(fi),fml_R(fi));
      }
    }
  }
  // remove redundant columns from fml_designvec(fi);
  //dvector& M=fml_designvec(ir,ip,it);
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    remove_extra_columns_from_fml_designvec(fi);
  }

  if (sum(data_fish_flags(2)))
  {
    new_do_build_part_for_projections2();
  }
}


void dvar_len_fish_stock_history::new_do_build_part(int fi,
  dmatrix& M,int& irow,int& icol)
{
  const MY_DOUBLE_TYPE twopi=2.0*3.1415926535897932;
  ivector ff27=column(fish_flags,27);
  ivector ff73=column(fish_flags,73);
  dvector leff=log_effort_by_fishery(fi);
  dvar_vector lfm=log(1.e-5+new_fm_level(fi));
  dvector vfmt;   //NMD_27jun2022
  int nft;
  if (sum(data_fish_flags(2)))
  {
    vfmt=fishing_incident_real_times(fi);
    nft=num_real_fish_times(fi);
  }
  else
  {
    vfmt=fishing_incident_times(fi);
    nft=num_fish_times(fi);
  }     //NMD_27jun2022
//  dvector vfmt=fishing_incident_times(fi);
  if (sum(data_fish_flags(2)))
  {
    if (!allocated(fml_Mbs_finyr(fi)))     //NMD_11jul2022
    {
      fml_Mbs_finyr(fi).allocate(1,fshtms_finyr(fi),
                                 1,implicit_fml_bounds(3,fi));
    }
    fml_Mbs_finyr(fi).initialize();     //NMD_11jul2022
  }
  int fin_cal_year;
  int cnt_inc=0;     //NMD_11jul2022
  MY_DOUBLE_TYPE fmtmax=max(vfmt);
  MY_DOUBLE_TYPE fmtmin=min(vfmt);

//  for (int nt=1;nt<=num_fish_times(fi);nt++)
  for (int nt=1;nt<=nft;nt++)
  {
    int ir=realization_region(fi,nt);
    int ip=realization_period(fi,nt);
    int it=realization_incident(fi,nt);  
    icol=0;
    
    if (missing_catch_for_incident_flag(ir,ip,it)==0
     &&   missing_effort_by_region_flag(ir,ip,it)==0)
    {
      MY_DOUBLE_TYPE fmt=vfmt(nt);
      irow++;
      icol++;
      if (icol<=M(irow).indexmax())
      {
        M(irow,icol)=1.0;
      }
      if (ff27(fi)>0)
      {
        icol++;
        if (icol<=M(irow).indexmax())
        {
          M(irow,icol)=
            sin(twopi*fmod(fmt,1));
        }
        icol++;
        if (icol<=M(irow).indexmax())
        {
          M(irow,icol)=
            cos(twopi*fmod(fmt,1));
        }
      }

      if (ff73(fi)>0)
      {
        icol++;
        MY_DOUBLE_TYPE u=(fmt-fmtmin)/(fmtmax-fmtmin);
        if (icol<=M(irow).indexmax())
        {
          M(irow,icol)=u;
        }
        for (int i=2;i<=ff73(fi);i++)
        {
          icol++;
          if (icol<=M(irow).indexmax())
          {
            M(irow,icol)=M(irow,icol-1)*u;
          }
        }
      }

      if (sum(data_fish_flags(2)))
      {
        if (age_flags(57) == 1)
        {
          fin_cal_year=year(ir,ip);
          if (eff_proj_fshry(fi) == 1 && catch_proj_fshry(fi)==1)
          {
            cerr << "This can't happen - a projection fishery being both"
                 << "catch- and effort-conditioned" << endl;
            ad_exit(1);
          }
          if ((eff_proj_fshry(fi) == 1 || catch_proj_fshry(fi)==1) &&
              fin_cal_year == last_real_year)
          {
            cnt_inc++;
            if (cnt_inc > fshtms_finyr(fi))
            {
              cerr << "Fishery time exceeds last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            if (ptr_fml_Mbs_mnth_finyr(fi,cnt_inc) != really_true_month(ir,ip))
            {
              cerr << "Finyr month not matching in last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            fml_Mbs_finyr(fi,cnt_inc)=M(irow);
          }
        }
        else
        {
          fin_cal_year=really_true_year(ir,ip);
          if (eff_proj_fshry(fi) == 1 && catch_proj_fshry(fi)==1)
          {
            cerr << "This can't happen - a projection fishery being both"
                 << "catch- and effort-conditioned" << endl;
            ad_exit(1);
          }
          if ((eff_proj_fshry(fi) == 1 || catch_proj_fshry(fi)==1) &&
              fin_cal_year == really_true_year(ir,last_real_year))
          {
            cnt_inc++;
            if (cnt_inc > fshtms_finyr(fi))
            {
              cerr << "Fishery time exceeds last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            if (ptr_fml_Mbs_mnth_finyr(fi,cnt_inc) != really_true_month(ir,ip))
            {
              cerr << "Finyr month not matching in last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            fml_Mbs_finyr(fi,cnt_inc)=M(irow);
          }
        }
      }
    }
    else if (missing_effort_by_region_flag(ir,ip,it)==0)
    {
      dvector& M=fml_designvec(ir,ip,it);
      if (!allocated(M)) M.allocate(1,implicit_fml_bounds(3,fi));
      //cout << "Missing obs number " << nt << endl;
      MY_DOUBLE_TYPE fmt=vfmt(nt);
      //x(irow)=nt;
      icol++;
      M(icol)=1.0;
      if (ff27(fi)>0)
      {
        icol++;
        M(icol)=
          sin(twopi*fmod(fmt,1));
        icol++;
        M(icol)=
          cos(twopi*fmod(fmt,1));
      }

      if (ff73(fi)>0)
      {
        icol++;
        MY_DOUBLE_TYPE u=(fmt-fmtmin)/(fmtmax-fmtmin);
        M(icol)=u;
        for (int i=2;i<=ff73(fi);i++)
        {
          icol++;
          M(icol)=M(icol-1)*u;
        }
      }
      
      if (sum(data_fish_flags(2)))
      {
        if (age_flags(57) == 1)
        {
          fin_cal_year=year(ir,ip);
          if (eff_proj_fshry(fi) == 1 && catch_proj_fshry(fi)==1)
          {
            cerr << "This can't happen - a projection fishery being both"
                 << "catch- and effort-conditioned" << endl;
            ad_exit(1);
          }
          if ((eff_proj_fshry(fi) == 1 || catch_proj_fshry(fi)==1) &&
              fin_cal_year == last_real_year)
          {
            cnt_inc++;
            if (cnt_inc > fshtms_finyr(fi))
            {
              cerr << "Fishery time exceeds last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            if (ptr_fml_Mbs_mnth_finyr(fi,cnt_inc) != really_true_month(ir,ip))
            {
              cerr << "Finyr month not matching in last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            fml_Mbs_finyr(fi,cnt_inc)=M;
          }
        }
        else
        {
          fin_cal_year=really_true_year(ir,ip);
          if (eff_proj_fshry(fi) == 1 && catch_proj_fshry(fi)==1)
          {
            cerr << "This can't happen - a projection fishery being both"
                 << "catch- and effort-conditioned" << endl;
            ad_exit(1);
          }
          if ((eff_proj_fshry(fi) == 1 || catch_proj_fshry(fi)==1) &&
              fin_cal_year == really_true_year(ir,last_real_year))
          {
            cnt_inc++;
            if (cnt_inc > fshtms_finyr(fi))
            {
              cerr << "Fishery time exceeds last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            if (ptr_fml_Mbs_mnth_finyr(fi,cnt_inc) != really_true_month(ir,ip))
            {
              cerr << "Finyr month not matching in last_real_year" << endl;
              cerr << "Fishery: " << fi << endl;
              ad_exit(1);
            }
            fml_Mbs_finyr(fi,cnt_inc)=M;
          }
        }
      }
    }
  }
}

dvariable dvar_len_fish_stock_history::
  old_build_implicit_catch_effort_design_matrix(void)
{
  ivector ff29=column(fish_flags,29);
  ivector ff60=column(fish_flags,60);
  ivector ff27=column(fish_flags,27);

  ivector group=ff29; 
  ivector num_ifmlrp=implicit_fml_bounds(2);
  dvariable fpen=0.0;
  dvector sigma;
  if (sum(group))
  {
    if (!pfml_group)
    {
      cerr << "pfml_group not allocated. fix this." << endl;
      ad_exit(1);
    }
    //group_flag_manager gfm(group);
    imatrix group_members=pfml_group->get_group_members();
    ivector num_in_group=pfml_group->get_num_in_group();
    int numgroups=pfml_group->get_numgroups();
    //cout << group_members << endl << endl;
    //cout << num_in_group << endl << endl;
    sigma.allocate(1,numgroups);
    for (int ig=1;ig<=numgroups;ig++)
    {
      int nrows=implicit_fml_bounds(1,ig);
      int ncols=implicit_fml_bounds(2,ig);
      if (ncols>0)
      {
        int irow=0;
        int icol=0;
        dvar_vector O(1,nrows);
        O.initialize();
        for (int in=1;in<=num_in_group(ig);in++)
        {
          int fi=group_members(ig,in);
          do_old_build_part(fi,O,irow);
        }
        // New ils_qr algorithm
        dvar_vector predicted;
        if (parest_flags(396))
        {
          dfils_manager dfilsm(fml_M(ig),O,2);
          //dvariable vr2=
          dfilsm.fit_data();
          predicted=fml_M(ig)*
            dfilsm.get_theta_hat()(1);
          implicit_fm_level_regression_pars
            (group_members(ig,1))(1,num_ifmlrp(ig))=
            dfilsm.get_theta_hat()(1);
        }
        else
        {
          dvar_vector ests=implicit_fm_level_regression_pars
            (group_members(ig,1))(1,num_ifmlrp(ig));
          predicted=fml_Q(ig)*ests;
        }
        
        dvar_vector diff=O-predicted;

        int icnt=1;  //NMD_20may2021
        for (int in=1;in<=num_in_group(ig);in++)
        {
          int fi=group_members(ig,in);
          dvar_vector lfm=log(1.e-20+new_fm_level(fi));
          for (int nt=1;nt<=num_fish_times(fi);nt++)
          {
            int ir=realization_region(fi,nt);
            int ip=realization_period(fi,nt);
            int it=realization_incident(fi,nt);
            if (missing_catch_for_incident_flag(ir,ip,it)==0
              &&   missing_effort_by_region_flag(ir,ip,it)==0)
            {
              log_pred_effort_by_fishery(fi,nt)=value(lfm(nt)-predicted(icnt));
              icnt++;
            }
          }
        }   //NMD_20may2021
  
        if (print_implicit_effort_flag && !af170q0)
        {
          dmatrix tmp(1,2,1,predicted.indexmax());
          adstring suffix=itoa(ig,10);
          tmp(1).fill_seqadd(1,1);
          tmp(2)=value(predicted);
          ofstream ofs1("pred"+suffix);
          ofs1 << trans(tmp) << endl;
          ofstream ofs2("obs"+suffix);
          tmp(2)=value(O);
          ofs2 << trans(tmp) << endl;
        } 
        int nobs=diff.indexmax();
        dvar_vector diff2=square(diff);
        const MY_DOUBLE_TYPE rsdcoff=sqrt(2.0/3.14159);
        switch (parest_flags(378))
        {
        case -1:
          // don't do anything
          break;
        case 0:
          {
            dvariable sighat=sum(sfabs(diff))/(nobs*rsdcoff);
            sigma(ig)=value(sighat);
            fpen+= 0.5*nobs*log(0.3+sighat);
          }
          break;
        case 1:
          {
            dvariable sighat=sum(sfabs(diff))/(nobs*rsdcoff);
            sigma(ig)=value(sighat);
            int fi=group_members(ig,1);
            dvariable a=exp(fish_pars(31,fi));
            fpen+=robust_normal_cauchy_mixture(a,sighat,diff2);
          }
          break;
        case 2:
          {
            dvariable vhat=sum(diff2)/nobs;
            dvariable sighat=sqrt(vhat);
            int fi=group_members(ig,1);
            dvariable a=exp(fish_pars(31,fi));
            if (parest_flags(362)>0)
            {
              sighat=1.0/parest_flags(362);
            }
            sigma(ig)=value(a*sighat);
            fpen+=robust_normal_mixture(a,sighat,diff2);
          }
          break;
        case 3:
          {
            int fi=group_members(ig,1);
            dvariable v=exp(fish_pars(32,fi));  // this will be 
                                          // the degrees of freedoom for
            dvariable sighat=0.025+sqrt(1.e-10+sum(diff2)/nobs);
            sigma(ig)=value(sighat);
            dvariable s=sqrt((v-2)/v)*sighat; 
            fpen+=sum(neg_log_student_density(v,s,diff));
          }   
          break;

        default:
          { 
            cout << "Illegal value for pf378 " << endl;
            ad_exit(1);
          }
        }
      }
    }
  }
  else
  {
    sigma.allocate(1,num_fisheries);
    for (int fi=1;fi<=num_fisheries;fi++)
    {
      adstring suffix=itoa(fi,10);
      int nrows=implicit_fml_bounds(1,fi);
      int ncols=implicit_fml_bounds(2,fi);
      if (ncols>0)
      {
        dvar_vector O(1,nrows);
        int irow=0;
        do_old_build_part(fi,O,irow);
      
        dvar_vector predicted;
        if (parest_flags(396))
        {
          dfils_manager dfilsm(fml_M(fi),O,2);
          dfilsm.fit_data();
          predicted=fml_M(fi)*
            dfilsm.get_theta_hat()(1);
          implicit_fm_level_regression_pars(fi)(1,num_ifmlrp(fi))=
            dfilsm.get_theta_hat()(1);
        }
        else
        {
          dvar_vector ests=implicit_fm_level_regression_pars
            (fi)(1,num_ifmlrp(fi));
          predicted=fml_Q(fi)*ests;
        }
        // ******************************************
        dvar_vector diff=O-predicted;

        dvar_vector lfm=log(1.e-05+new_fm_level(fi));  //NMD_20may2021
        int icnt=1;
        for (int nt=1;nt<=num_fish_times(fi);nt++)
        {
          int ir=realization_region(fi,nt);
          int ip=realization_period(fi,nt);
          int it=realization_incident(fi,nt);
          if (missing_catch_for_incident_flag(ir,ip,it)==0
            &&   missing_effort_by_region_flag(ir,ip,it)==0)
          {
            log_pred_effort_by_fishery(fi,nt)=value(lfm(nt)-predicted(icnt));
            icnt++;
          }
        }   //NMD_20may2021

        if (print_implicit_effort_flag && !af170q0)
        {
          dmatrix tmp(1,2,1,predicted.indexmax());
          adstring suffix=itoa(fi,10);
          tmp(1).fill_seqadd(1,1);
          tmp(2)=value(predicted);
          ofstream ofs1("pred"+suffix);
          ofs1 << trans(tmp) << endl;
          ofstream ofs2("obs"+suffix);
          tmp(2)=value(O);
          ofs2 << trans(tmp) << endl;
        } 
        // can get rid of this or put it in a reporting function
        int nobs=diff.indexmax();
        const MY_DOUBLE_TYPE rsdcoff=sqrt(2.0/3.14159);
        dvar_vector diff2=square(diff);

        switch (parest_flags(378))
        {
        case 0:
          {
            dvariable sighat=sum(sfabs(diff))/(nobs*rsdcoff);
            fpen+= 0.5*nobs*log(0.3+sighat);
            sigma(fi)=value(sighat);
          }
          break;
        case 1:
          {
            dvariable sighat=sum(sfabs(diff))/(nobs*rsdcoff);
            //fpen+= 0.5*nobs*log(0.3+sighat)+0.5*sum(diff2)/square(sighat);
            dvariable a=exp(fish_pars(31,fi));
            sigma(fi)=value(a*sighat);
            fpen+=robust_normal_cauchy_mixture(a,sighat,diff2);
          }
          break;
        case 2:
          {
            dvariable vhat=sum(diff2)/(nobs);
            dvariable sighat=sqrt(vhat);
            if (parest_flags(362)>0)
            {
              sighat=1.0/parest_flags(362);
            }
            dvariable a=exp(fish_pars(31,fi));
            sigma(fi)=value(a*sighat);
            fpen+=robust_normal_mixture(a,sighat,diff2);
          }
          break;
         
         // 888888888888888888888888888888888888888888888888888
         // 888888888888888888888888888888888888888888888888888
        case 3:
          {
            dvariable v=exp(fish_pars(32,fi));  // this will be 
                                          // the degrees of freedoom for
            dvariable sighat=0.025+sqrt(1.e-10+sum(diff2)/nobs);
            sigma(fi)=value(sighat);
            dvariable s=sqrt((v-2)/v)*sighat; 
            fpen+=sum(neg_log_student_density(v,s,diff));
          }   
          break;

         // 888888888888888888888888888888888888888888888888888
         // 888888888888888888888888888888888888888888888888888
        default:
          { 
            cout << "Illegal value for pf378 " << endl;
            ad_exit(1);
          }
        } 
      }
    }
  }
//  if (parest_flags(377))
  if (print_implicit_effort_flag && !af170q0)
  {
    dmatrix tmps(1,3,1,sigma.indexmax());
    tmps(1).fill_seqadd(1,1);
    tmps(2)=sigma;
    if (sum(group))
    {
      imatrix group_members=pfml_group->get_group_members();
      int numgroups=pfml_group->get_numgroups();
      for (int ig=1;ig<=numgroups;ig++)
      {
        int fi=group_members(ig,1);
        tmps(3,ig)=exp(value(fish_pars(31,fi)));
      }
    }
    else
    {
      tmps(3)=exp(value(fish_pars(31)));
    }	    
//    tmps(3)=exp(value(fish_pars(31)));
    ofstream ofss("sigma");
    ofss << trans(tmps) << endl;
  }
  return fpen;
}

dvariable robust_normal_cauchy_mixture(const prevariable&  a,
  const prevariable sighat,dvar_vector & r2)
{
  dvariable asig=a*sighat;
  dvariable avar=square(a*sighat);
  const MY_DOUBLE_TYPE c1=1.0/sqrt(2.0*3.141592653589793);
  const MY_DOUBLE_TYPE c2=1.0/3.141592653589793;
  int mmax=r2.indexmax();
  dvariable tmp=mmax*log(0.5+asig) -sum(log(
    0.95*c1*exp(-0.5*r2/avar)
   +0.05*c2/3.0/(1.0+r2/(9.0*avar))));
  return tmp;
}
dvariable robust_normal_mixture2(const prevariable&  a,
  const prevariable sighat,dvar_vector & r2)
{
  dvariable asig=a*sighat;
  dvariable avar=square(a*sighat);
  const MY_DOUBLE_TYPE c1=1.0/sqrt(2.0*3.141592653589793);
  int mmax=r2.indexmax();
  dvariable tmp=mmax*log(0.5+asig) -sum(log(1.e-10+
    0.95*c1*exp(-0.5*r2/avar)
   +0.05*c1/3.0*exp(-0.5*r2/(9.0*avar))));
  return tmp;
}

dvariable robust_normal_mixture(const prevariable&  a,
  const prevariable sighat,dvar_vector & r2)
{
  const MY_DOUBLE_TYPE c1=1.0/sqrt(2.0*3.141592653589793);
  dvariable sig=a*sighat;
  dvariable var=square(sig);
  int n=r2.indexmax();
  dvariable tmp=n*log(sig) -sum(log(1.e-10+
    0.95*c1*exp(-0.5*r2/var)
   +0.05*c1/3.0*exp(-0.5*r2/(9.0*var))));
  return tmp;
}


imatrix dvar_len_fish_stock_history::
  get_bounds_for_implicit_catch_effort_design_matrix(void)
{
  ivector ff29=column(fish_flags,29);
  ivector ff60=column(fish_flags,60);
  ivector ff27=column(fish_flags,27);

  ivector group=ff29; 
  imatrix bounds;

  if (sum(group))
  {
    group_flag_manager gfm(group);
    imatrix group_members=gfm.get_group_members();
    ivector num_in_group=gfm.get_num_in_group();
    int numgroups=gfm.get_numgroups();
    ivector ind(1,3);
    ind(1)=numgroups;
    ind(2)=numgroups;
    ind(3)=num_fisheries;
    bounds.allocate(1,3,1,ind);
    bounds.initialize();
    ivector& irow_vector=bounds(1);
    ivector& icol_vector=bounds(2);
    for (int ig=1;ig<=numgroups;ig++)
    {
      for (int in=1;in<=num_in_group(ig);in++)
      {
        int fi=group_members(ig,in);
        if (in==1) // only need to do this once per group
        {
          do_col_bounds_count(fi,icol_vector(ig));
        }
        do_row_bounds_count(fi,irow_vector(ig));
      }
      // can not have more parameters than observations so num cols can not be greater than num rows.
      bounds(2,ig)=min(bounds(2,ig),bounds(1,ig));
    }
    for (int i=1;i<=num_fisheries;i++)
    {
      bounds(3,i)=bounds(2,group(i));
    }
  }
  else
  {
    bounds.allocate(1,3,1,num_fisheries);
    bounds.initialize();
    ivector& irow_vector=bounds(1);
    ivector& icol_vector=bounds(2);
    for (int fi=1;fi<=num_fisheries;fi++)
    {
      do_col_bounds_count(fi,icol_vector(fi));
      do_row_bounds_count(fi,irow_vector(fi));
      // can not have more parameters than observations so num cols can not be greater than num rows.
      bounds(2,fi)=min(bounds(2,fi),bounds(1,fi));
    }
    bounds(3)=bounds(2);
  }
  // can not have more parameters than observations so num cols can not be greater than num rows.
  for (int i=1;i<=bounds(1).indexmax();i++)
  {
    //bounds(2,i)=min(bounds(2,i),bounds(1,i));
  }
  return bounds;
}

void dvar_len_fish_stock_history::
  new_catch_equations_calc_implicit_experiment_loop(dvar_vector& sv,
  ivector* pq_flag, dvariable & ffpen)
{
  if (af170q0 && loop_flag==0)
  {
    loop_flag=1;
    for (int i=1;i<=5;i++)
    {
      iloop=i;
      new_catch_equations_calc_implicit_experiment(sv,pq_flag,ffpen);
    }
  }
  new_catch_equations_calc_implicit_experiment(sv,pq_flag,ffpen);
}

void dvar_len_fish_stock_history::new_catch_equations_calc_implicit_experiment
  (dvar_vector& sv, ivector* pq_flag, dvariable & ffpen)
{
  //greport("beginning catch_equations_calc");
  ofstream ofs("q0eval_dbug.rpt",ios::app);
  if (generate_report && iloop==0)
  {
    ofs << " pq_flag iloop mn_csurv_chk  mn_csurv      mn_Xchk\
      mn_tot_mort  Obs_catch       Pred_catch" << endl;
  }
  tmprecr.initialize();
  tmpinitpop.initialize();
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  if (af170q0==0)
  {
    tot_mort.initialize();
  }
  else
  {
    tot_mort_q0.initialize();
  }
  // calculate the totalfishing mortality and survival rates for each
  // fishing period in each region
  //cout <<"newrshimp.cpp " <<  num_fish_periods << endl;
  int ir;
  if (num_regions>1 && !age_flags(114))
  {
    setup_diffusion();
  }
#if !defined(ZFEV)
     if (!parest_flags(143))
#endif
     {
       if (!parest_flags(175))
         mean_length_calc(0);
       else
         mean_length_calcx();
     }
#if !defined(ZFEV)
     else
     {
       // do this once at the beginning of the run
       if (parest_flags(144))
       {
         pcfsh->mean_length_calc();
         pcfsh->variance_calc();
       }
     }
#endif
  /*
  if (generate_report && iloop>0)
  {
    ofs << " " << endl;
    ofs << " before get inital population " << endl;
    q0debug_report(pq_flag,ofs);
  }
  */
  if (age_flags(177))
  {
    set_totpop();   //NMD_10feb2025
    if (pmsd)
    {
      pmsd->totpop=pmsd->totpop_coff;
    }
//    totpop=totpop_coff;
    if (pq_flag)
    {
      age_flags(192)=1;
    }
    ffpen+=get_initial_population(sv,0,pq_flag);
    age_flags(192)=0;
  }
  else
  {
    if (af170q0==0)
    {
//      get_population_multipliers(sv,pq_flag);
      if (!parest_flags(155)) get_population_multipliers(sv,pq_flag);
      get_initial_population(sv,0,0);
    }
    else if (!parest_flags(155))
    {
      if (age_flags(171)==0)
      {
        get_population_multipliers(sv,pq_flag);
      }
      // get the initial age structure
      if (!age_flags(94))
      {
        get_initial_age_structure(totpop,sv);
      }
      else
      {  
        xget_initial_age_structure_equilibrium();
      }
    }
    else
    {
      age_flags(192)=1;
      get_initial_population(sv,0,pq_flag);
      age_flags(192)=0;
    }
  }
  /*
  if (generate_report && iloop>0)
  {
    ofs << " " << endl;
    ofs << " after get_init_pop - running get numbers at age " << endl;
    q0debug_report(pq_flag,ofs);
  }
  */

  // do we have the variance for this?
  calculate_the_mean_weight();

  // set fm_level to crash if used and not set properly
  fm_level=1.e+60;
  log_pred_effort_by_fishery=1.e-10;   //NMD_20may2021 - so can take log

  new_get_numbers_at_age_implicit(sv,pq_flag,ffpen);

  get_annual_recruitment_flag=0;  //NMD_10jun2022

  ivector  mps = get_equilibrium_movements_per_season();
  int nya=get_numyears_average();

  //ivector * pq_flag=0;
  /*
  if (generate_report && iloop>0)
  {
    ofs << " " << endl;
    ofs << " before get initial survival " << endl;
    q0debug_report(pq_flag,ofs);
  }
  */

  if (nya>0)
  {
    dvar4_array mm = get_initial_equilibrium_survival
      (mps,nya,pq_flag);
  }
  else
  {
    cout << "need to check out nya=0 case " << endl;
  }
  //(ivector & mps,int nya,ivector * pq_flag)
  /*
  if (generate_report && iloop>0)
  {
    ofs << " " << endl;
    ofs << " after get initial survival " << endl;
    q0debug_report(pq_flag,ofs);
    ofs << " " << endl;
    ofs << " " << endl;
  }
  */

//  dvariable ftmp=implicit_catch_effort_relationship();
  dvariable ftmp=0.0;
  if (parest_flags(378)>0)  //NMD_13dec2021
  {
    ftmp=implicit_catch_effort_relationship();
  }
  cout << "VVVV!!! implicit_catch_effort_relationship"
         " regression penalty= " << ftmp << endl;
  ffpen+=ftmp;
  if (ppstf && !af170q0 && generate_report)  //NMD_5jun2024
  {
    ppstf->impl_fm_level_regr_pen=ftmp;
  }


//  get_implicit_catchability(*this);
  if (generate_report && !af170q0)   //NMD_18jun2024
  {
    if (do_fishery_projections_flag==1)
    {
      get_implicit_catchability_catch_conditioned(*this);
    }
    get_implicit_catchability(*this);
  }

  //do_fish_mort_intermediate_calcs();

  //NMD_10jun2022
  if (af170q0==1 && age_flags(171)==1)
  {
    dmatrix vtmp;      
//    if ((age_flags(94)==3  || age_flags(182)) && iloop==5)
    if ((age_flags(94)==3  || age_flags(182)) || iloop==5)
    {
      get_annual_recruitment_flag=1;
      vtmp=value(new_get_numbers_at_age_implicit(sv,pq_flag,ffpen));
      get_annual_recruitment_flag=0;
    }
  }
  //NMD_10jun2022


  if (age_flags(180))
  {   
    (*clogf)  <<  "num_fish" << endl;
    (*clogf)  <<  num_fish << endl;
  }
  if (age_flags(180))
  {   
    (*clogf)  <<  "catch" << endl;
    (*clogf)  <<  catch << endl;
  }
}

void dvar_len_fish_stock_history::q0debug_report(ivector* pq_flag,ofstream& ofs){
  int ipqf=0;
  if (pq_flag)
  {
    ipqf=1;
  }
  int jj=12;  //nage
  int ss=1;   //season
  int mv=1;   //movement
// Calculate tot_mort and catches
  int ir=1;
  int ip=1;
  int fi=2;
// tot_mort(1,num_regions,1,num_fish_periods,1,nage)
// - calculate the mean over all regions and ages for period: ip
  MY_DOUBLE_TYPE tmp=0;
  int icount=0;
  for (int rr=1;rr<=num_regions;rr++)
  {
    for (int j=1;j<=nage_by_region(rr);j++)
    {
      icount++;
      if (ipqf)
      {
        tmp+=value(tot_mort_q0(rr,ip,j));
      }
      else
      {
        tmp+=value(tot_mort(rr,ip,j));
      }	
    }
  }
  MY_DOUBLE_TYPE mn_tot_mort=tmp/icount;
// Catch(1,num_regions,1,num_fish_periods,1,nage)
// catch_q0(1,num_regions,1,num_fish_periods,1,nage)
  dvariable totcatch=0.0;
  dvector vtmp;
  dvariable sv27=get_sv_region(ir,27);
  dvariable sv28=get_sv_region(ir,28);

  if (pq_flag)
  {
    vtmp=exp(value(catch_q0(ir,ip,fi)));
  }
  else
  {
    vtmp=value(exp_catch(ir,ip,fi));
  }
  for (int j=1;j<=nage;j++)
  {  
    totcatch+=vtmp(j)*
      normal_length_to_weight(0.5,-3.5,
      mean_length(ir,ip,fi,j),sqrt(vars(ir,ip,fi,j)),
      value(sv27),value(sv28));
  }
  totcatch/=1000.;

  MY_DOUBLE_TYPE mn_csurv_chk=(sum(csurv_chk(iloop,ss,mv,jj)))/num_regions;
  MY_DOUBLE_TYPE mn_csurv=(sum(csurv(ss,mv,jj)))/num_regions;
  MY_DOUBLE_TYPE mn_Xchk=sum(Xchk(iloop))/num_regions;

  ofs << "  " << ipqf << "        " << iloop << "     " << setprecision(5)
    << mn_csurv_chk  << "        " << setprecision(5)
    << mn_csurv << "    " << fixed
    << mn_Xchk << "    " << setprecision(5)
    << mn_tot_mort << "      " << setprecision(5)
    << obs_tot_catch(ir,ip,fi) << "      " << setprecision(5)
    << totcatch << endl;
}


ivector gram_schmidt_remove_extra_columns(dmatrix& M)
{
  // This routine removes redundant columns in the design matrix and returns
  // a vector pointing to the remaining columns
  int m=M.indexmax();
  int n=M(m).indexmax();
  dmatrix TQ(1,n,1,m);
  TQ=trans(M);

  int ii=0;
  ivector columns(1,m);
  columns.initialize();
  for (int i=1;i<=n;i++)
  {
    MY_DOUBLE_TYPE a=norm(TQ(i));
    if (a<1.e-8)
    {
      // at least one column removec
    }
    else
    {
      columns(++ii)=i;
      TQ(i)/=a;
      for (int j=i+1;j<=n;j++)
      {
        MY_DOUBLE_TYPE a=TQ(i)*TQ(j);
        TQ(j)-=a*TQ(i);
      }
    }
  }
  return columns(1,ii);
}

void gram_schmidt_qr(dmatrix& M,dmatrix& Q,dmatrix& R)
{
  int m=M.indexmax();
  int n=M(m).indexmax();
  dmatrix TQ(1,n,1,m);
  TQ=trans(M);

  dmatrix ID(1,n,1,n);
  ID.initialize();
  for (int i=1;i<=n;i++)
  {
    ID(i,i)=1.0;
  }

  dmatrix TR(1,n,1,n);
  TR=ID;
  for (int i=1;i<=n;i++)
  {
    MY_DOUBLE_TYPE a=norm(TQ(i));
    TQ(i)/=a;
    TR(i)/=a;
    for (int j=i+1;j<=n;j++)
    {
      MY_DOUBLE_TYPE a=TQ(i)*TQ(j);
      TQ(j)-=a*TQ(i);
      TR(j)-=a*TR(i);
    }
  }

  Q=trans(TQ);
  R=inv(trans(TR));
}

dvar_vector log_student_density(const prevariable&  v,const dvar_vector& x)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  dvar_vector tmp(mmin,mmax);
  dvariable u=0.5*log(v)+betaln(0.5,v/2.0);
  dvariable  w=0.5*(v+1.0);
  for (int i=mmin;i<=mmax;i++)
  {
    tmp(i)=-w*log(1.0+square(x(i))/v)-u;
  }
  return tmp;
}

dvar_vector neg_log_student_density(const prevariable&  v,const prevariable& s,
  const dvar_vector& x)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  dvar_vector tmp(mmin,mmax);
  dvariable u=0.5*log(v)+betaln(0.5,v/2.0);
  dvariable  w=0.5*(v+1.0);
  //dvariable  sv=s*v;
  dvariable ls=log(s);
  dvar_vector e=x/s;
  for (int i=mmin;i<=mmax;i++)
  {
    tmp(i)=ls +w*log(1.0+square(e(i))/v)+u;
  }
  return tmp;
}

