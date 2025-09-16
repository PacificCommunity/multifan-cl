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
#if defined(_WIN32)
#  if defined(close)
#    undef close
#  endif
#endif

static void df_show_deriv(void)
{
  verify_identifier_string("DX");
  prevariable_position xpos=restore_prevariable_position();
  verify_identifier_string("YB");
  MY_DOUBLE_TYPE df=restore_prevariable_derivative(xpos);
  cout << "WWW The derivative = " << df << endl;
}

void  debug_show(const prevariable& x)
{
//    save_identifier_string("YB");
  const char * str1;
  str1="YB";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  x.save_prevariable_position();
//    save_identifier_string("DX");
  const char * str2;
  str2="DX";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  gradient_structure::GRAD_STACK1->set_gradient_stack(df_show_deriv);
}
 
dvariable  dvar_fish_stock_history::recinpop_orth(dvar3_array *pN,int override)
{
  if (!allocated(Ninit_orthogonal))
    {
    Ninit_orthogonal.allocate(1,num_regions,1,last_real_year);
  }
  Ninit_orthogonal.initialize();

  newfpen=0.0;
  int common_recruit_flag=0;
  int ns=1;
  ivector nc(numcomp.indexmin(),numcomp.indexmax());
  int nr=num_regions;
  ivector pf(parest_flags.indexmin(),parest_flags.indexmax());
  if (pmsd)
  {
    ivector sf2=column(pmsd->species_flags,2);
    int common_recruit_flag=sum(sf2);
    if (common_recruit_flag==0)
      ns=pmsd->num_species;
  }
  dvar_vector mn(1,ns);
  dvar_vector num(1,ns);
  mn.initialize();
  num.initialize();

  for (int is=1;is<=ns;is++)
  {
    if (!pmsd || is==1)
    {
      pf=parest_flags;
    }
    else
    {
      pf=pmsd->parest_flags(is);
    }
    if (pmsd)
    {
      nr=get_region_bounds(is)(2)-get_region_bounds(is)(1)+1;
    }
    int degree=pf(155);
    dvar_vector year_effect(1,num_real_years);
    dvar_matrix region_effect(1,nr,1,num_real_years);
    dvar_matrix season_effect(1,num_seasons,1,num_real_years);
    dvar3_array ses_reg_effect(1,num_seasons,1,nr,1,num_real_years);
    int ir;
    int dy,dr,ds,dsr;
    if (is==1)
    {
      nc=numcomp;
      dy=degree_yr;
      dr=degree_reg;
      ds=degree_ses;
      dsr=degree_ses_reg;
    }
    else
    {
      nc=pmsd->numcomp(is);
      is=pmsd->num_species;
      dy=pmsd->degree_yr(is);
      dr=pmsd->degree_reg(is);
      ds=pmsd->degree_ses(is);
      dsr=pmsd->degree_ses_reg(is);
    }
  
    dmatrix rrecr_polys_yr;
    dmatrix rrecr_polys_reg;
    dmatrix rrecr_polys_ses;
    dmatrix rrecr_polys_ses_reg;
    d3_array oorth_recr_basis;
    dvar_matrix oorth_recr_all; 
    imatrix sses_reg_recr_flags;
    dvar_matrix yyearly_recr_all;
  

    if (pmsd && is>1) 
    {
      rrecr_polys_yr=pmsd->recr_polys_yr(is);
      rrecr_polys_reg=pmsd->recr_polys_reg(is);
      rrecr_polys_ses=pmsd->recr_polys_ses(is);
      rrecr_polys_ses_reg=pmsd->recr_polys_ses_reg(is);
      oorth_recr_basis=pmsd->orth_recr_basis(is);
      oorth_recr_all=pmsd->orth_recr_all(is);
      sses_reg_recr_flags=pmsd->ses_reg_recr_flags(is);
      yyearly_recr_all=pmsd->yearly_recr_all(is);
    }
    else
    {
      rrecr_polys_yr=recr_polys_yr;
      rrecr_polys_reg=recr_polys_reg;
      rrecr_polys_ses=recr_polys_ses;
      rrecr_polys_ses_reg=recr_polys_ses_reg;
      oorth_recr_basis=orth_recr_basis;
      oorth_recr_all=orth_recr_all;
      sses_reg_recr_flags=ses_reg_recr_flags;
      yyearly_recr_all=yearly_recr_all;
    }
    int numbases=oorth_recr_basis.indexmax();
    dvar_matrix all_effects(1,numbases,1,num_real_years);
    int nry=(nyears-1)/num_seasons+1;
    int jj=0;
    int i,k;
    dvar3_array recr_by_yr(1,nry,1,num_seasons,1,nr);
    dvar3_array recr_by_seas(1,nry,1,num_seasons,1,nr);
    dvar3_array recr_by_reg(1,nry,1,num_seasons,1,nr);
    dvar3_array recr_by_seas_reg(1,nry,1,num_seasons,1,nr);
    dvar3_array recr_by_ses_reg(1,nry,1,num_seasons,1,nr);
    //dvar3_array recr_by_ses_reg(1,num_real_years,1,num_seasons,1,nr);
    recr_by_yr.initialize();
    recr_by_seas.initialize();
    recr_by_reg.initialize();
    recr_by_seas_reg.initialize();
    recr_by_ses_reg.initialize();
    {
      //ofstream ofs ("oorth_recr_all");
      //ofs << "oorth_recr_all" << endl;
      //ofs << setfixed() << setw (6) << setprecision(3) << oorth_recr_all 
      //  << endl << endl;
    }
    if (numcomp(1)!=1)
    {
      cerr << "This can not happen" << endl;
    }
    //if (degree_yr>0)
    //{
      for (i=1;i<=numcomp(1);i++)
      {
        jj++;
        if (degree_yr>=0)
        {
          all_effects(jj)=rrecr_polys_yr*oorth_recr_all(jj)(0,degree_yr-1);
          for (k=1;k<=num_real_years;k++)
          {
            recr_by_ses_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
            recr_by_yr(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
          }
        }
      }
    //}
    
    /*
    for (i=1;i<=numcomp(2);i++)
    {
      jj++;
      if (degree_ses>0)
      {
        all_effects(jj)=rrecr_polys_ses*oorth_recr_all(jj)(0,degree_ses-1);
        for (k=1;k<=num_real_years;k++)
        {
          recr_by_seas(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
          recr_by_ses_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
        }
      }
    }
    */
    for (i=1;i<=numcomp(2);i++)
    {
      jj++;
      if (degree_ses>0)
      {
        all_effects(jj)=rrecr_polys_ses*oorth_recr_all(jj)(0,degree_ses-1);
        for (k=1;k<=num_real_years;k++)
        {
          recr_by_seas(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
          recr_by_ses_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
        }
      }
    }
    for (i=1;i<=numcomp(3);i++)
    {
      jj++;
      if (degree_reg>0)
      {
        all_effects(jj)=rrecr_polys_reg*oorth_recr_all(jj)(0,degree_reg-1);
        for (k=1;k<=num_real_years;k++)
        {
          recr_by_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
          recr_by_ses_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
        }
      }
    }
    
    for (i=1;i<=numcomp(4);i++)
    {
      jj++;
      if (degree_ses_reg>0)
      {
        all_effects(jj)=rrecr_polys_ses_reg*oorth_recr_all(jj)(0,degree_ses_reg-1);
        for (k=1;k<=num_real_years;k++)
        {
          recr_by_seas_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
          recr_by_ses_reg(k)+=all_effects(jj,k)*oorth_recr_basis(jj);
        }
      }
    }

//    dvar_vector tmp(1,nyears);
    dvar_vector tmp(1,last_real_year);  //NMD_25jul2022
    // nr=get_region_bounds(is)(2)-get_region_bounds(is)(1)+1;
    if (common_recruit_flag)
    {
      if (is>1)
      {
        cerr << " Species can not be > 1 here" << endl;
        ad_exit(1);
      }
      int roffset=get_region_bounds(1)(2);
      for (ir=1;ir<=nr;ir++)
      {
//        for (int iy=1;iy<=nyears;iy++)
        for (int iy=1;iy<=last_real_year;iy++)  //NMD_25jul2022
        {
          int ty=(iy-1)/num_seasons+1;
          int ts=(iy-1)%num_seasons+1;
          if (sses_reg_recr_flags(ts,ir)==0)
          {
            Ninit_orthogonal(ir,iy)=-15.0;
            Ninit_orthogonal(ir+roffset,iy)=-15.0;
          }
          else
          {
            Ninit_orthogonal(ir,iy)=15.0+recr_by_ses_reg(ty,ts,ir);
            Ninit_orthogonal(ir+roffset,iy)=15.0+recr_by_ses_reg(ty,ts,ir);
          }
        }
      }
    }
    else
    {
      int roffset=0;
      if (is>1)
      {
        roffset=get_region_bounds(is-1)(2);
      }
//      dvar_matrix tmp(1,nr,1,nyears);  //NMD_25jul2022
      dvar_matrix tmp(1,nr,1,last_real_year);
      for (ir=1;ir<=nr;ir++)
      {
//        for (int iy=1;iy<=nyears;iy++)
        for (int iy=1;iy<=last_real_year;iy++)  //NMD_25jul2022
        {
          int ty=(iy-1)/num_seasons+1;
          int ts=(iy-1)%num_seasons+1;
          if (sses_reg_recr_flags(ts,ir)==0)
          {
            Ninit_orthogonal(ir+roffset,iy)=-15.0;
          }
          else
          {
            double minrec=10.0;
            if (parest_flags(389))
            {
              minrec=parest_flags(389);
            }
            Ninit_orthogonal(ir+roffset,iy)=minrec+recr_by_yr(ty,ts,ir);
            Ninit_orthogonal(ir+roffset,iy)=minrec+recr_by_yr(ty,ts,ir);
            Ninit_orthogonal(ir+roffset,iy)+=recr_by_seas(ty,ts,ir);
            Ninit_orthogonal(ir+roffset,iy)+=recr_by_reg(ty,ts,ir);
            Ninit_orthogonal(ir+roffset,iy)+=recr_by_seas_reg(ty,ts,ir);
            mn(is)+=Ninit_orthogonal(ir+roffset,iy);
            num(is)+=1;
          }
        }
      }
    }
  }
  /* will get rid of this penalty probably
  if (!pmsd)
  {
    dvariable mn2=mean(Ninit_orthogonal);
    cout << mn2-mn(1)/num(1) << endl;
    newfpen+=0.1*norm2(Ninit_orthogonal-mn2);
  }
  else
  {
    for (int is=1;is<=ns;is++)
    {
      int minr=get_region_bounds(is)(1);
      int maxr=get_region_bounds(is)(2);
      dvariable mn2=mean(Ninit_orthogonal.sub(minr,maxr));
      cout << mn2-mn(is) << endl;
      newfpen+=0.1*norm2(Ninit_orthogonal.sub(minr,maxr)-mn2);
    }
  }
  */
  copy_recruitment_to_population(Ninit_orthogonal,pN);


  dvariable fpen=0.0;
  return fpen;
}

int xinit_count=0;

void dvar_fish_stock_history::recinpop_orth_xinit(ofstream& xof,dvector& x,int& ii)
{
  int iisave=0;
  if (parest_flags(155)<0)
  {
    if (!pmsd || pmsd->num_species==1)
    {
      iisave=ii;
      set_value_inv(new_orth_recr,x,ii);
      xinit_message(xof,iisave,ii,"new_orth_recr_all");
    }
    else 
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (sum(sf2)==0)
      {
        iisave=ii;
        set_value_inv(new_orth_recr,x,ii);
        xinit_message(xof,iisave,ii,"new_orth_recr_all");
        for (int is=2;is<=pmsd->num_species;is++)
        {
          iisave=ii;
          set_value_inv(pmsd->new_orth_recr(is),x,ii);
          xinit_message(xof,iisave,ii,"new_orth_recr_all");
        }
      }
      else
      {
        if (sf2(1))
        {
          iisave=ii;
          set_value_inv(new_orth_recr,x,ii);
          xinit_message(xof,iisave,ii,"new_orth_recr_all");
        }
        for (int is=2;is<=pmsd->num_species;is++)
        {
          if (sf2(is))
          {
            iisave=ii;
            set_value_inv(pmsd->new_orth_recr(is),x,ii);
            xinit_message(xof,iisave,ii,"new_orth_recr_all");
          }
        }
      }
    }
  }
  else
  {
    int ns=1;
    int iisave=0;
    if (pmsd) ns=pmsd->num_species;
    for (int is=1;is<=ns;is++)
    {
      if (is==1)
      {
        int ir;
        int jj=0;
        for (int i=1;i<=numcomp(1);i++)
        {
          if (degree_yr>0)
          {
            iisave=ii;
            set_value_inv_partial(orth_recr_all(++jj),x,ii,degree_yr,100.);
            xinit_message(xof,iisave,ii,"orth_recr_all_year");
//            xinit_message(xof,iisave,ii,"orth_recr_all",orth_recr_all(jj)); //NMD
          }
          else
          {
            pmsd_error();
            iisave=ii;
            int by=1+ny_begin_yr-1;
            int ey=num_real_years-ny_end_yr+1;
            set_value_partial(yearly_recr_all(++jj),x,ii,ey-by+2,100.);
            xinit_message(xof,iisave,ii,"yearly_recr_all");
          }
        }
  
        for (int i=1;i<=numcomp(2);i++)
        {
          iisave=ii;
          set_value_inv_partial(orth_recr_all(++jj),x,ii,degree_ses,100.);
          xinit_message(xof,iisave,ii,"orth_recr_all_season");
//          xinit_message(xof,iisave,ii,"orth_recr_all",orth_recr_all(jj));  //NMD
        }
  
        for (int i=1;i<=numcomp(3);i++)
        {
          iisave=ii;
          set_value_inv_partial(orth_recr_all(++jj),x,ii,degree_reg,100.);
          xinit_message(xof,iisave,ii,"orth_recr_all_region");
//          xinit_message(xof,iisave,ii,"orth_recr_all",orth_recr_all(jj)); //NMD
        }
  
        for (int i=1;i<=numcomp(4);i++)
        {
          iisave=ii;
          set_value_inv_partial(orth_recr_all(++jj),x,ii,degree_ses_reg,100.);
          xinit_message(xof,iisave,ii,"orth_recr_all_season_region");
//          xinit_message(xof,iisave,ii,"orth_recr_all",orth_recr_all(jj)); //NMD
        }
      }
      else
      {
        ivector numcomp=pmsd->numcomp(is);
        ivector sf2=column(pmsd->species_flags,2);
        int degree_yr=pmsd->degree_yr(is);
        int degree_ses=pmsd->degree_ses(is);
        int degree_reg=pmsd->degree_reg(is);
        int degree_ses_reg=pmsd->degree_ses_reg(is);
        if (sum(sf2)==0)
        {
          int ir;
          int jj=0;
          for (int i=1;i<=numcomp(1);i++)
          {
            if (degree_yr>0)
            {
              iisave=ii;
              set_value_inv_partial(pmsd->orth_recr_all(is)(++jj),x,ii,
                degree_yr,10.);
              xinit_message(xof,iisave,ii,"orth_recr_all_year");
//              xinit_message(xof,iisave,ii,"orth_recr_all",
//                            pmsd->orth_recr_all(is)(jj)); //NMD
            }
            else
            {
              pmsd_error();
              int by=1+ny_begin_yr-1;
              int ey=num_real_years-ny_end_yr+1;
              set_value_partial(yearly_recr_all(++jj),x,ii,ey-by+2,10.);
            }
          }
    
          for (int i=1;i<=numcomp(2);i++)
          {
            iisave=ii;
            set_value_inv_partial(pmsd->orth_recr_all(is)(++jj),x,ii,
              degree_ses,10.);
            xinit_message(xof,iisave,ii,"orth_recr_all_season");
//            xinit_message(xof,iisave,ii,"orth_recr_all",
//                          pmsd->orth_recr_all(is)(jj-1)); //NMD
          }
    
          for (int i=1;i<=numcomp(3);i++)
          {
            iisave=ii;
            set_value_inv_partial(pmsd->orth_recr_all(is)(++jj),x,ii,
              degree_reg,10.);
            xinit_message(xof,iisave,ii,"orth_recr_all_region");
//            xinit_message(xof,iisave,ii,"orth_recr_all",
//                          pmsd->orth_recr_all(is)(jj)); //NMD
          }
    
          for (int i=1;i<=numcomp(4);i++)
          {
            iisave=ii;
            set_value_inv_partial(pmsd->orth_recr_all(is)(++jj),
              x,ii,degree_ses_reg,10.);
            xinit_message(xof,iisave,ii,"orth_recr_all_season_region");
//            xinit_message(xof,iisave,ii,"orth_recr_all",
//                          pmsd->orth_recr_all(is)(jj)); //NMD
          }
        }
        else
        {
          // do nothng here recruitments are assigned in another routine
          //pmsd_error();
        }
      }
    }
  }
  {
    /*
    ofstream ofs("orth-xinit"+str(parest_flags(392))+"_"+str(xinit_count++));
    ofs << " x.indexmin()    x.indexmax()    x " << endl;
    ofs << x.indexmin() << "  "  << x.indexmax() << "  "  << x << endl;
    ofs << "orth_recr_all" << endl;
    ofs << orth_recr_all << endl;
    ofs.close();
    //ad_exit(1);
    */
  }
}

MY_DOUBLE_TYPE dvar_fish_stock_history::orth_poly_fit(dvar_vector& vx)
{
  dvariable f;
  int ii=1;
  recinpop_orth_reset(vx,ii);
  recinpop_orth(0,3);
  f=norm2(Ninit_orthogonal-Ninit_standard);
  return value(f);
}

//dvector dvar_fish_stock_history::get_new_orthpolys(void)
void dvar_fish_stock_history::get_corresponding_orthogonal_coefficients(void)
{
  int nvar=recinpop_orth_size_count();
  independent_variables x(1,nvar);
  x.initialize();
  fmm fmc(nvar);
  MY_DOUBLE_TYPE f;
  dvector g(1,nvar);
  g.initialize();
  dvector xbest(1,nvar);
  MY_DOUBLE_TYPE fbest=1.e+50;;
  dvector gbest(1,nvar);
  gbest.fill_seqadd(1.e+50,0.);
  MY_DOUBLE_TYPE maxg;
  fmc.iprint=0;
  fmc.crit=1.e-6;
  gradient_structure::set_YES_DERIVATIVES();
  fmc.scroll_flag=0;
  fmc.imax=40;
  fmc.min_improve=0.;
  fmc.maxfn=100;
  fmc.ifn=0;
  {
    while (fmc.ireturn>=0)
    {
      int badflag=0;
      fmc.fmin(f,x,g);
      {
        maxg=max(g);
      }
      if (fmc.ireturn>0)
      {
        dvar_vector vx=dvar_vector(x);
        f=orth_poly_fit(vx);
        gradcalc(nvar,g);
      }
    }
    if (f < fbest)
    {
      fbest=f;
      gbest=g;
      xbest=x;
    }
  }
  gradient_structure::set_NO_DERIVATIVES();
  dvar_vector vxbest=dvar_vector(xbest);
  {
    //Npred.allocate(1,num_regions,1,last_real_year);
    d3_array tmp1(1,num_regions,1,last_real_year,1,2);
    d3_array tmp2(1,num_regions,1,last_real_year,1,2);    
  }
    
  int ii=1;
  recinpop_orth_reset(vxbest,ii);
}

int orth_reset_count=0;

void dvar_fish_stock_history::recinpop_orth_reset(dvar_vector& x,int& ii)
{
  if (parest_flags(155)<0)
  {
    if (!pmsd || pmsd->num_species==1)
    {
      set_value(new_orth_recr,x,ii);
    }
    else 
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (sum(sf2)==0)
      {
        set_value(new_orth_recr,x,ii);
        for (int is=2;is<=pmsd->num_species;is++)
        {
          set_value(pmsd->new_orth_recr(is),x,ii);
        }
      }
      else
      {
        pmsd_error();
        if (sf2(1))
        {
          set_value(new_orth_recr,x,ii);
        }
        for (int is=2;is<=pmsd->num_species;is++)
        {
          if (sf2(is))
            set_value(pmsd->new_orth_recr(is),x,ii);
        }
      }
    }
  }
  else
  {
    int ns=1;
    if (pmsd) ns=pmsd->num_species;
    for (int is=1;is<=ns;is++)
    {
      if (is==1)
      {
        int ir;
        int jj=0;
        for (int i=1;i<=numcomp(1);i++)
        {
          if (degree_yr>0)
          {
            set_value_partial(orth_recr_all(++jj),x,ii,degree_yr,100.);
          }
          else
          {
            int by=1+ny_begin_yr-1;
            int ey=num_real_years-ny_end_yr+1;
            set_value_partial(yearly_recr_all(++jj),x,ii,ey-by+2,100.);
          }
        }
  
        for (int i=1;i<=numcomp(2);i++)
        {
          set_value_partial(orth_recr_all(++jj),x,ii,degree_ses,100.);
        }
  
        for (int i=1;i<=numcomp(3);i++)
        {
          set_value_partial(orth_recr_all(++jj),x,ii,degree_reg,100.);
        }
  
        for (int i=1;i<=numcomp(4);i++)
        {
          set_value_partial(orth_recr_all(++jj),x,ii,degree_ses_reg,100.);
        }
      }
      else
      {
        ivector numcomp=pmsd->numcomp(is);
        ivector sf2=column(pmsd->species_flags,2);
        int degree_yr=pmsd->degree_yr(is);
        int degree_ses=pmsd->degree_ses(is);
        int degree_reg=pmsd->degree_reg(is);
        int degree_ses_reg=pmsd->degree_ses_reg(is);
        if (sum(sf2)==0)
        {
          int ir;
          int jj=0;
          for (int i=1;i<=numcomp(1);i++)
          {
            if (degree_yr>0)
            {
              set_value_partial(pmsd->orth_recr_all(is)(++jj),x,ii,
                degree_yr,10.);
            }
            else
            {
              pmsd_error();
              int by=1+ny_begin_yr-1;
              int ey=num_real_years-ny_end_yr+1;
              set_value_partial(yearly_recr_all(++jj),x,ii,ey-by+2,10.);
            }
          }
    
          for (int i=1;i<=numcomp(2);i++)
          {
            set_value_partial(pmsd->orth_recr_all(is)(++jj),x,ii,degree_ses,
              10);
          }
    
          for (int i=1;i<=numcomp(3);i++)
          {
            set_value_partial(pmsd->orth_recr_all(is)(++jj),x,ii,degree_reg,
              10);
          }
    
          for (int i=1;i<=numcomp(4);i++)
          {
            set_value_partial(pmsd->orth_recr_all(is)(++jj),
              x,ii,degree_ses_reg,10.);
          }
        }
        else
        {
          // do nothng here recruitments are assigned in another routine
          //pmsd_error();
        }
      }
    }
  }
}

dvariable dvar_fish_stock_history::recinpop_orth_penalties(void)
{
  dvariable fpen=0.0;
  dvariable like_level=0.0;  //NMD_05dec2022
  double penwght=0.01;
  if (parest_flags(155)>0)
  {
    int ns=1;
    if (pmsd) ns=pmsd->num_species;
    for (int is=1;is<=ns;is++)
    {
      if (is==1)
      {
        int ir;
        int jj=0;
        for (int i=1;i<=numcomp(1);i++)
        {
          if (degree_yr>0)
          {
            like_level+=penwght*norm2(orth_recr_all((jj+1)).sub(1,degree_yr)); //NMD_05dec2022
            fpen+=penwght*norm2(orth_recr_all(++jj).sub(1,degree_yr)); //NMD_04Oct2022
//            fpen+=penwght*norm2(orth_recr_all(++jj).sub(0,degree_yr));
          }
        }
        ppstf->orth_poly_penalty_by_level(1,1)=value(like_level);
        like_level=0.0;
  
        for (int i=1;i<=numcomp(2);i++)
        {
          //set_value_partial(orth_recr_all(++jj),x,ii,degree_ses,100.);
          like_level+=penwght*norm2(orth_recr_all((jj+1)).sub(0,degree_ses));
          fpen+=penwght*norm2(orth_recr_all(++jj).sub(0,degree_ses));
        }
        ppstf->orth_poly_penalty_by_level(1,2)=value(like_level);
        like_level=0.0;	  
  
        for (int i=1;i<=numcomp(3);i++)
        {
          //set_value_partial(orth_recr_all(++jj),x,ii,degree_reg,100.);
          like_level+=penwght*norm2(orth_recr_all((jj+1)).sub(0,degree_reg));
          fpen+=penwght*norm2(orth_recr_all(++jj).sub(0,degree_reg));
        }
        ppstf->orth_poly_penalty_by_level(1,3)=value(like_level);
        like_level=0.0;	  
  
        for (int i=1;i<=numcomp(4);i++)
        {
          //set_value_partial(orth_recr_all(++jj),x,ii,degree_ses_reg,100.);
          like_level=penwght*norm2(orth_recr_all((jj+1)).sub(0,degree_ses_reg));
          fpen+=penwght*norm2(orth_recr_all(++jj).sub(0,degree_ses_reg));
        }
        ppstf->orth_poly_penalty_by_level(1,4)=value(like_level);
        like_level=0.0;
      }
      else
      {
        ivector numcomp=pmsd->numcomp(is);
        ivector sf2=column(pmsd->species_flags,2);
        int degree_yr=pmsd->degree_yr(is);
        int degree_ses=pmsd->degree_ses(is);
        int degree_reg=pmsd->degree_reg(is);
        int degree_ses_reg=pmsd->degree_ses_reg(is);
        if (sum(sf2)==0)
        {
          int ir;
          int jj=0;
          for (int i=1;i<=numcomp(1);i++)
          {
            if (degree_yr>0)
            {
              //set_value_partial(pmsd->orth_recr_all(is)(++jj),x,ii,
              //  degree_yr,10.);
              like_level+=penwght*norm2(pmsd->orth_recr_all(is)((jj+1)).sub(1,degree_yr));  //NMD_04Oct2022
              fpen+=penwght*norm2(pmsd->orth_recr_all(is)(++jj).sub(1,degree_yr));  //NMD_04Oct2022
//              fpen+=penwght*norm2(pmsd->orth_recr_all(is)(++jj).sub(0,degree_yr));
            }
          }
          ppstf->orth_poly_penalty_by_level(is,1)=value(like_level);
          like_level=0.0;
    
          for (int i=1;i<=numcomp(2);i++)
          {
            //set_value_partial(pmsd->orth_recr_all(is)(++jj),x,ii,degree_ses,
            //  10);
            like_level+=penwght*norm2(pmsd->orth_recr_all(is)((jj+1)).sub(0,degree_ses));
            fpen+=penwght*norm2(pmsd->orth_recr_all(is)(++jj).sub(0,degree_ses));
          }
          ppstf->orth_poly_penalty_by_level(is,2)=value(like_level);
          like_level=0.0;
    
          for (int i=1;i<=numcomp(3);i++)
          {
            //set_value_partial(pmsd->orth_recr_all(is)(++jj),x,ii,degree_reg,
            //  10);
            like_level+=penwght*norm2(pmsd->orth_recr_all(is)((jj+1)).sub(0,degree_reg));
            fpen+=penwght*norm2(pmsd->orth_recr_all(is)(++jj).sub(0,degree_reg));
          }
          ppstf->orth_poly_penalty_by_level(is,3)=value(like_level);
          like_level=0.0;
    
          for (int i=1;i<=numcomp(4);i++)
          {
            //set_value_partial(pmsd->orth_recr_all(is)(++jj),
            //  x,ii,degree_ses_reg,10.);
            like_level+=penwght*norm2(pmsd->orth_recr_all(is)((jj+1)).sub(0,degree_ses_reg));
            fpen+=penwght*norm2(pmsd->orth_recr_all(is)(++jj).sub(0,degree_ses_reg));
          }
          ppstf->orth_poly_penalty_by_level(is,4)=value(like_level);
          like_level=0.0;
        }
        else
        {
          // do nothing here recruitments are assigned in another routine
          //pmsd_error();
          cout << " ERROR: multi-sex orth-poly recruitment error" << endl;
	  ad_exit(1);
        }
      }
    }
  }
  return fpen;
}

int dvar_fish_stock_history::recinpop_orth_size_count(void)
{
  int ii=0;
  if (parest_flags(155)<0)
  {
    if (!pmsd || pmsd->num_species==1)
    {
      ii+=size_count(new_orth_recr);
    }
    else 
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (sum(sf2)==0)
      {
        ii+=size_count(new_orth_recr);
        for (int is=2;is<=pmsd->num_species;is++)
        {
          ii+=size_count(pmsd->new_orth_recr(is));
        }
      }
      else
      {
        if (sf2(1))
        {
          ii+=size_count(new_orth_recr);
        }
        for (int is=2;is<=pmsd->num_species;is++)
        {
          if (sf2(is))
            ii+=size_count(pmsd->new_orth_recr(is));
        }
      }
    }
  }
  else
  {
    int dy,dr,ds,dsr;
    ivector nc;
    int ns=1;
    if (pmsd) ns=pmsd->num_species;
    if (ns>1)
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (sum(sf2)!=0)ns=1;
    }
    for (int is=1;is<=ns;is++)
    {
      if (is==1)
      {
        nc=numcomp;
        dy=degree_yr;
        dr=degree_reg;
        ds=degree_ses;
        dsr=degree_ses_reg;
      }
      else
      {
        nc=pmsd->numcomp(is);
        is=pmsd->num_species;
        dy=pmsd->degree_yr(is);
        dr=pmsd->degree_reg(is);
        ds=pmsd->degree_ses(is);
        dsr=pmsd->degree_ses_reg(is);
      }

      for (int i=1;i<=nc(1);i++)
      {
        if (dy>0)
        {
          ii+=dy;
        }
        else
        {
          pmsd_error();
        }
      }
      for (int i=1;i<=nc(2);i++)
      {
        ii+=ds;
      }

      for (int i=1;i<=nc(3);i++)
      {
        ii+=dr;
      }

      for (int i=1;i<=nc(4);i++)
      {
        ii+=dsr;
      }
    }
  }
  return ii;
}
