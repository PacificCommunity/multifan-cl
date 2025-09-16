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

  extern dvar_vector * psv;
#if !defined(ZFEV)
  extern dvar_len_fish_stock_history * pcfsh;
#endif


void  print_FF(dvar_fish_stock_history& fsh);
void  calc_FF_mean(dvar_fish_stock_history& fsh);

dvariable robmean(_CONST dvar_vector& v);

int check_pos(const dvar_vector& qq,int n)
{
  int pflag=1;
  for (int i=1;i<=n;i++)
  {
    if (value(qq(i))<=0.0)
    {
      pflag=0;
      break;
    }
  }
  return pflag;
}

 void dvar_len_fish_stock_history::catch_equations_calc_implicit(dvar_vector& sv,
   dvariable& ffpen)
 {
   dvariable fpen1=0.0;
   missing_catch_counter.initialize();
   ss2_flag=0;
   //greport("beginning catch_equations_calc");
   //totpop=implicit_totpop_coff;
   totpop=totpop_coff;
   tmprecr.initialize();
   tmpinitpop.initialize();
   tot_mort.initialize();
   // calculate the totalfishing mortality and survival rates for each
   // fishing period in each region
 
 #if !defined(ZFEV)
      if (!parest_flags(143))
 #endif
      {
        if (!parest_flags(175))
          mean_length_calc();
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
   //if (age_flags(69))
   if (num_regions>1 && !age_flags(114))
   {
     setup_diffusion();
   }
   // do we have the variance for this?
   calculate_the_mean_weight();
   ffpen+=get_initial_population(sv,0,0);
 
   int current_year=1;
   ivector rip(1,num_regions);
   rip=1;
   tot_mort.initialize();
   do
   {
     for (int ir=1;ir<=num_regions;ir++)
     {
       int& ip=rip(ir);
       if (year(ir,ip)==current_year)
       {
         num_fish(ir,ip)=N(ir,current_year);
         if (num_fish_incidents(ir,ip)>0)
         {
           do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
             missing_catch_counter,fm_level_devs,fpen1);
         }
         if (!pmsd)
         {
           tot_mort(ir,ip)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
         }
         else
         {
           tot_mort(ir,ip)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))
             +fraction(ir,ip));
         }
         survival(ir,ip)=mfexp(-tot_mort(ir,ip));
         do
         {
           if (ip>num_fish_periods(ir)) break;
           if (ip==num_fish_periods(ir)) 
           {
             if (num_fish_incidents(ir,ip)>0)
             {
               do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
                 missing_catch_counter,fm_level_devs,fpen1);
             }
             if (!pmsd)
             {
               tot_mort(ir,ip)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
             }
             else
             {
               tot_mort(ir,ip)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))
                 +fraction(ir,ip));
             }
             survival(ir,ip)=mfexp(-tot_mort(ir,ip));
             break;
           }
           if (year(ir,ip+1)==current_year)
           {
             num_fish(ir,ip+1)=num_fish(ir,ip)-tot_mort(ir,ip);
             ip++;
             if (num_fish_incidents(ir,ip)>0)
             {
               do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
                 missing_catch_counter,fm_level_devs,fpen1);
             }
             if (!pmsd)
             {
               tot_mort(ir,ip)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
             }
             else
             {
               tot_mort(ir,ip)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))
                 +fraction(ir,ip));
             }
             survival(ir,ip)=mfexp(-tot_mort(ir,ip));
           }
           else
           {
             num_fish(ir,ip+1,1)=N(ir,current_year+1,1);
 
             --num_fish(ir,ip+1)(2,nage)=
               num_fish(ir,ip)(1,nage-1)-tot_mort(ir,ip)(1,nage-1);
 
             num_fish(ir,ip+1,nage)=
               log(1.e-10 + mfexp(num_fish(ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
                 + mfexp(num_fish(ir,ip,nage)-tot_mort(ir,ip,nage)) );
   
             if (current_year<nyears)
             {
               N(ir,current_year+1)(2,nage)=num_fish(ir,ip+1)(2,nage);
             }
 
             ip++;
             
             if (ip<=num_fish_periods(ir)) 
             {
               if (num_fish_incidents(ir,ip)>0)
               {
                 do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
                   missing_catch_counter,fm_level_devs,fpen1);
               }
               if (!pmsd)
               {
                 tot_mort(ir,ip)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
               }
               else
               {
                 tot_mort(ir,ip)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))
                   +fraction(ir,ip));
               }
               survival(ir,ip)=mfexp(-tot_mort(ir,ip));
             }
             break;
           }
         }
         while (1);
       }
       else   // there were no fisheries in this region for this year
       {
         if (current_year<nyears)
         {
           if (!pmsd)
           {
             for (int j=1;j<nage;j++)      // Loop over age classes
             {
               N(ir,current_year+1,j+1)=N(ir,current_year,j)-
                 exp(nat_mort(current_year,j));
             }
           }
           else
           {
             for (int j=1;j<nage;j++)      // Loop over age classes
             {
               N(ir,current_year+1,j+1)=N(ir,current_year,j)-
                 exp(get_nat_mort_region(ir)(current_year,j));
             }
           }
         }
       }
     }
     if(current_year<nyears)
     {
       // Changed af(57) to af(53) J.H. 27 Nov 01
       if (age_flags(53))
       {
         if ( !((current_year+1)%age_flags(53)) )
         {
           if (num_regions>1) do_the_diffusion(current_year+1,sv,N);
         }
       }
       else
       {
         if (num_regions>1) do_the_diffusion(current_year+1,sv,N);
       }
     }
     current_year++;
     if (current_year == nyears-1)
     {
       //cout << "trap" << endl;
     }
   }
   while (current_year<=nyears);
   //greport("B catch_equations_calc");
 
   /*
   for (int ir=1;ir<=num_regions;ir++)
   {
     for (int ip=1;ip<=num_fish_periods(ir);ip++)  
     {
       for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
       {        
         catch(ir,ip,fi)=fish_mort(ir,ip,fi)-log(1.e-20+tot_mort(ir,ip))+
           log(1.e-20+(1.0-survival(ir,ip)))+num_fish(ir,ip);
         MY_DOUBLE_TYPE t=sum(exp(value(catch(ir,ip,fi))));
         MY_DOUBLE_TYPE otc;
         if (data_fish_flags(1,parent(ir,ip,fi))) 
         {
           otc=obs_tot_catch(ir,ip,fi);
         }
         else // convert from weigth to numbers
         {
           otc=obs_tot_catch(ir,ip,fi)/value(get_mean_weight(ir,ip,fi));
         }
         MY_DOUBLE_TYPE err=(t-otc)/(1.1e-5+t);
         if (fabs(err) > 1.e-5 && missing_catch_for_incident_flag(ir,ip,fi)==0)
         {
           cout << t << "  "   << obs_tot_catch(ir,ip,fi) 
              << "  " << err << endl;
         }
       }
     }
   }
   
   */
   int ir;
   if (ss2_flag==0)
   {
     for (ir=1;ir<=num_regions;ir++)
     {
       for (int ip=1;ip<=num_fish_periods(ir);ip++)  
       {
         for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
         {        
           dvar_vector fm=FF(ir,ip,fi)+incident_sel(ir,ip,fi);
           fish_mort(ir,ip,fi)=fm;
         }
       }
     }
   }
 
   do_fish_mort_intermediate_calcs();
 
   if (do_fishery_projections_flag==1 )
   {
     do_fish_mort_intermediate_projection_calcs();
   }
   dvar3_array * pnum_fish=0;
   dvar3_array * pN=0;
   dvar3_array * ptm=0;
   //if (!pq_flag)
   if (af170q0==0)
   {
     pnum_fish=&num_fish;
     ptm=&tot_mort;
     pN=&N;
   }
   else
   {
     ptm=&tot_mort_q0;
     pnum_fish=&num_fish_q0;
     pN=&N_q0;
   }
 
   dvar4_array * pc =0;
   dvar4_array * pfmc =0;
   if (af170q0==0)
   {
     pc = &catch;
     pfmc = &fish_mort_calcs;
   }
   else
   {
     pc = &catch_q0;
     pfmc = &fish_mort_calcs_q0;
   }
   
   for (ir=1;ir<=num_regions;ir++)
   {
     int tmp_nfp;
     if (do_fishery_projections_flag==0)
     {
       tmp_nfp=num_real_fish_periods(ir);  
     }
     else
     {
       tmp_nfp=num_fish_periods(ir);  
     }
    
     int ip;
     for (ip=1;ip<=tmp_nfp;ip++)  
     {
       for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
       {        
         (*pc)(ir,ip,fi)=(*pfmc)(ir,ip,fi)+(*pnum_fish)(ir,ip);
         int mmin=(*pc)(ir,ip,fi).indexmin();
         int mmax=(*pc)(ir,ip,fi).indexmax();
         for (int j=mmin;j<=mmax;j++) 
         {
           if (value((*pc)(ir,ip,fi,j))>1.e+15)
           {
              cout << (*pc)(ir,ip,fi) << " " << (*pfmc)(ir,ip,fi,j) 
                   << " " << (*pnum_fish)(ir,ip,j) << endl;
           }
         }
       }
     }
   }
   int ip;
   for (ip=1;ip<=3;ip++)  
   {
     cout << endl << "ip =" << ip << endl;
     for (ir=1;ir<=num_regions;ir++)
     {
       cout <<"rsh3imp.cpp " <<  (*pnum_fish)(ir,ip) << endl;
       for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
       {
         cout << (*pc)(ir,ip,fi) << " " << (*pfmc)(ir,ip,fi) << endl;
       }
     }
   }
   ad_exit(1);
   //greport("leaving catch_equations_calc");
   //print_FF(*this);
   get_implicit_catchability(*this);
   //calc_FF_mean_partial(*this);
 }

void  print_FF(dvar_fish_stock_history& fsh)
{
  ofstream ofs("FFreport");
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    dmatrix ttmp(1,fsh.num_fish_times(i),1,2);
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      ttmp(nt,1)=nt;
      fsh.implicit_catchability(i,nt)=
        mfexp(fsh.FF(rr,rp,ri)-fsh.effort(rr,rp,ri));
      ttmp(nt,2)= value(fsh.implicit_catchability(i,nt));
    }
    //ofs << "fishery " << i << endl;
    //ofs << ttmp << endl;
  }
}

void  get_implicit_catchability_catch_conditioned(dvar_fish_stock_history& fsh)
{
  ivector ff29=column(fsh.fish_flags,29);
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    dvector testc(1,fsh.num_real_fish_times(i));
    int nft;
    nft=fsh.num_fish_times(i);  // temp setting for development
    int nrft;
    nrft=fsh.num_real_fish_times(i);  // temp setting for development
    for (int nt=nrft+1;nt<=nft;nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      if (fsh.eff_proj_fshry(i)==1 || fsh.catch_proj_fshry(i)==1)
      {
        dvariable tmp2;
	tmp2.initialize();
        if (fsh.parest_flags(377)>0 && fsh.parest_flags(378)>0)
        {
          dvector dv=fsh.fml_designvec(rr,rp,ri);
          if (!allocated(dv))
          {
            cerr << "Projection incident has no matching incident in "
                 << "terminal year; fishery: " << i << endl;
            cerr << "Terminating execution" << endl;
            ad_exit(1);
          }
          ivector num_ifmlrp=fsh.implicit_fml_bounds(3);
          dvar_vector ests=fsh.implicit_fm_level_regression_pars(i)
            (1,num_ifmlrp(i));
          MY_DOUBLE_TYPE leff=fsh.log_effort_by_fishery(i,fsh.fish_times(rr,rp,ri));
          int gi=i;
          if (ff29(i))
          {
            gi=ff29(i);
          }
          tmp2 = dv * inv(fsh.fml_R(gi)) * ests;
        }
        else
        {
          int rmnth=fsh.really_true_month(rr,rp);
          for (int fnyri=1;fnyri<=fsh.fshtms_finyr(i);fnyri++)
          {
            int pmnth=fsh.ptr_fml_Mbs_mnth_finyr(i,fnyri);
            if (rmnth==pmnth)
            {
              tmp2=fsh.q_level_finyr(i,fnyri);
              break;
            }
          }
        }
        fsh.implicit_catchability_ccond(i,nt)=tmp2;
      }
    }
    /*
    if (sum(ff29)==0)   // no grouping for catchability
    {
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        fsh.q0(i)=
          mean(fsh.implicit_catchability(i)(1,fsh.num_real_fish_times(i)));
      }
    }
    else
    {
      int ngroups=fsh.gfish_index.indexmax();
      dvar_vector grouped_q0(1,ngroups);
      dvector numbers(1,ngroups);
      grouped_q0.initialize();
      numbers.initialize();
      int i;
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        grouped_q0(ff29(i))
          +=sum(fsh.implicit_catchability(i)(1,fsh.num_real_fish_times(i)));
        numbers(ff29(i))+=fsh.num_real_fish_times(i);
      }
      for (i=1;i<=ngroups;i++)
      {
        grouped_q0(i)/=numbers(i);
      }
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        fsh.q0(i)=grouped_q0(ff29(i));
      }
    }
    */
  }
}

  
void  get_implicit_catchability(dvar_fish_stock_history& fsh)
{
  // !! confusion ove whether this is ff60 or ff29
  //ivector ff29=column(fsh.fish_flags,29);
  ivector ff29=column(fsh.fish_flags,29);
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    dvector testc(1,fsh.num_real_fish_times(i));
//    if (fsh.missing_effort_flag(i)==0)
//    {
    for (int nt=1;nt<=fsh.num_real_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      if (fsh.missing_effort_by_region_flag(rr,rp,ri)==0)
      {
        fsh.implicit_catchability(i,nt)=
          log(fsh.fm_level(rr,rp,ri))
            -fsh.log_effort_by_fishery(i,nt);
      }
      else
      {
        fsh.implicit_catchability(i,nt)=
          log(fsh.fm_level(rr,rp,ri))
            -fsh.log_true_effort_by_fishery(i,nt);
      }
    }
    /*  
    }
    else
    {
      int nft;
      nft=fsh.num_real_fish_times(i)+1;
      for (int nt=1;nt<=fsh.num_real_fish_times(i);nt++)
//      for (int nt=1;nt<=nft;nt++)
      {
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        fsh.implicit_catchability(i,nt)=
          log(fsh.fm_level(rr,rp,ri))
            -fsh.log_true_effort_by_fishery(i,nt);
      }
    }
    */
    if (fsh.do_fishery_projections_flag==1)
    {
      int nft;
      nft=fsh.num_fish_times(i);
      int nrft;
      nrft=fsh.num_real_fish_times(i);
      for (int nt=nrft+1;nt<=nft;nt++)
      {
        if (fsh.eff_proj_fshry(i)==1 || fsh.catch_proj_fshry(i)==1)
        {
          fsh.implicit_catchability(i,nt)=fsh.implicit_catchability_ccond(i,nt);
        }
      }
    }
    /*    
    if (sum(ff29)==0)   // no grouping for catchability
    {
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        fsh.q0(i)=
          mean(fsh.implicit_catchability(i)(1,fsh.num_real_fish_times(i)));
      }
    }
    else
    {
      int ngroups=fsh.gfish_index.indexmax();
      dvar_vector grouped_q0(1,ngroups);
      dvector numbers(1,ngroups);
      grouped_q0.initialize();
      numbers.initialize();
      int i;
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        grouped_q0(ff29(i))
          +=sum(fsh.implicit_catchability(i)(1,fsh.num_real_fish_times(i)));
        numbers(ff29(i))+=fsh.num_real_fish_times(i);
      }
      for (i=1;i<=ngroups;i++)
      {
        grouped_q0(i)/=numbers(i);
      }
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        fsh.q0(i)=grouped_q0(ff29(i));
      }
    }
    */
  }
  cout << "HERE" << endl;
}

void  calc_FF_mean(dvar_fish_stock_history& fsh)
{
  //ofstream ofs("RMreport");
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    int nft=fsh.num_fish_times(i);
    dvar_vector lcat(1,fsh.num_fish_times(i));
    int mmin;
    int mmax;
    dvar_vector limpc=log(fsh.implicit_catchability(i));
    int nt;
    for (nt=2;nt<=nft+1;nt+=3)
    {
      if (nt<5)
      {
	mmin=1;
	mmax=7;
      }
      else if (nt>nft-6)
      {
	mmin=nft-6;
	mmax=nft;
      }
      else
      {
	mmin=nt-3;
	mmax=nt+3;
      }
      if (nt<=nft)
      {
        lcat(nt)=robmean(limpc(mmin,mmax));
        lcat(nt-1)=lcat(nt);
      }
      else	
      {
        lcat(nt-1)=robmean(limpc(mmin,mmax));
      }
      if (nt<nft)
      lcat(nt+1)=lcat(nt);
    }
    fsh.effort_dev_coffs(i)=limpc-lcat;
    //ofs << "fishery " << i << endl;
    for (nt=2;nt<=nft;nt++)
    {
      fsh.catch_dev_coffs(i,nt)=lcat(nt)-lcat(nt-1);
    }
    //ofs << ttmp << endl;
  }
}

void  calc_FF_mean_partial(dvar_fish_stock_history& fsh, int print_switch)
{
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    int nft=fsh.num_fish_times(i);
    dvar_vector lcat(1,fsh.num_fish_times(i));
    int mmin;
    int mmax;
    dvar_vector limpc=log(fsh.implicit_catchability(i));
    ivector& impb=fsh.imp_back(i);
    ivector& impf=fsh.imp_forward(i);
    int ntb=1;
    int nt;
    for (nt=1;nt<=nft;nt++)
    {
      if (nt==1)
      {
	int lb=impb(nt);
	int ub=impf(nt);
        lcat(nt)=robmean(limpc(lb,ub));
      }
      else if (fsh.between_times(i,nt))
      {
	int lb=impb(nt);
	int ub=impf(nt);
	int nmeans=(ub-lb+3)/3;
	int llb=lb-nmeans;
	dvar_vector tmp(llb,ub);
	tmp(llb,lb-1)=lcat(ntb);
	tmp(lb,ub)=limpc(lb,ub);
	ntb=nt;
        lcat(nt)=robmean(tmp);
      }
      else
      {
        lcat(nt)=lcat(nt-1);
      }
    }
    fsh.effort_dev_coffs(i)=limpc-lcat;
    {
      if (print_switch)
      {
        adstring number=str(i);
        ofstream ofs1("s."+number);
        ofstream ofs2("sm."+number);
        //ofs1 << "fishery " << i << endl;
        //ofs2 << "fishery " << i << endl;
        dvar_matrix ww(1,1,limpc.indexmin(),limpc.indexmax());
        ww(1)=exp(limpc);
        ofs1 << trans(ww) << endl;
        ww(1)=exp(lcat);
        ofs2 << trans(ww) << endl;
      }
      for (nt=2;nt<=nft;nt++)
      {
        fsh.catch_dev_coffs(i,nt)=lcat(nt)-lcat(nt-1);
      }
    }
  }
  //exit(1);
}


dvariable robfun(_CONST dvariable& m,_CONST dvar_vector& v)
{
  int n=v.indexmax()-v.indexmin()+1;
  dvariable x;
  dvar_vector d2=square(v-m);
  dvariable vhat=mean(d2);
#if !defined(NO_MY_DOUBLE_TYPE)
  x=0.5*n*log(vhat)-sum(log(exp(d2/(-1.0*vhat))+.05L));
#else
  x=0.5*n*log(vhat)-sum(log(exp(d2/(-1.0*vhat))+.05));
#endif
  return x;
}

dvariable robmp(_CONST dvariable m,_CONST dvar_vector& x)
{
  const MY_DOUBLE_TYPE a=1.0;
  const MY_DOUBLE_TYPE b=0.1;
  int n=x.indexmax()-x.indexmin()+1;
  dvar_vector r=x-m;
  dvar_vector r2=square(r);
  dvariable v=mean(r2)+1.e-10;
  dvariable psi=-2.*mean(r);
  dvar_vector phi = exp(-r2/(a*v));
  dvariable tmp;
  tmp=0.5*n/v*psi - elem_div(phi,phi+b)*r * (2.0+psi/v);
  return tmp;
}

/*
dvariable robmp(_CONST dvariable m,_CONST dvar_vector& v)
{
  dvariable x=(robfun(m+1.e-4,v)-robfun(m-1.e-4,v))/2.e-4;
  return x;
}
*/

dvariable robmpp(_CONST dvariable& m,_CONST dvar_vector& v)
{
  dvariable x=(robmp(m+2.e-4,v)-robmp(m-1.e-4,v))/3.e-4;
  return x;
}

dvariable robmean(_CONST dvar_vector& v)
{
  RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
  dvariable m=mean(v);
  //cout <<"rsh3imp.cpp " << m << endl;
  for (int i=1;i<=4;i++)
  {
    dvariable fp=robmp(m,v);
    dvariable fpp=robmpp(m,v);
    dvariable h=-robmp(m,v)/robmpp(m,v);
    m+=h; 
    //cout << m << " " << fp << " " << fpp << endl;
  }
  RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
  return m;
}

