/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define USE_DD_NOT
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

extern int global_rmax;
extern int dumflag;

int compare_switch=0;
  static void df_ddcatchimp(void);

void make_debug_info(dvar_matrix& M)
{
  dvar_matrix tmp;
  tmp.allocate(M);
  tmp=M;
  M=tmp;
}

static void newxxx(void)
{
}

void dvar_len_fish_stock_history::
  do_newton_raphson(int ir,int ip,dvariable& ffpen,dvar_vector& fmlev,
  ivector& mci,dvar_matrix& fmlevdevs,dvariable& ffpen1)
{
  //if (missing_catch_for_period_flag(ir,ip)==0  &&
  // zero_catch_for_period_flag(ir,ip)==0)
  switch (age_flags(92))
  {
  case 2:
    do_newton_raphson_with_totcatch(ir,ip,ffpen,fmlev,ffpen1);
  break;
  case 3:
    do_popes_approximation(ir,ip,ffpen,ffpen1);
   break;
  default:
    cerr << "illegal value for age_flags(92)" << endl;
    cerr << "exiting..." << endl;
    ad_exit(1);
  }
}

void dvar_len_fish_stock_history::
  do_newton_raphson(int ir,int ip,dvariable& ffpen,dvar_vector& fmlev,
  ivector& mci,dvar_matrix& fmlevdevs,ivector * pq_flag,
  dvariable& ffpen1)
{
  if (missing_catch_for_period_flag(ir,ip)==0  &&
   zero_catch_for_period_flag(ir,ip)==0)
  {
    do_newton_raphson_with_totcatch(ir,ip,ffpen,fmlev,ffpen1);
    int mmin=fmlev.indexmin();
    int mmax=fmlev.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      MY_DOUBLE_TYPE r=value(fmlev(i))/1.84077e-10;
      if ( (fabs(r-1.0) < 1.e-3) || ( ir==1 && ip==99) )
      {
        //cout << ir << "  "  << ip  
         //          << "  "  << fmlev(i) << endl;
      }
    }
  }
  else
  {
    do_newton_raphson_missing_totcatch(ir,ip,ffpen,fmlev,mci,
      fmlevdevs,pq_flag);
  }
}


vnrf::vnrf(const dvar_matrix& _fm,const dvar_vector& tm,
  const dvar_vector& surv,int _nf,int _ng,const dvar_matrix & _s,
  const dvar_vector& _M,const dvector& _Cobs,const dvar_vector& _N,
  const dvar_vector& _w,const ivector& _wtflag,MY_DOUBLE_TYPE _beta,
  const dvar_vector& _q, const dvar_fish_stock_history * _pfsh,int _ir,
  int _ip, MY_DOUBLE_TYPE _rmax,const ivector& _age_flags) : 
 F(_fm),nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
 q(_q),C(1,nf,_N.indexmin(),ng),Z(tm), S(surv),Chat(1,nf),
 wtflag(_wtflag),beta(_beta),phi(1.0),logq(1,nf),simflag(1),
 pfsh((dvar_fish_stock_history *) (_pfsh)),ir(_ir),ip(_ip), sN(1,_nf), swN(1,_nf),rmax(_rmax),
 have_weights(0),age_flags(_age_flags),pamin1(0.8),pamin2(0.9),
 omega(1,_ng),pwght1(100.0),pwght2(100.0)
{
  w=_w/1000.;
  fpen=0.0;
  fpen1=0.0;
  
  for (int fi=1;fi<=nf;fi++)
  {
    have_weights=1;
    switch (pfsh->age_flags(92))
    {
    case 0:     // add for tags  DF   May 08 08
    case 2:
      sN(fi)=elem_prod(s(fi),N);
      break;
    case 3:
    case 4:
      sN(fi)=elem_prod(s(fi),elem_prod(N,exp(-0.5*M))); // ss2 option
      break;
    default:
      cerr << "illegal age_flags(92) value" << endl;
      ad_exit(1);
    }
    if (wtflag(fi))
    {
      swN(fi)=elem_prod(w,sN(fi));
    }
  }
}

vnrf::vnrf(int _nf,int _ng,dvar_matrix _s,dvar_vector& _M,dvector _Cobs,
 dvar_vector _N,dvar_vector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
 dvar_vector _q, dvar_fish_stock_history * _pfsh,int _ir,int _ip,
 ivector _fishin) : 
 nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
 q(_q),C(1,nf,1,ng),F(1,nf,1,ng), Z(1,ng),S(1,ng),
 Chat(1,nf),
 wtflag(_wtflag),beta(_beta),phi(1.0),logq(1,nf),simflag(1),
 pfsh(_pfsh),ir(_ir),ip(_ip),sN(1,_nf),swN(1,_nf),fishin(_fishin),rmax(0.7) 
{
  w=_w/1000.;
  if (ip>pfsh->num_real_fish_periods(ir))
    rmax=3.0;

  for (int fi=1;fi<=nf;fi++)
  {
    sN(fi)=elem_prod(s(fi),N);
    if (wtflag(fi))
    {
      swN(fi)=elem_prod(w,sN(fi));
    }
  }
}

vnrfm::vnrfm(dvar_matrix& _fm,dvar_vector& tm,dvar_vector& surv,int _nf,
 int _ng,dvar_matrix & _s,dvar_vector& _M,
 dvector& _Cobs,dvar_vector& _N,dvar_vector& _w,ivector& _wtflag,MY_DOUBLE_TYPE _beta,
 dvar_vector& _q, dvar_fish_stock_history * _pfsh,int _ir,int _ip,
 MY_DOUBLE_TYPE _rmax, ivector& _age_flags,dvar_vector& mz ) : 
 vnrf(_fm,tm,surv,_nf,_ng,_s,_M,_Cobs,_N,_w,_wtflag,_beta,
 _q,_pfsh,_ir,_ip,_rmax,_age_flags) ,  missing_z(mz)
{}

void dvar_len_fish_stock_history::
  do_newton_raphson_with_totcatch(int ir,int ip,dvariable& ffpen,
  dvar_vector& q,dvariable& ffpen1)
{
  ivector ff55=column(fish_flags,55);
  ivector ff92=column(fish_flags,92);
  dvar_matrix * pfm=0;
  dvar_vector * ptm=0;
  dvar_vector * psurv=0;
  dvar_vector * pnumfish=0;
  if (af170q0==0)
  {
    pfm=&(fish_mort)(ir,ip);
    ptm=&(tot_mort)(ir,ip);
    psurv=&(survival)(ir,ip);
    pnumfish=&(num_fish)(ir,ip);
  }
  else
  {
    pfm=&(fish_mort_q0)(ir,ip);
    ptm=&(tot_mort_q0)(ir,ip);
    psurv=&(survival_q0)(ir,ip);
    pnumfish=&(num_fish_q0)(ir,ip);
  }
  MY_DOUBLE_TYPE sb=0.1;
  if (age_flags(159))
  {
    sb=age_flags(159)/10.;
  }
  MY_DOUBLE_TYPE sb1=0.1;
  dvar_vector enf=mfexp(*pnumfish);
  int nfi=num_fish_incidents(ir,ip);
  dvar_matrix sel(1,nfi,1,nage);
  dvar_matrix logsel(1,nfi,1,nage);
  dvar_vector temp_tm(1,nage);
  ivector wtflag(1,nfi);
  for (int fi=1;fi<=nfi;fi++)
  {
    //int i=parent(ir,ip,fi);
    //int rr=realization_region(i,1);
    //int rp=realization_period(i,1);
    //int ri=realization_incident(i,1);
    //logsel(fi)=incident_sel(rr,rp,ri);
    logsel(fi)=incident_sel(ir,ip,fi);
    sel(fi)=1.e-4+mfexp(logsel(fi));
    int pi=parent(ir,ip,fi);
    wtflag(fi)=data_fish_flags(1,pi);
  }

  dvar_matrix C(1,nfi,1,nage);
  MY_DOUBLE_TYPE beta=0.6;
  dvar_vector w(1,nage);
        
  if (age_flags(92)==3)
  {
    MY_DOUBLE_TYPE rmax=3.0;
    if (age_flags(116)!=0)
    {
      rmax=age_flags(116)/100.;
    }
     
    vnrf nrtester(*pfm,*ptm,*psurv,nfi,nage,sel,e_nat_mort_by_period(ir,ip),
      obs_tot_catch(ir,ip),enf,mean_weight(ir,ip,1),
      wtflag,beta,q,this,ir,ip,rmax,age_flags); 

    int printerswitch=0;
  
    ffpen+=nrtester.testnr();
  }
  else
  {
    int nfe1=get_non_zero_effort(ir,ip); //NMD_12mar2025
    int nfi1=get_num_for_nr(ir,ip);
    if (nfi1>0 || nfe1>0)
    {
      dvar_vector qtemp=ddnrf_interface(ir,ip,nfi1,wtflag);
      int ii1=0;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        int i=parent(ir,ip,fi);
        if (ff92(i)==0)   // not a survey sample
        {
          if (missing_catch_for_incident_flag(ir,ip,fi)==0)
             // && ( !af170q0 ||!ff55(i)) )
          {
            q(fi)=qtemp(++ii1);
          }
          else
          {
            q(fi)=exp(fm_level(ir,ip,fi));
          }
        }
        else
        {
           q(fi)=1.e-11;
        }
      }
    }
    else
    {
      q=1.e-11;
    }
    dvar_vector logq=log(q);
    ptm->initialize();
    for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
    {
      int i=parent(ir,ip,fi);
      int rr=realization_region(i,1);
      if (rr !=ir )
      {
        cout << "do_newton_raphson_with_totcatch - region mismatch " << endl;
        ad_exit(1);
      }
      // we don't have to do this every time !!!
      new_fm_level(i,fishery_realization_index(ir,ip,fi))=q(fi);
      (*pfm)(fi)=logq(fi)+incident_sel(ir,ip,fi);
      (*ptm)+=mfexp((*pfm)(fi));
    }
    (*ptm)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));

    temp_tm=*ptm;
    int jmin=(*ptm).indexmin();
    int jmax=(*ptm).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      if ( !(ip>num_real_fish_periods(ir) && do_fishery_projections_flag==1)
        && (*ptm)(j)>Zmax_fish)
      {
        dvariable u=(*ptm)(j)-Zmax_fish;
        dvariable v=Zmax_fish + u/(1.0 + 3.0*u/Zmax_fish);
        temp_tm(j)=v;
        Zmax_flag=1;
      }
      MY_DOUBLE_TYPE pen=100.0;
      double Zcut_mult=0.8;
      if (age_flags(189)>0)
      {
        Zcut_mult=age_flags(189)/100.0;
      }
      double Zcut=Zcut_mult*Zmax_fish;
      if ( !(ip>num_real_fish_periods(ir) && do_fishery_projections_flag==1)
        && (*ptm)(j)>Zcut)
      {
        MY_DOUBLE_TYPE pen=100.0;
        if (parest_flags(382))
        {
          pen=parest_flags(382);
        }
        ffpen1+=pen*square((*ptm)(j)-Zcut);
      }
    }
    (*ptm)=temp_tm;
    (*psurv)=mfexp(-(*ptm));
    if (0)
    {
      ofstream ofs("logfile2");
      ofs << "*psurv" << endl;
      ofs << *psurv << endl;
      ofs << "q" << endl;
      ofs << q << endl;
    }

    int i;
    for (i=1;i<=nfi;i++)
    {
      // check for non positive q
      if (q(i)<=0.0)
      {
        cout << "non positive q" << q(i) << endl;
        ad_exit(1);
      }
    }
  }
}


void dvar_len_fish_stock_history::
  do_newton_raphson_missing_totcatch(int ir,int ip,dvariable& ffpen,
  dvar_vector& fmlev,ivector& mci,dvar_matrix& fmlevdevs)
{
  dvar_matrix * pfm=0;
  dvar_vector * ptm=0;
  dvar_vector * psurv=0;
  dvar_vector * pnumfish=0;
  if (af170q0==0)
  {
    pfm=&(fish_mort)(ir,ip);
    ptm=&(tot_mort)(ir,ip);
    psurv=&(survival)(ir,ip);
    pnumfish=&(num_fish)(ir,ip);
  }
  else
  {
    pfm=&(fish_mort_q0)(ir,ip);
    ptm=&(tot_mort_q0)(ir,ip);
    psurv=&(survival_q0)(ir,ip);
    pnumfish=&(num_fish_q0)(ir,ip);
  }
  ptm->initialize();
  MY_DOUBLE_TYPE sb=0.1;
  if (age_flags(159))
  {
    sb=age_flags(159)/10.;
  }
  MY_DOUBLE_TYPE sb1=0.1;
  int nmiss=missing_catch_for_period_flag(ir,ip); 
  int nzero=zero_catch_for_period_flag(ir,ip); 
  int nfi=num_fish_incidents(ir,ip);
  int nfi1=nfi-nmiss-nzero;
  int fi;
  ivector fishin(1,nfi1);
  if (nfi1==0) 
  {
    for (fi=1;fi<=nfi;fi++)
    {
      //dvar_vector& tm=tot_mort(ir,ip);
      //dvar_matrix& fm1=fish_mort(ir,ip);

      if (zero_catch_for_incident_flag(ir,ip,fi))
      {
        fmlev(fi)=1.e-10;
      }
      else
      {
        int i=parent(ir,ip,fi);
        mci(i)++;
        dvariable tmp;
        if (af170q0==0 || q_flag(i)==0)
        {
          tmp=q0_miss(i)+effort(ir,ip,fi);
          if (mci(i)>1)
          {
            tmp+=fmlevdevs(i,mci(i));
          }
        }
        else
        {
          tmp=-20;
        }
        fmlev(fi)=exp(tmp);
        (*pfm)(fi)=tmp+incident_sel(ir,ip,fi);
        //cout << "tmp A " << tmp << endl;
        (*ptm)+=exp((*pfm)(fi));
      }
      if (!pmsd)
        (*ptm)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
      else
        (*ptm)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    }
  }
  else
  {
    dvar_vector enf=mfexp(num_fish(ir,ip));
    if (parest_flags(375))  // for two-species model att the numbers at age
    {                       // for both species
      enf+=mfexp(num_fish(ir+2,ip));  // need to get the corresponding 
    }                            // region for species 2
    dvar_vector& surv=survival(ir,ip);
    //dvariable tot_num_fish=sum(enf);
    dvar_vector& tm=*ptm;
    dvar_matrix& fm1=(*pfm);
    dvar_matrix sel(1,nfi1,1,nage);
    dvar_matrix logsel(1,nfi1,1,nage);
    dvar_matrix fm(1,nfi1);
    MY_DOUBLE_TYPE ortc=0;
    dvector otc(1,nfi1);
    int ii=1;
    //dvariable tnf=0.0;
    ivector wtflag(1,nfi1);
    dvar_vector mz(1,nage);
    mz.initialize();
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        int i=parent(ir,ip,fi);
//        int rr=realization_region(i,1);
//        int rp=realization_period(i,1);
//        int ri=realization_incident(i,1);
//        logsel(ii)=incident_sel(rr,rp,ri);
        logsel(ii)=incident_sel(ir,ip,fi);
        sel(ii)=mfexp(logsel(ii));
        //fm(ii)=fish_mort(ir,ip,fi);
        otc(ii)=obs_tot_catch(ir,ip,fi);
        //tnf+=enf(ii);
        ortc+=obs_region_tot_catch(ir,ip);
        wtflag(ii)=data_fish_flags(1,i);
        fishin(ii)=fi;
        ii++;
      }
      else if (missing_catch_for_incident_flag(ir,ip,fi)!=0)
      {
        int i=parent(ir,ip,fi);
        mci(i)++;
        int rr=realization_region(i,1);
        int rp=realization_period(i,1);
        int ri=realization_incident(i,1);

        dvariable tmp;
        if (af170q0==0 || q_flag(i)==0)
        {
          tmp=q0_miss(i)+effort(ir,ip,fi);
          if (mci(i)>1)
          {
            tmp+=fmlevdevs(i,mci(i));
          }
        }
        else
        {
          tmp=-20;
        }
        fmlev(fi)=exp(tmp);
        fm1(fi)=tmp+incident_sel(ir,ip,fi);
        mz+= exp(fm1(fi));
      }
      else if (zero_catch_for_incident_flag(ir,ip,fi)!=0)
      {
        fmlev(fi)=1.e-10;
      }
    }
  
    //dvar_matrix M(1,nfi1,1,nfi1);
    //M.initialize();
    dvar_matrix C(1,nfi1,1,nage);
    MY_DOUBLE_TYPE beta=0.6;
    dvar_vector w(1,nage);
    dvar_vector q1(1,nfi1);
    MY_DOUBLE_TYPE rmax=3.0;
    if (age_flags(116)!=0)
    {
      rmax=age_flags(116)/100.;
    }


    if (age_flags(92)==3)
    {
      vnrfm nrtester(fm,tm,surv,nfi1,nage,sel,e_nat_mort_by_period(ir,ip),otc,
        enf,mean_weight(ir,ip,1),wtflag,beta,q1,this,ir,ip,rmax,
        age_flags,mz); 
      ffpen+=nrtester.testnr();
    }
    else
    {

      // get q
      //cerr << "Not implemented for catch-conditioned version" << endl;
      //ad_exit(1);

    // ***********************************************************
    // ***********************************************************
    // ***********************************************************
      //q=ddnrf_interface(ir,ip,wtflag);
      {
        if (age_flags(116)>0)
        {
          rmax=age_flags(116)/100.;
        }
        else
        {
          rmax=3.0;
        }
        ddnrf ddnrtester(nfi1,nage,value(sel),
          value(e_nat_mort_by_period(ir,ip)),
          value(enf),otc,value(mean_weight(ir,ip,1)),
          wtflag,beta,this,ir,ip,rmax); 
    
        ddnrtester.testnr();
        dvector cq1(1,nfi1);
        for (fi=1;fi<=nfi1;fi++)
        {
    #if ((defined(USE_DD) && ((defined(__MSVC32__) &&  __MSVC32__ >=8) || defined(linux) || defined(__ADMING__) )) || defined(NO_MY_DOUBLE_TYPE))

          cq1(fi)=to_double(ddnrtester.get_q()(fi));
    #else
          cq1(fi)=ddnrtester.get_q()(fi);
    #endif
        }
        q1=nograd_assign(cq1);
    
    
        if (ddnrtester.icount>6)
        {
          cout << "A icount = " << ddnrtester.icount << endl;
        }
        //d3_array & dFdq=make_d3_array(ddnrtester.get_dFdq());
#if defined(USE_DD) || defined(NO_MY_DOUBLE_TYPE)
        const d3_array & dsdq=make_d3_array(ddnrtester.get_dsdq());
        const dmatrix & dNdq=make_dmatrix(ddnrtester.get_dNdq());
        const dmatrix & dMdq=make_dmatrix(ddnrtester.get_dMdq());
#else
        const d3_array & dsdq=ddnrtester.get_dsdq();
        const dmatrix & dNdq=ddnrtester.get_dNdq();
        const dmatrix & dMdq=ddnrtester.get_dMdq();
#endif
    
//          save_identifier_string("E18");
  const char * str1;
  str1="E18";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
        e_nat_mort_by_period(ir,ip).save_dvar_vector_position();
        
//          save_identifier_string("E17");
  const char * str2;
  str2="E17";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
        enf.save_dvar_vector_position();
        
//          save_identifier_string("E16");
  const char * str3;
  str3="E16";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
        sel.save_dvar_matrix_position();
        
//          save_identifier_string("E15");
  const char * str4;
  str4="E15";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
        dNdq.save_dmatrix_value();
        dNdq.save_dmatrix_position();
    
//          save_identifier_string("E14");
  const char * str5;
  str5="E14";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
        dMdq.save_dmatrix_value();
        dMdq.save_dmatrix_position();
    
//          save_identifier_string("E13");
  const char * str6;
  str6="E13";
  char* strx6=const_cast <char*> (str6);
  save_identifier_string(strx6);
        dsdq.save_d3_array_value();
        dsdq.save_d3_array_position();
    
        //save_identifier_string("E12");
        //dFdq.save_d3_array_value();
        //dFdq.save_d3_array_position();
    
//          save_identifier_string("E11");
  const char * str8;
  str8="E11";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
        q1.save_dvar_vector_position();
    
//          save_identifier_string("E10");
  const char * str9;
  str9="E10";
  char* strx9=const_cast <char*> (str9);
  save_identifier_string(strx9);
        gradient_structure::GRAD_STACK1->set_gradient_stack(df_ddcatchimp);
        //return v;
      }
    // ***********************************************************
    // ***********************************************************
    // ***********************************************************
  
  
      vnrfm nrtester(fm,tm,surv,nfi1,nage,sel,e_nat_mort_by_period(ir,ip),
        obs_tot_catch(ir,ip),enf,mean_weight(ir,ip,1),
        wtflag,beta,q1,this,ir,ip,rmax,age_flags,mz); 
  
      // get other parameters besides q
      ffpen+=nrtester.testnrx();
     
    }

    ii=1;
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        fm1(fi)=fm(ii);
        FF(ir,ip,fi)=log(q1(ii));
        fmlev(fi)=q1(ii);
        ii++;
      }
      else
      {
        if (zero_catch_for_incident_flag(ir,ip,fi)==0)
        {
        
          if (af170q0==0 || q_flag(parent(ir,ip,fi))==0)
          {
            FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
          }
          else
          {
            FF(ir,ip,fi)=-10;
          }
        }
        else
        {
          FF(ir,ip,fi)=-10.;
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::
  do_newton_raphson_missing_totcatch_fml(int ir,int ip,dvariable& ffpen,
  dvar_vector& fmlev,ivector& mci,dvar_matrix& fmlevdevs)
{
  dvar_matrix * pfm=0;
  dvar_vector * ptm=0;
  dvar_vector * psurv=0;
  dvar_vector * pnumfish=0;
  if (af170q0==0)
  {
    pfm=&(fish_mort)(ir,ip);
    ptm=&(tot_mort)(ir,ip);
    psurv=&(survival)(ir,ip);
    pnumfish=&(num_fish)(ir,ip);
  }
  else
  {
    pfm=&(fish_mort_q0)(ir,ip);
    ptm=&(tot_mort_q0)(ir,ip);
    psurv=&(survival_q0)(ir,ip);
    pnumfish=&(num_fish_q0)(ir,ip);
  }
  ptm->initialize();
  MY_DOUBLE_TYPE sb=0.1;
  if (age_flags(159))
  {
    sb=age_flags(159)/10.;
  }
  MY_DOUBLE_TYPE sb1=0.1;
  int nmiss=missing_catch_for_period_flag(ir,ip); 
  int nzero=zero_catch_for_period_flag(ir,ip); 
  int nfi=num_fish_incidents(ir,ip);
  int nfi1=nfi-nmiss-nzero;
  int fi;
  ivector fishin(1,nfi1);
  if (nfi1==0) 
  {
    for (fi=1;fi<=nfi;fi++)
    {
      //dvar_vector& tm=tot_mort(ir,ip);
      //dvar_matrix& fm1=fish_mort(ir,ip);

      if (zero_catch_for_incident_flag(ir,ip,fi))
      {
        fmlev(fi)=1.e-10;
      }
      else
      {
        int i=parent(ir,ip,fi);
        mci(i)++;
        dvariable tmp;
        if (af170q0==0 || q_flag(i)==0)
        {
          tmp=q0_miss(i)+effort(ir,ip,fi);
          if (mci(i)>1)
          {
            tmp+=fmlevdevs(i,mci(i));
          }
        }
        else
        {
          tmp=-20;
        }
        fmlev(fi)=exp(tmp);
        (*pfm)(fi)=tmp+incident_sel(ir,ip,fi);
        //cout << "tmp A " << tmp << endl;
        (*ptm)+=exp((*pfm)(fi));
      }
      if (!pmsd)
        (*ptm)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
      else
        (*ptm)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    }
  }
  else
  {
    dvar_vector enf=mfexp(num_fish(ir,ip));
    dvar_vector& surv=survival(ir,ip);
    //dvariable tot_num_fish=sum(enf);
    dvar_vector& tm=*ptm;
    dvar_matrix& fm1=(*pfm);
    dvar_matrix sel(1,nfi1,1,nage);
    dvar_matrix logsel(1,nfi1,1,nage);
    dvar_matrix fm(1,nfi1);
    MY_DOUBLE_TYPE ortc=0;
    dvector otc(1,nfi1);
    int ii=1;
    //dvariable tnf=0.0;
    ivector wtflag(1,nfi1);
    dvar_vector mz(1,nage);
    mz.initialize();
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        int i=parent(ir,ip,fi);
        logsel(ii)=incident_sel(ir,ip,fi);
        sel(ii)=mfexp(logsel(ii));
        //fm(ii)=fish_mort(ir,ip,fi);
        otc(ii)=obs_tot_catch(ir,ip,fi);
        //tnf+=enf(ii);
        ortc+=obs_region_tot_catch(ir,ip);
        wtflag(ii)=data_fish_flags(1,i);
        fishin(ii)=fi;
        ii++;
      }
      else if (missing_catch_for_incident_flag(ir,ip,fi)!=0)
      {
        int i=parent(ir,ip,fi);
        dvector dv=fml_designvec(ir,ip,fi);
        ivector num_ifmlrp=implicit_fml_bounds(3);
        dvar_vector ests=implicit_fm_level_regression_pars(i)(1,num_ifmlrp(i));
        MY_DOUBLE_TYPE leff=log_effort_by_fishery(i,fish_times(ir,ip,fi));
        dvariable tmp2 = dv * inv(fml_R(i)) * ests;
        dvariable log_fmlevel= tmp2+leff;
        cout << log_fmlevel << " " << exp(log_fmlevel)  << endl;
        dvar_vector logfm=log_fmlevel+incident_sel(ir,ip,fi);
        cout << "logfm" << endl << logfm << endl << " exp(logfm)" << endl   << exp(logfm) << endl;
        fm1(fi)=log_fmlevel+incident_sel(ir,ip,fi);
        mz+= exp(fm1(fi));
      }
    }
  
    //dvar_matrix M(1,nfi1,1,nfi1);
    //M.initialize();
    dvar_matrix C(1,nfi1,1,nage);
    MY_DOUBLE_TYPE beta=0.6;
    dvar_vector w(1,nage);
    dvar_vector q1(1,nfi1);
    MY_DOUBLE_TYPE rmax=3.0;
    if (age_flags(116)!=0)
    {
      rmax=age_flags(116)/100.;
    }

    if (age_flags(92)==3)
    {
      vnrfm nrtester(fm,tm,surv,nfi1,nage,sel,e_nat_mort_by_period(ir,ip),otc,
        enf,mean_weight(ir,ip,1),wtflag,beta,q1,this,ir,ip,rmax,
        age_flags,mz); 
      ffpen+=nrtester.testnr();
    }
    else
    {

      // get q
      //cerr << "Not implemented for catch-conditioned version" << endl;
      //ad_exit(1);

    // ***********************************************************
    // ***********************************************************
    // ***********************************************************
      //q=ddnrf_interface(ir,ip,wtflag);
      {
        if (age_flags(116)>0)
        {
          rmax=age_flags(116)/100.;
        }
        else
        {
          rmax=3.0;
        }
        ddnrf ddnrtester(nfi1,nage,value(sel),
          value(e_nat_mort_by_period(ir,ip)),
          value(enf),otc,value(mean_weight(ir,ip,1)),
          wtflag,beta,this,ir,ip,rmax); 
    
        ddnrtester.testnr();
        dvector cq1(1,nfi1);
        for (fi=1;fi<=nfi1;fi++)
        {
    #if ((defined(USE_DD) && ((defined(__MSVC32__) &&  __MSVC32__ >=8) || defined(linux) || defined(__ADMING__) )) || defined(NO_MY_DOUBLE_TYPE))
 
          cq1(fi)=to_double(ddnrtester.get_q()(fi));
    #else
          cq1(fi)=ddnrtester.get_q()(fi);
    #endif
        }
        q1=nograd_assign(cq1);
    
    
        if (ddnrtester.icount>6)
        {
          cout << "A icount = " << ddnrtester.icount << endl;
        }
        //d3_array & dFdq=make_d3_array(ddnrtester.get_dFdq());
#if defined(USE_DD) || defined(NO_MY_DOUBLE_TYPE)
        const d3_array & dsdq=make_d3_array(ddnrtester.get_dsdq());
        const dmatrix & dNdq=make_dmatrix(ddnrtester.get_dNdq());
        const dmatrix & dMdq=make_dmatrix(ddnrtester.get_dMdq());
#else
        const d3_array & dsdq=ddnrtester.get_dsdq();
        const dmatrix & dNdq=ddnrtester.get_dNdq();
        const dmatrix & dMdq=ddnrtester.get_dMdq();
#endif
    
//          save_identifier_string("E18");
  const char * str10;
  str10="E18";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
        e_nat_mort_by_period(ir,ip).save_dvar_vector_position();
        
//          save_identifier_string("E17");
  const char * str11;
  str11="E17";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
        enf.save_dvar_vector_position();
        
//          save_identifier_string("E16");
  const char * str12;
  str12="E16";
  char* strx12=const_cast <char*> (str12);
  save_identifier_string(strx12);
        sel.save_dvar_matrix_position();
        
//          save_identifier_string("E15");
  const char * str13;
  str13="E15";
  char* strx13=const_cast <char*> (str13);
  save_identifier_string(strx13);
        dNdq.save_dmatrix_value();
        dNdq.save_dmatrix_position();
    
//          save_identifier_string("E14");
  const char * str14;
  str14="E14";
  char* strx14=const_cast <char*> (str14);
  save_identifier_string(strx14);
        dMdq.save_dmatrix_value();
        dMdq.save_dmatrix_position();
    
//          save_identifier_string("E13");
  const char * str15;
  str15="E13";
  char* strx15=const_cast <char*> (str15);
  save_identifier_string(strx15);
        dsdq.save_d3_array_value();
        dsdq.save_d3_array_position();
    
        //save_identifier_string("E12");
        //dFdq.save_d3_array_value();
        //dFdq.save_d3_array_position();
    
//          save_identifier_string("E11");
  const char * str17;
  str17="E11";
  char* strx17=const_cast <char*> (str17);
  save_identifier_string(strx17);
        q1.save_dvar_vector_position();
    
//          save_identifier_string("E10");
  const char * str18;
  str18="E10";
  char* strx18=const_cast <char*> (str18);
  save_identifier_string(strx18);
        gradient_structure::GRAD_STACK1->set_gradient_stack(df_ddcatchimp);
        //return v;
      }
    // ***********************************************************
    // ***********************************************************
    // ***********************************************************
  
  
      vnrfm nrtester(fm,tm,surv,nfi1,nage,sel,e_nat_mort_by_period(ir,ip),
        obs_tot_catch(ir,ip),enf,mean_weight(ir,ip,1),
        wtflag,beta,q1,this,ir,ip,rmax,age_flags,mz); 
  
      // get other parameters besides q
      ffpen+=nrtester.testnrx();
     
    }

    ii=1;
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        fm1(fi)=fm(ii);
        FF(ir,ip,fi)=log(q1(ii));
        fmlev(fi)=q1(ii);
        ii++;
      }
      else
      {
        if (zero_catch_for_incident_flag(ir,ip,fi)==0)
        {
        
          if (af170q0==0 || q_flag(parent(ir,ip,fi))==0)
          {
            FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
          }
          else
          {
            FF(ir,ip,fi)=-10;
          }
        }
        else
        {
          FF(ir,ip,fi)=-10.;
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::
  do_newton_raphson_missing_totcatch(int ir,int ip,dvariable& ffpen,
  dvar_vector& fmlev,ivector& mci,dvar_matrix& fmlevdevs,ivector * pq_flag)
{
  dvar_matrix * pfm=0;
  dvar_vector * ptm=0;
  dvar_vector * psurv=0;
  dvar_vector * pnumfish=0;
  if (af170q0==0)
  {
    pfm=&(fish_mort)(ir,ip);
    ptm=&(tot_mort)(ir,ip);
    psurv=&(survival)(ir,ip);
    pnumfish=&(num_fish)(ir,ip);
  }
  else
  {
    pfm=&(fish_mort_q0)(ir,ip);
    ptm=&(tot_mort_q0)(ir,ip);
    psurv=&(survival_q0)(ir,ip);
    pnumfish=&(num_fish_q0)(ir,ip);
  }
  ptm->initialize();
  MY_DOUBLE_TYPE sb=0.1;
  if (age_flags(159))
  {
    sb=age_flags(159)/10.;
  }
  MY_DOUBLE_TYPE sb1=0.1;
  int nmiss=missing_catch_for_period_flag(ir,ip); 
  int nzero=zero_catch_for_period_flag(ir,ip); 
  int nfi=num_fish_incidents(ir,ip);
  int nfi1=nfi-nmiss-nzero;
  int fi;
  ivector fishin(1,nfi1);
  if (nfi1==0) 
  {
    for (fi=1;fi<=nfi;fi++)
    {
      //dvar_vector& tm=tot_mort(ir,ip);
      //dvar_matrix& fm1=fish_mort(ir,ip);

      if (zero_catch_for_incident_flag(ir,ip,fi))
      {
        fmlev(fi)=1.e-10;
      }
      else
      {
        int i=parent(ir,ip,fi);
        mci(i)++;
        dvariable tmp;
        if (!(*pq_flag)(i))
        {
          tmp=q0_miss(i)+effort(ir,ip,fi);
        }
        else
        {
          tmp=-20.0;
        }
        if (mci(i)>1)
        {
          tmp+=fmlevdevs(i,mci(i));
        }
        fmlev(fi)=exp(tmp);
        (*pfm)(fi)=tmp+incident_sel(ir,ip,fi);
        //cout << "tmp B " << tmp << endl;
        (*ptm)+=exp((*pfm)(fi));
      }
      if (!pmsd)
        (*ptm)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
      else
        (*ptm)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    }
  }
  else
  {
    const dvar_vector& _enf=mfexp(num_fish(ir,ip));
	dvar_vector& enf = (dvar_vector&) (_enf);
    dvar_vector& surv=survival(ir,ip);
    //dvariable tot_num_fish=sum(enf);
    dvar_vector& tm=*ptm;
    dvar_matrix& fm1=(*pfm);
    dvar_matrix sel(1,nfi1,1,nage);
    dvar_matrix logsel(1,nfi1,1,nage);
    dvar_matrix fm(1,nfi1);
    MY_DOUBLE_TYPE ortc=0;
    dvector otc(1,nfi1);
    int ii=1;
    //dvariable tnf=0.0;
    ivector wtflag(1,nfi1);
    dvar_vector mz(1,nage);
    mz.initialize();
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        int i=parent(ir,ip,fi);
        logsel(ii)=incident_sel(ir,ip,fi);
        sel(ii)=mfexp(logsel(ii));
        //fm(ii)=fish_mort(ir,ip,fi);
        otc(ii)=obs_tot_catch(ir,ip,fi);
        //tnf+=enf(ii);
        ortc+=obs_region_tot_catch(ir,ip);
        wtflag(ii)=data_fish_flags(1,i);
        fishin(ii)=fi;
        ii++;
      }
      else if (missing_catch_for_incident_flag(ir,ip,fi)!=0)
      {
        int i=parent(ir,ip,fi);
        mci(i)++;

        dvariable tmp;
        if (!(*pq_flag)(i))
        {
          tmp=q0_miss(i)+effort(ir,ip,fi);
        }
        else
        {
          tmp=-20.0+effort(ir,ip,fi);
        }
        if (mci(i)>1)
        {
          tmp+=fmlevdevs(i,mci(i));
        }
        fmlev(fi)=exp(tmp);
        fm1(fi)=tmp+incident_sel(ir,ip,fi);
        mz+= exp(fm1(fi));
      }
      else if (zero_catch_for_incident_flag(ir,ip,fi)!=0)
      {
        fmlev(fi)=1.e-10;
      }
    }
  
    //dvar_matrix M(1,nfi1,1,nfi1);
    //M.initialize();
    dvar_matrix C(1,nfi1,1,nage);
    MY_DOUBLE_TYPE beta=0.6;
    dvar_vector w(1,nage);
    dvar_vector q1(1,nfi1);
    MY_DOUBLE_TYPE rmax=3.0;
    if (age_flags(116)!=0)
    {
      rmax=age_flags(116)/100.;
    }

    if (age_flags(92)==3)
    {
      vnrfm nrtester(fm,tm,surv,nfi1,nage,sel,e_nat_mort_by_period(ir,ip),otc,
        enf,mean_weight(ir,ip,1),wtflag,beta,q1,this,ir,ip,rmax,
        age_flags,mz); 
      ffpen+=nrtester.testnr();
    }
    else
    {

      // get q
      cerr << "Not imlemented for catch-conditioned version" << endl;
      //ad_exit(1);

    // ***********************************************************
    // ***********************************************************
    // ***********************************************************
      //q=ddnrf_interface(ir,ip,wtflag);
      {
        ddnrf ddnrtester(nfi1,nage,value(sel),
          value(e_nat_mort_by_period(ir,ip)),
          value(enf),otc,value(mean_weight(ir,ip,1)),
          wtflag,beta,this,ir,ip,rmax); 
    
        ddnrtester.testnr();
        dvector cq1(1,nfi1);
        for (fi=1;fi<=nfi1;fi++)
        {
    #if ((defined(USE_DD) && ((defined(__MSVC32__) &&  __MSVC32__ >=8) || defined(linux) || defined(__ADMING__) )) || defined(NO_MY_DOUBLE_TYPE))
          cq1(fi)=to_double(ddnrtester.get_q()(fi));
    #else
          cq1(fi)=ddnrtester.get_q()(fi);
    #endif
        }
        q1=nograd_assign(cq1);
    
    
        if (ddnrtester.icount>6)
        {
          cout << "B icount = " << ddnrtester.icount << endl;
        }
        //d3_array & dFdq=make_d3_array(ddnrtester.get_dFdq());
#if defined(USE_DD) || defined(NO_MY_DOUBLE_TYPE)
        const d3_array & dsdq=make_d3_array(ddnrtester.get_dsdq());
        const dmatrix & dNdq=make_dmatrix(ddnrtester.get_dNdq());
        const dmatrix & dMdq=make_dmatrix(ddnrtester.get_dMdq());
#else
        const d3_array & dsdq=ddnrtester.get_dsdq();
        const dmatrix & dNdq=ddnrtester.get_dNdq();
        const dmatrix & dMdq=ddnrtester.get_dMdq();
#endif
    
//          save_identifier_string("E18");
  const char * str19;
  str19="E18";
  char* strx19=const_cast <char*> (str19);
  save_identifier_string(strx19);
        e_nat_mort_by_period(ir,ip).save_dvar_vector_position();
        
//          save_identifier_string("E17");
  const char * str20;
  str20="E17";
  char* strx20=const_cast <char*> (str20);
  save_identifier_string(strx20);
        enf.save_dvar_vector_position();
        
//          save_identifier_string("E16");
  const char * str21;
  str21="E16";
  char* strx21=const_cast <char*> (str21);
  save_identifier_string(strx21);
        sel.save_dvar_matrix_position();
        
//          save_identifier_string("E15");
  const char * str22;
  str22="E15";
  char* strx22=const_cast <char*> (str22);
  save_identifier_string(strx22);
        dNdq.save_dmatrix_value();
        dNdq.save_dmatrix_position();
    
//          save_identifier_string("E14");
  const char * str23;
  str23="E14";
  char* strx23=const_cast <char*> (str23);
  save_identifier_string(strx23);
        dMdq.save_dmatrix_value();
        dMdq.save_dmatrix_position();
    
//          save_identifier_string("E13");
  const char * str24;
  str24="E13";
  char* strx24=const_cast <char*> (str24);
  save_identifier_string(strx24);
        dsdq.save_d3_array_value();
        dsdq.save_d3_array_position();
    
        //save_identifier_string("E12");
        //dFdq.save_d3_array_value();
        //dFdq.save_d3_array_position();
    
//          save_identifier_string("E11");
  const char * str26;
  str26="E11";
  char* strx26=const_cast <char*> (str26);
  save_identifier_string(strx26);
        q1.save_dvar_vector_position();
    
//          save_identifier_string("E10");
  const char * str27;
  str27="E10";
  char* strx27=const_cast <char*> (str27);
  save_identifier_string(strx27);
        gradient_structure::GRAD_STACK1->set_gradient_stack(df_ddcatchimp);
        //return v;
      }
  
      vnrfm nrtester(fm,tm,surv,nfi1,nage,sel,e_nat_mort_by_period(ir,ip),
        obs_tot_catch(ir,ip),enf,mean_weight(ir,ip,1),
        wtflag,beta,q1,this,ir,ip,rmax,age_flags,mz); 
  
      // get other parameters besides q
      ffpen+=nrtester.testnrx();
     
    }

    ii=1;
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        fm1(fi)=fm(ii);
        FF(ir,ip,fi)=log(q1(ii));
        fmlev(fi)=q1(ii);
        ii++;
      }
      else
      {
        if (zero_catch_for_incident_flag(ir,ip,fi)==0)
        {
          FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
        }
        else
        {
          FF(ir,ip,fi)=-10.;
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::
  do_newton_raphson_missing_totcatch(int ir,int ip,dvariable& ffpen,
    dvar_vector& fmlev)
{
  dvar3_array * pnum_fish=0;
  dvar3_array * pN=0;
  dvar3_array * ptm=0;
  dvar4_array * pfm=0;
  dvar_vector * psurv=0;
  //if (!pq_flag)
  if (af170q0==0)
  {
    pnum_fish=&num_fish;
    ptm=&tot_mort;
    pN=&N;
    pfm=&fish_mort;
    psurv=&(survival)(ir,ip);
  }
  else
  {
    ptm=&tot_mort_q0;
    pnum_fish=&num_fish_q0;
    pN=&N_q0;
    pfm=&fish_mort_q0;
    psurv=&(survival_q0)(ir,ip);
  }
  MY_DOUBLE_TYPE sb=0.1;
  if (age_flags(159))
  {
    sb=age_flags(159)/10.;
  }
  MY_DOUBLE_TYPE sb1=0.1;
  int nmiss=missing_catch_for_period_flag(ir,ip); 
  int nzero=zero_catch_for_period_flag(ir,ip); 
  int nfi=num_fish_incidents(ir,ip);
  int nfi1=nfi-nmiss-nzero;
  int fi;
  ivector fishin(1,nfi1);
  if (nfi1==0) 
  {
    for (fi=1;fi<=nfi;fi++)
    {
      if (zero_catch_for_incident_flag(ir,ip,fi))
      {
        FF(ir,ip,fi)=-10.;
      }
      else
      {
        FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
      }
    }
  }
  else
  {
    const dvar_vector& enf=mfexp((*pnum_fish)(ir,ip));
    dvariable tot_num_fish=sum(enf);
    dvar_matrix sel(1,nfi1,1,nage);
    dvar_matrix logsel(1,nfi1,1,nage);
    dvar_matrix fm(1,nfi1);
    MY_DOUBLE_TYPE ortc=0;
    dvector otc(1,nfi1);
    int ii=1;
    //dvariable tnf=0.0;
    ivector wtflag(1,nfi1);
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        int i=parent(ir,ip,fi);
        logsel(ii)=incident_sel(ir,ip,fi);
        sel(ii)=mfexp(logsel(ii));
        fm(ii)=fish_mort(ir,ip,fi);
        otc(ii)=obs_tot_catch(ir,ip,fi);
        //tnf+=enf(ii);
        ortc+=obs_region_tot_catch(ir,ip);
        wtflag(ii)=data_fish_flags(1,i);
        fishin(ii)=fi;
        ii++;
      }
    }
  
    //dvar_matrix M(1,nfi1,1,nfi1);
    //M.initialize();
    dvar_matrix C(1,nfi1,1,nage);
    MY_DOUBLE_TYPE beta=0.6;
    //dvar_vector w(1,nage);
    dvar_vector q(1,nfi1);

    dvar_vector enmp;
    if (!pmsd)
        enmp=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
      else
        enmp=mfexp(get_nat_mort_region(ir)(year(ir,ip))
           +fraction(ir,ip));
    for (fi=1;fi<=nfi;fi++)
    {
      if( !(missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0) )
      {
        enmp+=exp((*pfm)(ir,ip,fi));
        //cout << "removed for test" << endl;
      }
    }
    vnrf nrtester(nfi1,nage,sel,enmp,otc,
      enf,mean_weight(ir,ip,1),wtflag,beta,q,this,ir,ip,fishin); 

    nrtester.get_initial_q();
    nrtester.calculate_F(); 
    nrtester.calculate_Z(); 
    nrtester.calculate_S(); 


    nrtester.testnr();

    ii=1;
    (*ptm)(ir,ip)=enmp;  // add fishing mortality for missing catch fisheries 
    //cout << "ptm A " << (*ptm)(ir,ip) << endl;
                // and nat mort for this period
    for (fi=1;fi<=nfi;fi++)
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0
        && zero_catch_for_incident_flag(ir,ip,fi)==0)
      {
        (*pfm)(ir,ip,fi)=nrtester.get_F()(ii);
        
        //cout << " pfm " << (*pfm)(ir,ip,fi) << endl;
        (*ptm)(ir,ip)+=exp((*pfm)(ir,ip,fi));  // add fish mort for this 
                                      // fishery to tmort
        FF(ir,ip,fi)=log(q(ii));
        fmlev(fi)=q(ii);  //NMD_1jul2020
        ii++;
      }
      else
      {
        if (zero_catch_for_incident_flag(ir,ip,fi)==0)
        {
          FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
        }
        else
        {
          FF(ir,ip,fi)=-10.;
        }
      }
    }
    (*psurv)=exp(-(*ptm)(ir,ip));
  }
}

extern MY_DOUBLE_TYPE rmax_penalty_weight; 
extern MY_DOUBLE_TYPE rmax_multiplier;
  
 dvariable vnrf::testnrx(void)
 {
   dvariable fpen=0.0;
   int fi;
   int log_flag=0;
   calculate_F(); 
   calculate_Z(); 
   calculate_S(); 
   int mmin=Z.indexmin();
   int mmax=Z.indexmax();
   for (int j=mmin;j<=mmax;j++)
   {
     if (value(Z(j))>rmax_multiplier*rmax)
     {
       fpen+=rmax_penalty_weight*square(square(log(Z(j)/(rmax_multiplier*rmax))));
     }
     
     if (value(Z(j))>rmax)
     {
       age_flags(181)=2;
       //fpen+=100.*square(Z(j)-rmax);
       dvariable dd=Z(j)-rmax;
       dvariable ppen=0.0;
       dvariable Zstar=rmax-posfun(0.2-dd,0.1,ppen)+0.2;

       dvariable lambda;
       if (pfsh->age_flags(92)!=3 && pfsh->age_flags(92)!=4)
       {
         if (Zstar<=M(j))
         {
           cerr << "Error in newton-raphson missing totcatch"
                << endl << " Zstar is less than M(j))"
                << endl << " You need to increase Zstar " << endl
                << " at present Zstar = " << Zstar  << " and M(j) = "
                << M(j) << endl;
           ad_exit(1);
         }
         lambda=(Zstar-M(j))/(Z(j)-M(j));
         dvariable ppen2=0.0;
         lambda=posfun(lambda,0.01,ppen2);
         ppen+=ppen2;
       }
       else
         lambda=Zstar/Z(j);

       for (fi=1;fi<=nf;fi++)
       {
         F(fi,j)*=lambda;

         if (value(F(fi,j))<0.0) 
           cout << "neg F " << value(F(fi,j)) << endl;
       }
       Z(j)=Zstar;
     }
   }
   if (pfsh->age_flags(92)!=3)
   {
     if (pfsh->age_flags(92)==4) // for Pope need to add M to get true Z
     {
       Z+=M;
     }
     S=exp(-Z);
   }
   else
   {
     Z=M-log(1.0-Z);
     S=exp(-Z);
   }
   // check for <= F
   {
     int mmin=F.indexmin();
     int mmax=F.indexmax();
     for (int i=mmin;i<=mmax;i++)
     {
       int cmin=F(i).indexmin();
       int cmax=F(i).indexmax();
       for (int j=cmin;j<=cmax;j++)
       {
         if (F(i,j)<0.0)
         {
           cout << "F(" << i << "," << j << ") < 0 in nrcatch3.cpp " << F(i,j) << endl;
         }
       }
     }

   }
   
   F=log(1.e-20+F);
   return fpen;
 }
  
 dvariable vnrf::testnr(void)
 {
   dvariable fpen=0.0;
   int fi;
   int log_flag=0;
   
   MY_DOUBLE_TYPE ndvalues[42];
   if (!log_flag)
   {
     icount=1;
     if(pfsh->age_flags(92) !=3)  // check for ss2 option
     {
       do
       {
         // newton raphson loop
         dvar_vector d=calculate_d();
 
         if (norm(value(q))> 1.e+4)
         {
           cerr << "nractch3.cpp Not enough fish for projection in region "
                << endl << ir <<  "  period  " << ip << endl;
           //break;
         }
         dvariable nd=norm(d);
         ndvalues[icount]=value(nd);
         if (nd<1.e-10) 
         {
           //cout << ir << " " << ip << " " << icount  << endl;
           break;
         }
         if (icount == 8 ) 
         {
           //break;
         }
 
         if (icount > 40 ) 
         {
           cout << ir << " " << ip << " " << icount  << " " 
                << setscientific() << setprecision(4) << nd << endl;
           break;
         }
 
         
         //if (pfsh->age_flags(92)==4)
         //  J=calculate_J();
         //else
         dvar_matrix J=calculate_J_vanal();
         
         if (compare_switch)
         {
           dvar_matrix J1=calculate_J_vanal();
           cout <<"nrcatch3.cpp " << setprecision(12) << J << endl;
           cout <<"nrcatch3.cpp " << setprecision(12) << J1 << endl;
           cout <<"nrcatch3.cpp " << setprecision(12) << J-J1 << endl;
         }
        
         //double detJ=det(value(J));
         //if (fabs(detJ<=1.e-14))
         //  cout << " detJ = " << detJ << endl;
         dvar_vector h=solve(J,d);
         MY_DOUBLE_TYPE cd=norm(value(J)*value(h)-value(d));
         if (cd>1.e-12)
         {
           int icount=0;
           do 
           {
             dvar_vector delta=J*h-d;
             h-=solve(J,delta);
             cd=norm(value(J)*value(h)-value(d));
             //cout <<  " cd = " << setscientific() << setprecision(5) 
             //   << cd << endl;
             icount++;
           }
           while(cd>1.e-12 && icount<3);
         }
         int mflag=0;
         do
         {
           q-=h;
           
           //cout << " q " << q << endl;

           for (fi=1;fi<=nf;fi++)
           {
             mflag=0;
             if (q(fi)<0) 
             {
               q+=h;
               h/=2.0;
               mflag=1;
               break;
             }
           }
         }
         while(mflag);
         
         icount++;
       }
       while(1);
     } 
     else
     {
       get_initial_q();
       calculate_F(); 
       calculate_Z(); 
       calculate_S(); 
     }
   }
   else
   {
      get_initial_q();
      return fpen;
      get_initial_q();
      icount=1;
      do
      {
        // newton raphson loop
        dvar_vector d=calculate_ld();
        dvariable nd=norm(d);
        if (nd<1.e-8) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 5)
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 15 && nd < 1.e-6) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 25 && nd < 1.e-4) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 40 ) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        dvar_matrix J=calculate_lJ();

        dvar_vector h=solve(trans(J),d);
        logq-=h;
        
        icount++;
      }
      while(1);
    }

    int mmin=Z.indexmin();
    int mmax=Z.indexmax();
    for (int j=mmin;j<=mmax;j++)
    {
      if (value(Z(j))>rmax_multiplier*rmax)
      {
        fpen+=rmax_penalty_weight*square(log(Z(j)/(rmax_multiplier*rmax)));
      }
      if (value(Z(j))>rmax)
      {
        //fpen+=100.*square(Z(j)-rmax);
        dvariable dd=Z(j)-rmax;
        dvariable ppen=0.0;
        dvariable Zstar=rmax-posfun(0.2-dd,0.1,ppen)+0.2;

        dvariable lambda;
        if (pfsh->age_flags(92)!=3 && pfsh->age_flags(92)!=4)
        {
          if (Zstar<=M(j))
          {
            cerr << "Error in newton-raphson missing totcatch"
                 << endl << " Zstar is less than M(j))"
                 << endl << " You need to increase rmax " << endl
                 << " at present Zstar = " << Zstar  << " and M(j) = "
                 << M(j) << endl;
            ad_exit(1);
          }
          lambda=(Zstar-M(j))/(Z(j)-M(j));
          if (value(lambda)<0.0)
          {
            cout << " lambda = " << lambda << " "
                 << " Zstar = " << Zstar << " "
                 << " M(j) = " << M(j) << endl;
          }
        }
        else
        {
          lambda=Zstar/Z(j);
          if (value(lambda)<0.0)
          {
            cout << " Zstar = " << Zstar << " "
                 << " Z(j) = " << Z(j) << endl;
          }
        }


        for (fi=1;fi<=nf;fi++)
        {
          F(fi,j)*=lambda;
        }
        if (Zstar > 0.99999)
        {
          cout << "CC " << Zstar << endl;
        }
        Z(j)=Zstar;
      }
    }
    if (pfsh->age_flags(92)!=3)
    {
      if (pfsh->age_flags(92)==4) // for Pope need to add M to get true Z
      {
        Z+=M;
      }
      S=exp(-Z);
    }
    else
    {
      Z=M-log(1.0-Z);
      S=exp(-Z);
    }
    F=log(1.e-20+F);
    return fpen;
  }
  
 dvariable vnrf::testnry(void)
 {
   dvariable fpen=0.0;
   int fi;
   int log_flag=0;
   
   MY_DOUBLE_TYPE ndvalues[42];
   icount=1;
   get_initial_q_ss2();
   do
   {
     // newton raphson loop
     dvar_vector d=calculate_d();
 
     if (norm(value(q))> 1.e+4)
     {
       cerr << "not enough fish for projection in region "
            << endl << ir <<  "  period  " << ip << endl;
       //break;
     }
     dvariable nd=norm(d);
     ndvalues[icount]=value(nd);
     if (nd<1.e-10) 
     {
       //cout << ir << " " << ip << " " << icount  << endl;
       break;
     }
     if (icount == 8 ) 
     {
       //break;
     }
 
     if (icount > 40 ) 
     {
       cout << ir << " " << ip << " " << icount  << " " 
            << setscientific() << setprecision(4) << nd << endl;
       break;
     }
 
     
     //if (pfsh->age_flags(92)==4)
     //  J=calculate_J();
     //else
     dvar_matrix J=calculate_J_vanal();
     
     if (compare_switch)
     {
       dvar_matrix J1=calculate_J_vanal();
       cout <<"nrcatch3.cpp " << setprecision(12) << J << endl;
       cout <<"nrcatch3.cpp " << setprecision(12) << J1 << endl;
       cout <<"nrcatch3.cpp " << setprecision(12) << J-J1 << endl;
     }
    
     //double detJ=det(value(J));
     //if (fabs(detJ<=1.e-14))
     //  cout << " detJ = " << detJ << endl;
     dvar_vector h=solve(J,d);
     MY_DOUBLE_TYPE cd=norm(value(J)*value(h)-value(d));
     if (cd>1.e-12)
     {
       int icount=0;
       do 
       {
         dvar_vector delta=J*h-d;
         h-=solve(J,delta);
         cd=norm(value(J)*value(h)-value(d));
         //cout <<  " cd = " << setscientific() << setprecision(5) 
         //   << cd << endl;
         icount++;
       }
       while(cd>1.e-12 && icount<3);
     }
     int mflag=0;
     do
     {
       q-=h;
       
       //cout << " q " << q << endl;

       for (fi=1;fi<=nf;fi++)
       {
         mflag=0;
         if (q(fi)<0) 
         {
           q+=h;
           h/=2.0;
           mflag=1;
           break;
         }
       }
     }
     while(mflag);
     
     icount++;
   }
   while(1);

    int mmin=Z.indexmin();
    int mmax=Z.indexmax();
    for (int j=mmin;j<=mmax;j++)
    {
      if (value(Z(j))>rmax_multiplier*rmax)
      {
        fpen+=rmax_penalty_weight*square(log(Z(j)/(rmax_multiplier*rmax)));
      }
      if (value(Z(j))>rmax)
      {
        //fpen+=100.*square(Z(j)-rmax);
        dvariable dd=Z(j)-rmax;
        dvariable ppen=0.0;
        dvariable Zstar=rmax-posfun(0.2-dd,0.1,ppen)+0.2;

        dvariable lambda;
        if (pfsh->age_flags(92)!=3 && pfsh->age_flags(92)!=4)
        {
          if (Zstar<=M(j))
          {
            cerr << "Error in newton-raphson missing totcatch"
                 << endl << " Zstar is less than M(j))"
                 << endl << " You need to increase rmax " << endl
                 << " at present Zstar = " << Zstar  << " and M(j) = "
                 << M(j) << endl;
            ad_exit(1);
          }
          lambda=(Zstar-M(j))/(Z(j)-M(j));
          if (value(lambda)<0.0)
          {
            cout << " lambda = " << lambda << " "
                 << " Zstar = " << Zstar << " "
                 << " M(j) = " << M(j) << endl;
          }
        }
        else
        {
          lambda=Zstar/Z(j);
          if (value(lambda)<0.0)
          {
            cout << " Zstar = " << Zstar << " "
                 << " Z(j) = " << Z(j) << endl;
          }
        }


        for (fi=1;fi<=nf;fi++)
        {
          F(fi,j)*=lambda;
        }
        if (Zstar > 0.99999)
        {
          cout << "CC " << Zstar << endl;
        }
        Z(j)=Zstar;
      }
    }
    if (pfsh->age_flags(92)!=3)
    {
      if (pfsh->age_flags(92)==4) // for Pope need to add M to get true Z
      {
        Z+=M;
      }
      S=exp(-Z);
    }
    else
    {
      Z=M-log(1.0-Z);
      S=exp(-Z);
    }
    F=log(1.e-20+F);
    return fpen;
  }
  
 void vnrf::get_initial_q(void)
 {
   int fi;
   int j;
    // get initial q values
    for (fi=1;fi<=nf;fi++)
    {
      dvariable ssum=0.0;
      if (!wtflag(fi))
      {
        //sum=s(fi)*N;
        ssum=sum(sN(fi));
        q(fi)=Cobs(fi)/ssum;
        if (q(fi)<0.0)
        {
          cout << "0" << endl;
        }
      }
      else
      {
        //sum=s(fi)*elem_prod(w,N);
        ssum=sum(swN(fi));
        q(fi)=Cobs(fi)/ssum;
        if (q(fi)<0.0)
        {
          cout << "0" << endl;
        }
      }
    }
    logq=log(1.e-15+q);
  }

 void vnrf::get_initial_q_ss2(void)
 {
   int fi;
    // get initial q values
    dvar_vector N1=elem_prod(N,exp(-0.5*M));
    dvar_vector wN1=elem_prod(w,N1);
    for (fi=1;fi<=nf;fi++)
    {
      dvariable ssum=0.0;
      if (!wtflag(fi))
      {
        ssum=s(fi)*N1;
        q(fi)=Cobs(fi)/ssum;
        if (q(fi)<0.0)
        {
          cout << "0" << endl;
        }
      }
      else
      {
        ssum=s(fi)*elem_prod(w,N1);
        q(fi)=Cobs(fi)/ssum;
        if (q(fi)<0.0)
        {
          cout << "0" << endl;
        }
      }
    }
    logq=log(1.e-15+q);
  }

  int getopt(int argc,char * argv[],const char * s)
  {
    int i;
    int sflag=0;
    for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],s)==0)
      {
        sflag=1;
        break;
      }
    }
    return sflag;
  }
  
  dvar_vector vnrf::calculate_ld()
  {
    q=exp(logq);
    return calculate_d();
  }
  
  dvar_vector vnrf::calculate_d()
  {
   calculate_F(); 
   calculate_Z(); 
   calculate_S(); 
   //calculate_T(); 
   //calculate_That(); 
   calculate_Chat(); 
   //dvar_vector d=Chat-phi*Cobs;
   dvar_vector d=Chat-Cobs;
   return d;
  }
 // 
  dvar_matrix vnrf::calculate_J()
  {
  
    int fi;
    MY_DOUBLE_TYPE diff1=1.e-4;
    MY_DOUBLE_TYPE diff2=3.e-4;
    dvar_matrix J(1,nf,1,nf);
    //cout << " enter diff " << endl;
    //cin >> diff; 
    for (fi=1;fi<=nf;fi++)
    {
      dvariable sd1=(diff1-q(fi))+q(fi);
      dvariable tmp1=q(fi);
      q(fi)+=sd1;
      dvar_vector fp1=calculate_d();
      q(fi)=tmp1;
      q(fi)-=sd1;
      dvar_vector fm1=calculate_d();
      q(fi)=tmp1;
      dvar_vector t1=(fp1-fm1)/(2.0*sd1);
      //J(fi)=(fp1-fm1)/(2.0*sd1);

     
      dvariable sd2=(diff2-q(fi))+q(fi);
      dvariable tmp2=q(fi);
      q(fi)+=sd2;
      dvar_vector fp2=calculate_d();
      q(fi)=tmp2;
      q(fi)-=sd2;
      dvar_vector fm2=calculate_d();
      q(fi)=tmp2;
      dvar_vector t2=(fp2-fm2)/(2.0*sd2);
      dvar_vector tmp=calculate_d();

      J(fi)=(9.0*t1-t2)/8.0;
      
    }
    return trans(J);
  }

  d3_array vnrf::calculate_alpha()
  {
    int fi,j;
    MY_DOUBLE_TYPE diff=1.e-6;
    d3_array u_alpha(1,nf,1,ng,1,nf);
    dvar_matrix J(1,nf,1,nf);
    //cout << " enter diff " << endl;
    //cin >> diff; 
    for (fi=1;fi<=nf;fi++)
    {
      for (j=1;j<=ng;j++)
      {
        dvariable sd=(diff-s(fi,j))+s(fi,j);
        dvariable tmp=s(fi,j);
        s(fi,j)+=sd;
        dvar_vector fp=calculate_d();
        //cout << " fp = " << fp << endl;
        s(fi,j)=tmp;
        s(fi,j)-=sd;
        dvar_vector fm=calculate_d();
        s(fi,j)=tmp;
        u_alpha(fi,j)=value((fp-fm))/(2.0*value(sd));
      }
    }
    return u_alpha;
  }

  dvar_matrix vnrf::calculate_lJ()
  {
  
    int fi;
    dvar_matrix J(1,nf,1,nf);
    MY_DOUBLE_TYPE diff=1.e-5;
    for (fi=1;fi<=nf;fi++)
    {
      dvariable sd=(diff-logq(fi))+logq(fi);
      dvariable tmp=logq(fi);
      logq(fi)+=sd;
      dvar_vector fp=calculate_ld();
      logq(fi)=tmp;
      logq(fi)-=sd;
      dvar_vector fm=calculate_ld();
      logq(fi)=tmp;
      J(fi)=(fp-fm)/(2.0*sd);
    }
    return J;
  }
 // 
 // 
 //     
  void vnrf::calculate_F(void)
  {
    int fi;
    int j;
    int mmin=F(1).indexmin();
    int mmax=F(1).indexmax();

    for (fi=1;fi<=nf;fi++)
    {
      F(fi)=q(fi)*s(fi);
      for (j=mmin;j<=mmax;j++)
      {
        if (value(F(fi,j))<0.0) 
        {
          cout << "neg F " << value(F(fi,j)) << endl;
          ad_exit(1);
        }
      }
    }
  }
  void vnrfm::calculate_Z(void)
  {
    vnrf::calculate_Z();
    get_Z()+=missing_z;
  }

  void vnrf::calculate_Z(void)
  {
    int fi;
    int j;
    Z.initialize();
    int mmin=Z.indexmin();
    for (fi=1;fi<=nf;fi++)
    {
      for (j=mmin;j<=ng;j++)
      {
        Z(j)+=F(fi,j);
      }
    }
    if (af92!=3 && af92!=4)
    {
      Z+=M;
    }
  }


  void vnrf::calculate_S()
  {
    if (af92 !=3)
    {
      for (int j=1;j<=ng;j++)
      { 
        if (value(Z(j))<=rmax)
        {
            S(j)=exp(-Z(j));
        }
        else
        {
          dvariable dd=Z(j)-rmax;
          S(j)=exp(-rmax)-exp(-rmax)*dd;
          dvariable ppen=0.0;
          //Z(j)=rmax-posfun(0.2-dd,0.1,ppen)+0.2;
        }
      }
    }
    else
    {
      for (int j=1;j<=ng;j++)
      { 
        if (Z(j) < pamin1)
        {
          S(j)=1-Z(j);
        }
        else if (Z(j) < pamin2)
        {
          S(j)=1-Z(j);
          fpen+=pwght1*square(Z(j)-pamin1);
        }
        else
        {
          fpen+=pwght1*square(Z(j)-pamin1);
          dvariable ptmp=0.0;
          dvariable tmp=pamin2-posfun(pamin2-Z(j),0.02,ptmp);
          fpen1+=pwght2*ptmp;
          omega(j)=tmp/Z(j);
          Z(j)=tmp;
          for (int fi=1;fi<=nf;fi++)
          {
            F(fi,j)*=omega(j);
          }
        }
      } 
    }
  }

  void vnrf::calculate_Chat(void)
  {
    int fi,j;
    Chat.initialize();
    int mmin=Z.indexmin();
    for (j=mmin;j<=ng;j++)
    {
      dvariable tmp=1.0/Z(j)*(1.0-S(j))*N(j);
      if (pfsh->age_flags(92)==4)  // Popes aprox
      {
        tmp*=exp(-0.5*M(j));
      }
      for (fi=1;fi<=nf;fi++)
      {
        C(fi,j)=F(fi,j)*tmp;
      }
    }
    for (fi=1;fi<=nf;fi++)
    {
      if (!wtflag(fi))
      {
        Chat(fi)=sum(C(fi));
      }
      else
      {
        Chat(fi)=w*C(fi);
      }
    }
  }

        
  dvar_matrix vnrf::calculate_J_vanal(void)
  {
    int i,j,k;
    dvar_matrix Jac(1,nf,1,nf);
    Jac.initialize();
    for (i=1;i<=nf;i++)
    {
      for (k=1;k<=nf;k++)
      {
        for (j=1;j<=ng;j++)
        {
          if (i==k)  
          {
            if (Z(j)<=rmax)
            {
              dvariable tmp=((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*S(j))*s(k,j)/Z(j)*N(j);
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
            else
            {
              dvariable tmp=((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*exp(-rmax))*s(k,j)/Z(j)*N(j);
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
          }
          else
          {
            if (Z(j)<=rmax)
            {
              dvariable tmp=( -(1.0-S(j))/Z(j)+S(j))*s(k,j)*F(i,j)/Z(j)*N(j);
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
            else
            {
              dvariable tmp=
                ( -(1.0-S(j))/Z(j)+exp(-rmax))*s(k,j)*F(i,j)/Z(j)*N(j);
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
          }
        }
      }
    }
    return Jac;
  }
  dmatrix vnrf::calculate_J_anal(void)
  {
    int i,j,k;
    dvar_matrix Jac(1,nf,1,nf);
    Jac.initialize();
    for (i=1;i<=nf;i++)
    {
      for (k=1;k<=nf;k++)
      {
        for (j=1;j<=ng;j++)
        {
          if (i==k)  
          {
            if (Z(j)<=rmax)
            {
              Jac(i,k)+=
                ((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*S(j))*s(k,j)/Z(j)*N(j);
            }
            else
            {
              Jac(i,k)+=
                ((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*exp(-rmax))*s(k,j)/Z(j)*N(j);
            }
          }
          else
          {
            if (Z(j)<=rmax)
            {
              Jac(i,k)+=
                ( -(1.0-S(j))/Z(j)+S(j))*s(k,j)*F(i,j)/Z(j)*N(j);
            }
            else
            {
              Jac(i,k)+=
                ( -(1.0-S(j))/Z(j)+exp(-rmax))*s(k,j)*F(i,j)/Z(j)*N(j);
            }
          }
        }
      }
    }
    return value(Jac);
  }
        
  d3_array vnrf::calculate_alpha_ders(void)
  {
    int i,j,k;
    dvar3_array u_alpha(1,nf,1,ng,1,nf);
    for (i=1;i<=nf;i++)
    {
      for (k=1;k<=nf;k++)
      {
        for (j=1;j<=ng;j++)
        {
          if (i==k)  
          {
            u_alpha(k,j,i)=
              ((1.0-S(j))*(1.0-F(i,j)/Z(j))+F(i,j)*S(j))*q(k)/Z(j)*N(j);
          }
          else
          {
            u_alpha(k,j,i)=
              ( -(1.0-S(j))/Z(j)+S(j))*q(k)*F(i,j)/Z(j)*N(j);
              //(-S(j)*(1.0-F(i,j)/Z(j))+F(i,j)*S(j))*q(k)/Z(j)*N(j);
          }
        }
      }
    }
    return value(u_alpha);
  }
        
  dvector vnrf::calculate_N_ders(void)
  {
    int i,j,k;
    dvar_vector N_alpha(1,ng);
    N_alpha.initialize();
    for (i=1;i<=nf;i++)
    {
      for (j=1;j<=ng;j++)
      {
        N_alpha(j)+=F(i,j)/Z(j)*(1.0-S(j));
      }
    }
    return value(N_alpha);
  }

   static void df_ddcatchimp(void)
   {
     verify_identifier_string("E10");
     dvar_vector_position vpos=restore_dvar_vector_position();
     verify_identifier_string("E11");
 
     //d3_array_position dFdqpos=restore_d3_array_position();
     //d3_array dFdq= restore_d3_array_value(dFdqpos);
     //verify_identifier_string("E12");
 
     d3_array_position dsdqpos=restore_d3_array_position();
     d3_array dsdq= restore_d3_array_value(dsdqpos);
     verify_identifier_string("E13");
 
     dmatrix_position dMdqpos=restore_dmatrix_position();
     dmatrix dMdq= restore_dmatrix_value(dMdqpos);
     verify_identifier_string("E14");
 
     dmatrix_position dNdqpos=restore_dmatrix_position();
     dmatrix dNdq= restore_dmatrix_value(dNdqpos);
     verify_identifier_string("E15");
 
     dvar_matrix_position selpos=restore_dvar_matrix_position();
     verify_identifier_string("E16");
 
     dvar_vector_position Npos=restore_dvar_vector_position();
     verify_identifier_string("E17");
 
     dvar_vector_position Mpos=restore_dvar_vector_position();
     verify_identifier_string("E18");
 
     dvector dfu=restore_dvar_vector_derivatives(vpos);
 
     //dmatrix dfF=dFdq*dfu;
     dmatrix dfs=dsdq*dfu;
     dvector dfN=dNdq*dfu;
     dvector dfM=dMdq*dfu;
 
     if (norm2(dfs)> 1.e+30 || norm2(dfN) > 1.e+30 || norm2(dfM) > 1.e+30)
     {
       cout << "A1 possible overflow" << endl;
     } 
     dfs.save_dmatrix_derivatives(selpos);
     dfN.save_dvector_derivatives(Npos);
     dfM.save_dvector_derivatives(Mpos);
   }
 // 
 //   dvar_vector dvar_fish_stock_history::ddnrf_interface(int ir,int ip,
 //     ivector& wtflag)
 //   {
 //     int fi;
 //     MY_DOUBLE_TYPE beta=0.6;
 //     dvar_vector enf=mfexp(num_fish(ir,ip));
 //     int nfi=num_fish_incidents(ir,ip);
 // 
 //     // these are dependent variables that we don't need ?
 //     dvar_matrix sel(1,nfi,1,nage);
 //     for (fi=1;fi<=nfi;fi++)
 //     {
 //       int i=parent(ir,ip,fi);
 //       int rr=realization_region(i,1);
 //       int rp=realization_period(i,1);
 //       int ri=realization_incident(i,1);
 //       // we don't have to do this every time !!!
 //       sel(fi)=mfexp(incident_sel(rr,rp,ri));
 //     }
 //     MY_DOUBLE_TYPE rmax=1.7;
 //     if (age_flags(92)==2) rmax*=1.4;
 // 
 //     ddnrf ddnrtester(nfi,nage,value(sel),value(e_nat_mort_by_period(ir,ip)),
 //       value(enf),obs_tot_catch(ir,ip),value(mean_weight_yr(ir,ip)),
 //       wtflag,beta,this,ir,ip,rmax); 
 // 
 //     ddnrtester.testnr();
 // 
 //     if (0)
 //     {
 //        ofstream ofs("logfile1");
 //    
 //        ofs << "ddnrtester.S" << endl;
 //        ofs << ddnrtester.get_S() << endl;
 //        ofs << "ddnrtester.q" << endl;
 //        ofs << ddnrtester.get_q() << endl;
 //     }
 //    
 // 
 //     dvector cv(1,nfi);
 //     for (fi=1;fi<=nfi;fi++)
 //     {
 // #if ((defined(__MSVC32__) &&  __MSVC32__ >=8) || defined(linux) || defined(__ADMING__))
 //       cv(fi)=to_double(ddnrtester.get_q()(fi));
 // #else
 //       cv(fi)=ddnrtester.get_q()(fi);
 // #endif
 //     }
 //     dvar_vector v=nograd_assign(cv);
 // 
 // 
 //     if (ddnrtester.icount>8)
 //     {
 //       cout << "D icount = " << ddnrtester.icount << endl;
 //     }
 //     //d3_array & dFdq=make_d3_array(ddnrtester.get_dFdq());
 //     const d3_array & dsdq=make_d3_array(ddnrtester.get_dsdq());
 //     const dmatrix & dNdq=make_dmatrix(ddnrtester.get_dNdq());
 //     const dmatrix & dMdq=make_dmatrix(ddnrtester.get_dMdq());
 // 
 //     save_identifier_string("E18");
 //     e_nat_mort_by_period(ir,ip).save_dvar_vector_position();
 //     
 //     save_identifier_string("E17");
 //     enf.save_dvar_vector_position();
 //     
 //     save_identifier_string("E16");
 //     sel.save_dvar_matrix_position();
 //     
 //     save_identifier_string("E15");
 //     dNdq.save_dmatrix_value();
 //     dNdq.save_dmatrix_position();
 // 
 //     save_identifier_string("E14");
 //     dMdq.save_dmatrix_value();
 //     dMdq.save_dmatrix_position();
 // 
 //     save_identifier_string("E13");
 //     dsdq.save_d3_array_value();
 //     dsdq.save_d3_array_position();
 // 
 //     //save_identifier_string("E12");
 //     //dFdq.save_d3_array_value();
 //     //dFdq.save_d3_array_position();
 // 
 //     save_identifier_string("E11");
 //     v.save_dvar_vector_position();
 // 
 //     save_identifier_string("E10");
 //     gradient_structure::GRAD_STACK1->set_gradient_stack(df_ddcatchimp);
 //     return v;
 //   } 
