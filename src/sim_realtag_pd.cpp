/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
//# define COUNT_FISH 
extern mf_pvm_manager * mf_pvm;
MY_DOUBLE_TYPE multiplier=1.e+12;
void myxxx(int &per){;}

dvariable max(const dvar_vector & v,int & nn);
#ifdef __MSVC32__
  dvariable age_at_length_calcxx(MY_DOUBLE_TYPE& v,dvar_vector& gml,int nslots);
#else
  dvariable age_at_length_calcxx(const dvariable& v,dvar_vector& gml,int nslots);
#endif
      
int do_only_projection_periods_flag=0;
int compare(dmatrix data,int i);
dmatrix subsort(int level,ivector c,dmatrix M);

int see_if_caught(dvector & cprob,MY_DOUBLE_TYPE sr,random_number_generator & rng)
{
  
  MY_DOUBLE_TYPE psc=sum(cprob)/sr;
  MY_DOUBLE_TYPE psnc=1.0-sum(cprob)/sr;
  MY_DOUBLE_TYPE e=randu(rng);
  if (psc>1.0 || psc < 0.0)
  {
    cerr << "psc too large or too small " << psc  << endl;
    ad_exit(1);
  }
  if (e<psc)
  {
    //cout << e << "<" << psc <<" so caught in simulation period" << endl;
    MY_DOUBLE_TYPE e1=randu(rng);
    dvector v=cprob/sr;
    v/=sum(v);
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    MY_DOUBLE_TYPE sum=0;
    int i;
    for (i=mmin;i<=mmax;i++)
    {
      sum+=v(i);
      if (e1<sum)
      {
        //cout << "caught in place " << i << endl;
        break;
      }
    }
    return i;
  }
  else
  {
    //cout << e << ">" << psnc <<" so not caught during simulation period" << endl;
    return 0;
  }
}

MY_DOUBLE_TYPE dvar_len_fish_stock_history::get_simulated_real_tag_fish_age(int il,
 random_number_generator& rng)
{  
  dvector itra(1,nage);
  
  itra.initialize();
  
#if !defined(NO_MY_DOUBLE_TYPE)
  MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
#else
  MY_DOUBLE_TYPE len=tag_shlen+(il-0.5)*tag_filen;
#endif
  dvariable age;
  if (!parest_flags(175) && !parest_flags(174))
  {
    age=age_at_length_calc(len,vb_coff,nage,parest_flags);
  }
  else
  {
    age=age_at_length_calcxx(len,gml,nage+1);
  }
  MY_DOUBLE_TYPE cage=value(age);
  if (cage<=1)
  {
    itra(1)+=1;
  }
  else if (cage>=nage)
  {
    itra(nage)+=1;
  }
  else
  {
    MY_DOUBLE_TYPE sf;
    sf=value(daves_kludge1(age));
    int jj=int(cage);
    MY_DOUBLE_TYPE tp= sf*1;

    itra(jj)+= 1;
    itra(jj)-=tp;

    itra(jj+1)+=tp;
  } 

  // use itra to get the one or two ages
  int ii=0;
  dvector p(1,2);
  ivector ag(1,2);
  p.initialize();
  ag.initialize();
  for (int i=1;i<=nage;i++)
  {
    if (itra(i)>1.e-12)
    {
      p(++ii)=itra(i);
      ag(ii)=i;
    }
  }
  int iage=0;
#if !defined(NO_MY_DOUBLE_TYPE)
  if (abs(sum(p)-1.0L)> 1.e-6 || ii==0 || ii>2)
#else
  if (abs(sum(p)-1.0)> 1.e-6 || ii==0 || ii>2)
#endif
  {
    cout << "Simulated tag fish age proportion is out of bounds" << endl;
  }
  if (ii==1) 
  {
    iage=ag(1);
  }
  else
  {
    MY_DOUBLE_TYPE e=randu(rng);
    if (e<p(1))     // age is age(1)
    {
      iage=ag(1);
    }
    else
    {
      iage=ag(2);
    }
  }
  return iage;
} 
ofstream * pofsvvv = 0;

void dvar_len_fish_stock_history::
  check_if_tag_survived_and_was_caught(int it,int iage,
  random_number_generator& rng,int il,const dmatrix & _data1,int & ii)
{
  ADUNCONST(dmatrix,data1)
  //cout << column(alltagsurv(iage),1) << endl;
  MY_DOUBLE_TYPE sr=sum(column(alltagsurv(iage),1));
  //cout << sr<< endl;
  MY_DOUBLE_TYPE e1=randu(rng);
  int ss=0;
  if (e1>sr)  
  {
    //cout << e1 << ">" << sr <<  " so tag died before projection period "
    //     << endl;
  }
  else
  {
    //cout << e1 << "<" << sr <<  " so tag survived to projection period "
    //     << endl;
    dvector cprob=column(allprob(iage),1);
    //cout << cprob << endl;
    MY_DOUBLE_TYPE psc=sum(cprob)/sr;
    if (psc>1.0 || psc < 0.0)
    {
      cerr << "psc too large or too small " << psc  << endl;
      ad_exit(1);
    }
    //cout << "prob of survivor being caught " << psc << endl;
    //cout << "prob of survivor not being caught " << 1.0-psc << endl;
    int ic=see_if_caught(cprob,sr,rng);
    if (ic) // fish was caught
    {
      int ir=int(allprob(iage,ic,2));
      int ip=int(allprob(iage,ic,3));
      int fi=int(allprob(iage,ic,4));
      int fishery=parent(ir,ip,fi);
      int ty=really_true_year(ir,ip)+year1-1;
      int tm=really_true_month(ir,ip);
      ii++;
      data1(ii,1)=it;
      data1(ii,2)=iage;
      data1(ii,3)=(tag_shlen-tag_filen)+(il*tag_filen);      //NMD_22May2018
      data1(ii,4)=fishery;
      data1(ii,5)=ty;
      data1(ii,6)=tm;
    }
  }
}
  
void dvar_len_fish_stock_history::
  make_fake_recapture(int it, dmatrix & _data1,int & ii)
{
  ADUNCONST(dmatrix,data1)
  ii++;
  data1(ii,1)=it;
  data1(ii,2)=0;
  data1(ii,3)=0;
  data1(ii,4)=0;
  data1(ii,5)=0;
  data1(ii,6)=0;
}
  
int dvar_len_fish_stock_history::check_if_tag_exists(MY_DOUBLE_TYPE d,
  random_number_generator &rng)
{
  MY_DOUBLE_TYPE e=randu(rng);
  if (e<d)    // this fish exists
  {
    //cout << "tag exists " << endl;
    return 1;
  }
  else
  {
    //cout << "tag did not exist" << endl;
    return 0;
  }
}
int get_random_seed(void)
{
  int iseed=1101;
  int iseed1=0;
  {
    ifstream ifs("simseed"); 
    {
      ifs >> iseed1;
    }
    if (!ifs)
    {
//      ifs.close();
      #if defined(close)
      #  undef close
            ifs.close();
      #  define close _close
      #else
            ifs.close();
      #endif

      ofstream ofs("simseed"); 
      {
        ofs << iseed+2;
      }
      return iseed;
    }
    else
    {
//      ifs.close();
      #if defined(close)
      #  undef close
            ifs.close();
      #  define close _close
      #else
            ifs.close();
      #endif
      
      ofstream ofs("simseed"); 
      {
        ofs << iseed1+2;
      }
      return iseed1;
    }
  }
}


void dvar_len_fish_stock_history::generate_real_simulated_tagsx(void)
{
  int iseed=get_random_seed();
  random_number_generator rng(iseed);
  pofsvvv = new ofstream("simrealtag.dat");
  dmatrix data1(1,150000,1,6);
  data1.initialize();
  int idata=0;
  for (int it=1;it<=num_tag_releases;it++)
  {
    (*pofsvvv) << "Tag release " << it << endl;

    make_fake_recapture(it,data1,idata);

    get_real_simtag_probs(sv,it);
  
    dvector& itrl=initial_tag_release_by_length(it);
    dvector itra(1,nage);
    int itags=0;
    int intnf=0;
    MY_DOUBLE_TYPE fracnf=0.0;
    for (int il=1;il<=tag_nlint;il++)
    {
#if !defined(NO_MY_DOUBLE_TYPE)
      for (MY_DOUBLE_TYPE d=itrl(il);d>0.0;d-=1.0L)
#else
      for (MY_DOUBLE_TYPE d=itrl(il);d>0.0;d-=1.0)
#endif
      {
        if (!check_if_tag_exists(d,rng)) break;
  
        int iage=static_cast<int>(get_simulated_real_tag_fish_age(il,rng));
        check_if_tag_survived_and_was_caught(it,iage,rng,il,data1,idata);
      }
    }
  }
  make_real_tag_data_report(data1,idata,"report.realtag_");

  delete pofsvvv;
  pofsvvv=0;
  //ad_exit(1);
}

void dvar_fish_stock_history::get_real_simtag_probs
  (dvar_vector& sv,int it)
 {
   //dmatrix prob(1,1000,1,5);
   ofstream xpofs("getrsp_movestuff");

   int break_flag=0;
   epooled_tagnum_fish_recr.initialize();
   pooledtagN=-20;
   if (allocated(allprob))
     allprob.deallocate();
   //allprob.allocate(1,nage,1,1000,1,5);
   allprob.allocate(1,nage);
   for (int iage=1;iage<=nage;iage++)
   {
     if (initial_tag_release_by_age(it,iage)>0)
     {
       allprob(iage).allocate(1,10000,1,5);
     }
   }
   allprob.initialize();
   if (allocated(alltagsurv))
     alltagsurv.deallocate();
   alltagsurv.allocate(1,nage,1,num_regions,1,2);
   alltagsurv.initialize();
   for (int iage=1;iage<=nage;iage++)
   {
     if (initial_tag_release_by_age(it,iage)>0)
     {
       dmatrix prob=allprob(iage);
       dmatrix tagsurv=alltagsurv(iage);
       int current_year=0;
       ivector rip(1,num_regions);
       current_year=tag_year(it);
       rip=initial_tag_period(it);
       int irr=tag_region(it);
       int finished_flag=0;
       ivector tmp_mp(1,num_regions);
       ivector tmp_yr(1,num_regions);
       ivector tmp_mn(1,num_regions);
       ivector tmp_ip(1,num_regions);
       tmp_mp.initialize();
       tmp_yr.initialize();
       tmp_mn.initialize();
       tmp_ip.initialize();
       int eflag=0;
       int ir;
       for (ir=1;ir<=num_regions;ir++)
       {
         tagnum_fish(it,ir)=-20.0;
       }
       tagnum_fish(it,irr,rip(irr),iage)=log(multiplier);
   
       do
       {
         finished_flag=1;
         eflag=0;
         for (int ir=1;ir<=num_regions;ir++)
         {
           break_flag=0;
           int& ip=rip(ir);
           int ttp;
           ttp=terminal_tag_period(it,ir);
           do
           {
             if ( ip>=num_fish_periods(ir) || ip>ttp ) 
             {
               eflag=1;
               break;
             }
             finished_flag=0;
             if (ip < ttp)
             {
               if (year(ir,ip+1)==year(ir,ip))
               { 
                 if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) 
                    || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
                 {
                     tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-tot_mort(ir,ip);
                 }
                 else
                 {
                   tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-nrtm(it,ir,ip);
                 }
                 tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-tot_mort(ir,ip);
               }
               else
               {
                 if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) 
                   || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
                 {
                   dvar_vector  tnf;
                   dvar_vector  tnfp;
                   tnf=tagnum_fish(it,ir,ip);
                   tnfp=tagnum_fish(it,ir,ip+1);
                   int jmin=tnf.indexmin();
                   if (nage>2 && jmin<nage-1)
                     --tnfp(jmin+1,nage-1)
                       =tnf(jmin,nage-2)
                       -tot_mort(ir,ip)(jmin,nage-2);
          
                   if (jmin<nage)
                   {
                     tnfp(nage) 
                       = log(1.e-20
                         +mfexp(tnf(nage)-tot_mort(ir,ip,nage))
                         + mfexp(tnf(nage-1)-tot_mort(ir,ip,nage-1)));
                   }
                   else
                   {
                     tnfp(nage) 
                       = tnf(nage)-tot_mort(ir,ip,nage);
                   }
                 }
                 else
                 {
                   int minage=nrtm(it,ir,ip).indexmin();
                   if (nage>minage+1)
                     --tagnum_fish(it,ir,ip+1)(minage+1,nage-1)
                      =tagnum_fish(it,ir,ip)(minage,nage-2)
                      -nrtm(it,ir,ip)(minage,nage-2);
          
                   tagnum_fish(it,ir,ip+1,nage) 
                     = log( 1.e-20
                        +(mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
                         )
                     + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1)));
                 }
               }
             }
             else
             {
               eflag=1;
               if (!tag_flags(it,1) || ip >= initial_tag_period(it,ir)
                 +tag_flags(it,1))
               {
                 int jmin;
                 jmin=tagnum_fish(it,ir,ip).indexmin();
                 
                 if (nage>2 && jmin<nage-1)
                   --epooled_tagnum_fish_recr(ir,ip+1)(jmin+1,nage-1)
                     +=mfexp(tagnum_fish(it,ir,ip)(jmin,nage-2)
                     -tot_mort(ir,ip)(jmin,nage-2));
       
                 if (nage>2 && jmin<nage)
                   epooled_tagnum_fish_recr(ir,ip+1,nage) 
                     += mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage))
                     +mfexp(tagnum_fish(it,ir,ip,nage-1)-tot_mort(ir,ip,nage-1));
                 else
                   epooled_tagnum_fish_recr(ir,ip+1,nage) 
                     += mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage));
               }
               else
               {
                 if (nage>2)
                   --epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1)
                    +=mfexp(tagnum_fish(it,ir,ip)(1,nage-2)
                    -nrtm(it,ir,ip)(1,nage-2));
      
                 epooled_tagnum_fish_recr(ir,ip+1,nage) 
                   += mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
                   + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1));
               }
             }
             if (move_flags(ir,ip+1))
             {
               tmp_ip(ir)=ip+1;
               tmp_yr(ir)=year(ir,ip+1);
               tmp_mn(ir)=month(ir,ip+1);
               tmp_mp(ir)=move_index(ir,ip+1);
          
               if (ir==num_regions)
               {
                 print_movement_stuff(xpofs,tmp_mp,tmp_yr,tmp_mn,tmp_ip,
                   num_regions);
                 xpofs << endl;

                 int diff_flag=0;
                 for (int i=2;i<=num_regions;i++)
                 {
                   if (tmp_mp(i) !=tmp_mp(1))
                     diff_flag=1;
                 }

                 if (diff_flag)
                 {
                   cout << "sanity error" << endl;
                   cout << "map = " << tmp_mp(1,ir) << endl;
                   cout << "year = " << tmp_yr(1,ir) << endl;
                   cout << "month = " << tmp_mn(1,ir) << endl;
                   ad_exit(1);
                 }
               }
               break_flag=1;
             } 
             ip++;
           }
           while (!break_flag); 
         }
         // move the tags  -- need to deal with pooled ones
         {
           if (num_regions>1 && !eflag)
           {
             check_sanity(tmp_mp);
             dvar_matrix tmp=
               fast_diffusion_calcs(nage,num_regions,tagnum_fish(it),
               Dad(tmp_mp(1)),rip);
             for (int ir=1;ir<=num_regions;ir++)
             {
               tagnum_fish(it,ir,rip(ir))=log(5.e-10+tmp(ir));
             }
           }
         }
       } // need to decide when to quit
       while (!finished_flag);
   
       const MY_DOUBLE_TYPE one_plus=1.e0+1.e-12;
     
       int id=0;
       imatrix * pitp;
       imatrix * pttp;
       pitp=&initial_tag_period;
       pttp=&terminal_tag_period;
   
       id=0;
       for (int ir=1;ir<=num_regions;ir++)
       {
         int per=(*pttp)(it,ir);
         int per1=num_real_fish_periods(ir);
         if (per1>=per)
         {
           cerr << "hosed" << endl;
           ad_exit(1);
         }
         myxxx(per);
         myxxx(per1);
         int jmin=tagnum_fish(it,ir,(*pttp)(it,ir)).indexmin();
         if (jmin<nage)
         {
           cerr << "Have not dealt with this yet" << endl;
           ad_exit(1);
         }
  
         /*
         int imin=tagnum_fish(it,ir,per1).indexmin();
         int imax=tagnum_fish(it,ir,per1).indexmax();
         cout << ir << "  " 
              << tagnum_fish(it,ir,per).indexmin() << "  "
              << tagnum_fish(it,ir,per).indexmax() << "  "
              << max(exp(tagnum_fish(it,ir,per)))/multiplier << "  "
              << max(exp(tagnum_fish(it,ir,per1+1)))/multiplier << "  "
              << max(exp(tagnum_fish(it,ir,per1)-tot_mort(ir,per1)(imin,imax)))/
                   multiplier << endl;
         */
       }
       for (int ir=1;ir<=num_regions;ir++)
       {
         int nn=0;
         int per1=initial_tag_period(it,ir);
         if (do_only_projection_periods_flag)
         {
           per1=num_real_fish_periods(ir)+1;
         }
         tagsurv(ir,1)=exp(value(max(tagnum_fish(it,ir,per1),nn)))/
            multiplier;
         tagsurv(ir,2)=nn;
         for (int ip=per1;ip<=(*pttp)(it,ir);ip++)
         {
           if (!tag_flags(it,1) || 
             ip >= initial_tag_period(it,ir)+tag_flags(it,1))
           {
             for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
             {        
               int jmin=tagcatch(it,ir,ip,fi).indexmin();
               dvariable tc=
                exp(max(fish_mort_calcs(ir,ip,fi)(jmin,nage)
                  +tagnum_fish(it,ir,ip),nn));
               prob(++id,1)=value(tc)/multiplier;
               // correct for reporting rate
               MY_DOUBLE_TYPE rr=get_simulation_reporting_rate(it,ir,ip,fi,1);
               prob(id,1)*=rr;
               prob(id,2)=ir;
               prob(id,3)=ip;
               prob(id,4)=fi;
               prob(id,5)=nn;
             }
           }
           else
           {
             for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
             {        
               int jmin=tagcatch(it,ir,ip,fi).indexmin();
               int nn=0;
               {
                 nrfm(it,ir,ip,fi);
                 nrtm(it,ir,ip);
                 nrsurv(it,ir,ip)(jmin,nage);
                 tagnum_fish(it,ir,ip);
               }
               dvariable tc=exp(max(nrfm(it,ir,ip,fi)-log(1.e-10+nrtm(it,ir,ip))
                 +log(one_plus-nrsurv(it,ir,ip)(jmin,nage))+
                 tagnum_fish(it,ir,ip),nn));
               prob(++id,1)=value(tc)/multiplier;
               // correct for reporting rate
               MY_DOUBLE_TYPE rr=get_simulation_reporting_rate(it,ir,ip,fi,1);
               prob(id,1)*=rr;
               prob(id,2)=ir;
               prob(id,3)=ip;
               prob(id,4)=fi;
               prob(id,5)=nn;
             }
           }
         }
       }
       //cout <<  endl;
     }
   }
   //cout << "quit for test" << endl;
   //ad_exit(1);
 }
void dvar_len_fish_stock_history::make_real_tag_data_report(dmatrix& data1,
  int ii,const char * s)
//  int ii,char s[])
{
  int psi=projection_sim_index;
//  ofstream ofs3(adstring("tmp_")+str(psi));
//  ofs3 << setw(6) << data1.sub(1,ii) << endl;
  dvector tag_recoveries(1,num_tag_releases);
  tag_recoveries.initialize();
  int nsort=5;
  ivector sort_columns(1,nsort);
  sort_columns(1)=1;
  sort_columns(2)=5;
  sort_columns(3)=6;
  sort_columns(4)=4;
  sort_columns(5)=3;
  dmatrix sdata=sort(data1.sub(1,ii),sort_columns(1));
    
  for (int i=2;i<=nsort;i++)
  {
    sdata=subsort(2,sort_columns,sdata);
  }

  ivector nrecaps(1,num_tag_releases);
  int maxlen=int(tag_shlen+tag_nlint*tag_filen);
  nrecaps.initialize();

  for (int it=1;it<=num_tag_releases;it++)
  {
    int number=1;
    for (int i=2;i<=ii;i++)
    {
      if (sdata(i-1,1) == it && sdata(i,1) == it)
      {
        if (!compare(sdata,i))
        {
//    - rows have different stratification, increment row numbers
          number++;
        }
        else
        {
//    - rows have same stratification, frequency is >1, do nothing
        }
      }
    }
    nrecaps(it)=number;
  }

  
  //ofstream ofs10(adstring("report.tag_")+str(psi));
  ofstream ofs10(adstring(s)+str(psi));
  
  ofs10 << "# RELEASE GROUPS    STARTING LENGTH    NUMBER INTERVALS    "
              "INTERVAL LENGTH" << endl;
  ofs10 << "      " << setw(4) << num_tag_releases << "              "  << tag_shlen 
        << "                  " << tag_nlint << "                   " << tag_filen << endl;
  
  ofs10 << "# TAG RECOVERIES" << endl;

  ofs10 << setw(5) << nrecaps-1 << endl;

  ofs10 << "# XX " << endl;
  int nr=1;
  int it= static_cast<int>(sdata(1,1));
  int need_print=1;
  int itold=it;
  int minage=1;
  while (need_print<it)
  {
    ofs10 << "#  " << need_print << "  RELEASE REGION    YEAR    MONTH   "
          " Tag_program Simulation " << endl;

    int relregion=tag_region(need_print);
    int trelmonth=true_tag_month(need_print);
    int trelyear=true_tag_year(need_print);
    int sitp=initial_tag_period(need_print,relregion);
    //int sti=tag_incident(need_print);
    int sti=1;
    dvector ml=value(mean_length(relregion,sitp,sti)(minage,nage));
    dvector sd=value(sdevs(relregion,sitp,sti)(minage,nage));
    ofs10 << "         " << relregion << "           " << trelyear 
          << "         " << trelmonth << endl;
    ofs10 <<  initial_tag_release_by_length(need_print)  << endl;
    need_print++;
  }

  need_print=it+1;
  ofs10 << "#  " << it << "  RELEASE REGION    YEAR    MONTH   "
          " Tag_program Simulation " << endl;
  int relregion=tag_region(it);
  int trelmonth=true_tag_month(it);
  int trelyear=true_tag_year(it);
  int sitp=initial_tag_period(it,relregion);
  //int sti=tag_incident(it);
  int sti=1;

  dvector ml=value(mean_length(relregion,sitp,sti)(minage,nage));
  dvector sd=value(sdevs(relregion,sitp,sti)(minage,nage));
  //dvector lendist(1,tag_nlint);
  /*
  for (int il=1;il<=tag_nlint;il++)
  {
#if !defined(NO_MY_DOUBLE_TYPE)
    MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
#else
    MY_DOUBLE_TYPE len=tag_shlen+(il-0.5)*tag_filen;
#endif
    lendist(il)=cinitial_tag_release_by_age(it)(minage,nage)
      *exp(-0.5*square(elem_div(len-ml,sd)));
  }
  */

  ofs10 << "         " << relregion << "           " << trelyear 
        << "         " << trelmonth << endl;
  ofs10 <<  initial_tag_release_by_length(it)  << endl;
  //ofs10 <<  num_tags_at_length(it)  << endl;
  ofs10 << "# LENGTH RELEASE    FISHERY    RECAP YEAR    RECAP MONTH    "
            "NUMBER" << endl;
           /*
            data1(ii,1)=it;
            data1(ii,2)=j;
            data1(ii,3)=int(length);
            data1(ii,4)=fishery;
            data1(ii,5)=ty;
            data1(ii,6)=tm;
           */

  int number=1;
  for (int i=2;i<=ii;i++)
  {
    int i1=i-1;
    if (!compare(sdata,i))
    {
      if (sdata(i1,3)>0)
      {
        ofs10 << "      " << setw(4) << sdata(i1,3) 
            <<  "           "  << setw(3) <<  sdata(i1,4)
            << "          " << setw(4) << sdata(i1,5)   
            << "          " << setw(4) << sdata(i1,6)   
            << "      " << setw(4) << number << endl;
      }
      number=1;
    }
    else
    {
      number++;
    }
    if (sdata(i,1)>itold)
    {
      it=static_cast<int>(sdata(i,1));
      while (need_print<it)
      {
        ofs10 << "#  " << need_print << "  RELEASE REGION    YEAR    MONTH   "
              " Tag_program Simulation " << endl;
    
        int relregion=tag_region(need_print);
        int trelmonth=true_tag_month(need_print);
        int trelyear=true_tag_year(need_print);
        int sitp=initial_tag_period(need_print,relregion);
        //int sti=tag_incident(need_print);
        int sti=1;
        dvector ml=value(mean_length(relregion,sitp,sti)(minage,nage));
        dvector sd=value(sdevs(relregion,sitp,sti)(minage,nage));
        ofs10 << "         " << relregion << "           " << trelyear 
              << "         " << trelmonth << endl;
        ofs10 <<  initial_tag_release_by_length(need_print)  << endl;
        need_print++;
      }

      relregion=tag_region(it);
      trelmonth=true_tag_month(it);
      trelyear=true_tag_year(it);
      //int sti=tag_incident(it);
      int sti=1;
      need_print=it+1;
      ofs10 << "#  " << it << "  RELEASE REGION    YEAR    MONTH   "
          " Tag_program Simulation " << endl;
      ofs10 << "            " << relregion << "           " << trelyear 
            << "      " << trelmonth << endl;

      dvector ml=value(mean_length(relregion,sitp,sti)(minage,nage));
      dvector sd=value(sdevs(relregion,sitp,sti)(minage,nage));
     /*
      dvector lendist(1,tag_nlint);
      for (int il=1;il<=tag_nlint;il++)
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
#else
        MY_DOUBLE_TYPE len=tag_shlen+(il-0.5)*tag_filen;
#endif
        lendist(il)=cinitial_tag_release_by_age(it)(minage,nage)
          *exp(-0.5*square(elem_div(len-ml,sd)));
      }
      lendist/=sum(lendist);
      */
     // ofs10 <<  num_tags_at_length(it)  << endl;
      ofs10 <<  initial_tag_release_by_length(it)  << endl;

      ofs10 << "# LENGTH RELEASE    FISHERY    RECAP YEAR    RECAP MONTH    "
            "NUMBER" << endl;
      itold=it;
    }
  }
  int i1=ii;
  if (sdata(i1,3)>0)
  {
    ofs10 << "      " << setw(4) << sdata(i1,3) 
        <<  "           "  << setw(3) <<  sdata(i1,4)
        << "          " << setw(4) << sdata(i1,5)   
        << "          " << setw(4) << sdata(i1,6)   
        << "      " << setw(4) << number << endl;
  }
  cout << "finished" << endl;
}
