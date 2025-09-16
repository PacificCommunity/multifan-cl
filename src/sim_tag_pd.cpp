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
dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int maxage,MY_DOUBLE_TYPE delta);

dvariable dvar_len_fish_stock_history::get_simtag_pd(int it, int iage, int irr,
  dmatrix & prob,int & id,const char * s)
{
  MY_DOUBLE_TYPE multiplier=1000000;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  dvariable ffpen=0.0;
  epooled_tagnum_fish_recr.initialize();
  sim_epooled_tagnum_fish_recr.initialize();
  pooledtagN=-20;
  ofstream xpofs("getsp_movestuff");

  if (simulation_seeds(5) && projection_sim_index > 0)
  {
    sim_get_initial_tag_population();
  }
  int xsim_tag_flag=1;
  int sim_tag_flag=1;
  sim_tagnum_fish.initialize();
  {
    int sit=0;

    int current_year=0;
    ivector rip(1,num_regions);
    current_year=sim_tag_year(it);
    rip=sim_initial_tag_period(it);
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
      sim_tagnum_fish(it,ir)=-20.0;
    }
    sim_tagnum_fish(it,irr,rip(irr),iage)=log(multiplier);
    do
    {
      finished_flag=1;
      eflag=0;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        int ttp;
        ttp=num_fish_periods(ir);
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
              if (!num_fish_incidents(ir,ip) 
                 || !tag_flags(it,1) 
                 || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
              {
                sim_tagnum_fish(it,ir,ip+1)=
                  sim_tagnum_fish(it,ir,ip)-tot_mort(ir,ip);
              }
              sim_tagnum_fish(it,ir,ip+1)=sim_tagnum_fish(it,ir,ip)
                -tot_mort(ir,ip);
            }
            else
            {
              if (!num_fish_incidents(ir,ip) 
                || !tag_flags(it,1) 
                || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
              {
                dvar_vector  tnf;
                dvar_vector  tnfp;
                tnf=sim_tagnum_fish(it,ir,ip);
                tnfp=sim_tagnum_fish(it,ir,ip+1);
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
            }
          }
          else
          {
            eflag=1;
            if (!tag_flags(it,1) 
              || ip >= initial_tag_period(it,ir)
              +tag_flags(it,1))
            {
              int jmin;
              jmin=sim_tagnum_fish(it,ir,ip).indexmin();
              
              if (nage>2 && jmin<nage-1)
                --sim_epooled_tagnum_fish_recr(ir,ip+1)(jmin+1,nage-1)
                  +=mfexp(sim_tagnum_fish(it,ir,ip)(jmin,nage-2)
                  -tot_mort(ir,ip)(jmin,nage-2));
     
              if (nage>2 && jmin<nage)
                sim_epooled_tagnum_fish_recr(ir,ip+1,nage) 
                  += mfexp(sim_tagnum_fish(it,ir,ip,nage)
                   -tot_mort(ir,ip,nage))
                  +mfexp(sim_tagnum_fish(it,ir,ip,nage-1)
                      -tot_mort(ir,ip,nage-1));
              else
                sim_epooled_tagnum_fish_recr(ir,ip+1,nage) 
                  +=mfexp(sim_tagnum_fish(it,ir,ip,nage)
                   -tot_mort(ir,ip,nage));
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
            fast_diffusion_calcs(nage,num_regions,sim_tagnum_fish(it),
            Dad(tmp_mp(1)),rip,1.e-12);
          for (int ir=1;ir<=num_regions;ir++)
          {
            sim_tagnum_fish(it,ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
  }

  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-12;

  int sit=0;
  imatrix pitp;
  imatrix pttp;
  pitp=sim_initial_tag_period;
  int mmin=pitp.indexmin();
  int mmax=pitp.indexmax();
  pttp.allocate(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    pttp(i)=num_fish_periods;
  } 
  int rr=sim_tag_region(it);
  int savety=really_true_year(rr,sim_initial_tag_period(it,rr))+year1-1;
  int savetm=really_true_month(rr,sim_initial_tag_period(it,rr));

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=pitp(it,ir);ip<=pttp(it,ir);ip++)
    {
      if (!tag_flags(it,1) || 
        ip >= initial_tag_period(it,ir)+tag_flags(it,1))
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          dvar_vector& tc=sim_tagcatch(it,ir,ip,fi);
          int jmin=tc.indexmin();
          tc=
           exp(fish_mort_calcs(ir,ip,fi)(jmin,nage)+sim_tagnum_fish(it,ir,ip));
          MY_DOUBLE_TYPE pp=sum(value(tc))/multiplier;
          if (pp>1.e-6)
          {
            int nn=0;
            MY_DOUBLE_TYPE maxval=-10;
            for (int i=jmin;i<=nage;i++)
            {
              if (tc(i)>maxval)
              {
                maxval=value(tc(i));
                nn=i;
              }
            }
            if (id>=prob.indexmax())
            {
              cerr << "Need to increase dimension of prob in tag simulation"
                   << endl;
            }
            prob(++id,1)=value(tc(nn))/multiplier;
            // correct for reporting rate
            MY_DOUBLE_TYPE rr=get_simulation_reporting_rate(it,ir,ip,fi,0);
            prob(id,1)*=rr;
            prob(id,2)=ir;
            prob(id,3)=ip;
            prob(id,4)=fi;
            prob(id,5)=nn;
          }
        }
      }
    }
  }
/*
# if defined(COUNT_FISH)
  print_tag_accounting_info();
  exit(1);
#endif
*/
  return ffpen;
}


MY_DOUBLE_TYPE dvar_fish_stock_history::get_simulation_reporting_rate(int it,
  int ir,int ip,int fi,int real_flag)
{
  MY_DOUBLE_TYPE tmp=1.0;
  if (real_flag==1)
  {
    if (age_flags(198))
      tmp=value(tag_rep_rate(it,ir,ip,fi));
    else
      tmp=value(rep_rate(ir,ip,fi));
  }
  else
  {
    //tmp=0.5;
    int p=parent(ir,ip,fi);
    tmp=sim_tag_fish_rep(p,it); 
  }
  return tmp;
}

void dvar_len_fish_stock_history::generate_simulated_tags(void)
{
  if (!allocated(csim_initial_tag_release_by_age))
  {
    csim_initial_tag_release_by_age.allocate(1,sim_num_tag_releases);
  }
  int psi=projection_sim_index;
  int minage=1;
//  ofstream ofs(adstring("simulated_tag_data_")+str(psi));
  dmatrix data1(1,150000,1,6);
  int ii=0;
  data1.initialize();
  MY_DOUBLE_TYPE maxlen=tag_shlen+tag_nlint*tag_filen;
  sim_tagcatch.initialize();  //NMD_6dec2017
  int iseed=simulation_seeds(5);    //NMD_19dec2019
  random_number_generator rng(iseed);    //Df_19dec2019
  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    sim_pooled_tagcatch.initialize();
    int ir=sim_tag_region(it);
    int ip=sim_initial_tag_period(it,ir);
    int fi=sim_tag_incident(it);
    int yr1=sim_tag_year(it);

    dvector ecatch=exp(value(catch(ir,ip,fi)));
    csim_initial_tag_release_by_age(it)=
      ecatch/(sum(ecatch)+1.0e-15);

    make_fake_recapture(it,data1,ii);
    // generate the number of tags in each age class
//    random_number_generator rng(103);
 //   random_number_generator rng(iseed);    //NMD_19dec2019
    ivector agesample(1,sim_tag_numbers_released(it));
    agesample.fill_multinomial(rng,csim_initial_tag_release_by_age(it));
    ivector sim_num_tags_at_age(minage,nage);
    sim_num_tags_at_age.initialize();
    for (int i=1;i<=sim_tag_numbers_released(it);i++)
    {
      sim_num_tags_at_age(agesample(i))++;
    }

    sim_num_tags_at_length(it)=get_sim_num_tags_at_length(it,
      sim_num_tags_at_age,rng);

    sim_obs_tagcatch.initialize();
    d3_array prob1(minage,nage);
    dmatrix prob(1,1000,1,5);
    for (int j=minage;j<=nage;j++)
    {
      int sample_size=sim_num_tags_at_age(j);
      if (sample_size>0)
      {
        prob.initialize();
        int id=0;
        get_simtag_pd(it,j,ir,prob,id,"file_1");
        int id1=id;
        int id2=id;
//        ofstream ofscc("file_2b");
//        for (int i=1;i<=id;i++)
//        {
//          ofscc  << " " << setw(10) << setprecision(5) 
//                 << setscientific() << prob(i,1) 
//                 << " " << setw(4) << setprecision(5) 
//                 << setfixed() << ivector(prob(i))(2,5) 
//                 << endl;
//        }
      
        dvector p(1,id+1);
        p(2,id+1).shift(1)=column(prob,1)(1,id);
        p(1)=1.0-sum(p(2,id+1));
        prob1(j).allocate(1,id,1,5);
        prob1(j)=prob.sub(1,id);
        ivector sample(1,sample_size);
        sample.fill_multinomial(rng,p);
        dvector eps(1,sample_size);
        eps.fill_randn(rng);
        for (int i=1;i<=sample_size;i++)
        {
          if (sample(i)>1)
          {
            int i1=sample(i)-1;
            int trelmonth=sim_true_tag_month(it);
            int trelyear=sim_true_tag_year(it);
            int relregion=sim_tag_region(it);
            int sitp=sim_initial_tag_period(it,relregion);
            int relyear=sim_true_tag_year(it);
            int sti=sim_tag_incident(it);
            MY_DOUBLE_TYPE sd=value(sdevs(relregion,sitp,sti,j));
            MY_DOUBLE_TYPE ml=value(mean_length(relregion,sitp,sti,j));
//            MY_DOUBLE_TYPE length=min(max(shlen,ml+sd*eps(i)),maxlen);
            MY_DOUBLE_TYPE length=min(max(tag_shlen,ml+sd*eps(i)),maxlen);
//          NMD_8dec2022 - to use the tagging data specified shortest length
            int num_tag=static_cast<int>(prob(sample(i),1));
            int ir=static_cast<int>(prob(i1,2));
            int ip=static_cast<int>(prob(i1,3));
            int fi=static_cast<int>(prob(i1,4));
            int age=static_cast<int>(prob(i1,5));
            int fishery=parent(ir,ip,fi);
            int ty=really_true_year(ir,ip)+year1-1;
            int tm=really_true_month(ir,ip);
            ii++;
            data1(ii,1)=it;
            data1(ii,2)=j;
            data1(ii,3)=int(length);
            data1(ii,4)=fishery;
            data1(ii,5)=ty;
            data1(ii,6)=tm;
          }
        }
        for (int i=id1+1;i<=id2;i++)
        {
          int num_tag=static_cast<int>(prob(i,1));
          int ir=static_cast<int>(prob(i,2));
          int ip=static_cast<int>(prob(i,3));
          int fi=static_cast<int>(prob(i,4));
          int age=static_cast<int>(prob(i,5));
          int fishery=parent(ir,ip,fi);
          int ty=really_true_year(ir,ip)+year1-1;
          int tm=really_true_month(ir,ip);

        }
      }
    }
   /*
    ofstream ofsz(adstring("allprob_")+str(psi) + "_" +str(it));
    for (int j=minage;j<=nage;j++)
    {
      ofsz << "Age " << j << endl;
      dvector p=column(prob1(j),1);
      ofsz << "Probability  not returned "  << 1.0-sum(p) << endl;
      ofsz << prob1(j) << endl;
    }
   ad_exit(1);
   */
  }
  //ofstream ofs10(adstring("report.tag_")+str(psi));
  make_tag_data_report(data1,ii,"report.simtag_");
  //ofstream ofs3(adstring("tmp_")+str(psi));
  //ofs3 << data << endl;
  //ad_exit(1);
}

dmatrix subsort(int level,ivector c,dmatrix M)
{
  int maxlevel=c.indexmax();
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int c1=c(level-1);
  int c2=c(level);
  MY_DOUBLE_TYPE x=M(mmin,c1);
  int lb=mmin;
  for (int i=mmin+1;i<=mmax;i++)
  {
    if (x<M(i,c1))
    {
      M.sub(lb,i-1)=sort(M.sub(lb,i-1),c2);
      if (level<maxlevel)
      {
        M.sub(lb,i-1)=subsort(level+1,c,M.sub(lb,i-1));
      }
      lb=i;
      x=M(i,c1);
    }
  }
  if (lb<mmax)
  {
    M.sub(lb,mmax)=sort(M.sub(lb,mmax),c2);
    if (level<maxlevel)
    {
      M.sub(lb,mmax)=subsort(level+1,c,M.sub(lb,mmax));
    }
  }
  return M;
}
      
int compare(dmatrix data,int i)
{
 /*
  data1(ii,1)=it;
  data1(ii,2)=j;
  data1(ii,3)=int(length);
  data1(ii,4)=fishery;
  data1(ii,5)=ty;
  data1(ii,6)=tm;
  */
  if (data(i,1) !=data(i-1,1) ||
     data(i,3) !=data(i-1,3) ||
     data(i,4) !=data(i-1,4) ||
     data(i,5) !=data(i-1,5) ||
     data(i,6) !=data(i-1,6) )
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

void dvar_len_fish_stock_history::make_tag_data_report(dmatrix& data1,
  int ii,const char * s)
//  int ii,char s[])
{
  int psi=projection_sim_index;
//  ofstream ofs3(adstring("tmp_")+str(psi));
//  ofs3 << setw(6) << data1.sub(1,ii) << endl;
  dvector tag_recoveries(1,sim_num_tag_releases);
  tag_recoveries.initialize();
 /*
  data1(ii,1)=it;
  data1(ii,2)=j;
  data1(ii,3)=int(length);
  data1(ii,4)=fishery;
  data1(ii,5)=ty;
  data1(ii,6)=tm;
  */
  int nsort=5;
  ivector sort_columns(1,nsort);
  sort_columns(1)=1;
  sort_columns(2)=5;
  sort_columns(3)=6;
  sort_columns(4)=4;
  sort_columns(5)=3;
  dmatrix sdata=sort(data1.sub(1,ii),sort_columns(1));
//  ofstream ofs5(adstring("stmp2_")+str(psi));
//  ofs5 << setw(6) << sdata << endl;
    
  for (int i=2;i<=nsort;i++)
  {
    sdata=subsort(2,sort_columns,sdata);
  }
//  ofstream ofs4(adstring("stmp_")+str(psi));
//  ofs4 << setw(6) << sdata << endl;

  ivector nrecaps(1,sim_num_tag_releases);
  int maxlen=int(tag_shlen+tag_nlint*tag_filen);
  nrecaps.initialize();

  for (int it=1;it<=sim_num_tag_releases;it++)
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
  ofs10 << "      " << setw(4) << sim_num_tag_releases << "              "  << tag_shlen 
        << "                  " << tag_nlint << "                   " << tag_filen << endl;
  
  ofs10 << "# TAG RECOVERIES" << endl;

  ofs10 << setw(5) << nrecaps-1 << endl;

  ofs10 << "# XX " << endl;
  int nr=1;
  int it=static_cast<int>(sdata(1,1));
  int itold=it;
  int minage=1;

  ofs10 << "#  " << it << "  RELEASE REGION    YEAR    MONTH   "
          " Tag_program Simulation " << endl;
  int relregion=sim_tag_region(it);
  int trelmonth=sim_true_tag_month(it);
  int trelyear=sim_true_tag_year(it);
  int sitp=sim_initial_tag_period(it,relregion);
  int sti=sim_tag_incident(it);

  dvector ml=value(mean_length(relregion,sitp,sti)(minage,nage));
  dvector sd=value(sdevs(relregion,sitp,sti)(minage,nage));
  dvector lendist(1,tag_nlint);
  for (int il=1;il<=tag_nlint;il++)
  {
    MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
    lendist(il)=csim_initial_tag_release_by_age(it)(minage,nage)
      *exp(-0.5*square(elem_div(len-ml,sd)));
  }

  ofs10 << "         " << relregion << "           " << trelyear 
        << "         " << trelmonth << endl;
  ofs10 <<  sim_num_tags_at_length(it)  << endl;
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
      relregion=sim_tag_region(it);
      trelmonth=sim_true_tag_month(it);
      trelyear=sim_true_tag_year(it);
      int sti=sim_tag_incident(it);
      ofs10 << "#  " << it << "  RELEASE REGION    YEAR    MONTH   "
          " Tag_program Simulation " << endl;
      ofs10 << "            " << relregion << "           " << trelyear 
            << "      " << trelmonth << endl;

      dvector ml=value(mean_length(relregion,sitp,sti)(minage,nage));
      dvector sd=value(sdevs(relregion,sitp,sti)(minage,nage));
      dvector lendist(1,tag_nlint);
      for (int il=1;il<=tag_nlint;il++)
      {
        MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
        lendist(il)=csim_initial_tag_release_by_age(it)(minage,nage)
          *exp(-0.5*square(elem_div(len-ml,sd)));
      }
      lendist/=sum(lendist);
      ofs10 <<  sim_num_tags_at_length(it)  << endl;

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

ivector dvar_len_fish_stock_history::get_sim_num_tags_at_length(int it,
  ivector sim_num_tags_at_age,random_number_generator& rng)   // DF_19dec2019
{
 // int iseed=simulation_seeds(5);    //NMD_19dec2019
//  random_number_generator rng(469);
//  random_number_generator rng(iseed);        //NMD_19dec2019
  int mmin=sim_num_tags_at_age.indexmin();
  int mmax=sim_num_tags_at_age.indexmax();
  int relregion=sim_tag_region(it);
  int sitp=sim_initial_tag_period(it,relregion);
  int sti=sim_tag_incident(it);
  dvector sd=value(sdevs(relregion,sitp,sti));
  dvector ml=value(mean_length(relregion,sitp,sti));
  MY_DOUBLE_TYPE maxlen=tag_shlen+tag_nlint*tag_filen;
  ivector tmp(1,tag_nlint);
  tmp.initialize();
  for (int i=mmin;i<=mmax;i++)
  {
    int ns=sim_num_tags_at_age(i);
    if (ns>0)
    {
      dvector eps(1,ns);
      eps.fill_randn(rng);
      for (int j=1;j<=ns;j++)
      {
        MY_DOUBLE_TYPE len=min(max(shlen,ml(i)+sd(i)*eps(j)),maxlen);
        int ind=min(tag_nlint,int((len-tag_shlen)/tag_filen)+1); //NMD13Dec2017
//        int ind=int((len-tag_shlen)/tag_filen)+1;	
        tmp(ind)+=1;
      }
    }
  }
  return tmp;
}
