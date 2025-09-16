/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvar_matrix dvar_fish_stock_history::get_equilibrium_age_structure
  ( int navg,int nomit,MY_DOUBLE_TYPE lambda,dvar_matrix& tmpcatch)
{
  pmsd_error();
  int ir;

  MY_DOUBLE_TYPE log_lambda=log(lambda);
  int current_year=nyears-navg+1;
  ivector rip(1,num_regions);
  ivector beginning_rip(1,num_regions);
  ivector ending_rip(1,num_regions);
  // set rip to point to the first fishing period in the n-r'th year 
  // for each region
  for (ir=1;ir<=num_regions;ir++)
  {
    int ip=num_fish_periods(ir);
    do
    { 
      ip--;
    }
    while(year(ir,ip)>=nyears-nomit+1);
    ending_rip(ir)=ip+1;
  } 
  for (ir=1;ir<=num_regions;ir++)
  {
    int ip=num_fish_periods(ir);
    do
    { 
      ip--;
    }
    while(year(ir,ip)>=nyears-navg+1);
    beginning_rip(ir)=ip+1;
  } 
  dvar_matrix tmpfish(1,num_regions,1,nage);
  dvar_matrix equilibrium_population(1,num_regions,1,nage);


  dvar3_array local_tot_mort(1,num_regions,
    beginning_rip,ending_rip,1,nage);

  dvar3_array local_survival(1,num_regions,
    beginning_rip,ending_rip,1,nage);

  dvar_vector X(1,num_regions);
  dvar_matrix M(1,num_regions,1,num_regions);

  local_tot_mort.initialize();
  tmpcatch.initialize();
  tmpfish=-100.0;

  imatrix local_nfi(1,num_regions,beginning_rip,ending_rip);
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      local_nfi(ir,ip)=num_fish_incidents(ir,ip);
    }
  }

  dvar4_array local_fish_mort(1,num_regions,
    beginning_rip,ending_rip,1,local_nfi,1,nage);

  dvar4_array local_fish_mort_calcs(1,num_regions,
    beginning_rip,ending_rip,1,local_nfi,1,nage);

  dvar4_array local_catch(1,num_regions,
    beginning_rip,ending_rip,1,local_nfi,1,nage);

      //ofstream ofs("morts"); 
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      if (allocated(local_fish_mort(ir,ip)))
      {
        local_fish_mort(ir,ip)=log_lambda+fish_mort(ir,ip);
      }
      dvar_vector& tm=local_tot_mort(ir,ip);
      dvar_matrix& fm=local_fish_mort(ir,ip);
      
   
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        tm+=mfexp(fm(fi));
      }
      //ofs << tm << endl << endl;
      
      tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    
      local_survival(ir,ip)=mfexp(-tm);
    }
  }

  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      const dvar_vector& tmp1=log(1.e-10+local_tot_mort(ir,ip));
      const dvar_vector& tmp2=log(one_plus-local_survival(ir,ip));
      const dvar_vector& tmp3=tmp2-tmp1;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        local_fish_mort_calcs(ir,ip,fi)=local_fish_mort(ir,ip,fi)+tmp3;
      }
    }
  } 

     
  int finished_flag=1;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  ivector tmp_ip(1,num_regions);
  dvar3_array enum_fish(1,num_regions,beginning_rip,ending_rip,1,nage);
  enum_fish=-100;
  for (ir=1;ir<=num_regions;ir++)
  {
    enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
  }


  // get equilibrium age structure for the first nage-1
  // age classes. loopflag=1 until we have got to last age class.
  int loopflag=1;
  ivector years_covered(1,num_regions);
  years_covered=1;
  int tmpfish_flag=0;
  do
  {
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          // different from ordinary catch equations
          // because we need to get the number of fish after 
          // the last period
          //if (ip>=ending_rip(ir)) break;
          if (ip>ending_rip(ir)) 
          {
            break;
          }
          else if (ip == ending_rip(ir))
          {
            // tmpfish holds the population after each iteration
            tmpfish(ir)(1)=pop_delta(ir);
            --tmpfish(ir)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
            //cout <<"testmsy2.cpp " << setprecision(6) << exp(tmpfish) << endl;
            // this assumes that the last fishing period is at the end 
            // of a  year it may be necessary to exclude a 
            // fishing period or two at the end to ensure this. 
            // there is a switch.
            tmpfish_flag=1;
            years_covered(ir)++;
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break;
          }

          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
          }
          else
          {
            years_covered(ir)++;
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                -local_tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
  
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          if (tmpfish_flag)
          {
            tmpfish_flag=0;
            tmpfish=log(5.e-10+fast_diffusion_calcs(nage,num_regions,tmpfish,
              Dad(tmp_mp(1))));
            for (ir=1;ir<=num_regions;ir++)
            {
              enum_fish(ir,rip(ir))=tmpfish(ir);
            }
          }
          else
          {
            dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
              Dad(tmp_mp(1)),rip);
            for (ir=1;ir<=num_regions;ir++)
            {
              enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
            }
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    
    for (ir=1;ir<=num_regions;ir++)
    {
      enum_fish(ir)=-100;
      enum_fish(ir,beginning_rip(ir))=tmpfish(ir);
      enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
    }
    for (ir=1;ir<=num_regions;ir++)
    {
      if (years_covered(ir)>=nage) 
      {
        check_sanity(years_covered);
        loopflag=0;
        break;
      }
    }
  } // for icount
  while(loopflag);
  // get the contibution to the last age class from other
  // age classes at equilibrium
  equilibrium_population=exp(tmpfish);
  for (ir=1;ir<=num_regions;ir++)
  {
    // zero out last age class
    enum_fish(ir,beginning_rip(ir),nage)=-100.0;
    equilibrium_population(ir,nage)=0.0;
  }

 
  // I don't think we need this because contribution to the last age class
  // from younger age classes is already  in equilibrium_population
  //do
  
 
  rip=beginning_rip;
  tmpfish_flag=0;
  do
  {
    finished_flag=1;
    for (int ir=1;ir<=num_regions;ir++)
    {
      break_flag=0;
      int& ip=rip(ir);
      do
      {
        // different from ordinary catch equations
        // because we need to get the number of fish after 
        // the last period
        //if (ip>=ending_rip(ir)) break;
        if (ip>ending_rip(ir)) 
        {
          break;
        }
        else if (ip == ending_rip(ir))
        {
          // tmpfish holds the population after each iteration
          tmpfish(ir)(1)=pop_delta(ir);
          --tmpfish(ir)(2,nage-1)=
            enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
          tmpfish(ir,nage)=
            log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
              -local_tot_mort(ir,ip,nage-1))
              + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
          tmpfish_flag=1;
          years_covered(ir)++;
          tmp_ip(ir)=ip;
          tmp_mp(ir)=move_index(ir,ip);
          //cout <<"testmsy2.cpp " << setprecision(6) << exp(tmpfish) << endl;
          // this assumes that the last fishing period is at the end 
          // of a  year it may be necessary to exclude a 
          // fishing period or two at the end to ensure this. 
          // there is a switch.
          break;
        }

        finished_flag=0;
        if (year(ir,ip+1)==year(ir,ip))
        {
          enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
        }
        else
        {
          // age the fish
          enum_fish(ir,ip+1,1)=pop_delta(ir);

          if (nage>2)
          --enum_fish(ir,ip+1)(2,nage-1)=
            enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);

          enum_fish(ir,ip+1,nage)=
            log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
              -local_tot_mort(ir,ip,nage-1))
              + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );

        }
        if (move_flags(ir,ip))
        {
          tmp_ip(ir)=ip;
          tmp_mp(ir)=move_index(ir,ip);
          break_flag=1;
        } 
        ip++;
      }
      while (!break_flag); 
    }
    // move the fish
    {
      if (num_regions>1)
      {
        int ir;
        check_sanity(tmp_mp);
        if (tmpfish_flag)
        {
          tmpfish_flag=0;
          tmpfish=log(5.e-10+fast_diffusion_calcs(nage,num_regions,tmpfish,
            Dad(tmp_mp(1))));
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=tmpfish(ir);
          }
        }
        else
        {
          dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
            Dad(tmp_mp(1)),rip);
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    }
  } // need to decide when to quit
  while (!finished_flag);
  
  
  // this is the constribution to the last age class made by younger
  // fish at equilibrium
  for (ir=1;ir<=num_regions;ir++)
  {
    X(ir)=exp(tmpfish(ir)(nage));
  }

  // now get the contribution of one nage class fish in each region
  // we should only need to do the catch equations for the last age class 
  // here

  dmatrix Id=identity_matrix(1,num_regions);
  for (int ir1=1;ir1<=num_regions;ir1++)
  {
    tmpfish_flag=0;
    for (ir=1;ir<=num_regions;ir++)
    {
      enum_fish(ir)=-100;
    }
    enum_fish(ir1,beginning_rip(ir1),nage)=0.0;
    //do
    {
      rip=beginning_rip;
      do
      {
        finished_flag=1;
        for (int ir=1;ir<=num_regions;ir++)
        {
          break_flag=0;
          int& ip=rip(ir);
          do
          {
            // different from ordinary catch equations
            // because we need to get the number of fish after 
            // the last period
            //if (ip>=ending_rip(ir)) break;
            if (ip>ending_rip(ir)) 
            {
              break;
            }
            else if (ip == ending_rip(ir))
            {
              // tmpfish holds the population after each iteration
              //tmpfish(ir)(1)=pop_delta(ir);
              //--tmpfish(ir)(2,nage-1)=
              //  enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
              tmpfish(ir,nage)=enum_fish(ir,ip,nage)
                -local_tot_mort(ir,ip,nage);
              tmpfish_flag=1;
              tmp_ip(ir)=ip;
              tmp_mp(ir)=move_index(ir,ip);
              //cout <<"testmsy2.cpp " << setprecision(6) << exp(tmpfish) << endl;
              // this assumes that the last fishing period is at the end 
              // of a  year it may be necessary to exclude a 
              // fishing period or two at the end to ensure this. 
              // there is a switch.
              //years_covered(ir)++;
              break;
            }

            finished_flag=0;
            if (year(ir,ip+1)==year(ir,ip))
            {
              enum_fish(ir,ip+1,nage)=enum_fish(ir,ip,nage)
                -local_tot_mort(ir,ip,nage);
            }
            else
            {
              // age the fish
              // no recruitment
              //enum_fish(ir,ip+1,1)=pop_delta(ir);
    
              //if (nage>2)
              //--enum_fish(ir,ip+1)(2,nage-1)=
              //  enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
    
              enum_fish(ir,ip+1,nage)=enum_fish(ir,ip,nage)
                -local_tot_mort(ir,ip,nage);
                //log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                //  -local_tot_mort(ir,ip,nage-1))
                //  + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
    
            }
            if (move_flags(ir,ip))
            {
              tmp_ip(ir)=ip;
              tmp_mp(ir)=move_index(ir,ip);
              break_flag=1;
            } 
            ip++;
          }
          while (!break_flag); 
        }
        // move the fish
        {
          if (num_regions>1)
          {
            int ir;
            check_sanity(tmp_mp);
            if (tmpfish_flag)
            {
              tmpfish_flag=0;
              tmpfish=log(5.e-10+fast_diffusion_calcs(nage,num_regions,tmpfish,
                Dad(tmp_mp(1))));
              for (ir=1;ir<=num_regions;ir++)
              {
                enum_fish(ir,rip(ir))=tmpfish(ir);
              }
            }
            else
            {
              dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
                Dad(tmp_mp(1)),rip);
              for (ir=1;ir<=num_regions;ir++)
              {
                enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
              }
            }
          }
        }
      } // need to decide when to quit
      while (!finished_flag);
      
      for (ir=1;ir<=num_regions;ir++)
      {
        M(ir,ir1)= exp(tmpfish(ir)(nage));
        //M(ir,ir1)= exp(enum_fish(ir,ending_rip(ir),nage));
      }
    } // for icount
    //while(1);
  }
  //cout <<"testmsy2.cpp " << M << endl;
  // solve for equilbrium terminal age class numbers
  dvar_vector tn=inv(Id-M)*X; 
  //cout <<"testmsy2.cpp " << setprecision(6) << setscientific() << inv(Id-M)*X << endl; 
  for (ir=1;ir<=num_regions;ir++)
  {
    // set last age class
    equilibrium_population(ir,nage)=tn(ir);
  }
  tmpfish=log(equilibrium_population);
  
  // Now can test to see if equilibrium population is
  // really at equilibrium
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        tmpcatch(ir)+=mfexp(local_fish_mort_calcs(ir,ip,fi)
          +tmpfish(ir));
      }
    }
  }
  tmpcatch/=double(navg);

  int iswitch=0;
  if (iswitch==0)
  {
    return equilibrium_population;
  }

  for (ir=1;ir<=num_regions;ir++)
  {
    enum_fish(ir,beginning_rip(ir))=
      log(equilibrium_population(ir));
  }
  //do
  {
    tmpfish_flag=0;
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          //if (ip>=ending_rip(ir)) break;
          if (ip>ending_rip(ir)) 
          {
            break;
          }
          else if (ip == ending_rip(ir))
          {
            // tmpfish holds the population after each iteration
            tmpfish(ir)(1)=pop_delta(ir);
            --tmpfish(ir)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
            //cout <<"testmsy2.cpp " << setprecision(6) << exp(tmpfish) << endl;
            // this assumes that the last fishing period is at the end 
            // of a  year it may be necessary to exclude a 
            // fishing period or two at the end to ensure this. 
            // there is a switch.
            tmpfish_flag=1;
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break;
          }
          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
          }
          else
          {
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                -local_tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          if (tmpfish_flag)
          {
            tmpfish_flag=0;
            tmpfish=log(5.e-10+fast_diffusion_calcs(nage,num_regions,tmpfish,
              Dad(tmp_mp(1))));
            for (ir=1;ir<=num_regions;ir++)
            {
              enum_fish(ir,rip(ir))=tmpfish(ir);
            }
          }
          else
          {
            dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
              Dad(tmp_mp(1)),rip);
            for (ir=1;ir<=num_regions;ir++)
            {
              enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
            }
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    
   
    dvar_matrix test_population(1,num_regions,1,nage);
    
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpfish(ir)=enum_fish(ir,ending_rip(ir));
    }
    test_population=exp(tmpfish);
    cout << "XXX" << setscientific() << setprecision(16)  << norm2(equilibrium_population-test_population)<< endl;
    cout << "XXX" << setscientific() << setprecision(16)  << norm2(equilibrium_population)<< endl;
    cout << "XXX" << setscientific() << setprecision(16)  << norm2(test_population)<< endl;
    ofstream ofs("eqpop");
    for (ir=1;ir<=num_regions;ir++)
    {
      ofs << setprecision(5) << setscientific() 
          << equilibrium_population(ir,nage) << " " 
          << test_population(ir,nage) << endl;
    }
   
  } // for icount
  return equilibrium_population;
}

dvar_matrix dvar_fish_stock_history::get_equilibrium_age_structure
  ( int navg,int nomit,const prevariable& lambda,
    dvar_matrix& tmpcatch)
{
  int ir;

  dvariable log_lambda=log(lambda);
  int current_year=nyears-navg+1;
  ivector rip(1,num_regions);
  ivector beginning_rip(1,num_regions);
  // set rip to point to the first fishing period in the n-r'th year 
  // for each region
  ivector ending_rip(1,num_regions);
  for (ir=1;ir<=num_regions;ir++)
  {
    int ip=num_fish_periods(ir);
    do
    { 
      ip--;
    }
    while(year(ir,ip)>=nyears-nomit+1);
    ending_rip(ir)=ip+1;
  } 
  for (ir=1;ir<=num_regions;ir++)
  {
    int ip=ending_rip(ir);
    do
    { 
      ip--;
    }
    while(year(ir,ip)>=nyears-navg+1);
    beginning_rip(ir)=ip+1;
  } 
  dvar_matrix tmpfish(1,num_regions,1,nage);
  dvar_matrix equilibrium_population(1,num_regions,1,nage);


  dvar3_array local_tot_mort(1,num_regions,
    beginning_rip,ending_rip,1,nage);

  dvar3_array local_survival(1,num_regions,
    beginning_rip,ending_rip,1,nage);

  dvar_vector X(1,num_regions);
  dvar_matrix M(1,num_regions,1,num_regions);

  local_tot_mort.initialize();
  tmpcatch.initialize();

  imatrix local_nfi(1,num_regions,beginning_rip,ending_rip);
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      local_nfi(ir,ip)=num_fish_incidents(ir,ip);
    }
  }

  dvar4_array local_fish_mort(1,num_regions,
    beginning_rip,ending_rip,1,local_nfi,1,nage);

  dvar4_array local_fish_mort_calcs(1,num_regions,
    beginning_rip,ending_rip,1,local_nfi,1,nage);

  dvar4_array local_catch(1,num_regions,
    beginning_rip,ending_rip,1,local_nfi,1,nage);

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      if (allocated(local_fish_mort(ir,ip)))
      {
        local_fish_mort(ir,ip)=log_lambda+fish_mort(ir,ip);
      }
      dvar_vector& tm=local_tot_mort(ir,ip);
      dvar_matrix& fm=local_fish_mort(ir,ip);
      
   
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        tm+=mfexp(fm(fi));
      }
      
      tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    
      local_survival(ir,ip)=mfexp(-tm);
    }
  }

  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      const dvar_vector& tmp1=log(1.e-10+local_tot_mort(ir,ip));
      const dvar_vector& tmp2=log(one_plus-local_survival(ir,ip));
      const dvar_vector& tmp3=tmp2-tmp1;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        local_fish_mort_calcs(ir,ip,fi)=local_fish_mort(ir,ip,fi)+tmp3;
      }
    }
  } 

     
  int finished_flag=1;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  ivector tmp_ip(1,num_regions);
  dvar3_array enum_fish(1,num_regions,beginning_rip,ending_rip,1,nage);
  enum_fish=-100;
  for (ir=1;ir<=num_regions;ir++)
  {
    enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
  }


  // get equilibrium age structure for the first nage-1
  // age classes
  int loopflag=1;
  do
  {
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          if (ip>=ending_rip(ir)) break;
          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
          }
          else
          {
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                -local_tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
  
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
            Dad(tmp_mp(1)),rip);
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpfish(ir)(1)=pop_delta(ir);
      tmpfish(ir)(2,nage)=++enum_fish(ir,ending_rip(ir))(1,nage-1);
      tmpfish(ir)(nage)=log(exp(tmpfish(ir)(nage))
        +exp(enum_fish(ir,ending_rip(ir))(nage)));
      enum_fish(ir)=-100;
      enum_fish(ir,beginning_rip(ir))=tmpfish(ir);
      enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
    }
    for (ir=1;ir<=num_regions;ir++)
    {
      if (exp(tmpfish(ir)(nage))>1.e-6) 
      {
        loopflag=0;
        break;
      }
    }
  } // for icount
  while(loopflag);
  // get the contibution to the last age class from other
  // age classes at equilibrium
  equilibrium_population=exp(tmpfish);
  for (ir=1;ir<=num_regions;ir++)
  {
    // zero out last age class
    enum_fish(ir,beginning_rip(ir),nage)=-100.0;
    equilibrium_population(ir,nage)=0.0;
  }

  //loopflag=1;
  //do
  // just do this once I think
  {
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          if (ip>=ending_rip(ir)) break;
          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
          }
          else
          {
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                -local_tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
  
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
            Dad(tmp_mp(1)),rip);
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    
    
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpfish(ir)(1)=pop_delta(ir);
      tmpfish(ir)(2,nage)=++enum_fish(ir,ending_rip(ir))(1,nage-1);
      tmpfish(ir)(nage)=log(exp(tmpfish(ir)(nage))
        +exp(enum_fish(ir,ending_rip(ir))(nage)));
      enum_fish(ir)=-100;
      enum_fish(ir,beginning_rip(ir))=tmpfish(ir);
      enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
    }
    /*
    for (ir=1;ir<=num_regions;ir++)
    {
      if (exp(tmpfish(ir)(nage))>1.e-6) 
      {
        loopflag=0;
        break;
      }
    }
    */
  } // for icount
  //while(loopflag);
  for (ir=1;ir<=num_regions;ir++)
  {
    X(ir)=exp(tmpfish(ir)(nage));
  }

  // now get the contribution of one nage class fish in each region

  dmatrix Id=identity_matrix(1,num_regions);
  for (int ir1=1;ir1<=num_regions;ir1++)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      enum_fish(ir)=-100;
    }
    enum_fish(ir1,beginning_rip(ir1),nage)=0.0;
    //do
    {
      rip=beginning_rip;
      do
      {
        finished_flag=1;
        for (int ir=1;ir<=num_regions;ir++)
        {
          break_flag=0;
          int& ip=rip(ir);
          do
          {
            if (ip>=ending_rip(ir)) break;
            finished_flag=0;
            if (year(ir,ip+1)==year(ir,ip))
            {
              enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
            }
            else
            {
              // age the fish
              // no recruitment
              //enum_fish(ir,ip+1,1)=pop_delta(ir);
    
              if (nage>2)
              --enum_fish(ir,ip+1)(2,nage-1)=
                enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
    
              enum_fish(ir,ip+1,nage)=
                log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                  -local_tot_mort(ir,ip,nage-1))
                  + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
    
            }
            if (move_flags(ir,ip))
            {
              tmp_ip(ir)=ip;
              tmp_mp(ir)=move_index(ir,ip);
              break_flag=1;
            } 
            ip++;
          }
          while (!break_flag); 
        }
        // move the fish
        {
          if (num_regions>1)
          {
            int ir;
            check_sanity(tmp_mp);
            dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
              Dad(tmp_mp(1)),rip);
            for (ir=1;ir<=num_regions;ir++)
            {
              enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
            }
          }
        }
      } // need to decide when to quit
      while (!finished_flag);
      
      for (ir=1;ir<=num_regions;ir++)
      {
        //M(ir,ir1)= exp(tmpfish(ir)(nage));
        M(ir,ir1)= exp(enum_fish(ir,ending_rip(ir),nage));
      }
    } // for icount
    //while(1);
  }
  // solve for equilbrium terminal age class numbers
  dvar_vector tn=inv(Id-M)*X; 
  cout <<"testmsy2.cpp " << setprecision(6) << setscientific() << inv(Id-M)*X << endl; 
  for (ir=1;ir<=num_regions;ir++)
  {
    // set last age class
    equilibrium_population(ir,nage)=tn(ir);
  }
  tmpfish=log(equilibrium_population);
  
  // Now can test to see if equilibrium population is
  // really at equilibrium
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=ending_rip(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        tmpcatch(ir)+=mfexp(local_fish_mort_calcs(ir,ip,fi)
          +tmpfish(ir));
      }
    }
  }
  tmpcatch/=double(navg);

  return equilibrium_population;

  for (ir=1;ir<=num_regions;ir++)
  {
    enum_fish(ir,beginning_rip(ir))=
      log(equilibrium_population(ir));
  }
  //do
  {
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          if (ip>=ending_rip(ir)) break;
          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
          }
          else
          {
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                -local_tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
            Dad(tmp_mp(1)),rip);
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    
   
    dvar_matrix test_population(1,num_regions,1,nage);
    
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpfish(ir)(1)=pop_delta(ir);
      tmpfish(ir)(2,nage)=++enum_fish(ir,ending_rip(ir))(1,nage-1);
      tmpfish(ir)(nage)=log(exp(tmpfish(ir)(nage))
        +exp(enum_fish(ir,ending_rip(ir))(nage)));
      enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
    }
    test_population=exp(tmpfish);
    cout <<"testmsy2.cpp " << norm2(equilibrium_population-test_population)<< endl;
    ofstream ofs("eqpop");
    for (ir=1;ir<=num_regions;ir++)
    {
      ofs << setprecision(5) << setscientific() 
          << equilibrium_population(ir,nage) << " " 
          << test_population(ir,nage) << endl;
    }
   
  } // for icount
}

