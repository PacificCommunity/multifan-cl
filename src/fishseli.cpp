/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"
#include "scbd.hpp"

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg);
//void set_value_inv(_CONST dvar_vector& _x,_CONST ivector& y, dvector& v,
//  int& ii, _CONST MY_DOUBLE_TYPE& fmin, _CONST MY_DOUBLE_TYPE& fmax,_CONST MY_DOUBLE_TYPE& s);
//void set_value(dvar_vector& x,ivector& y,dvar_vector& v, int& ii, 
//  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, dvariable& fpen,_CONST MY_DOUBLE_TYPE& s);
void print_big_effort(dvar_vector& edc,ivector& zcf,ofstream & ofs);

  int size_count(const prevariable& s) {return 1;}


int dvar_fish_stock_history::fishing_selectivity_interface_nvar(void)
{
  int nv=0;
  ivector ff57=column(fish_flags,57);
  int mmin=ff57.indexmin();
  int mmax=ff57.indexmax();
  int seloption=sum(ff57);
  ivector ff64=column(fish_flags,64);
  int ageseloption=sum(ff64);
  if (ageseloption)
  {
    nv+=num_active_partial(ageselcoff,column(fish_flags,64),
      column(fish_flags,24),column(fish_flags,63));
  }
  if (seloption==0)
  {
    //nv+=num_active_partial(selcoff,column(fish_flags,48),
    //  column(fish_flags,24),column(fish_flags,3));
    return nv;
  }
  else
  {
    ivector range=column(fish_flags,3);
    ivector group=column(fish_flags,24);
    ivector flags=column(fish_flags,48);
    ivector option=column(fish_flags,57);
    ivector spline_degree=column(fish_flags,61);
    if (mmin != flags.indexmin() || mmax != flags.indexmax()
      ||  mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "A incompatible array bounds in num_active() " << endl;
      ad_exit(1);
    }
    int maxg=max(group);
    if (maxg)
    {
      int maxf=max(flags);
      if (maxf)
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)==0)
          {
            cerr << "At present for common (grouped) parameters"
                 " all active flags must be on or off " << endl;
          }
        }  
      }
      
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      int flag_value=flags(kin);
      int option_value=option(kin);
      switch(ff57(kin))
      {
      case 0:
       // if (flags(kin)) nv+=size_count_partial(selcoff(kin),range(kin));
        break;
      case 1:  // logistic ?
        if (flags(kin)) nv+=2;
        break;
      case 2:  // double_normal ?
        if (flags(kin)) nv+=3;
        break;
      case 3:  // cubic spline
        if (spline_degree(kin)==0 || spline_degree(kin)>=nage-1)
        {
          cerr << "Illegal value for spline degree in selectivity"
            " parameterization -- value is " << spline_degree(kin) <<
            " it must lie between 1 and " << nage-1 << endl;
          ad_exit(1);
        }
        else
        {
          if (flags(kin))
          {
            if (spline_degree(kin)>=0)
            {
             // nv+=size_count_partial(selcoff(kin),spline_degree(kin));
            }
            else
            {
              switch(spline_degree(kin))
              {
              case -1:
                nv+=2;
                break;
              case -2:
                nv+=3;
                break;
              default:
                cerr << "Illegal value for spline degree in selectivity"
                  " parameterization -- value is " << spline_degree(kin) 
                  << endl;
              }
            }
          }
        }
        break;
      default:
        cerr << "Illegal flag value for fish_flags(" << kin
             << ",57)" << endl;
        ad_exit(1);
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
          if (option(key(i))!=option_value)
          {
            cerr << "Error -- grouped initial parameters have unequal options"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          option_value=option(kin);
          switch(ff57(kin))
          {
          case 0:
           // if (flags(kin)) nv+=size_count_partial(selcoff(kin),range(kin));
            break;
          case 1:  // logistic ?
            if (flags(kin)) nv+=2;
            break;
          case 2:  // logistic ?
            if (flags(kin)) nv+=3;
            break;
          case 3:  // cubic spline
            //if (flags(kin)) nv+=size_count_partial(selcoff(kin),
             //  spline_degree(kin));
            break;
          default:
            cerr << "Illegal flag value for fish_flags(" << kin
                 << ",57)" << endl;
            ad_exit(1);
          }
        } 
      }
    }
    else
    {
      int i=0;
      for (i=mmin;i<=mmax;i++)
      {
        switch(ff57(i))
        {
        case 0:
          //if (flags(i)) nv+=size_count_partial(selcoff(i),range(i));
          break;
        case 1:  // logistic ?
          if (flags(i)) 
          {
            nv+=size_count(fish_pars(9,i));
            nv+=size_count(fish_pars(10,i));
          }
          break;
        case 2:  // double normal
          if (flags(i)) 
          {
            nv+=size_count(fish_pars(9,i));
            nv+=size_count(fish_pars(10,i));
            nv+=size_count(fish_pars(11,i));
          }
          break;
        case 3:
         // if (flags(i)) nv+=size_count_partial(selcoff(i),spline_degree(i));
          break;
        default:
          cerr << "Illegal flag value for fish_flags(" << i
               << ",57)" << endl;
          ad_exit(1);
        }
      }
    }
    return nv;
  }
}
void dvar_fish_stock_history::fishing_selectivity_interface
  (dvar_vector& x, int& ii,const prevariable& _pen)
{
  {
    ivector group=column(fish_flags,24);
    ofstream ofs("group");
    ofs << group << endl;
  }
  ofstream ofs("fishrep");
  prevariable & pen = (prevariable&) _pen;
  // check if new selectivity options are used
  ivector ff57=column(fish_flags,57);
  int mmin=ff57.indexmin();
  int mmax=ff57.indexmax();
  int seloption=sum(ff57);
  ivector ff64=column(fish_flags,64);
  int ageseloption=sum(ff64);
  if (ageseloption)
  {
    set_value_partial(ageselcoff,x,ii,column(fish_flags,63),-10.,10.,pen,
      column(fish_flags,64),column(fish_flags,24),1000.);
  }
  if (seloption==0)
  {
    //set_value_partial(selcoff,x,ii,column(fish_flags,3),-20.,7.,pen,
    //  column(fish_flags,48),column(fish_flags,24),1000.);
  }
  else
  {
    // check options by group
    ivector range=column(fish_flags,3);
    ivector group=column(fish_flags,24);
    ivector flags=column(fish_flags,48);
    ivector option=column(fish_flags,57);
    ivector spline_degree=column(fish_flags,61);
    MY_DOUBLE_TYPE fmin=-20.;
    MY_DOUBLE_TYPE fmax=7.;
    MY_DOUBLE_TYPE scale=1000.;
    
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      int option_value=option(kin);
      switch(ff57(kin))
      {
      case 0:
        //if (flags(kin)) set_value_partial(selcoff(kin),x,ii,range(kin),
        //  fmin,fmax,pen,scale);
        break;
      case 1:  // logistic ?
        if (flags(kin)) 
        {
          set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
          set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
        }
        break;
      case 2:  // double normal
        if (flags(kin)) 
        {
          set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
          set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
          set_value(fish_pars(11,kin),x,ii,-5.,5.,pen,scale);
        }
        break;
      case 3:  // cubic spline
        if (flags(kin)) 
        {
          if (spline_degree(kin)>=0)
          {
            ofs << "fishery " << kin << endl
                << selcoff(kin) << endl
                << x(ii,ii+spline_degree(kin)) << endl;
            //set_value_partial(selcoff(kin),x,ii,spline_degree(kin),
            //  fmin,fmax,pen,scale);
            ofs <<  selcoff(kin) << endl;
          }
          else
          {
            switch(spline_degree(kin))
            {
            case -1:
              set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
              set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
              break;
            case -2:
              set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
              set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
              set_value(fish_pars(11,kin),x,ii,-5.0,5.0,pen,scale);
              break;
            default:
              cerr << "Illegal value for spline degree in selectivity"
                " parameterization -- value is " << spline_degree(kin) 
                << endl;
            }
          }
        }
        break;
      default:
        cerr << "Illegal flag value for fish_flags(" << kin
             << ",57)" << endl;
        ad_exit(1);
      }
      ofstream * pofs=0;
      //ofstream * pofs=new ofstream("fishsel");
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          switch(ff57(kin))
          {
          case 0:
          case 3:
            if (spline_degree(kin)>=0)
            {
              //selcoff(key(i))=selcoff(key(i-1));
//              if (key(i)==5)
//                cout << "A" << endl;
            }
            else
            {
//              if (key(i)==5)
//                cout << "B" << endl;
              switch(spline_degree(kin))
              {
              case -1:
                fish_pars(9,key(i))=fish_pars(9,key(i-1));
                fish_pars(10,key(i))=fish_pars(10,key(i-1));
                break;
              case -2:
                fish_pars(9,key(i))=fish_pars(9,key(i-1));
                fish_pars(10,key(i))=fish_pars(10,key(i-1));
                fish_pars(11,key(i))=fish_pars(11,key(i-1));
                break;
              default:
                cerr << "Illegal value for spline degree in selectivity"
                  " parameterization -- value is " << spline_degree(kin) 
                  << endl;
              }
            }
            break;
          case 1:
//              if (key(i)==5)
//                cout << "C" << endl;
            fish_pars(9,key(i))=fish_pars(9,key(i-1));
            fish_pars(10,key(i))=fish_pars(10,key(i-1));
            break;
          case 2:
//              if (key(i)==5)
//                cout << "D" << endl;
            fish_pars(9,key(i))=fish_pars(9,key(i-1));
            fish_pars(10,key(i))=fish_pars(10,key(i-1));
            fish_pars(11,key(i))=fish_pars(11,key(i-1));
            break;
          default:
            cerr << "Illegal flag value for fish_flags(" << kin
               << ",57)" << endl;
          ad_exit(1);
          }

          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
          if (option(key(i))!=option_value)
          {
            cerr << "Error -- grouped initial parameters have unequal options"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          option_value=option(kin);
          switch(ff57(kin))
          {
          case 0:
//              if (key(i)==5)
//                cout << "D" << endl;
           // if (flags(kin)) set_value_partial(selcoff(kin),x,ii,range(kin),
           //   fmin,fmax,pen,scale);
            break;
          case 1:  // logistic ?
//              if (key(i)==5)
//                cout << "E" << endl;
            if (flags(kin)) 
            {
              set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
              set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
            }
            break;
          case 2:  // double normal
//              if (key(i)==5)
//                cout << "F" << endl;
            if (flags(kin)) 
            {
              set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
              set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
              set_value(fish_pars(11,kin),x,ii,-5.,5.,pen,scale);
            }
            break;
          case 3:
//              if (key(i)==5)
//                cout << "G" << endl;
            if (flags(kin)) 
            {
              if (spline_degree(kin)>=0)
              {
               // set_value_partial(selcoff(kin),x,ii,spline_degree(kin),
                //  fmin,fmax,pen,scale);
              }
              else
              {
                switch(spline_degree(kin))
                {
                case -1:
//              if (key(i)==5)
//                cout << "H" << endl;
                  set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
                  set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
                  break;
                case -2:
//              if (key(i)==5)
//                cout << "I" << endl;
                  set_value(fish_pars(9,kin),x,ii,-5.0,5.0,pen,scale);
                  set_value(fish_pars(10,kin),x,ii,-.9,10.,pen,scale);
                  set_value(fish_pars(11,kin),x,ii,-5.,5.,pen,scale);
                  break;
                default:
                  cerr << "Illegal value for spline degree in selectivity"
                    " parameterization -- value is " << spline_degree(kin) 
                    << endl;
                }
              }
            }
            break;
          default:
            cerr << "Illegal flag value for fish_flags(" << kin
                 << ",57)" << endl;
            ad_exit(1);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        switch(ff57(i))
        {
        case 0:
          //if (flags(i)) set_value_partial(selcoff(i),x,ii,range(i),
          //  fmin,fmax,pen,scale);
          break;
        case 1:  // logistic ?
          if (flags(i)) 
          {
            set_value(fish_pars(9,i),x,ii,-5.0,5.0,pen,scale);
            set_value(fish_pars(10,i),x,ii,-.9,10.,pen,scale);
          }
          break;
          case 2:  // double normal
            if (flags(i)) 
            {
              set_value(fish_pars(9,i),x,ii,-5.0,5.0,pen,scale);
              set_value(fish_pars(10,i),x,ii,-.9,10.,pen,scale);
              set_value(fish_pars(11,i),x,ii,-5.,5.,pen,scale);
            }
            break;
        case 3:
          if (flags(i)) 
          {
            if (spline_degree(i)>=0)
            {
             // set_value_partial(selcoff(i),x,ii,spline_degree(i),
              //  fmin,fmax,pen,scale);
            }
            else
            {
              switch(spline_degree(i))
              {
              case -1:
                set_value(fish_pars(9,i),x,ii,-5.0,5.0,pen,scale);
                set_value(fish_pars(10,i),x,ii,-.9,10.,pen,scale);
                break;
              case -2:
                set_value(fish_pars(9,i),x,ii,-5.0,5.0,pen,scale);
                set_value(fish_pars(10,i),x,ii,-.9,10.,pen,scale);
                set_value(fish_pars(11,i),x,ii,-5.,5.,pen,scale);
                break;
              default:
                cerr << "Illegal value for spline degree in selectivity"
                  " parameterization -- value is " << spline_degree(i) 
                  << endl;
              }
            }
          }
          break;
        default:
          cerr << "Illegal flag value for fish_flags(" << i
               << ",57)" << endl;
          ad_exit(1);
        }
      }
    }
  }
}

void dvar_fish_stock_history::fishing_selectivity_interface_inv
  (dvector& x, int& ii)
{
  ivector ff57=column(fish_flags,57);
  int mmin=ff57.indexmin();
  int mmax=ff57.indexmax();
  int seloption=sum(ff57);
  ivector ff64=column(fish_flags,64);
  int ageseloption=sum(ff64);
  if (ageseloption)
  {
    set_value_inv_partial(value(ageselcoff),x,ii,column(fish_flags,63),
      -10.,10.,column(fish_flags,64),column(fish_flags,24),1000.);
  }
  if (seloption==0)
  {
   // set_value_inv_partial(value(selcoff),x,ii,column(fish_flags,3),-20.,7.,
   //   column(fish_flags,48),column(fish_flags,24),1000.);
  }
  else
  {
    // check options by group
    ivector range=column(fish_flags,3);
    ivector group=column(fish_flags,24);
    ivector flags=column(fish_flags,48);
    ivector option=column(fish_flags,57);
    MY_DOUBLE_TYPE fmin=-20.;
    MY_DOUBLE_TYPE fmax=7.;
    MY_DOUBLE_TYPE scale=1000.;
    ivector spline_degree=column(fish_flags,61);
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      int option_value=option(kin);

      switch(ff57(kin))
      {
      case 0:
       // if (flags(kin)) set_value_inv_partial(value(selcoff(kin)),x,ii,
        //  range(kin),fmin,fmax,scale);
        break;
      case 1:  // logistic ?
        if (flags(kin)) 
        {
          set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
          set_value_inv(fish_pars(10,kin),x,ii,-0.9,10.,scale);
        }
        break;
      case 2:  // double normal
        if (flags(kin)) 
        {
          set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
          set_value_inv(fish_pars(10,kin),x,ii,-0.9,10.,scale);
          set_value_inv(fish_pars(11,kin),x,ii,-5.0,5.0,scale);
        }
        break;
      case 3:
        if (flags(kin)) 
        {
          if (spline_degree(kin)>=0)
          {
            //set_value_inv_partial(value(selcoff(kin)),x,ii,
            //  spline_degree(kin),fmin,fmax,scale);
          }
          else
          {
            switch(spline_degree(kin))
            {
            case -1:
              set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
              set_value_inv(fish_pars(10,kin),x,ii,-.9,10,scale);
              break;
            case -2:
              set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
              set_value_inv(fish_pars(10,kin),x,ii,-.9,10,scale);
              set_value_inv(fish_pars(11,kin),x,ii,-5.0,5.0,scale);
              break;
            default:
              cerr << "Illegal value for spline degree in selectivity"
                " parameterization -- value is " << spline_degree(kin) 
                << endl;
            }
          }
        }
        break;
      default:
        cerr << "Illegal flag value for fish_flags(" << kin
             << ",57)" << endl;
        ad_exit(1);
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
          if (option(key(i))!=option_value)
          {
            cerr << "Error -- grouped initial parameters have unequal options"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          option_value=option(kin);
          switch(ff57(kin))
          {
          case 0:
            //if (flags(kin)) set_value_inv_partial(value(selcoff(kin)),x,ii,
            //  range(kin),fmin,fmax,scale);
            break;
          case 1:  // logistic ?
            if (flags(kin)) 
            {
              set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
              set_value_inv(fish_pars(10,kin),x,ii,-.9,10,scale);
            }
            break;
          case 2:  // double normal
            if (flags(kin)) 
            {
              set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
              set_value_inv(fish_pars(10,kin),x,ii,-0.9,10.,scale);
              set_value_inv(fish_pars(11,kin),x,ii,-5.0,5.0,scale);
            }
            break;
          case 3:
            if (flags(kin)) 
            {
              if (spline_degree(kin)>=0)
              {
                //set_value_inv_partial(value(selcoff(kin)),x,ii,
                //  spline_degree(kin),fmin,fmax,scale);
              }
              else
              {
                switch(spline_degree(kin))
                {
                case -1:
                  set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
                  set_value_inv(fish_pars(10,kin),x,ii,-.9,10,scale);
                  break;
                case -2:
                  set_value_inv(fish_pars(9,kin),x,ii,-5.0,5.0,scale);
                  set_value_inv(fish_pars(10,kin),x,ii,-.9,10,scale);
                  set_value_inv(fish_pars(11,kin),x,ii,-5.0,5.0,scale);
                  break;
                default:
                  cerr << "Illegal value for spline degree in selectivity"
                    " parameterization -- value is " << spline_degree(kin) 
                    << endl;
                  ad_exit(1);
                }
              }
            }
            break;
          default:

            cerr << "Illegal flag value for fish_flags(" << kin
                 << ",57)" << endl;
            ad_exit(1);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        switch(ff57(i))
        {
        case 0:
          //if (flags(i)) set_value_inv_partial(value(selcoff(i)),x,ii,
          //  range(i),fmin,fmax,scale);
          break;
        case 1:  // logistic ?
          if (flags(i)) 
          {
            set_value_inv(fish_pars(9,i),x,ii,-5.0,5.0,scale);
            set_value_inv(fish_pars(10,i),x,ii,-.9,10,scale);
          }
          break;
        case 2:  // double normal
          if (flags(i)) 
          {
            set_value_inv(fish_pars(9,i),x,ii,-5.0,5.0,scale);
            set_value_inv(fish_pars(10,i),x,ii,-.9,10,scale);
            set_value_inv(fish_pars(11,i),x,ii,-5.0,5.0,scale);
          }
          break;
        case 3:
          //if (flags(i)) set_value_inv_partial(value(selcoff(i)),x,ii,
          //  spline_degree(i),fmin,fmax,scale);
          break;
        default:
          cerr << "Illegal flag value for fish_flags(" << i
               << ",57)" << endl;
          ad_exit(1);
        }
      }
    }
  }
}
#undef HOME_VERSION
