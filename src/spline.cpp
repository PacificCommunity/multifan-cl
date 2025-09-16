/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
  
void block_error();
extern dvar_vector * const_vb;
extern dvar_vector * const_var;

dvariable lenbound(const dvar_vector& vb,int age,int nage,dvar_vector ap3,
		   ivector parest_flags)
{
  dvariable rho=exp(-vb(3));
  dvariable tmp=vb(1)+(vb(2)-vb(1))*
    (1.0-pow(rho,age-1))/(1.0-pow(rho,nage-1));
  if (parest_flags(173) && (age>1 && age<=parest_flags(173)))
  {
    tmp+= ap3(age-1);
  }
  return tmp;
}

void dvar_len_fish_stock_history::get_scaled_two_species_lengths(void)  
{
  // this routine is needed on the condition
  //      if (age_flags(193) && pmsd)
  //
  ivector nd=column(fish_flags,3);
  ivector ff26=column(fish_flags,26);
  ivector ff57=column(fish_flags,57);
  ivector ff74=column(fish_flags,74);
  ivector ff75=column(fish_flags,75);
  ivector ff71=column(fish_flags,71);
  if (!allocated(tlength))
    tlength.allocate(1,num_fisheries,ff75+1,nd);
  int i,j;
  for (i=1;i<=num_fisheries;i++)
  {
    if (nd(i)==0) nd(i)=nage;
  }

  dvar_vector ap3;
  ap3=age_pars(3);
  dvar_vector sap3;
  if(pmsd) sap3=pmsd->age_pars(2,3);

  dvariable rho=0.0;
  dvariable l1=0.0;
  dvariable ln=0.0;
  dvariable l2=0.0; // Length at a=ff75+1
  dvariable l3=0.0; // length at a=nd=ff3
  dvariable smaxln=0;
  dvariable sminl1=0;
  ofstream * zzofs=0;
  for (i=1;i<=num_fisheries;i++)
  {
    int ir=fishery_regions(i);
    if (!pmsd  || pmsd->region_species_pointer(ir)==1)
    {
      rho=exp(-vb_coff(3));
      l1=vb_coff(1);
      ln=vb_coff(2);
    }
    else
    {
      rho=exp(-pmsd->vb_coff(pmsd->region_species_pointer(ir),3));
      l1=pmsd->vb_coff(pmsd->region_species_pointer(ir),1);
      ln=pmsd->vb_coff(pmsd->region_species_pointer(ir),2);
    }
    if (!parest_flags(173))
    {
      l2=l1+(ln-l1)*(1.0-pow(rho,ff75(i)+1-1))/
        (1.0-pow(rho,nage-1)); // YT 2017-3-22
      l3=l1+(ln-l1)*(1.0-pow(rho,nd(i)-1))/
        (1.0-pow(rho,nage-1));   // YT 2017-3-22
    }
    else
    {
      if (!pmsd  || pmsd->region_species_pointer(ir)==1)
      {
        if(ff75(i))
        {
          l2=l1+(ln-l1)*(1.0-pow(rho,ff75(i)+1-1))/
            (1.0-pow(rho,nage-1))+ap3(ff75(i));
        }
        else
        {
          l2=l1+(ln-l1)*(1.0-pow(rho,ff75(i)+1-1))/
            (1.0-pow(rho,nage-1));
        }
        l3=l1+(ln-l1)*(1.0-pow(rho,nd(i)-1))/
          (1.0-pow(rho,nage-1))+ap3((nd(i)-2));
      }
      else
      {
        if(ff75(i))
        {
          l2=l1+(ln-l1)*(1.0-pow(rho,ff75(i)+1-1))/
            (1.0-pow(rho,nage-1))+sap3(ff75(i));
        }
        else
        {
          l2=l1+(ln-l1)*(1.0-pow(rho,ff75(i)+1-1))/
            (1.0-pow(rho,nage-1));
        }
        l3=l1+(ln-l1)*(1.0-pow(rho,nd(i)-1))/
          (1.0-pow(rho,nage-1))+sap3((nd(i)-2));
      }
    }
    int mmin=tlength(i).indexmin();
    if (age_flags(193) && pmsd)
    {
      dvariable div=1.0/(1.0-pow(rho,nd(i)-mmin));
      if (pmsd->num_species !=2)
      {
        cerr << "option for af192 only implemented for two species/sexes"
             << endl;
        ad_exit(1);
      }

//      smaxln=smax(lenbound(vb_coff,nd(i),nage),
//        lenbound(pmsd->vb_coff(2),nd(i),nage));
//    This is for the case with different ff75 across sex YT 2017-4-12 
      int nrf=pmsd->num_real_fisheries;
      int i1=i;
      int i2=i;
      if(i<=nrf){
        i1=i;
        i2=i+nrf; // This should be fishery of 2ned species/gender
      }else{
        i1=i-nrf;
        i2=i;   // This should be fishery of 2ned species/gender
      }      
      // 1st nd(i) was replaced with nd(i1), 2nd nd(i) was replaced with nd(i2)
//      smaxln=smax(lenbound(vb_coff,nd(i1),nage,age_pars,parest_flags),   
//        lenbound(pmsd->vb_coff(2),nd(i2),nage,age_pars,parest_flags));   

//      sminl1=smin(lenbound(vb_coff,ff75(i1)+1,nage,age_pars,parest_flags),   
//          lenbound(pmsd->vb_coff(2),ff75(i2)+1,nage,age_pars,parest_flags));  
      dvar_vector ap3=age_pars(3);
      dvar_vector sap3=pmsd->age_pars(2,3);

      dvariable bnd1=0;
      bnd1=lenbound(vb_coff,nd(i1),nage,ap3,parest_flags);
      dvariable bnd2=0;
      bnd2=lenbound(pmsd->vb_coff(2),nd(i2),nage,sap3,parest_flags);
      smaxln=smax(bnd1,bnd2);

      bnd1=lenbound(vb_coff,ff75(i1)+1,nage,ap3,parest_flags);
      bnd2=lenbound(pmsd->vb_coff(2),ff75(i1)+1,nage,sap3,parest_flags);
      sminl1=smin(bnd1,bnd2);
	
      pmsd->cminlength=sminl1;
      pmsd->cmaxlength=smaxln;
      dvariable ll=(l2-sminl1)/(smaxln-sminl1);
      dvariable uu=(l3-sminl1)/(smaxln-sminl1);
      dvar_vector groffset;

      if (parest_flags(173))
      {
        groffset.allocate(mmin+1,parest_flags(173));
        // check for case where mmin<2 ff75=0
        if (!pmsd  || pmsd->region_species_pointer(ir)==1)
        {
          for (j=mmin+1;j<=parest_flags(173);j++)
          {
            groffset(j)=ap3(j-1);
            groffset(j)=groffset(j)/(smaxln-sminl1);
          }
        }
        else
        {
          for (j=mmin+1;j<=parest_flags(173);j++)
          {
            groffset(j)=sap3(j-1);
            groffset(j)=groffset(j)/(smaxln-sminl1);
          }
        }
      }
      
      tlength(i,mmin)=ll;
      for (j=mmin+1;j<=nd(i);j++)
      {
        tlength(i,j)=ll+(uu-ll)*(1.0-pow(rho,j-mmin))*div;
      }
      for (j=mmin+1;j<=nd(i);j++)
      {
        if (parest_flags(173) && j<=parest_flags(173))
        {
          tlength(i,j)+=groffset(j);
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::get_scaled_lengths(void)  
{
  int i,j;
  ivector nd=column(fish_flags,3);
  ivector ff26=column(fish_flags,26);
  ivector ff57=column(fish_flags,57);
  ivector ff74=column(fish_flags,74);
  ivector ff75=column(fish_flags,75);
  ivector ff71=column(fish_flags,71);
  if (!allocated(tlength))
    tlength.allocate(1,num_fisheries,ff75+1,nd);
  for (i=1;i<=num_fisheries;i++)
  {
    int ir=fishery_regions(i);
    dvariable rho;
    if (!pmsd  || pmsd->region_species_pointer(ir)==1)
    {
      rho=exp(-vb_coff(3));
    }
    else
    {
      rho=exp(-pmsd->vb_coff(pmsd->region_species_pointer(ir),3));
    }
    int mmin=tlength(i).indexmin();
   
    dvariable div=1.0/(1.0-pow(rho,nd(i)-1));
    tlength(i,mmin)=0.0;
    for (j=mmin+1;j<=nd(i);j++)
    {
      tlength(i,j)=(1.0-pow(rho,j-mmin))*div;
    }
  }
}
void dvar_len_fish_stock_history::
  get_scaled_single_or_two_species_lengths(void)  
{
  if (age_flags(193) && pmsd)
  {
    if (pmsd->num_species !=2)
    {
      cerr << "option for af192 only implemented for two species/sexes"
           << endl;
      ad_exit(1);
    }
    get_scaled_two_species_lengths();
  }
  else
  {
    get_scaled_lengths();
  }
}

void dvar_len_fish_stock_history::splines(void)  
{
  ivector nd=column(fish_flags,3);
  ivector ff26=column(fish_flags,26);
  ivector ff57=column(fish_flags,57);
  ivector ff74=column(fish_flags,74);
  ivector ff75=column(fish_flags,75);
  ivector ff71=column(fish_flags,71);
  ofstream * zzofs=0;
  if (generate_report)
  {
   // zzofs=new ofstream("spline_selectivity");
  }
  int i,j;
  for (i=1;i<=num_fisheries;i++)
  {
    if (nd(i)==0) nd(i)=nage;
  }

  /*
  cout << tlength(1) << endl;
  cout << 2.0*tlength(1)-1 << endl;
  cout << tlength(15) << endl;
  cout << 2.0*tlength(15)-1 << endl;
  ad_exit(1);
  if (zzofs) 
  {
    (*zzofs) << "scaled lengths at age" << endl;
    (*zzofs) << tlength << endl;
  }
  */
  ivector tmplength(1,num_fisheries);
  for (i=1;i<=num_fisheries;i++)
  {
    if (ff26(i)==3)
      tmplength(i)=nlint;
    else
      tmplength(i)=nd(i);
  }
  
  dvar_matrix tempsel(1,num_fisheries,1,tmplength);
  for (i=1;i<=num_fisheries;i++)
  {
    if (fish_flags(i,57) ==3 )
    {
      int sd=fish_flags(i,61);
      if (sd>=0)
      {
        dvector x(1,sd);
        dvector xx(1,nd(i));
        MY_DOUBLE_TYPE longlen=fmid(nlint)+0.5*filen;
    
        if (ff26(i)==3)
        {
          x.fill_seqadd(shlen,(longlen-shlen)/(sd-1));
        }
        else
        {
          x.fill_seqadd(0,1.0/(sd-1));
          xx.fill_seqadd(0,1.0/(nd(i)-1));
        }
        for (int is=1;is<=ff74(i);is++)
        {
          for (int ib=1;ib<=ff71(i)+1;ib++)
          {
            dvar_vector tsel(1,sd);
            tsel=bs_selcoff(i,is,ib)(1,sd);
            bs_selmean(i,is,ib)=log(mean(exp(tsel)));
            tsel-=bs_selmean(i,is,ib);
            vcubic_spline_function csf(x,tsel);
            int np=100;
            dvector points(1,np);
#if !defined(NO_MY_DOUBLE_TYPE)
            points.fill_seqadd(0,1.0/(np-1.0L));
#else
            points.fill_seqadd(0,1.0/(np-1.0));
#endif
            switch (ff26(i))
            {
            case 0:
              break;
            case 1:
            case 2:
                {
                  bstempsel(i,is,ib)=csf(tlength(i));
                  /*
                  if (zzofs) 
                  {
                    dvar_vector sizes=csf(points);
                    dvar_matrix tmp(1,2,1,np);
                    tmp(1)=llvector(i)+(uuvector(i)-llvector(i))*points;
                    tmp(2)=exp(sizes);
                    tmp(2)/=max(tmp(2));
                    (*zzofs) << "fishery " << i << " " << trans(tmp) << endl;
                  }
                  */
                }
              break;
            case 3:
              bstempsel(i,is,ib)=csf(fmid);
              //cout << "BBB TEMPORARY DAVE" << endl;
              //bsplinesel(i,is,ib)=mfexp(bstempsel(i,is,ib));
        
              if (nwint>0)
              {
                bswtempsel(i,is,ib)=csf(wmid);  //NMD_20jan2023
                //bswsplinesel(i,is,ib)=mfexp(csf(wmid));
              }
              break;
            default:
              cerr << "illegal value for fish_flags(" << i <<",61)"
                   << endl;
              ad_exit(1);
            }
          }
        }
        if (generate_report)
        {
          for (int is=1;is<=ff74(i);is++)
          {
            for (int ib=1;ib<=ff71(i)+1;ib++)
            {
              dvar_vector tsel(1,sd);
              tsel=bs_selcoff(i,is,ib)(1,sd);
              bs_selmean(i,is,ib)=log(mean(exp(tsel)));
              tsel-=bs_selmean(i,is,ib);
              vcubic_spline_function csf(x,tsel);
              switch (ff26(i))
              {
              case 0:
                bstempsel(i,is,ib)=csf(xx);     
                break;
              case 1:
              case 2:
                  bstempsel(i,is,ib)=csf(tlength(i));
                  if (zzofs) 
                  {
                    (*zzofs) << "fishery " << i << " " << exp(bstempsel(i,is,ib)) << endl;
                  }
                break;
              case 3:
                bstempsel(i,is,ib)=csf(fmid);
                //cout << "BBB TEMPORARY DAVE" << endl;
                //bsplinesel(i,is,ib)=mfexp(bstempsel(i,is,ib));
          
                if (nwint>0)
                {
                  bswtempsel(i,is,ib)=csf(wmid);  //NMD_20jan2023
                  //bswsplinesel(i,is,ib)=mfexp(csf(wmid));
                }
                break;
              default:
                cerr << "illegal value for fish_flags(" << i <<",61)"
                     << endl;
                ad_exit(1);
              }
            }
          }
        }
      }
      else
      {
        switch(sd)
        {
        case -1:
          {
            const prevariable & fp9=fish_pars(9,i);
            dvariable fp10=1.0+fish_pars(10,i);
            dvariable nineteen=19.0;
            {
              MY_DOUBLE_TYPE jmid=0.5*(1+nlint);
              MY_DOUBLE_TYPE jmid1=1.0-jmid;
              
              splinesel(i)= 1.0/(1+pow(nineteen,((fmid-jmid)/jmid1-fp9)/fp10));
            }
            {
              MY_DOUBLE_TYPE jmid=0.5*(1+nwint);
              MY_DOUBLE_TYPE jmid1=1.0-jmid;
              wsplinesel(i)= 1.0/(1+pow(nineteen,((wmid-jmid)/jmid1-fp9)/fp10));
            }
          }
          break;
        case -2:
          {
            const prevariable & fp9=fish_pars(9,i);
            dvariable fp10=1.0+fish_pars(10,i);
            const dvariable & efp11=exp(fish_pars(11,i));
            MY_DOUBLE_TYPE jmid=0.5*(1+nlint);
            MY_DOUBLE_TYPE jmid1=1.0-jmid;
            int j;
            if (!allocated(splinesel(i)))
              splinesel(i).allocate(1,nlint);
              
            for (j=1;j<=nlint;j++)
            {
              MY_DOUBLE_TYPE fj=(j-jmid)/jmid1;
              if (fj<=value(fp9))
              {
                splinesel(i,j)= pow(2,-square(fj/(fp10*efp11)));
              }
              else
              {
                splinesel(i,j)= pow(2,-square(fj*efp11/fp10));
              }
            }
            MY_DOUBLE_TYPE wjmid=0.5*(1+nwint);
            MY_DOUBLE_TYPE wjmid1=1.0-jmid;
            if (!allocated(wsplinesel(i)))
              wsplinesel(i).allocate(1,nwint);
            for (j=1;j<=nwint;j++)
            {
              MY_DOUBLE_TYPE fj=(j-wjmid)/wjmid1;
              if (fj<=value(fp9))
              {
                wsplinesel(i,j)= pow(2,-square(fj/(fp10*efp11)));
              }
              else
              {
                wsplinesel(i,j)= pow(2,-square(fj*efp11/fp10));
              }
            }
          }    
          break;
        default:
          cerr << "illegal value for fish_flags(" << i <<",61)"
               << endl;
          ad_exit(1);
        }
      }
    }
  }
}
