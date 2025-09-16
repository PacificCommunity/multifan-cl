/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#if !defined(__MSVC32__)
//#  pragma interface
#endif
#ifndef NEWMULT_HPP
#define NEWMULT_HPP

#include <cifstrem.h>
#include <adstring.hpp>
//#define catch Catch
#include <stdlib.h>
#include <memory>
#define CLOSE close

#if (defined(MY_REAL_DOUBLE))
#  include "extras.h"
#endif

void check_group_flag_sanity(const ivector& flag,const ivector& group);
void check_group_for_holes(const ivector& v);
int getdim(int ff57,int ff61,int nage);
#include "old_set.hpp"

extern adtimer * profile_timer_1;
extern adtimer * profile_timer_2;
extern ofstream * clogf;
extern int ad_argc;
extern int _data_file_version_no; 
extern int fish_nonmv_flag;   //NMD_14Sep2018 
/*
extern int pvm_switch;
extern int pvm_save_switch;
*/
extern char ** ad_argv;
extern adstring ad_root;

#if (!defined(__MSVC32__) || __MSVC32__<=7)
class mf_pvm_manager;
extern mf_pvm_manager * mf_pvm;
#endif

using namespace std;    //NMD_24Mar2015

void mfcl_error(void);
// The class fish_stock_history contains all the relevant data and
// parameters which determine the exploitation of the stock
// the data are ordered by fishing period. A fishing period is
// a time (month and week of month for length frequency purposes)
// during which at least one fishing incident occurs. A fishing
// incident is an instance of a fishery.  A historical fishery
// consists of a time series of fishing incidents. Fishing incdients
// are grouped into fisheries because fisheries are supposed to reflect
// regularities in the fishing process which means that different
// fishing incidents within a fishery share parameters -- the idea is
// to reduce the number of free parameters to be estimated to a minimum.
// A realization of a fishery is a fishing incident for that fishery
// For each realization of a fishery the realization_period is the
// fishing period during which that particular realization occurred,
// while the realization incident is the fishing incident to which that
// particular realization corresponded.

   void check_index(int index_min,int index_max,int i,char *);
   class par_cifstream;
   class par_ifstream;
   class par_ofstream;

  class dvar_fish_stock_history;
  class fishery_catch_at_age_record_array;
  class fishery_freq_record_array;
  class fishery_header_record_array;
  class fishery_header_record;

 class cubic_spline_array;
  typedef cubic_spline_array * pcubic_spline_array;

  class group_manager_1;

class dvar_len_fish_stock_history;

class xgroup_manager  
{
  int ngroups;   // number of groups things are grouped into
  ivector ng;     // the number in each group
  ivector group_ptr; // picks an ungrouped vector 
  ivector inv_group_ptr; // picks an ungrouped vector 
public:  
  int allocated(void) 
  { 
    if (ngroups>0) 
      return 1;
    else
      return 0;
  }
  xgroup_manager(void) : ngroups(0) {}  
  void safe_allocate(const ivector& flags,const ivector& groups);
  int get_ngroups(void){ return ngroups;}
  ivector& get_group_ptr(void){ return group_ptr;}
  ivector& get_inv_group_ptr(void){ return inv_group_ptr;}
};


class newstyle_frq_data
{
public:
  int frq_file_version;
  ivector new_frq_flags;  // flags read from top; of new style frq file
  adstring_array fishery_names; // the user supplied names for each fishery
  imatrix new_fish_flags;  // flags read from top; of new style frq file
  //newstyle_frq_data& operator = (newstyle_frq_data& t);
};

  class header_record
  {
  public:
    int year;
    int movement_period;
    int month;
    int week;
    header_record& operator = (header_record& hr);
    friend MY_DOUBLE_TYPE decimal_time(header_record& hr);
    //friend int operator > (header_record& fr1,header_record& fr);
    //friend int operator < (header_record& fr1,header_record& fr);
    friend istream& operator >> (istream& ifs,header_record& fhr);
    friend cifstream& operator >> (cifstream& ifs,header_record& fhr);
    friend ostream& operator << (ostream& ofs,header_record& fhr);
    header_record(){ year=0;month=0;week=0;}
    header_record(header_record& hr)
    {
      year=hr.year;
      movement_period=hr.movement_period;
      month=hr.month;
      week=hr.week;
    }

    int& get_year(void) { return year;}
    int& get_month(void) { return month;}
    int& get_week(void) { return week;}

    //friend ostream& operator << (ostream& ofs,header_record& fhr);


  }; // end of class header_record

  class freq_record : public header_record // encapsulates the data for a
					   // length frequency record
  {
  protected:
    dvector freq;

  public:
    freq_record() : freq() {}

    freq_record(int nlfreq) : freq(1,nlfreq) {}

    void allocate(int nlint)
    {
      freq.allocate(1,nlint);
      freq.initialize();
    }
    freq_record& operator = (freq_record& fr1);

    friend istream& operator >> (istream& infile, freq_record& fr);

    friend cifstream& operator >> (cifstream& infile, freq_record& fr);

    virtual int operator > (freq_record& fr)
    {
      if (year<fr.year)
      {
	return 0;
      }
      if (year>fr.year)
      {
	return 1;
      }
      if (month<fr.month)
      {
	return 0;
      }
      if (month>fr.month)
      {
	return 1;
      }
      if (week<fr.week)
      {
	return 0;
      }
      if (week>fr.week)
      {
	return 1;
      }
      return 0;
    }
  }; // end of class freq_record


  class fishery_header_record : public header_record
  {
  //protected:
  public:
    int fishery;
    MY_DOUBLE_TYPE obs_tot_catch;
    MY_DOUBLE_TYPE effort;
    MY_DOUBLE_TYPE effort_weight;
    ivector species;
    ivector species1;
    ivector species2;
    ivector species3;
    friend class fishery_header_record_array;
    friend ostream& operator << (ostream& ifs,fishery_header_record& fhr);
    friend istream& operator >> (istream& ifs,fishery_header_record& fhr);
    friend cifstream& operator >> (cifstream& ifs,fishery_header_record& fhr);
    friend int operator > (fishery_header_record& fr1,
					    fishery_header_record& fr);
    friend int operator < (fishery_header_record& fr1,
					    fishery_header_record& fr);

    //friend ostream& operator << (ostream& ofs,fishery_header_record& fhr);

  public:
    fishery_header_record() {fishery=0;}

    fishery_header_record(fishery_header_record& fhr) : header_record(fhr),
      species(fhr.species),species1(fhr.species1)
    {
      fishery    = fhr.fishery;
      obs_tot_catch  = fhr.obs_tot_catch;
      effort = fhr.effort;
      effort_weight = fhr.effort_weight;
    }

    int get_species(void)
    {
      int ssum=sum(species);
      if (ssum!=1)
      {
        cerr << "Error in fishery_header_record -- ssum = " << ssum << endl;
        ad_exit(1);
      }
      for (int i=species.indexmin();i<=species.indexmax();i++)
      {
        if (species(i)==1)
          return i;
      }
      cerr << "Error in fishery_header_record  -- no species"
           << endl << species << endl;
      ad_exit(1);
    }

    fishery_header_record& operator = (const fishery_header_record& fhr)
     // : header_record(fhr)
    {
      year = fhr.year;
      movement_period=fhr.movement_period;
      month = fhr.month;
      week = fhr.week;
      fishery = fhr.fishery;
      obs_tot_catch = fhr.obs_tot_catch;
      effort = fhr.effort;
      effort_weight = fhr.effort_weight;
      if (allocated(fhr.species))
      {
        species.allocate(fhr.species.indexmin(),fhr.species.indexmax());
        species=fhr.species;
      }
      if (allocated(fhr.species1))
      {
        species1.allocate(fhr.species1.indexmin(),fhr.species1.indexmax());
        species1=fhr.species1;
      }
      if (allocated(fhr.species2))
      {
        species2.allocate(fhr.species2.indexmin(),fhr.species2.indexmax());
        species2=fhr.species2;
      }
      if (allocated(fhr.species3))
      {
        species3.allocate(fhr.species3.indexmin(),fhr.species3.indexmax());
        species3=fhr.species3;
      }
      return *this;
    }

    int& get_fishery(void) { return fishery;}
    MY_DOUBLE_TYPE& get_obs_tot_catch(void) { return obs_tot_catch;}
    MY_DOUBLE_TYPE& get_effort(void) { return effort; }

  }; // end of class fishery_header_record

  class fishery_data_record : public fishery_header_record
  {
  protected:
    virtual void fit(void)=0;
    fishery_data_record(void) : fishery_header_record() {}
  };


  class fishery_data_record_ptr
  {
    fishery_data_record * ptr;
    friend class fishery_data_record_array;
  };

  class fishery_data_record_array
  {
    int index_min;
    int index_max;
    int indexmin(){return index_min;}
    int indexmax(){return index_max;}
    int size(){return index_max-index_min+1;}
    fishery_data_record_ptr * ptr;
  public:

    fishery_data_record_array(int min,int max)
    {
      index_min=min;
      index_max=max;
      int sz=size();
      ptr=new fishery_data_record_ptr[sz];
      if (ptr==NULL)
      {
	cerr << "Error allocating memory in"
		" fishery_data_record_array(int min,int max)" <<endl;
	exit(1);
      }
    }

    fishery_data_record& operator [] (int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i,
	   "fishery_data_record::operator [] (int i)");
      #endif
      return *(ptr[i].ptr);
    }

    fishery_data_record& elem(int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i,
	   "fishery_data_record::operator [] (int i)");
      #endif
      return *(ptr[i].ptr);
    }

  };


  class fishery_catch_at_age_record : public fishery_data_record
			    // encapsulates the catch at age data
			    // associated with a fishery
  {
    friend class fishery_catch_at_age_record_array;
    friend class fishery_header_record;
  protected:
    //fishery_catch_at_age_record(void) : fishery_header_record() {}
    fishery_catch_at_age_record(void){};

    dvector props;
    void allocate(int nage)
    {
      props.allocate(1,nage);
      props.initialize();
    }

    void allocate(fishery_header_record& fhr,int ng)
    {
      fishery    = fhr.fishery;
      obs_tot_catch  = fhr.obs_tot_catch;
      effort = fhr.effort;
      effort_weight = fhr.effort_weight;
      props.allocate(1,ng);
      props.initialize();
    }

    friend ostream& operator <<
	   (ostream& ofs,fishery_catch_at_age_record& fhr);
    friend istream& operator >>
	   (istream& ifs,fishery_catch_at_age_record& fhr);
    friend cifstream& operator >>
	   (cifstream& ifs,fishery_catch_at_age_record& fhr);
    void fit(void){ cout << "fit in catch_at_age_record"<<endl;}
  }; // end class fishery_catch_at_age_record



  class fishery_freq_record : public fishery_data_record
			    // encapsulates the data
			   // for a length frequency associated with a fishery
  {
  //protected:
  public:
    static int real_nwint;
    static int age_nage;
    static int age_age1;
    static int wgroup_factor;
    static int real_nlint;
    static int group_factor;
    dvector freq;
    dvector wfreq;
    dvector age_freq;

    friend class fishery_freq_record_array;
    friend class fishery_header_record;

    void allocate(int nlint,int nwint);

    fishery_freq_record& operator = (fishery_freq_record& fr1);

    friend istream& operator >> (istream& infile, fishery_freq_record& fr);

    friend ostream& operator << (ostream& ofs,fishery_freq_record& fhr);

    friend cifstream& operator >> (cifstream& infile, fishery_freq_record& fr);

    fishery_freq_record(int nlfreq) : freq(1,nlfreq){}
    fishery_freq_record(void){}
    void fit(void){ cout << "fit in freq_record"<<endl;}

    MY_DOUBLE_TYPE& get_freq(int i) { return freq(i); }

  }; //  end class fishery_freq_record


  class freq_record_array
  // this class is intended to hold an array of length frequency
  // records in such a form that they can be read in sequentially from
  // a user's FRQ file. then they will probably be sorted to be used in
  // a fish_stock_history structure
  {
  public:
    int            nlfreq;
    int            index_min;
    int            index_max;
    int            nlint;  // number of length intervals
    MY_DOUBLE_TYPE         shlen;  // the shortest length
    MY_DOUBLE_TYPE         filen;  // the width of the length intervals
  private:
    freq_record *  ptr;

  public:
     int& nli(void) { return nlint;}
     MY_DOUBLE_TYPE& wli(void) { return filen; }
     MY_DOUBLE_TYPE& lmin(void) { return shlen; }

    int indexmin(){return index_min;}
    int indexmax(){return index_max;}
    int size(){return index_max-index_min+1;}
    freq_record& operator [] (int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i," freq_record& operator [] (int i)");
      #endif
      return ptr[i];
    }

    freq_record& operator () (int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i," freq_record& operator [] (int i)");
      #endif
      return ptr[i];
    }
  // !!!!!!!!!!!!!!!  temporary patch 
  //  void allocate(void);
  
    void allocate(void)
    {
      for (int i=indexmin();i<=indexmax();i++)
      {
	((*this)(i)).freq_record::allocate(nlint);
      }
    }
  

    freq_record_array(int min,int max,int nfrq,int _nlint,
      int _shlen,int _filen)
    {
      nlfreq=nfrq;
      index_min=min;
      index_max=max;
      nlint=_nlint;
      shlen=_shlen;
      filen=_filen;
      int sz=size();
      ptr=new freq_record [sz];
      if (ptr==NULL)
      {
	cerr << "Error allocating memory in freq_record_array"<<endl;
	exit(21);
      }
      ptr-=indexmin();
      allocate();
    }
    ~freq_record_array()
    {
      if (ptr)
      {
        ptr+=indexmin();
        delete [] ptr;
        ptr=NULL;
      }
    }
  };  // class freq_record_array
  class multi_species_data;
  typedef multi_species_data * pmulti_species_data;

  class fishery_freq_record_array
  {
  // this class is intended to hold an array of FISHERY length frequency
  // records in such a form that they can be read in sequentially from
  // a user's FRQ file. Then they will probably be sorted to be used in
  // a fish_stock_history structure.
  public:
    int            nlfreq; // number of len frq data sets
    int            nlint;  // number of length intervals
    MY_DOUBLE_TYPE         shlen;  // the shortest length
    MY_DOUBLE_TYPE         filen;  // the width of the length intervals
    int            nwint;  // number of length intervals
    MY_DOUBLE_TYPE         wshlen;  // the shortest length
    MY_DOUBLE_TYPE         wfilen;  // the width of the length intervals
    int            index_min;
    int            index_max;
    ivector        fishery_regions;
    //ivector        regmin;
    //ivector        regmax;
    //ivector        num_fish_times;
    //imatrix        realization_period;
    //imatrix        num_fish_periods;
    //imatrix realization_incident;
    //imatrix realization_region;
    //imatrix header_record_index;
    void get_header_array_index(ivector& nft, 
      imatrix& realization_period,imatrix& realization_incident,
      imatrix& header_record_index);
    fishery_freq_record * ptr;
    fishery_freq_record * ptr1;   // add these for femail data
  private:
    void allocate();
  public:

    int& nli(void) { return nlint;}
    MY_DOUBLE_TYPE& wli(void) { return filen; }
    MY_DOUBLE_TYPE& lmin(void) { return shlen; }
    int& nsamp(void) { return nlfreq; }

    int indexmin(){return index_min;}
    int indexmax(){return index_max;}
    int size(){return index_max-index_min+1;}
    int numfreq(){return nlfreq;}
    fishery_freq_record& operator [] (int i);

    fishery_freq_record& operator () (int i);

    fishery_freq_record& elem(int i)
    {
      return ptr[i];
    }
    fishery_freq_record_array(int min,int max,int _nlint,MY_DOUBLE_TYPE _shlen,
      MY_DOUBLE_TYPE _filen, int _nwint,MY_DOUBLE_TYPE _wshlen,MY_DOUBLE_TYPE _wfilen);
    //fishery_freq_record_array(len_fish_stock_history& fsh);

    ~fishery_freq_record_array();

  
    dvar_len_fish_stock_history get_history_data(int ntg,int num_regions,
      int nage,ivector& parest_flags,ivector&,ivector&,ivector&,
      imatrix&,dvector&,int,int _mfactor,int& _first_time,ivector&,
      int,imatrix& ses_reg);

  
    dvar_len_fish_stock_history get_history_data(int ntg,int num_regions,
      int nage,ivector& parest_flags,ivector&,ivector&,ivector&,
      imatrix&,dvector&,int,int _mfactor,int& _first_time,ivector&,
      int,imatrix& ses_reg,pmulti_species_data & pmsd);
  

    friend ostream& operator << (ostream& ofs,
		      fishery_freq_record_array& fcara);

    friend ofstream& operator << (ofstream& ofs,
		      fishery_freq_record_array& fcara);

    void sort(void);

    int& get_fishery(int i) { return elem(i).get_fishery();}
    MY_DOUBLE_TYPE& get_obs_tot_catch(int i) { return elem(i).get_obs_tot_catch();}
    MY_DOUBLE_TYPE& get_effort(int i) { return elem(i).get_effort(); }
    int& get_year(int i) { return elem(i).get_year();}
    int& get_month(int i) { return elem(i).get_month();}
    int& get_week(int i) { return elem(i).get_week();}

    MY_DOUBLE_TYPE& get_freq(int i, int j) { return elem(i).get_freq(j); }

  };  // end class fishery_freq_record_array

  class fishery_catch_at_age_record_array
  {
  // this class is intended to hold an array of FISHERY catch at age
  // records in such a form that they can be read in sequentially from
  // a user's FRQ file. Then they will probably be sorted to be used in
  // a fish_stock_history structure.

    int            nage;
    int            index_min;
    int            index_max;
    fishery_catch_at_age_record * ptr;
    void allocate(int nage)
    {
      for (int i=indexmin();i<=indexmax();i++)
      {
	(*this)(i).allocate(nage);
      }
    }

    void allocate(fishery_header_record_array& fhra,int ng);
  public:
    int indexmin(){return index_min;}
    int indexmax(){return index_max;}
    int size(){return index_max-index_min+1;}
    int numage(){return nage;}
    fishery_catch_at_age_record& operator [] (int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i," fishery_catch_at_age_record& operator [] (int i)");
      #endif
      return ptr[i];
    }

    fishery_catch_at_age_record& elem(int i)
    {
      #ifdef SAFE_ARRAYS
        check_index(index_min,index_max,i," fishery_catch_at_age_record& elem(int i)");
      #endif
      return ptr[i];
    }

    fishery_catch_at_age_record& operator () (int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i," freq_record& operator [] (int i)");
      #endif
      return ptr[i];
    }

    fishery_catch_at_age_record_array(int min,int max,int ng)
    {
      nage=ng;
      index_min=min;
      index_max=max;
      int sz=size();
      ptr=new fishery_catch_at_age_record [sz];
      if (ptr==NULL)
      {
        cerr << "Error allocating memory in fishery_freq_record_array"<<endl;
        exit(21);
      }
      ptr-=indexmin();
      allocate(nage);
    }

    friend ostream& operator <<
             (ostream& ofs, fishery_catch_at_age_record_array& fcara);


    fishery_catch_at_age_record_array(fishery_header_record_array& fhra,
      int ng);

    ~fishery_catch_at_age_record_array()
    {
      ptr+=indexmin();
      delete [] ptr;
      ptr=NULL;
    }

    void sort(void);
  };  // end class fishery_catch_at_age_record_array

  fishery_freq_record_array read_fisheries_data(int total_fishing_incidents,
    int nlfreq,adstring& frq_file_name);

  class void_pointer
  {
    void * ptr;   // can change void to whatever for "safe" programing
  public:
    void_pointer& operator = (void * vptr) { ptr=vptr; return (*this); }
  };

  class void_pointer_array
  {
  // this class is intended to hold an array of pointers to
  // anything
    int            index_min;
    int            index_max;
    void_pointer * ptr;
  public:
    int indexmin(){return index_min;}
    int indexmax(){return index_max;}
    int size(){return index_max-index_min+1;}

    void_pointer& operator [] (int i)
    {
      #ifdef SAFE_ARRAYS
        check_index(index_min,index_max,i," fishery_freq_record& operator [] (int i)");
      #endif
      return ptr[i];
    }
    void_pointer_array(int min,int max)
    {
      index_min=min;
      index_max=max;
      int sz=size();
      ptr=new void_pointer [sz];
      if (ptr==NULL)
      {
        cerr << "Error allocating memory in void_pointer_array"<<endl;
        exit(21);
      }
      ptr-=indexmin();
    }

    ~void_pointer_array()
    {
      ptr+=indexmin();
      delete [] ptr;
      ptr=NULL;
    }
  }; // end class void_pointer_array

  class fishery_freq_record_pointer
  {
    fishery_freq_record * ptr;   // can change void to whatever for "safe" programing
  public:
    fishery_freq_record_pointer& operator = (fishery_freq_record * vptr)
    {
      ptr=vptr;
      return (*this);
    }
    fishery_freq_record_pointer& operator = (fishery_freq_record_pointer&  v)
    {
      ptr=v.ptr;
      return (*this);
    }
    fishery_freq_record& operator * (void)
    {
      return *ptr;
    }

  };

  class fishery_freq_record_pointer_array
  {
  // this class is intended to hold an array of pointers to
  // anything
    int            index_min;
    int            index_max;
    fishery_freq_record_pointer * ptr;
  public:
    int indexmin(){return index_min;}
    int indexmax(){return index_max;}
    int size(){return index_max-index_min+1;}

    fishery_freq_record_pointer& operator [] (int i)
    {
      #ifdef SAFE_ARRAYS
        check_index(index_min,index_max,i," fishery_freq_record& operator [] (int i)");
      #endif
      return ptr[i];
    }
    fishery_freq_record_pointer& operator () (int i)
    {
      #ifdef SAFE_ARRAYS
        check_index(index_min,index_max,i," fishery_freq_record& operator [] (int i)");
      #endif
      return ptr[i];
    }
    fishery_freq_record_pointer_array(int min,int max)
    {
      index_min=min;
      index_max=max;
      int sz=size();
      ptr=new fishery_freq_record_pointer [sz];
      if (ptr==NULL)
      {
        cerr << "Error allocating memory in fishery_freq_record_pointer_array"<<endl;
        exit(21);
      }
      ptr-=indexmin();
    }

    ~fishery_freq_record_pointer_array()
    {
      ptr+=indexmin();
      delete [] ptr;
      ptr=NULL;
    }
  }; // end class fishery_freq_record_pointer_array

  class multi_species_data
  {
  public:
    int tag_index;
    imatrix tag_region_bounds;
    ivector nage_by_region;
    ivector degree_yr;
    ivector degree_reg;
    ivector degree_ses;
    ivector degree_ses_reg;
    d4_array orth_recr_basis;
    dvar3_array      orth_recr_all; 
    int  get_num_real_regions(void) { return num_real_regions;}
    int combined_tags_flag;
    ivector num_tag_release_by_species;
    ivector sim_num_tag_release_by_species;
    i4_array fisc; // fish_incident_species_code;
    i4_array fisc_lf; // fish_incident_species_code;
    i4_array fisc_wf; // fish_incident_species_code;
    i3_array fisn; // fish_incident_species_code;
    i3_array fisn_lf; // fish_incident_species_code;
    i3_array fisn_wf; // fish_incident_species_code;
    imatrix  tag_species_index;
    i4_array sp_in_catch;
    i4_array reg_in_catch;
    i4_array sp_in_lf;
    i4_array reg_in_lf;
    i4_array sp_in_wf;
    i4_array reg_in_wf;
    dvar_matrix rec_delta;
    dvar_vector recmean;
    dvar_vector totpop_coff;
    dvar_vector totpop;
    dvar_vector rec_init_diff;
    dvar_matrix recr;
    dvar_matrix recr1;
    dvar_matrix bh_recr_devs;
    dvar_matrix bh_predicted_recruits;   //NMD_12Apr2021
    imatrix parest_flags;
    imatrix age_flags;
    imatrix historical_age_flags;  //NMD_23Apr2015
    imatrix historical_parest_flags;  //NMD_23Apr2015
    ivector nage;
    dvar_matrix sv;
    int current_species;
    dvar_matrix pmature;  // percent maturity for species 2 to num_species
    dvar_matrix cpmature_at_length;
    // percent maturity at length for species 2 to num_species NMD_aug28_2018
    dvar3_array age_pars; //NMD 4Nov2011
	imatrix region_bounds;
    int num_species;
    int num_real_regions;
    int num_real_fisheries;
    ivector num_tag_groups;
    ivector tag_species_pointer;
    imatrix numcomp;
    i3_array ses_reg_recr_flags;
    imatrix tag_species_flag;
    imatrix sim_tag_species_flag;
    imatrix species_by_region;
    ivector region_species_pointer;
    ivector num_region_by_species;
    ivector rowoffset;
    ivector kludged_fishery_regions;
    int num_kludged_regions;
    void allocate(int num_species,int num_regions);
    multi_species_data(int num_species,int num_regions);
    multi_species_data(void);
    dvar_vector nat_mort_coff;  //NMD 02Nov2011
    dvar3_array nat_mort;       // The annual mortality rate for each species, year and age class //NMD 02Nov2011
    dvar_matrix vb_coff;
    dvector fmin1;
    dvector fmax1;
    dvector fminl;
    dvector fmaxl;
    dvector rhomin;
    dvector rhomax;
    dvar_matrix var_coff;
    dvar_matrix global_vars;
    dvector vmin1;
    dvector vmax1;
    dvector vminl;
    dvector vmaxl;
    dvector len_wt_coff;
    dvar_vector phi;
    dvar_vector alpha;
    dvar_vector beta;
    dvar_vector steepness;
    dmatrix wmid;
    dmatrix wm2;
    dvar_matrix region_rec_diff_sums;
    dvar3_array sdevs_yr;
    dvar3_array length_month_yr;
    dvar3_array F_by_age_by_year;
    imatrix species_flags;
    dvariable cminlength;   // min and max of two sex growth curves used
    dvariable cmaxlength;   // for common rescaling to [0,1]
    ivector recr_degree_yr;
    ivector recr_degree_reg;
    ivector recr_degree_ses;
    ivector num_new_weights;
    dvar_matrix new_orth_recr;
    d3_array OR;
    d3_array recr_polys_yr;
    d3_array recr_polys_reg;
    d3_array recr_polys_ses;
    d3_array recr_polys_ses_reg;
    ivector  ny_begin_yr;
    ivector  ny_begin_reg;
    ivector  ny_begin_ses;
    ivector  ny_begin_ses_reg;
    ivector  nd_begin_yr;
    ivector  nd_begin_reg;
    ivector  nd_begin_ses;
    ivector  nd_begin_ses_reg;
    ivector  ny_end_yr;
    ivector  ny_end_reg;
    ivector  ny_end_ses;
    ivector  ny_end_ses_reg;
    ivector  nd_end_yr;
    ivector  nd_end_reg;
    ivector  nd_end_ses;
    ivector  nd_end_ses_reg;
    dvar3_array yearly_recr_all;  

  };  

  typedef multi_species_data * pmulti_species_data;

  class fishery_header_record_array
  {
    int index_min;
    int index_max;
    fishery_header_record * ptr;
  public:

    int indexmin(void){return index_min;}
    int indexmax(void){return index_max;}
    int size(void){return index_max-index_min+1;}
    fishery_header_record_array(int min,int max)
    {
      index_min=min;
      index_max=max;
      int sz=size();
      ptr=new fishery_header_record [sz];
      if (ptr==NULL)
      {
        cerr << "Error allocating memory in fishery_header_record_array"<<endl;
        exit(21);
      }
      ptr-=min;
    }
    ~fishery_header_record_array()
    {
      ptr+=indexmin();
      delete [] ptr;
      ptr=NULL;
    }
    fishery_header_record& operator () (int i)
    {
      #ifdef SAFE_ARRAYS
        check_index(index_min,index_max,i," fishery_header_record_array::"
                 "operator () (int i)");
      #endif
      return ptr[i];
    }
    fishery_header_record& operator [] (int i)
    {
      #ifdef SAFE_ARRAYS
	check_index(index_min,index_max,i," fishery_header_record_array::"
		 "operator [] (int i)");
      #endif
      return ptr[i];
    }
    fishery_header_record& elem (int i)
    {
      #ifdef SAFE_ARRAYS
        check_index(index_min,index_max,i," fishery_header_record_array::"
                 "operator () (int i)");
      #endif
      return ptr[i];
    }

  
  dvar_fish_stock_history get_history_data(int ntg, int nage,int num_regions,
      ivector& parest_flags,ivector&,ivector&,ivector&,imatrix&,int,
      int _mfactor,int& _first_time,ivector&,
      int,imatrix& ses_reg);
 
    dvar_fish_stock_history get_history_data(int ntg,
      int num_regions,int nage,ivector& _parest_flags,ivector& regmin,
      ivector& regmax,ivector& dataswitch,imatrix& Dflags,int month_1,
      int _mfactor,int& _first_time,ivector& mw,int _direction_flag,
      imatrix& _season_region_flags,pmulti_species_data & pmsd);
    

  }; // end class fishery_header_record_array


#ifdef __ZTC__
  class par_ifstream : public ifstream
  {
  public:
    par_ifstream(adstring& s) : ifstream( (char *) s ) , ios(&buffer){}
    par_ifstream(char * s) : ifstream(s) , ios(&buffer) {}
    //friend par_ifstream& operator >> (par_ifstream& pif,
					 fish_stock_history& fsh);
  }; // end class par_ifstream

  class par_cifstream : public cifstream
  {
  public:
    par_cifstream(adstring& s) : cifstream( (char *) s) , ios(&buffer){}
    par_cifstream(char * s) : cifstream(s)  , ios(&buffer){}
    friend par_cifstream& operator >> (par_cifstream& pif, dvar_fish_stock_history& fsh);
  }; // end class par_cifstream


  class par_ofstream : public ofstream
  {
  public:
    par_ofstream(adstring& s) : ofstream( (char *) s), ios(&buffer) {}
    par_ofstream(char * s) : ofstream(s), ios(&buffer) {}
   // friend par_ofstream& operator << (par_ofstream& pof,
	//				      fish_stock_history& fsh);
    friend par_ofstream& operator << (par_ofstream& pof,
					      dvar_fish_stock_history& fsh);
  }; // end class par_ofstream

#else
  class par_ifstream : public ifstream
  {
  public:
    par_ifstream(adstring& s) : ifstream( (char *) s ) {}
    par_ifstream(char * s) : ifstream(s) {}
    //friend par_ifstream& operator >> (par_ifstream& pif,
//					 fish_stock_history& fsh);
  }; // end class par_ifstream


  class par_cifstream : public cifstream
  {
  public:
    par_cifstream(adstring& s) : cifstream( (char *) s) {}
    par_cifstream(char * s) : cifstream(s) {}
    friend par_cifstream& operator >> (par_cifstream& pif, dvar_fish_stock_history& fsh);
  }; // end class par_cifstream


  class par_ofstream : public ofstream
  {
  public:
    par_ofstream(adstring& s) : ofstream( (char *) s) {}
    par_ofstream(const char * s) : ofstream(s) {}    //NMD_24mar2023
//    par_ofstream(char * s) : ofstream(s) {}
   // friend par_ofstream& operator << (par_ofstream& pof,
	//				      fish_stock_history& fsh);
    friend par_ofstream& operator << (par_ofstream& pof,
					      dvar_fish_stock_history& fsh);
  }; // end class par_ofstream

#endif

par_ofstream& operator << (par_ofstream& pof,int);
par_ofstream& operator << (par_ofstream& pof,const char *);
par_ofstream& operator << (par_ofstream& pof,const prevariable&);
par_ofstream& operator << (par_ofstream& pof,const dvar_vector&);
par_ofstream& operator << (par_ofstream& pof,const dvar_matrix&);
par_ofstream& operator << (par_ofstream& pof,const dvar3_array&);
par_ofstream& operator << (par_ofstream& pof,MY_DOUBLE_TYPE);
par_ofstream& operator << (par_ofstream& pof,const dvector&);
par_ofstream& operator << (par_ofstream& pof,const ivector&);
par_ofstream& operator << (par_ofstream& pof,const dmatrix&);
par_ofstream& operator << (par_ofstream& pof,const d3_array&);

//ofstream& operator << (ofstream& pof,fish_stock_history& fsh);
//ofstream& operator << (ofstream& pof,len_fish_stock_history& fsh);

int max(int x,int y);
int min(int x,int y);


ifstream& operator >> (ifstream& ifs,fishery_header_record& fhr);
ifstream& operator >> (ifstream& ifs,fishery_header_record_array& fhra);
cifstream& operator >> (cifstream& ifs,fishery_header_record& fhr);
cifstream& operator >> (cifstream& ifs,fishery_header_record_array& fhra);
//cifstream& operator >> (cifstream& cif, fish_stock_history& fsh);
cifstream& operator >> (cifstream& ifs,
		      fishery_catch_at_age_record_array& fcara);
istream& operator >> (istream& ifs,
		      fishery_catch_at_age_record_array& fcara);
istream& operator >>
           (istream& ifs,fishery_catch_at_age_record& fhr);
cifstream& operator >>
	   (cifstream& ifs,fishery_catch_at_age_record& fhr);
cifstream& operator >> (cifstream& ifs,header_record& fhr);
cifstream& operator >> (cifstream& ofs,fishery_header_record& fhr);
//ostream& operator << (ostream& ofs, fish_stock_history& fsh);

cifstream& operator >> (cifstream& ifs,
		      fishery_freq_record_array& fcara);
istream& operator >> (istream& ifs,
		      fishery_freq_record_array& fcara);
ostream& operator <<
	   (ostream& ifs,fishery_freq_record_array& fhr);

ostream& operator <<
	   (ostream& ifs,fishery_freq_record& fhr);

void read_again(dvar_len_fish_stock_history& fsh,adstring& file_root_par);

//uostream& operator << (uostream& us, len_fish_stock_history& fsh);
//uostream& operator << (uostream& us, fish_stock_history& fsh);
//uistream& operator >> (uistream& us, fish_stock_history& fsh);
//uistream& operator >> (uistream& us, len_fish_stock_history& fsh);
uostream& operator << (uostream& us, dvar_len_fish_stock_history& fsh);
uostream& operator << (uostream& us, dvar_fish_stock_history& fsh);
int is_ascii_file(const char* file_name);

class index_weight
{
  dvector wt;
public:
  index_weight(int max);
  MY_DOUBLE_TYPE operator () (int n){ return wt(n);}
  MY_DOUBLE_TYPE operator [] (int n){ return wt(n);}
};

int read_newstyle_file(adstring& full_datafile_path);

void read_newstyle_data(cifstream& infile,newstyle_frq_data& tmpfrq);

void do_between_time_calculations(dvar_fish_stock_history& fsh);
void greport(const char * s);

void length_basedsel_calc(dvar_matrix& lengthbsel,dvar_fish_stock_history& fsh,
  dvar_vector& vb_coff,dvar_vector& var_coff,int nlint, dvector& fmid);

dvariable setup_alpha(int j,dvar_vector& sv,int nage,ivector& age_flags);

MY_DOUBLE_TYPE setup_alpha(int j,dvector& sv,int nage,ivector& age_flags);

void read_environmental_data(adstring& full_datafile_path,
  dvar_len_fish_stock_history& fsh);

dvariable exploitation_penalty(dvar_len_fish_stock_history& fsh);

void lbsel_calc(dvar_matrix& lbsel,dvar_fish_stock_history& fsh,
  dvar_vector& vb_coff,dvar_vector& var_coff);

dvar_vector cutoff_sel(dvar_fish_stock_history& fsh, dvar_vector& vb_coff,
  dvar_vector& var_coff);

dvariable age_at_length_calc(const prevariable& v,const dvar_vector& vb_coff,int nage);

dvariable age_at_length_calc(const prevariable& v,const dvar_vector& vb_coff,int nage,
  const dvar_vector& sigma);
  
dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvar_vector& vb_coff,int nage,
  const dvar_vector& sigma);

dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvar_vector& vb_coff,int nage,
  const ivector & pf);

dvar_matrix cutoff_sel_report(dvar_fish_stock_history& fsh, dvar_vector& vb_coff,
  dvar_vector& var_coff);


int option_match(int argc,char * argv[], char * string);
int option_match(int argc,char * argv[], char * string,const int& nopt);

par_cifstream& operator >> (par_cifstream& pif, dvar_len_fish_stock_history& fsh);

void fitting_procedure(dvar_len_fish_stock_history& fsh);
dvariable objective_function(dvar_fish_stock_history& fsh);
void sort(fishery_freq_record_array& flfra,ivector& fishery_regions,
  fishery_freq_record_array& flfra1,ivector& mmin,ivector& mmax);
MY_DOUBLE_TYPE fcomp(BOR_CONST dvar_len_fish_stock_history& fsh,BOR_CONST dvar_vector& x,
  int print_switch);
void write_report(ostream& ofs,dvar_fish_stock_history& fsh);
imatrix normalize_year(fishery_freq_record_array& fra,int _month_1,int& df,
  int& m1,int& first_data_month);
int sub_main(adstring& directory_path,adstring& full_parfile_path,
  adstring& full_datafile_path,adstring& full_output_parfile_path);
int read_newstyle_file(adstring& full_datafile_path);
void do_switch_changes(imatrix& cls,dvar_len_fish_stock_history& fsh);
void set_command_line_switches(imatrix& cl_switches,int argc,
  char * argv[],int opos);
void getroot(adstring& adr,adstring& fdp);

int greater_than(fishery_freq_record **tmp,int j, ivector& fr);
    
void deal_with_command_line(int argc,char * argv[],
  adstring& full_datafile_path,adstring& full_input_parfile_path,
  adstring& full_output_parfile_path);

void set_gtradient_structure_globals(void);

void set_region_incidence_matrix(int nswitch,cifstream& infile,
  ivector& dataswitch,imatrix& Dflags,int& num_fisheries,int& num_regions,
  ivector& fishery_regions,int& tag_group_flag,dvector& _region_area,int & _data_file_version_no);

void setup_Dyf_age_dep(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags);


#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(BOR_CONST prevariable& );
#endif

//void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
//  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
//  dvar_vector& diff_coffs3,ivector& age_flags);
void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags,
  dvar_vector & diff_coff_region_means,
  dvar_vector & diff_coff_region_sums,
  dvariable& diff_coff_region_means_sum,
  dvariable& diff_coff_region_tot,
  ivector & parest_flags);


dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,int year);

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar_matrix& N,
  dvar3_array& Dad,int magage=0);

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int magage=0);

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int maxage,pmulti_species_data & pmsd);

  MY_DOUBLE_TYPE age_at_length_calc(MY_DOUBLE_TYPE v,const dvector& vb_coff,int nage);


imatrix group_pointers(const ivector& v);

  void normalize_effort_data(dvar_fish_stock_history& fsh);

void set_true_months(dvar_len_fish_stock_history& fsh);


void get_effort_data_by_fishery(dvar_fish_stock_history& fsh);

dvariable ratio_first_last_biocalc(dvar_len_fish_stock_history& fsh);
dvariable total_catch_or_weight_fit_ms(dvar_len_fish_stock_history& fsh,
  int print_switch,int);
dvariable total_catch_or_weight_fit(dvar_len_fish_stock_history& fsh,
  int print_switch,int);
dvariable stock_recruit_bh(dvar_len_fish_stock_history& fsh,
  ofstream *,dvar_matrix&,ivector * pq_flag=0 );

void get_msy_pt(dvariable& Bmsy,dvariable & Msy,dvariable& k,
  const dvariable& m, dvariable& r);

//class function_minimizer_exception
//{
//  public:
 //   function_minimizer_exception(void){}
//};

void length_calc(dvariable& tt,const prevariable& vb1,
    const prevariable& vb2,const prevariable& t1,const prevariable& t2);


void yield_per_recruit_analysis(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvar_matrix& Fay);



dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvariable& rho,
  const dvariable& vbdiff,const dvariable& vbtmp,const dvar_vector& vb_coff,int nage);

prevariable& daves_kludge1(const prevariable&);

/*
extern struct pvmhostinfo *ad_hostp;
extern int ad_nhost;
extern ivector * ad_stid;
*/
int size_count(const dvar_matrix& w,int flags,const ivector& group);

void set_value_inv(const dmatrix& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,
    int flags,const ivector& group);

  void set_value_inv(const dvar_matrix& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,int flags,const ivector& group);

void set_value(const dvar_matrix&,const dvar_vector&,const int&, MY_DOUBLE_TYPE,
  MY_DOUBLE_TYPE,const  dvariable&,const  imatrix&,const  imatrix&);
int num_active(const dvar_matrix&,const  dvar_matrix&, dvar_matrix&);
void set_value_inv(const dvar_matrix&,const  dvector&, int&, MY_DOUBLE_TYPE, 
  MY_DOUBLE_TYPE, dvar_matrix&, dvar_matrix&);
void set_value(const dvar_vector& x,const dvector& y,const dvar_vector& v, const int& ii);

void set_value(const dvar_vector& x,const dvector& y,const dvar_vector& v, const int& ii, 
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& fpen);

void set_value_inv(const dvar_vector& x,const dvector& y ,const dvector& v, const int& ii);

void set_value_inv(const dvar_vector& x,const dvector& y,const dvector& v, const int& ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax);


void set_value_inv(const dvar_vector& x,const dvector& y,const dvector& v, const int& ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);

void set_value(const dvar_vector& x,const dvector& y,const dvar_vector& v, const int& ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,MY_DOUBLE_TYPE s);

void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,const ivector& group);


void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags);


int num_active(const dvar_vector& w,const ivector& flags,const ivector& group);

int num_active(const dvar_vector& w,const ivector& flags);

void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group);

void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags);


 void set_value(const dvar_matrix& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const MY_DOUBLE_TYPE s,int flags,const ivector& group);

void set_value(const dvar_matrix& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group);
int num_active(const dvar_matrix& w,const ivector& flags,const ivector& group);
void set_value_inv(const dvar_matrix& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& flags,const ivector& group);

void set_value_partial(const dvar4_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE scale);

void set_value_partial(const dvar4_array& _w,const dvar_vector& x,
  const int& ii,int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,int flags,MY_DOUBLE_TYPE scale);

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE scale);

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,int flags,MY_DOUBLE_TYPE scale);
int num_active_partial(const dvar4_array& w,int range);
int num_active_partial(const dvar4_array& w,int flags,
  int range);
int num_active_partial(const dvar3_array& w,int flags,
  int range);
void set_value_inv_partial(const dvar_matrix& w,const dvector& x,const int& ii,
  int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale);
void set_value_partial(const dvar4_array& _w,const dvector& x,
  const int& ii,const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  int flags,MY_DOUBLE_TYPE scale);

void set_value_partial(const dvar4_array& _w,const dvar_vector& x,
  const int& ii,const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,MY_DOUBLE_TYPE scale);


void set_value_inv_partial(const dvar4_array& _w,const dvector& x,
  const int& ii,int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  int flags,MY_DOUBLE_TYPE scale);

 // set_value_inv_partial(kludged_equilib_coffs,
 //    x, ii, parest_flags(374),-3.0,3.0,parest_flags(373),1.0L);



void set_value_partial(const dvar_matrix& w,const dvar_vector& x,const int& ii,const ivector& range,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,const ivector& group);

void set_value_inv_partial(const dvar_matrix& w,const dvector& x,const int& ii,const ivector& range,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group);

int num_active_partial(const dvar_matrix& w,const ivector& flags,const ivector& group,const ivector& range);

void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale);

void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale);

void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,const ivector& onsw,
  const ivector& gpsw,MY_DOUBLE_TYPE scale);

void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& flags,MY_DOUBLE_TYPE scale);

void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE scale);

void set_value_inv_partial(const dmatrix& w,const dvector& x,const int& ii,const ivector& range,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale);

void set_value_inv_partial(const dvar4_array& w,const dvector& x,const int& ii,
  const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale);
void set_value_inv(const dvar_vector& _x,const dvar_vector& v,
  const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);


//void set_value(const dvar_vector& _x,const dvar_vector& v,
//  const int& _ii,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,MY_DOUBLE_TYPE s);

void put_in_indep_vars(MY_DOUBLE_TYPE idv, MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax, dvector indep_var, dvector indep_var_lo,
  dvector indep_var_hi,int& jj);
void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj);
void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj, int num);
void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj, imatrix fshgptr);
void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj, ivector fflgs, ivector fgrps);
void put_in_indep_vars(dmatrix idv, MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax, dvector indep_var, dvector indep_var_lo,
  dvector indep_var_hi,int& jj);
void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj, ivector ffgroup);
void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj,ivector n,ivector fgrps);
void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj, imatrix fflgs, imatrix fgrps);
void put_in_indep_vars(d4_array idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj);
void put_in_indep_vars(d4_array idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jj,MY_DOUBLE_TYPE n);



void  get_seasonal_catchability_pars_index(dvar_len_fish_stock_history& fsh);

int set_month1_flag(ivector& dataswitch);

d3_array total_wght_calc(dvar_len_fish_stock_history& fsh,int print_switch);

dvariable square_fit_wght(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq);

dvariable square_fit_t(dvar_len_fish_stock_history& fsh, d3_array& total_freq);
dvariable square_fit_t_wght(dvar_len_fish_stock_history& fsh, d3_array& total_wght);

#endif //define NEWMULT_HPP

dvariable catch_at_weight_fit(dvar_len_fish_stock_history& fsh,
  int print_switch);

dvariable square_fita_wght(dvar_len_fish_stock_history& fsh,
  d3_array& total_wght);

dvariable normal_length_to_weight(MY_DOUBLE_TYPE delta,MY_DOUBLE_TYPE xmin,
  const dvariable& _mean,const dvariable & _sigma,const dvariable& _a,
  const dvariable & _b);

MY_DOUBLE_TYPE get_node_date_for_region_and_period(int rr,int rp,
  const dvar_fish_stock_history& fsh);

int max_index(const dvar_vector & v);
int max_index(const dvector & v);
int min_index(const dvar_vector & v);
int min_index(const dvector & v);

void get_global_fish_periods(dvar_fish_stock_history& fsh);
void month_doubling_kludge(fishery_header_record_array& fra,
  int mult,ivector& regmin,int& first_time,pmulti_species_data & pmsd);

  class movement_info
  {
  public:
    void operator ++ (void); 
    int year;
    ivector move_weeks;
    int current;
    void allocate(const ivector & mw, int year); 
    int month(void);
    int week(void);
    int num_periods(void);
    //int (void);
    void initialize(const fishery_header_record& t1);
    void initialize_terminal(const fishery_header_record& t1);
  };

void check_sanity(ivector &t);

void check_sanity(ivector &t,dvar3_array& tg,ivector& rip,int it,
  dvar_fish_stock_history& fsh);

void read_move_weeks(ivector& move_weeks,ivector& dataswitch,
  cifstream & infile,int month_1);

int search_for_string(cifstream & in,
		      const string& s,const string& es);  //NMD_24Mar2015

dvariable boundp(const prevariable& xx, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const MY_DOUBLE_TYPE& s);

void set_value(const dvar_vector& x,_CONST dvar_vector& v, const int& _ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const MY_DOUBLE_TYPE s);
  void set_value_inv(const dvar_matrix& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,
    int flags,const ivector& group);
  void set_value_inv(const dmatrix& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    int flags,const ivector& group);
  void set_value(const dvar_matrix& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,int flags,const ivector& group);
  int num_active(const dmatrix& w,int flags,const ivector& group);
  void set_value(const dvar_matrix& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,MY_DOUBLE_TYPE s,int flags,const ivector& group);

  int num_active(const dvar_matrix& w,int flags,const ivector& group);
  dvector get_generic_mean_lengths(dvar_len_fish_stock_history& fsh);

dvar_vector fast_diffusion_calcs_for_plusgroup_biomass(int nage,
  int num_regions,dvar_matrix& B, dvar3_array& Dad,ivector& rip);

void get_taggroup_assignments(dvar_len_fish_stock_history& fsh);
void set_tag_group_assignments(dvar_len_fish_stock_history& fsh);
void mfsend_int_to_slaves(int n);
int mfget_int_from_master(void);
void setup_pvm_borland_debug(int pvm_switch,int& ad_nhost,int narch,
  struct pvmhostinfo * ad_hostp );
void setup_pvm(int pvm_switch,int& ad_nhost,int narch);
  //struct pvmhostinfo * ad_hostp );
int check_working_directory_option(int argc,char *argv[]);

 void F0biomass_calcs(adstring& filename,dvar_len_fish_stock_history& fsh,
      dvar_matrix& sel);
void send_taggroup_assignments(const ivector& min_group,
  const ivector& max_group,dvar_len_fish_stock_history& fsh);

dvariable yield_analysis_bh_daves_folly(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvariable alpha,dvariable beta,
 MY_DOUBLE_TYPE lwc,const dvar_vector& pmature);

void mfsend_bytes_to_slaves(char * buffer,int n);
void mfget_bytes_from_master(char * buffer,int n);

#if !defined(__MINGW__) && !defined(__MSVC32__) &&  (__BORLANDC__  <= 0x0550) 
struct timeb {
      long  time;
      short millitm;
      short timezone;
      short dstflag;
     };
#else
#include <sys/timeb.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#if !defined(__MINGW__) && defined(linux)	
  int  ftime(struct timeb *);
#else
  void  ftime(struct timeb *);
#endif

#ifdef __cplusplus
}
#endif

class mf_timer
{
  timeb t0;
  timeb t1;
public:
  void initialize(void);
  int report(void);
};
void print_identifier_stuff(ofstream& ofs);
void print_yield_stuff(ofstream* pof,dvar_len_fish_stock_history& fsh,
   int lastFmlt,  dvar_vector& F_agg,dvar_vector& eq_yield, dvar_vector& Beq,
    dvar_vector& sBeq, int lastFmult);

void print_yield_stuff(ofstream* pof,dvar_len_fish_stock_history& fsh,
   int lastFmlt,  dvector& F_agg,dvector& eq_yield, dvector& Beq,
    dvector& sBeq, int lastFmult);


void get_biomass_ratio_by_time(dvar_len_fish_stock_history& fsh,
  int numflag,const dvar_vector& pmature,dvar_vector& tmp2,dvar_vector& tmp3,
  const dvector& tmp4);
void get_biomass_ratio_by_time(dvar_len_fish_stock_history& fsh,
  int numflag,const dvar_vector& pmature,dvar_vector& tmp2,dvar_vector& tmp3,
  const dvar_vector& tmp4);

void get_biomass_ratio_by_time(dvar_len_fish_stock_history& fsh,
  int numflag,const dvar_vector& pmature,dvector& tmp2,dvector& tmp3,
  dvector& tmp4);
dvariable selmean_penalty(dvar_fish_stock_history& fsh,
  int print_switch);


  class nrf
  {
    int nf;
    int ng;
    dvector q;
    dvector logq;
    dmatrix s;
    dvector N;
    dmatrix F;
    dmatrix C;
    dvector Z;
    dvector S;
    dvector wtflag;
    dvector w;
    MY_DOUBLE_TYPE M;
    dvector Cobs;
    dvector Chat;
    MY_DOUBLE_TYPE T;
    MY_DOUBLE_TYPE That;
    MY_DOUBLE_TYPE R;
    MY_DOUBLE_TYPE beta;
    MY_DOUBLE_TYPE phi;
    int simflag;
  public:
    nrf(void);

    nrf(int _nf,int _ng,dmatrix _s,MY_DOUBLE_TYPE _M,dvector _Cobs,
     dvector _N,dvector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
     dvector _q) : 
     nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
     q(_q),C(1,nf,1,ng),F(1,nf,1,ng), Z(1,ng), S(1,ng),Chat(1,nf),w(_w),
     wtflag(_wtflag),beta(_beta),phi(1.0L),logq(1,nf),simflag(1) {}

    void get_initial_q(void);
    dvector calculate_d(void);
    void calculate_F(void); 
    void calculate_Z(void); 
    void calculate_S(void); 
    void calculate_T(void); 
    void calculate_That(void); 
    void calculate_Chat(void); 
    void testnr(void);
    void simulate(void);
    dmatrix calculate_J(void);
    dmatrix calculate_J_anal(void);

  };


  class vnrf
  {
  public:
    int af92;
    int icount;
    int have_weights;
  private:
    MY_DOUBLE_TYPE rmax;
    ivector fishin;
  public:
    int nf;
  private:
    int ir;
    int ip;
    int ng;
    dvar_fish_stock_history * pfsh;
    dvar_vector q;
    ivector age_flags;
    dvar_vector logq;
    dvar_matrix s;
    dvar_vector N;
    dvar_matrix sN;
    dvar_matrix swN;
    dvar_matrix F;
    dvar_matrix C;
    dvar_vector Z;
    dvar_vector omega;
    dvar_vector S;
    ivector wtflag;
    dvar_vector w;
    dvar_vector M;
    dvector Cobs;
    dvar_vector Chat;
    dvariable T;
    dvariable That;
    dvariable R;
    dvariable beta;
    dvariable phi;
    int simflag;
    MY_DOUBLE_TYPE pamin1;
    MY_DOUBLE_TYPE pamin2;
    MY_DOUBLE_TYPE pwght1;
    MY_DOUBLE_TYPE pwght2;
    dvariable fpen;
    dvariable fpen1;
    
  public:
    vnrf(void);

   vnrf(const dvar_matrix& fm,const dvar_vector& tm,const dvar_vector& surv,
     int _nf,int _ng,const dvar_matrix & _s,
     const dvar_vector& _M,const dvector& _Cobs,const dvar_vector& _N,
     const dvar_vector& _w,const ivector& _wtflag,MY_DOUBLE_TYPE _beta,
     const dvar_vector& _q, const dvar_fish_stock_history * _pfsh,int _ir,
     int _ip,MY_DOUBLE_TYPE _rmax,const ivector& age_flags);
  /*
   vnrf(dvar_matrix& fm,dvar_vector& tm,dvar_vector& surv,
     int _nf,int _ng,dvar_matrix & _s,
     dvar_vector& _M,dvector& _Cobs,dvar_vector& _N,dvar_vector& _w,
     ivector& _wtflag,MY_DOUBLE_TYPE _beta,dvar_vector& _q, 
     dvar_fish_stock_history * _pfsh,int _ir,int _ip,MY_DOUBLE_TYPE _rmax,
     ivector& age_flags);
   */

    vnrf(int _nf,int _ng,dvar_matrix _s,dvar_vector& _M,dvector _Cobs,
     dvar_vector _N,dvar_vector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
     dvar_vector _q, dvar_fish_stock_history * _pfsh,int _ir,int _ip,
     ivector _fishin);

    void get_initial_q(void);
    void get_initial_q_ss2(void);
    dvar_vector calculate_d(void);
    dvar_vector calculate_ld(void);
    void calculate_F(void); 
    virtual void calculate_Z(void); 
    void calculate_S(void); 
    //void calculate_T(void); 
    //void calculate_That(void); 
    void calculate_Chat(void); 
    dvariable testnr(void);
    dvariable testnrx(void);
    dvariable testnry(void);
    void simulate(void);
    dvar_matrix calculate_J(void);
    dvar_matrix calculate_lJ(void);
    dvar_matrix & get_F(void) { return F;}

    dmatrix calculate_J_anal(void);
    dvar_matrix calculate_J_vanal(void);
    d3_array calculate_alpha_ders(void);
    dvector calculate_N_ders(void);
    d3_array calculate_alpha(void);
    dvar_vector& get_Z(void){ return Z;}
  };

  class vnrfm : public vnrf
  {
    dvar_vector missing_z;
  public:
   vnrfm(dvar_matrix& fm,dvar_vector& tm,dvar_vector& surv,
     int _nf,int _ng,dvar_matrix & _s,
     dvar_vector& _M,dvector& _Cobs,dvar_vector& _N,dvar_vector& _w,
     ivector& _wtflag,MY_DOUBLE_TYPE _beta,dvar_vector& _q, 
     dvar_fish_stock_history * _pfsh,int _ir,int _ip,MY_DOUBLE_TYPE _rmax,
     ivector& age_flags,dvar_vector& _missing_z);
    virtual void calculate_Z(void); 
  };


  void normalize_projection_year(int month_1,int direction_flag,int min_year,
    imatrix& fdat,dvar_fish_stock_history& fsh);

void set_value(const dvar_vector& _x,const dvar_vector& v,const int& _ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& fpen,MY_DOUBLE_TYPE s,const ivector& flags);
void set_value(const dvar_vector& _x,const dvar_vector& v,const int& _ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& fpen,MY_DOUBLE_TYPE s,const ivector& flags,
  MY_DOUBLE_TYPE off_value);


void set_value_inv(const dvar_vector& x,const dvector& _v,const int& _ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,const ivector& flags);

dvariable yield_analysis_bh_penalty(dvar_len_fish_stock_history& fsh,
 dvariable alpha,dvariable beta,const dvar_vector& pmature);
 //dvariable alpha,dvariable beta,dvector& pmature,dvar_matrix& Fay);
//dvariable alpha,dvariable beta,MY_DOUBLE_TYPE lwc,dvector& pmature,dvar_matrix& Fay);

dvariable yield_analysis_bh(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvariable alpha,dvariable beta,
 /*double lwc,*/const dvar_vector& pmature);

dvariable stock_recruit_bh_steep(dvar_len_fish_stock_history& fsh,
 ofstream * pof,ivector * pq_flag);

dvariable calculate_the_biomass(int i,dvar_fish_stock_history& fsh,
  dvar3_array & N,  int numflag,const dvar_vector& pmature);

void get_bh_alpha_and_beta(dvar_len_fish_stock_history& fsh,
  dvariable& alpha,dvariable& beta);

dvariable get_bh_recruitment_for_projections(int i,
  dvar_len_fish_stock_history& fsh,dvariable& alpha,dvariable& beta);

dmatrix orthpoly_constant_begin_end(int n,int deg,int begin_degree,
  int nconst_begin,int end_degree,int nconst_end);

imatrix read_datafile_flaginfo(cifstream& infile,int num_fish,
  int data_file_version_no,imatrix &,pmulti_species_data& pmsd);

dvar_vector get_average_recruitment_for_projections(dvar3_array& N,
  int num_regions,int last_real_year,int imin,int imax);

//dvariable get_average_log_recruitment_for_projections(dvar3_array& N,
//  int num_regions, int last_real_year, int imin, int imax);   //NMD13Dec2012

dvariable get_average_log_recruitment_for_projections(dvar3_array& N,
  int num_regions, int last_real_year, int imin, int imax, int af57,
  int pf232, int pf233, int af182);   //NMD_1Sep2015

dvar_vector get_average_recruitment(dvar3_array& N,int num_regions,
  int last_real_year);

void  get_implicit_catchability(dvar_fish_stock_history& fsh);
void  get_implicit_catchability_catch_conditioned(dvar_fish_stock_history& fsh);
dvar_vector vget_generic_mean_lengths(dvar_len_fish_stock_history& fsh);

#if defined(USE_PTHREADS)
  void pthreads_master_send_signal_to_slave(void);
  extern pthread_key_t admb_pthread_key;
#endif

  void my_do_nothing(void*);

  class vnrfss2
  {
  public:
    int icount;
    dvar_vector& get_N(void) { return N; }
    dvar_vector& get_Nafter(void) { return Nafter; }
  private:
    const MY_DOUBLE_TYPE rmax;
    ivector fishin;
    int nf;
    int ir;
    int ip;
    int ng;
    dvar_len_fish_stock_history * pfsh;
    dvar_vector q;
    dvar_vector logq;
    dvar_matrix s;
    dvar_vector N;
    dvar_vector Nbefore;
    dvar_vector Nafter;
    dvar_matrix sN;
    dvar_matrix swN;
    dvar_matrix F;
    dvar_matrix Fraw;
    dvar_matrix C;
    dvar_vector Z;
    dvar_vector Zraw;
    dvar_vector S;
    dvar_vector Sraw;
    ivector wtflag;
    dvar_vector w;
    dvar_vector M;
    dvector Cobs;
    dvar_vector Chat;
    dvariable T;
    dvariable That;
    dvariable R;
    dvariable beta;
    dvariable phi;
    int simflag;
  public:
    vnrfss2(void);

   vnrfss2(dvar_matrix& _fm,int _nf,int _ng,dvar_matrix & _s,dvar_vector& _M,dvector& _Cobs,
     dvar_vector& _N,dvar_vector& _w,ivector& _wtflag,MY_DOUBLE_TYPE _beta,
     dvar_vector& _q, dvar_len_fish_stock_history * _pfsh,int _ir,int _ip,
     MY_DOUBLE_TYPE _rmax);
  /*
     : nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
     q(_q),C(1,nf,1,ng),F(1,nf,1,ng), Z(1,ng), S(1,ng),Chat(1,nf),w(_w),
     wtflag(_wtflag),beta(_beta),phi(1.0L),logq(1,nf),simflag(1),
     pfsh(_pfsh),ir(_ir),ip(_ip) {}
   */

    vnrfss2(int _nf,int _ng,dvar_matrix _s,dvar_vector& _M,dvector _Cobs,
     dvar_vector _N,dvar_vector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
     dvar_vector _q, dvar_len_fish_stock_history * _pfsh,int _ir,int _ip,
     ivector _fishin);

  /*
    vnrfss2(int _nf,int _ng,dvar_matrix _s,dvariable _M,dvector _Cobs,
     dvar_vector _N,dvar_vector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
     dvar_vector _q, dvar_len_fish_stock_history * _pfsh,int _ir,int _ip,
     ivector _fishin) : 
     nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
     q(_q),C(1,nf,1,ng),F(1,nf,1,ng), Z(1,ng), S(1,ng),Chat(1,nf),w(_w),
     wtflag(_wtflag),beta(_beta),phi(1.0L),logq(1,nf),simflag(1),
     pfsh(_pfsh),ir(_ir),ip(_ip),sN(1,_nf),swN(1,_nf),fishin(_fishin) 
    {}
   */

    void get_initial_q(void);
    dvar_vector calculate_d(void);
    dvar_vector calculate_ld(void);
    void calculate_F(void); 
    void calculate_Fraw(void); 
    void calculate_Z(void); 
    void calculate_Zraw(void); 
    dvariable calculate_S(void); 
    void calculate_T(void); 
    void calculate_That(void); 
    void calculate_Chat(void); 
    dvariable testnr(void);
    void simulate(void);
    dvar_matrix calculate_J(void);
    dvar_matrix calculate_lJ(void);

    dmatrix calculate_J_anal(void);
    dvar_matrix calculate_J_vanal(void);
    d3_array calculate_alpha_ders(void);
    dvector calculate_N_ders(void);
    d3_array calculate_alpha(void);
  };


#if (!defined(__MSVC32__) &&  !defined(__BORLANDC__) ) 
# define _export
#endif
//_export

//extern "C" _export  void dd_newton_raphson(int n,MY_DOUBLE_TYPE * v,MY_DOUBLE_TYPE * diag,
 //   MY_DOUBLE_TYPE * ldiag, MY_DOUBLE_TYPE * yy);
  void set_value_inv_exp(const dvar_vector& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax, const ivector& flags,
    const ivector& group,MY_DOUBLE_TYPE s);

  void set_value_exp(const dvar_vector& _w,const dvar_vector& x,const int& _ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,
    const ivector& group,MY_DOUBLE_TYPE s);

  void set_value_inv_exp(const prevariable& x,const dvector& _v,const int& _ii,
    MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);

  void set_value_exp(const prevariable& _x,const dvar_vector& v,const int& _ii, 
    MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& fpen,MY_DOUBLE_TYPE s);

#if defined(USE_DD) || defined(USE_DD_NOT)
//#include <ddfvar.hpp>
#if !defined(NO_MY_DOUBLE_TYPE) 
  class ddnrf
  {
  public:
    int icount;
    int have_weights;
  private:
    const MY_DOUBLE_TYPE rmax;
    ivector fishin;
  public:
    int nf;
  private:
    int ir;
    int ip;
    int ng;
    dvar_fish_stock_history * pfsh;
    dvector q;
    dvector logq;
    dmatrix s;
    dvector N;
    dmatrix sN;
    dmatrix swN;
    dmatrix F;

    d3_array dfF;
    d3_array dfsN;
    d3_array dfswN;
    d3_array dfs;
    d3_array dfC;
    dmatrix dfS;
    dmatrix dfq;
    dmatrix dfZ;
    dmatrix dfChat;
    dmatrix dfw;
    dmatrix dfM;
    dmatrix dfN;

    d3_array dFdq;
    d3_array dsdq;
    dmatrix dSdq;
    dmatrix dqdq;
    dmatrix dZdq;
    dmatrix dwdq;
    dmatrix dMdq;
    dmatrix dNdq;

    dmatrix C;
    dvector Z;
    dvector S;
    ivector wtflag;
    dvector w;
    dvector M;
    dvector Cobs;
    dvector Chat;
    MY_DOUBLE_TYPE T;
    MY_DOUBLE_TYPE That;
    MY_DOUBLE_TYPE R;
    MY_DOUBLE_TYPE beta;
    MY_DOUBLE_TYPE phi;
    int simflag;
  public:
    ddnrf(void);
    dvector & get_q(void) { return q;}
    dvector & get_S(void) { return S;}
    d3_array & get_dFdq(void) { return dFdq;}
    d3_array & get_dsdq(void) { return dsdq;}
    dmatrix & get_dNdq(void) { return dNdq;}
    dmatrix & get_dZdq(void) { return dZdq;}
    dmatrix & get_dMdq(void) { return dMdq;}
    dmatrix & get_dwdq(void) { return dwdq;}

   ddnrf(int _nf,int _ng,const dmatrix & _s,
     const dvector& _M,const dvector& _N,dvector& _Cobs,const dvector& _w,
     ivector& _wtflag,MY_DOUBLE_TYPE _beta,
     dvar_fish_stock_history * _pfsh,int _ir,int _ip,MY_DOUBLE_TYPE _rmax);

    ddnrf(int _nf,int _ng,dmatrix _s,dvector& _M,dvector _Cobs,
     dvector _N,dvector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
     dvector _q, dvar_fish_stock_history * _pfsh,int _ir,int _ip,
     ivector _fishin,MY_DOUBLE_TYPE _rmax);

    void get_initial_q(void);
    dvector calculate_d(void);
    dvector calculate_ld(void);
    void calculate_sN(void); 
    void calculate_F(void); 
    void calculate_Z(void); 
    void calculate_S(void); 
    //void calculate_T(void); 
    //void calculate_That(void); 
    void calculate_Chat(void); 
    MY_DOUBLE_TYPE testnr(void);
    void simulate(void);
    dmatrix calculate_J(void);
    dmatrix calculate_lJ(void);
    void df_calculate_sN(void); 
    void df_calculate_F(void); 
    void df_calculate_Z(void); 
    void df_calculate_S(void); 
    void df_calculate_Chat(void); 
    void calculate_derivatives();
    void implict_derivative_calculations(void);
    void newton_raphson_loop(void);
    void get_initial_q_ss2(void);

    dmatrix calculate_J_anal(void);
    dmatrix calculate_J_vanal(void);
    d3_array calculate_alpha_ders(void);
    dvector calculate_N_ders(void);
    d3_array calculate_alpha(void);
  };
  void test_fpu(void);
  dvector mfget_dvector_from_slave(int i);
  dvector solve(const banded_symmetric_dmatrix& _S,const dmatrix& u,
    const dmatrix& v,const dvector& x);

#else
#include <ddfvar.hpp>
  class ddnrf
  {
  public:
    int icount;
    int have_weights;
  private:
    const dd_real rmax;
    ivector fishin;
  public:
    int nf;
  private:
    int ir;
    int ip;
    int ng;
    dvar_fish_stock_history * pfsh;
    ddvector q;
    ddvector logq;
    ddmatrix s;
    ddvector N;
    ddmatrix sN;
    ddmatrix swN;
    ddmatrix F;

    dd3_array dfF;
    dd3_array dfsN;
    dd3_array dfswN;
    dd3_array dfs;
    dd3_array dfC;
    ddmatrix dfS;
    ddmatrix dfq;
    ddmatrix dfZ;
    ddmatrix dfChat;
    ddmatrix dfw;
    ddmatrix dfM;
    ddmatrix dfN;

    dd3_array dFdq;
    dd3_array dsdq;
    ddmatrix dSdq;
    ddmatrix dqdq;
    ddmatrix dZdq;
    ddmatrix dwdq;
    ddmatrix dMdq;
    ddmatrix dNdq;

    ddmatrix C;
    ddvector Z;
    ddvector S;
    ivector wtflag;
    ddvector w;
    ddvector M;
    ddvector Cobs;
    ddvector Chat;
    dd_real T;
    dd_real That;
    dd_real R;
    dd_real beta;
    dd_real phi;
    int simflag;
  public:
    ddnrf(void);
    ddvector & get_q(void) { return q;}
    ddvector & get_S(void) { return S;}
    dd3_array & get_dFdq(void) { return dFdq;}
    dd3_array & get_dsdq(void) { return dsdq;}
    ddmatrix & get_dNdq(void) { return dNdq;}
    ddmatrix & get_dZdq(void) { return dZdq;}
    ddmatrix & get_dMdq(void) { return dMdq;}
    ddmatrix & get_dwdq(void) { return dwdq;}

   ddnrf(int _nf,int _ng,const ddmatrix & _s,
     const ddvector& _M,const ddvector& _N,dvector& _Cobs,const ddvector& _w,
     ivector& _wtflag,dd_real _beta,
     dvar_fish_stock_history * _pfsh,int _ir,int _ip,dd_real _rmax);

    ddnrf(int _nf,int _ng,ddmatrix _s,ddvector& _M,dvector _Cobs,
     ddvector _N,ddvector _w,ivector _wtflag,dd_real _beta,
     ddvector _q, dvar_fish_stock_history * _pfsh,int _ir,int _ip,
     ivector _fishin,double _rmax);

    void get_initial_q(void);
    ddvector calculate_d(void);
    ddvector calculate_ld(void);
    void calculate_sN(void); 
    void calculate_F(void); 
    void calculate_Z(void); 
    void calculate_S(void); 
    //void calculate_T(void); 
    //void calculate_That(void); 
    void calculate_Chat(void); 
    dd_real testnr(void);
    void simulate(void);
    ddmatrix calculate_J(void);
    ddmatrix calculate_lJ(void);
    void df_calculate_sN(void); 
    void df_calculate_F(void); 
    void df_calculate_Z(void); 
    void df_calculate_S(void); 
    void df_calculate_Chat(void); 
    void calculate_derivatives();
    void implict_derivative_calculations(void);
    void newton_raphson_loop(void);
    void get_initial_q_ss2(void);

    ddmatrix calculate_J_anal(void);
    ddmatrix calculate_J_vanal(void);
    dd3_array calculate_alpha_ders(void);
    ddvector calculate_N_ders(void);
    dd3_array calculate_alpha(void);
  };
  void test_fpu(void);
  dvector mfget_dvector_from_slave(int i);
  ddvector solve(const banded_symmetric_ddmatrix& _S,const ddmatrix& u,
    const ddmatrix& v,const ddvector& x);
#endif  // #if !defined(NO_MY_DOUBLE_TYPE) 

#endif //#if defined(USE_DDOUBLE)

dvariable biomass_target_calcs(dvar_len_fish_stock_history& fsh,
   dvariable & Bmsy, dvariable & sBmsy,const dvar_vector& p);


  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvariable& rho,
    const dvariable& vbdiff,const dvariable& vbtmp,const dvar_vector& vb_coff,int nage,const ivector& pf);

  MY_DOUBLE_TYPE get_flag(int i,MY_DOUBLE_TYPE deflt);
//double get_fmsy_pen_wt(int i);
int check(const ivector& v,int n);
#include "makebig2.hpp"
void set_value_exp(const prevariable& _x,const dvar_vector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& _fpen,MY_DOUBLE_TYPE s);
void set_value_exp(const prevariable& _x,const dvar_vector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& _fpen);
dvar_vector rowstack(const dvar_matrix& M);
void rowunstack(const dvar_vector& _w,const dvar_matrix& _M);
dvector rowstack(const dmatrix& M);
ivector rowstack(const imatrix& M);

void set_value_inv_exp(const dvar_matrix& _mw,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax, const imatrix& mflags,const imatrix& mgroup,
    MY_DOUBLE_TYPE s);
void set_value_inv(const dvar_matrix& _mw,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const imatrix& mflags,const imatrix& mgroup);
int num_active(const dvar_matrix& _mw,const imatrix& mflags,const imatrix& mgroup);
void set_value_exp(const dvar_matrix& _mw,const dvar_vector& x,
    const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,
    const imatrix& mflags,const imatrix& mgroup);
void set_value_exp(const dvar_matrix& _mw,const dvar_vector& x,
      const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,
      const imatrix& mflags,const imatrix& mgroup,MY_DOUBLE_TYPE s);
  void read_selectivity_blocks(dvar_len_fish_stock_history& fsh,
    ivector& ff71);
dvar_matrix colsub(const dvar_matrix& _M,int lb,int ub);
dmatrix colsub(const dmatrix& _M,int lb,int ub);
void read_age_length_records(cifstream & cif,int nlint,int nage,
  dvar_len_fish_stock_history& fsh);
dvar_vector choleski_solve(_CONST dvar_matrix& MM,const dvar_vector& vv,
   const prevariable& det,const int& sgn);

  void set_value(const dvar3_array& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax, const dvariable& pen,MY_DOUBLE_TYPE s,ivector& flags,const ivector& group);
  int num_active(const dvar3_array& w,ivector& flags,const ivector& group);
  void set_value_inv(const dvar3_array& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group);
  int size_count(const dvar3_array& w,ivector&  flags,const ivector& group);
void set_value(const dvar4_array& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,int flags,const ivector& group);
  void set_value(const dvar4_array& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group);
  void set_value_inv(const dvar4_array& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group);
  int size_count(const dvar4_array& w,ivector&  flags,const ivector& group);
void write_bs_selcoff(par_ofstream& pof,dvar4_array& bsc,int nage);
void read_bs_selcoff(cifstream& pof,dvar4_array& bsc,int nage);
void read_bs_selcoff(cifstream& pof,d4_array& bsc,int nage);
dvariable bs_selmean_penalty(dvar_fish_stock_history& fsh,
  int print_switch);
void logistic_sel_calc(dvar_len_fish_stock_history& fsh);
dvariable bs_sel_pen(dvar_fish_stock_history& fsh);
dvariable normalized_bs_sel_pen(dvar_fish_stock_history& fsh);
void grouped_flag_sanity_check(const ivector gp,const ivector& flags);
dvariable lognormal_multinomial_log_likelihood
   (const dvariable& vN,const dvector& q,const dvar_vector& vp,
   const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix& S);
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvariable&  _vrho,const dvariable&  _vvar,const dvector& _eps,int ic,
  ccubic_spline_array * pccsa);


dvariable old_wght_self_scaling_multinomial_nore
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch);
dvariable new_wght_self_scaling_multinomial_nore
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch);
dvariable new_len_self_scaling_multinomial_nore
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch);
dvariable old_len_self_scaling_multinomial_nore
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch);
dvariable wt_self_scaling_multinomial_re(dvar_len_fish_stock_history& fsh);
dvariable len_self_scaling_multinomial_re(dvar_len_fish_stock_history& fsh);
dvariable wght_self_scaling_multinomial_nore(dvar_len_fish_stock_history& fsh,
  d3_array & wtl,int print_switch);
dvariable len_self_scaling_multinomial_nore(dvar_len_fish_stock_history& fsh,
  d3_array & tl,int print_switch);
dvariable wght_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array & tw,int print_switch);
dvariable xlen_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array & tf,int print_switch);
dvariable len_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array & tf,int print_switch);
dvariable new_len_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array & tf,int print_switch);
dvariable new_wght_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array & tf,int print_switch);

 void error_msg2(const char *s);
dvar_matrix get_initial_equilibrium_plus_group(int nseas,int nreg,int nage);

void set_value_inv(const dvar_vector& _x,const dvar_vector& v,int flag,
  const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);
void set_value(const dvar_vector& _x,const dvar_vector& v,int flag,
  const int& _ii,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,MY_DOUBLE_TYPE s);
#if defined(AD_FPEXCEPTION)
#include<float.h>
class ad_float_exception
{
public:
  int errtype;
  ad_float_exception(int i) : errtype(i){}
};

extern "C" void adfloat_except(int k);
/*
{
  throw ad_float_exception(1);
}
*/
#include <exception>
//using namespace std;

class myexception: public std::exception
{
  virtual const char* what() const throw()
  {
    return "My exception happened";
  }
};
#include "tridiagonal_dmatrix.h"
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix & S,
  const banded_lower_triangular_dmatrix & _chinv,
  const symmetric_tridiagonal_dmatrix& lensinv);

dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  //const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix & S,
  const dvariable& vrho,const dvariable& vvar,
   const dvector& eps,int ic=0);

dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvariable& vrho,const dvariable& vvar,int ic=0);

dvar_matrix get_average_recruitment_by_season(dvar3_array& N,
  int num_regions,int last_real_year,int af57);

dvar_matrix get_average_recruitment_by_season_for_region
  (dvar_matrix& average_recruitment_by_season,int num_regions,int af57);

dvar_matrix get_scalar_recruitment_by_season_for_region
  (dvar_matrix& average_recruitment_by_season_for_regions,
   int num_regions,int af57);

dmatrix  get_learner_code_coffs_low();

dvariable len_dm_nore(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch);
dvariable wght_dm_nore(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch);
dvariable len_dm_nore_notc(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch);
dvariable wght_dm_nore_notc(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch);

void len_tail_compress(dvar_len_fish_stock_history& fsh);
void wght_tail_compress(dvar_len_fish_stock_history& fsh);

dvar_vector chinv_multiply(const prevariable rho,const dvar_vector& v);
dvar_vector vbiocalc(dvar_len_fish_stock_history& fsh,int pqflag=0);
dvar_vector adult_vbiocalc(dvar_len_fish_stock_history& fsh,int pqflag=0);
dvar_matrix reg_vbiocalc(dvar_len_fish_stock_history& fsh,int q0flag=0);
dvar_matrix unnormalized_reg_vbiocalc(dvar_len_fish_stock_history& fsh,
  int q0flag=0);
dvar_matrix unnormalized_adult_reg_vbiocalc(dvar_len_fish_stock_history& fsh,int q0flag=0);

dvariable xlognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvariable&  _vrho,const dvariable&  _vvar,const dvector& eps,int ic);
void set_option_flag(const char * s,int& ss,int argc,char * argv[]);
void normalize_tag_dates(ivector & tag_year, ivector & tag_month, ivector & true_tag_year, int month_1, int year1, int direction_flag, int month_factor, int first_time,int it);

void sanity_check_2(ivector& initial_tag_year,imatrix& year,
  imatrix& initial_tag_period);
dmatrix make_dmatrix(const banded_lower_triangular_dmatrix& B);
  dvector operator * (const dvector &,const banded_lower_triangular_dmatrix& B);
banded_lower_triangular_dmatrix get_ltdchinv(MY_DOUBLE_TYPE rho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax);
banded_lower_triangular_dmatrix get_ltdchinv(MY_DOUBLE_TYPE rho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax);


  void check_number(cifstream& cif,int n,ivector pf);
  imatrix months_analyzer(imatrix& month_detector);
ivector setup_group_selectivity_pointer(ivector& months_used,int month_1,
  int ff74);
  void set_value_inv(const dvar_matrix & w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,
    const ivector& group);
  void set_value_inv(const dvar_matrix & w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group);
  int num_active(const dvar_matrix& w,ivector& flags,imatrix& mflags,
    const ivector& group);
  void set_value(const dvar_matrix & w,const dvar_vector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,
    MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group);
void set_value_inv(const dvar_vector& x,const dvector& _v,const int& _ii,
  ivector mflags,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);
  void zero_effdev_sanity_check(dvar_matrix& w,ivector& flags,imatrix& mflags,
    const ivector& group);

dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv,
  const banded_symmetric_dmatrix& _stdsinv,const dvector& _ee);
MY_DOUBLE_TYPE lognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
  const dvector& q,const dvector& p,
  const dvector& eta,const banded_symmetric_dmatrix& stdsinv);
dvariable lognormal_multinomial_log_likelihood_cubic_spline
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvector& _eps,int ic,
  banded_symmetric_dvar_matrix& vbsd2,const dvector& eta_hat);
dvector xlt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w);
void xdflt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w,
  const dvector& _dfw,const dvector& _dfx,
  const banded_lower_triangular_dmatrix & _dfL);
dvector get_cross_derivatives(const dmatrix& M,const dvector& q,
  const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
  const dvector dfeta,const dvector& p);
dvector get_cross_derivatives(
  const banded_lower_triangular_dmatrix L1,
  const banded_symmetric_dmatrix& BH,
  const dvector& w,
  const dvector& q,
  const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
  const dvector dfeta,const dvector& p);
dmatrix xlt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& ww,
  const dmatrix& N);
dmatrix inv(const banded_lower_triangular_dmatrix & L,
  const dvector& ww);
dvariable len_multiplier_calcs (dvar_len_fish_stock_history& fsh,
  int ir,int ip,int fi,int pp,int gp,int sz,const dvector& average_sample_size,
  const dvector& OP);
dvariable wght_multiplier_calcs (dvar_len_fish_stock_history& fsh,
  int ir,int ip,int fi,int pp,int gp,int sz,const dvector& average_sample_size,
  const dvector& OP);

void  add_test_par(const prevariable & _y,const prevariable& _x);
void  add_test_par1(const prevariable & _y,const prevariable& _x);
void  check_test_par(const dvar_vector & _y);
void  check_test_par(const dvar_vector & _y);
void  check_test_par(const prevariable & _y);
ivector get_inv_group_ptr(ivector v);
void print_correlation_matrices(const banded_symmetric_dvar_matrix& vbsd2,
  ofstream & ofcv,int parent,const prevariable& rho);
dvar_vector age_at_length_calc(int nlint,MY_DOUBLE_TYPE shlen,MY_DOUBLE_TYPE filen,
  const dvar_vector& vb_coff,int nage,const ivector& pf,int mult);
dvar_vector maturity_length_to_age(dvar_vector& alc,
  dvector& cpmature_at_length,int nage,int nlint,int mult);

void get_yield_at_effort_daves_folly(dvar_len_fish_stock_history& fsh,
  prevariable& lambda,prevariable & ABn,prevariable & TBn,prevariable & Cn,
  dvar_vector& Cn_by_region,
  const dvar_vector& pmature,prevariable& alpha,prevariable& beta);
void get_yield_at_effort_daves_folly(dvar_len_fish_stock_history& fsh,
  MY_DOUBLE_TYPE lambda,prevariable & ABn,prevariable & TBn,prevariable & Cn,
  dvar_vector& Cn_by_region,
  dvar_vector& ABn_by_region,
  dvar_vector& TBn_by_region,
  const dvar_vector& pmature,prevariable& alpha,prevariable& beta);

dvariable get_yield_at_multiplier2(const dvariable& lambda, dvar_vector& F,
  int nyears,dvar_vector& n1,const dvar_vector& pmature,dvariable& alpha,
  dvariable& beta,dvar_len_fish_stock_history& fsh,
  dvariable& Bmsy,dvariable& sBmsy);

dvariable get_yield_at_multiplier2(const dvariable& lambda, 
  dvar_vector& _F,int nyears,
  dvar_vector& n1,const dvar_vector& _pmature,dvariable& alpha,
  dvariable& beta,
  dvar_len_fish_stock_history& fsh,dvariable& Bmsy,dvariable& sBmsy,
  dvar_vector& _FS,const dvar_vector& _psmature);

dvariable smax1(const dvar_vector & y);
void print_movement_stuff(ofstream& xpofs, ivector& tmp_mp, ivector& tmp_yr,
  ivector& tmp_mn,ivector& tmp_ip,int num_regions);
dvector setm11(const dvector & v);
dvar_vector setm11(const dvar_vector & v);

 int num_active(const dvar_vector& w,int flags);
 int num_active(const dvar_matrix& w,int flags);

void set_value_inv(const dvar_vector& _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,int flags,MY_DOUBLE_TYPE& scale);


 void correlation_report(adstring& s,dvar_fish_stock_history& fsh);

void set_value(const prevariable& _x,const dvar_vector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& _fpen,MY_DOUBLE_TYPE s);

void set_value_inv(const dvar_matrix& _x,const dvector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,int flag,MY_DOUBLE_TYPE s);

void set_value_inv(const prevariable& _x,const dvector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,int flag,MY_DOUBLE_TYPE s);

void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags,
  dvar_matrix & diff_coff_region_means,
  dvar_matrix & diff_coff_region_sums,
  dvar_vector& diff_coff_region_means_sum,
  dvar_vector& diff_coff_region_tot,
  ivector & parest_flags,pmulti_species_data & pmsd);

void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags,
  dvar_matrix & diff_coff_region_means,
  dvar_matrix & diff_coff_region_sums,
  dvar_vector& diff_coff_region_means_sum,
  dvar_vector& diff_coff_region_tot,
  ivector & parest_flags);
void  set_value(dvar5_array & v,MY_DOUBLE_TYPE x);
dmatrix get_orthogonal_diffusion_matrix(int n);
void lapack_symmetric_eigen(int n,dvector & vM,dvector& a,
  dvector& eigenvalues,int & ierr);
dmatrix make_dmatrix(int,int,int,int,dvector &);
dvariable my_gamma_density(const prevariable&  x,const prevariable & mu,
  const prevariable & tau);
dvariable xcumd_gamma(const prevariable&  x,const prevariable&  a);
dvariable xgammln(const prevariable& xx);
df1_two_variable xgammln(const df1_two_variable& xx);

//dvar_vector censored_gamma(const prevariable & tau,dvar_vector& o,
//  const dvar_vector& mu,MY_DOUBLE_TYPE e);
  int num_active(const dvar_matrix& w,ivector& flags,imatrix& mflags,
    const ivector& group,std::shared_ptr<group_manager_1>& pgroup_manager_1,
    const ivector& old_flags);

  //int num_active(const dvar_matrix& w,ivector& flags,imatrix& mflags,
  //  const ivector& group,group_manager_1 * pgroup_manager_1,
  //  const ivector& old_flags);

  void set_value(const dvar_matrix & _w,const dvar_vector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,
    MY_DOUBLE_TYPE s,std::shared_ptr<group_manager_1>& pgroup_manager_1);

  /*
  void set_value(const dvar_matrix & _w,const dvar_vector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,
    MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group,std::shared_ptr<group_manager_1>& pgroup_manager_1);
  */

  void set_value(const dvar_matrix & _w,const dvar_vector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,
    MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group,group_manager_1 * pgroup_manager_1);

  void set_value_inv(const dvar_matrix & _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group,std::shared_ptr<group_manager_1>& pgroup_manager_1);

  /*
  void set_value_inv(const dvar_matrix & _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group,group_manager_1 * pgroup_manager_1);
  */

  void set_value_inv(const dvar_matrix & _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,
    std::shared_ptr<group_manager_1>& pgroup_manager_1);

  void gram_schmidt_qr(dmatrix& M,dmatrix& Q,dmatrix& R);

  void check_group_for_holes(const ivector& v);
  dvariable robust_normal_cauchy_mixture(const prevariable& a,
    const prevariable sighat,dvar_vector & r2);
  dvariable robust_normal_mixture(const prevariable& a,
    const prevariable sighat,dvar_vector & r2);

  class group_flag_manager
  {
    int numgroups;
    imatrix group_members;
    ivector num_in_group;
  public:
    int& get_numgroups(void) { return numgroups;}
    imatrix& get_group_members(void) { return group_members;}
    ivector& get_num_in_group(void) { return num_in_group;}
    group_flag_manager(ivector& group);
  };

  void mytestcin(void);

  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const dvar_vector& M );

  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const dvar_vector& M, int fi);

  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const dvar_matrix& M );

  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const dvar_matrix& M ,ivector& jmax);

  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const imatrix& mflags,const imatrix& mgroup,const dvar_matrix& M );

  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const dvar3_array& M );
  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const dvar4_array& M );
  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    int vmax,const dvar4_array& M );
  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
    const imatrix& fshgptr, const dvar4_array& M );
  dvar_vector lineup(dvar_matrix& m);
  dvar_matrix stack(dvar_vector& v,int nr);  // break up vector into nc pieces
  adstring get_string_marker_30(void);
  void insert_identifier_string_30(const char * _s);
  void set_grad_profiler_30(const char * s);
  int save_identifier_string_dbg_and_opt(char* str);
  void xinit_message(ofstream& xof,int iisave,int ii,const char * msg);
  int test_orthpoly(void);


  class  blocked_orthogonal_design
  {
    int nyears;
    int year_degree;
    int row_degree;
    int col_degree;
    int row_col_degree;
    int max_degree;
    int nrows;
    int ncols;
    int nc;
    int nr;
    int nrc;
    dmatrix Xnodes;
    int ngrand;    // this is always equal to 1 ?
    int nrow_effect; // this is nrows-1 or less depending on row_column flags
    int ncol_effect;// this is ncols-1 or less depending on row_column flags
    int nrowcol_effect; // this is the rest
    dvar_matrix orthpoly;
    d3_array Blocks;
    imatrix row_column_flags;
    void calculate_xnodes(void);
    void calculate_Blocks(void);
  public:
    int get_nr(void)
    {
      return nr;
    }
  
    int get_nc(void)
    {
      return nc;
    }
  
    int get_nrc(void)
    {
      return nrc;
    }
  
    int get_year_degree(void)
    {
      return year_degree;
    }
  
    int get_row_degree(void)
    {
      return row_degree;
    }
  
    int get_col_degree(void)
    {
      return col_degree;
    }
  
    int get_row_col_degree(void)
    {
      return row_col_degree;
    }
  
    dmatrix get_Xnodes(void)
    {
      return Xnodes;
    }
    d3_array get_Blocks(void)
    {
      return Blocks;
    }
    blocked_orthogonal_design(int _nyears,int _nrows,int _ncols,
      imatrix _row_column_flags,int year_degree,int row_degree,
      int col_degree,int row_col_degree);
  };

class  hemholtz_contrast    // Hemholtz contrast, I think;
{

protected:
  dmatrix C;

public:

  dmatrix & get_C(void) { return C; }

  dvector operator * (dvector& v)
  {
    return get_C() * v;
  }

  dvar_vector operator * (dvar_vector& v)
  {
    return get_C() * v;
  }
  hemholtz_contrast(int n);  // Hemholtz contrast, I think;
  hemholtz_contrast(void){}  // Hemholtz contrast, I think;
  void allocate(int n);
  int nvcal(void);
  void xinit(dvector& x,int & ii);
  void reset(dvar_vector& x,int & ii);
};

class  df_orthogonal_contrast    
{

protected:
  dmatrix C;
  dvar_vector v;

public:

  dmatrix & get_C(void) { return C; }
  dvar_vector & get_v(void) { return v; }


  dvector operator * (dvector& w)
  {
    return get_C() * w;
  }

  dvar_vector operator * (dvar_vector& v)
  {
    return get_C() * v;
  }
  df_orthogonal_contrast(int n);  
  df_orthogonal_contrast(void){}  
  void allocate(int n);
  int nvcal(void);
  void xinit(dvector& x,int & ii);
  virtual void reset(dvar_vector& x,int & ii);
};
ostream& operator << (ostream& s,const df_orthogonal_contrast& H);
cifstream& operator >> (cifstream& s,const df_orthogonal_contrast& H);

ostream& operator << (ostream& s,const hemholtz_contrast& H);
void hem_test(int n);

class new_group_manager  
{
  int ngroups;   // number of groups things are grouped into
  ivector ng;     // the number in each group
  imatrix groups;  // shape 1,ngroups,1,ng(i)
  ivector grouped_vector_ptr; // picks an ungrouped vector 
                             // value for each group

  
public:
  imatrix&  get_groups(void) { return groups; }
  ivector& get_grouped_vector_ptr(void) { return grouped_vector_ptr;}
  ivector& get_ng(void) { return ng;}
  int get_ngroups(void) { return ngroups; }
  new_group_manager(void): ngroups(0) { }
  int allocated(void) 
  { 
    if (ngroups>0) 
      return 1;
    else
      return 0;
  }

  new_group_manager(ivector& v)
  {
    allocate(v);
  }
  void allocate(ivector& v);
};

 
class  grouped_df_orthogonal_contrast : public df_orthogonal_contrast 
{
  new_group_manager  gmgr;
protected:
  dvar_vector ungrouped_v;   // vector of all ungrouped parameters
  //df_orthogonal_contrast C;
public:
  grouped_df_orthogonal_contrast(void); 
  grouped_df_orthogonal_contrast(int n); 
  void allocate();
  grouped_df_orthogonal_contrast(ivector & v); 
  dvar_vector & get_ungrouped_v(void) { return ungrouped_v; }
  int get_ngroups(void) { return gmgr.get_ngroups(); }
  void xinit(dvector& x,int & ii,ivector& grouper);
  int allocated(void) { return gmgr.allocated(); }
  virtual void reset(dvar_vector& x,int & ii);
//  int grouped_df_orthogonal_contrast::nvcal(ivector& groups);
  int nvcal(ivector& groups);
};

void grouped_flag_sanity_check(const ivector flags,const ivector& gp,int index,int group_index);
void grouped_flag_sanity_check(const ivector flags,const ivector& gp,const ivector& ovrd,
       int index,int group_index,int override_index);
ivector gram_schmidt_remove_extra_columns(dmatrix& M);

dmatrix lapack_luinv(const dmatrix& M);
dmatrix lapack_matrix_multiplcation(const dmatrix& A,const dmatrix& B);
dvar_vector neg_log_student_density(const prevariable&  v,const prevariable& s,
  const dvar_vector& x);

dvariable betaln(double a,const prevariable& b );

dvariable beta(const double a,const prevariable& b );

dvar_vector betaln(const double a,const dvar_vector& b );

#endif
