/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"

#ifdef __BCPLUSPLUS__
  #include <iomanip.h>
#endif
#ifdef __ZTC__
  #include <iomanip.hpp>
  #include <fstream.hpp>
  #include <strstrea.hpp>
#endif
#ifdef __NDPX__
  #include <fstream.h>
  #include <sstream.h>
  extern "C" {
  #include <dos.h>
  }
#endif

    static int stupid_flag=0;

    istream& operator >> (istream& infile, freq_record& fr)
    {
      infile >> * (header_record *) (&fr);
      MY_DOUBLE_TYPE tmp;
      infile >> tmp;
      if (tmp < -1.e-10)
      {
        //fr.freq.initialize();
        fr.freq.deallocate();
      }
      else
      {
        dvector tmp_freq(1,fishery_freq_record::real_nlint);
        tmp_freq(1)=tmp; 
        int i;
        for (i=2;i<=tmp_freq.indexmax();i++)
        {
          infile >> tmp_freq(i); 
        }
        if (fishery_freq_record::group_factor==1)
        {
          fr.freq=tmp_freq;
        }     
        else
        {
          fr.freq.initialize();
          // group the length frequencies
          for (i=1;i<=tmp_freq.indexmax();i++)
          {
            int iind=(i-1)/fishery_freq_record::group_factor+1;
            fr.freq(iind)+=tmp_freq(i);
          }
        }
      }
      if (!infile)
      {
        cerr << "Error reading freq record" << endl;
        exit(1);
      }
      return infile;
    }

    cifstream& operator >> (cifstream& infile, freq_record& fr)
    {
      infile >> * (header_record *) (&fr);
      infile >> fr.freq;
      if (!infile)
      {
        cerr << "Error reading freq record" << endl;
        exit(1);
      }
      return infile;
    }


ostream& operator << (ostream& ofs,
                      fishery_freq_record_array& fcara)
{
  ofs << "# the number of  the number of  the shortest  interval   grouping" 
      << endl;
  ofs << "#  data sets    length intervals     length        width   factor"
      << endl;
  ofs << "   "  << fcara.nlfreq
       << "          "<< fcara.nlint << "     "  << fcara.shlen
       << "        " << fcara.filen << "   1   "   << endl;
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ofs << fcara(i);
  }
  return ofs;
}

ofstream& operator << (ofstream& ofs,
                      fishery_freq_record_array& fcara)
{
  ofs << "# the number of  the number of  the shortest  interval   grouping" 
      << endl;
  ofs << "#  data sets    length intervals     length        width   factor"
      << endl;
  ofs << "   "  << fcara.nlfreq
       << "          "<< fcara.nlint << "     "  << fcara.shlen
       << "        " << fcara.filen << "   1   "   << endl;
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ofs << fcara(i);
  }
  return ofs;
}

istream& operator >> (istream& ifs,
                      fishery_freq_record_array& fcara)
{
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ifs >> fcara(i);
  }
  //fcara.sort();
  //ofstream ofs("sorted.frq");
  //ofs << fcara;
  return ifs;
}

cifstream& operator >> (cifstream& ifs,
                      fishery_freq_record_array& fcara)
{
  char Errmsg[]= " Error reading fishery record ";
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ifs >> fcara(i);
    if (!ifs)
    {
      if (ifs.eof())
      {  
        cerr << endl << "  Read past end of file reading frq data" << endl;
      }
      cerr << Errmsg << i << endl;
      cerr << "previous record was  year: " << fcara(i-1).year
        << " month: " << fcara(i-1).month << " week: "
        << fcara(i-1).week <<  " fishery: "
        << fcara(i-1).fishery << endl << " freqs :" << endl
        << fcara(i-1).freq << endl;
      if (allocated(fcara(i-1).species) && _data_file_version_no>6)
        cerr << "species: " << fcara(i-1).species << " ";
      if (allocated(fcara(i-1).species1) && _data_file_version_no>7)
        cerr << "species1: " << fcara(i-1).species1;
      if (allocated(fcara(i-1).species2) && _data_file_version_no>8)
        cerr << "species2: " << fcara(i-1).species2;
      if (allocated(fcara(i-1).species3) && _data_file_version_no>8)
        cerr << "species3: " << fcara(i-1).species3;
       cerr << endl;
      cerr << "current record contains  year: " << fcara(i).year
        << " month: " << fcara(i).month << " week: "
        << fcara(i).week <<  " fishery: "
        << fcara(i).fishery << endl << " freqs :" << endl
        << fcara(i).freq << endl;
      if (allocated(fcara(i).species) && _data_file_version_no>6)
        cerr << "species: " << fcara(i).species << " ";
      if (allocated(fcara(i).species1) && _data_file_version_no>7)
        cerr << "species1: " << fcara(i).species1;
       cerr << endl;
      ad_exit(1);
      //ifs >> fhr.obs_tot_catch >> fhr.effort;
    }
   /*
    if ( fcara(i).fishery==4 && stupid_flag==0)
    {
      if (stupid_flag<=8)
      {
        fcara(i).obs_tot_catch=-1;
        stupid_flag++;
      }
    }
   */
    if (  (1 > fcara(i).week || 4 < fcara(i).week) ||
      ( 1 > fcara(i).month || 12 < fcara(i).month ) )
    {
      cerr << Errmsg << i << endl;
      cerr << " year was " <<  fcara(i).year << endl;
      cerr << " month was " <<  fcara(i).month << endl;
      cerr << " week was " <<  fcara(i).week << endl;
      if (i>1)
      {
        cerr << "previous record was  year: " << fcara(i-1).year
          << " month: " << fcara(i-1).month << " week: "
          << fcara(i-1).week <<  " fishery: "
          << fcara(i-1).fishery << endl << " freqs :" << endl;
          if (allocated(fcara(i-1).freq))
            cerr << "length freqs" << endl << fcara(i-1).freq << endl;
          if (allocated(fcara(i-1).wfreq))
            cerr << "weight freqs" << endl << fcara(i-1).wfreq << endl;
          if (allocated(fcara(i-1).age_freq)) 
            cerr << "age freqs" << endl << fcara(i-1).age_freq << endl;
      }
      ad_exit(1);
    }

    if (1 > fcara(i).month || 12 < fcara(i).month)
    {
      cerr << Errmsg << i << endl;
      cerr << " month was " <<  fcara(i).month << endl;
      ad_exit(1);
    }


    //if (fcara(i).year > 1995)
    //{
    //  cout <<"lmult_io.cpp " << fcara(i).year  << endl;
    //}
    //if (fcara(i).year < 1940 || fcara(i).year > 2020)
    //{
    //  cerr << "error reading record" << i << endl;
    //  if (i>1)
    //  {
    //    cerr << "previous record was  year: " << fcara(i-1).year
    //      << " month: " << fcara(i-1).month << " week: "
    //  << fcara(i-1).week << endl << " freqs :" << endl
    //          << fcara(i-1).freq << endl;
    //  }
    //}
  }
  return ifs;
}

    ostream& operator <<
           (ostream& ofs,fishery_freq_record& fhr)
    {
      ofs << * (fishery_header_record *) (&fhr);
      ofs << "#\n#  The length proportions in the catch\n#" << endl;
      ofs <<  setfixed()<< ivector(fhr.freq) << endl;
      return ofs;
    }

    cifstream& operator >> (cifstream& infile, fishery_freq_record& fr)
    {
      infile >> * (fishery_header_record *) (&fr);
      if (!infile)
      {
        cerr << "Error reading freq record" << endl;
        //exit(1);
      }
      else
      {
        // read length data  !! Dave F 
        if (fishery_freq_record::real_nlint)
        {
          MY_DOUBLE_TYPE tmp;
          infile >> tmp;
          if (tmp < -1.e-10)
          {
            //fr.freq.initialize();
            if (allocated(fr.freq))
              fr.freq.deallocate();
          }
          else
          {
            dvector tmp_freq(1,fishery_freq_record::real_nlint);
            tmp_freq.initialize(); 
            tmp_freq(1)=tmp; 
            int i;
            for (i=2;i<=tmp_freq.indexmax();i++)
            {
              infile >> tmp_freq(i); 
            }
            if (sum(tmp_freq)<1.e-50)
            {
               cerr << "You entered a zero length frequency set" << endl
                    << "This should have a -1 in the first entry"
                    << endl;
               
               cerr << * (fishery_header_record *) (&fr) << endl;
               ad_exit(1);
            }
            if (fishery_freq_record::group_factor==1)
            {
              if (allocated(fr.freq))
                fr.freq.deallocate();
              fr.freq.allocate(1,fishery_freq_record::real_nlint);
              fr.freq=tmp_freq;
            }     
            else
            {
              int grouped_nlint=tmp_freq.indexmax()/
                fishery_freq_record::group_factor;
              if (tmp_freq.indexmax()%fishery_freq_record::group_factor)
                grouped_nlint+=1;
              if (allocated(fr.freq))
                fr.freq.deallocate();
              fr.freq.allocate(1,grouped_nlint);
 
              fr.freq.initialize();
              // group the length frequencies
              for (i=1;i<=tmp_freq.indexmax();i++)
              {
                int iind=(i-1)/fishery_freq_record::group_factor+1;
                fr.freq(iind)+=tmp_freq(i);
              }
            }
          }
        }
        if (!infile)
        {
          cerr << "Error reading freq record" << endl;
          //exit(1);
        }
        // read weight data  !! Dave F 
        if (fishery_freq_record::real_nwint)
        {
          MY_DOUBLE_TYPE tmp;
          infile >> tmp;
          if (tmp < -1.e-10)
          {
            //fr.wfreq.initialize();
            fr.wfreq.deallocate();
          }
          else
          {
            dvector tmp_freq(1,fishery_freq_record::real_nwint);
            tmp_freq(1)=tmp; 
            int i;
            for (i=2;i<=tmp_freq.indexmax();i++)
            {
              infile >> tmp_freq(i); 
            }
            if (sum(tmp_freq)<1.e-50)
            {
               cerr << "You entered a zero weight frequency set" << endl
                    << "This should have a -1 in the first entry"
                    << endl;
               cerr << * (fishery_header_record *) (&fr) << endl;
               ad_exit(1);
            }
            if (fishery_freq_record::wgroup_factor==1)
            {
              fr.wfreq=tmp_freq;
            }     
            else
            {
              if (allocated(fr.wfreq))
              {
                fr.wfreq.deallocate();
              }
              int ntmp=(fishery_freq_record::real_nwint-1)/
                fishery_freq_record::wgroup_factor+1;
              fr.wfreq.allocate(1,ntmp);
              fr.wfreq.initialize();
              // group the length frequencies
              for (i=1;i<=tmp_freq.indexmax();i++)
              {
                int iind=(i-1)/fishery_freq_record::wgroup_factor+1;
                fr.wfreq(iind)+=tmp_freq(i);
              }
            }
          }
        }
        if (fishery_freq_record::age_nage)
        {
          MY_DOUBLE_TYPE tmp;
          infile >> tmp;
          if (tmp < -1.e-10)
          {
            //fr.wfreq.initialize();
            fr.age_freq.deallocate();
          }
          else
          {
            dvector tmp_age_freq(1,fishery_freq_record::age_nage);
            tmp_age_freq(1)=tmp; 
            for (int i=2;i<=tmp_age_freq.indexmax();i++)
            {
              infile >> tmp_age_freq(i); 
            }
            fr.age_freq=tmp_age_freq;
            // no grouping for age data
           /*
            if (fishery_freq_record::group_factor==1)
            {
              fr.wfreq=tmp_freq;
            }     
            else
            {
              fr.wfreq.initialize();
              // group the length frequencies
              for (i=1;i<=tmp_freq.indexmax();i++)
              {
                int iind=(i-1)/fishery_freq_record::wgroup_factor+1;
                fr.wfreq(iind)+=tmp_freq(i);
              }
            }
           */
          }
        }
        if (!infile)
        {
          cerr << "Error reading freq record" << endl;
          //exit(1);
        }
      }
      return infile;
    }

    istream& operator >> (istream& infile, fishery_freq_record& fr)
    {
      infile >> * (fishery_header_record *) (&fr);
      if (!infile)
      {
        cerr << "Error reading freq record" << endl;
        exit(1);
      }
      MY_DOUBLE_TYPE tmp;
      infile >> tmp;
      if (tmp < -1.e-10)
      {
        fr.freq.initialize();
      }
      else
      {
        dvector tmp_freq(1,fishery_freq_record::real_nlint);
        tmp_freq(1)=tmp; 
        int i;
        for (i=2;i<=tmp_freq.indexmax();i++)
        {
          infile >> tmp_freq(i); 
        }
        if (fishery_freq_record::group_factor==1)
        {
          fr.freq=tmp_freq;
        }     
        else
        {
          fr.freq.initialize();
          // group the length frequencies
          for (i=1;i<=tmp_freq.indexmax();i++)
          {
            int iind=(i-1)/fishery_freq_record::group_factor+1;
            fr.freq(iind)+=tmp_freq(i);
          }
        }
      }
      if (!infile)
      {
        cerr << "Error reading freq record" << endl;
        exit(1);
      }
      return infile;
    }

#undef HOME_VERSION
