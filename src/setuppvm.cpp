/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#pragma implementation "newmult.hpp"
#pragma implementation "variable.hpp"
#include "all.hpp"
#include "safe_mem.h"
//#include <ctrl87.h>
#include <stdio.h>

#if defined(WIN32)
#  include <windows.h>
#endif
#if defined(USE_ADPVM)
extern "C" {
#include "pvm3.h"
}
#include <adpvm2.h>
#endif //#if defined(USE_ADPVM)

ivector get_slave_information(adstring_array & hostnames,
  adstring_array & prognames,adstring_array & exedirs,
  adstring_array & workdirs,ivector& mult);
ivector match_slave_names(int nhost,pvmhostinfo * hinf,adstring_array& hnames);
extern adstring full_datafile_path;
extern adstring full_input_parfile_path;
extern adstring full_output_parfile_path;

extern int ad_argc;
extern char ** ad_argv;
int check_number_on_master(void);

char exe_directory[2][25]={"g:\\bin","/bin"};
char working_directory[2][25]={"g:\\temp","/tmp"};
char exe_name[2][25]={"gmult1","gmult"};

const char * mfpvm_execdir(int i) { return exe_directory[i]; }
const char * mfpvm_wkdir(int i) { return working_directory[i]; }
extern imatrix cl_switches;
const int MAXMYARGV=500;

mf_pvm_manager::mf_pvm_manager(void)
{
  ad_hostp=0;
  ad_nhost=0;
  //ad_stid=0;
  pvm_switch=0;
  pvm_save_switch=0;
  minspawn=0;
  maxspawn=-1;
  narch=-1;
}

extern mf_pvm_manager * mf_pvm;

#if defined(USE_ADPVM)
void setupargs(adstring& progname,adstring& exedir,adstring& workdir,
  char ** myargv)
{
  int ii=0;
  int nopt,nnopt;

  strcpy (myargv[ii++],"-exe");
  
  strcpy (myargv[ii++],progname);
  
  strcpy (myargv[ii++],"-ep");
  
  strcpy (myargv[ii++],(char*)exedir);
    
  strcpy (myargv[ii++],"-exeargs");
  
  strcpy (myargv[ii++],full_datafile_path);
  
  strcpy (myargv[ii++],full_input_parfile_path);
  
  strcpy (myargv[ii++],full_output_parfile_path);
  
  strcpy (myargv[ii++],"-slave");

  strcpy (myargv[ii++],"-wd");

  strcpy (myargv[ii++],(char*) workdir);

 /*
  // send this stuff after slaves are spawned
  if (allocated(cl_switches))
  {
    if (ii+3*cl_switches.indexmax()+3>=MAXMYARGV)
    {
      cerr << "Need to increase MAXMYARGV" << endl;
    }
    strcpy (myargv[ii++],"-switches");
    strcpy (myargv[ii++],str(cl_switches.indexmax()));
    for (int i=1;i<=cl_switches.indexmax();i++)
    {
      for (int j=1;j<=3;j++)
      {
        strcpy (myargv[ii++],str(cl_switches(i,j)));
      }
    }
  }
        
    
  if ((nopt=option_match(ad_argc,ad_argv,"-file",nnopt))>0)
  {
    strcpy (myargv[ii++],"-file");
    strcpy (myargv[ii++],ad_argv[nopt+1]);
  }
 */

  myargv[ii++]=0;

}
  


void mf_pvm_manager::setup_pvm(void)
  //struct pvmhostinfo * ad_hostp )  /* get configuration */
{
  int i,check;
  int ii=0;
  if (pvm_switch==1)
  {
    typedef char * charptr;
    char ** myargv=0;
    myargv = new charptr[MAXMYARGV];

    
    for (i=0;i<MAXMYARGV;i++) myargv[i] = new char[100];

   /*
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"-dbg");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"td32");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"-exe");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"gmult1");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"-ep");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"g:\\multifan_binary");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"-exeargs");
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],full_datafile_path);
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],full_input_parfile_path);
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],full_output_parfile_path);
  
    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"-slave");

    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"-wd");

    myargv[ii] = new char[100];
    strcpy (myargv[ii++],"g:\\temp");

    myargv[ii++] =0;
   */

    pvm_setopt(PvmRoute, PvmRouteDirect);  /* channel for communication */
    
    /* get and display configuration of the parallel machine */
    pvm_config( &ad_nhost, &narch, &ad_hostp );  
        /* get configuration */
    //for (i = 1; i < ad_nhost; i++)
    for (i = 0; i < ad_nhost; i++)
    { 
       printf("\t%s\n", ad_hostp[i].hi_name);
    }
    


    // read slave information can't have spaces in path names 
    // microsoft bite the dust.

    adstring_array hostnames;
    adstring_array exedirs; 
    adstring_array prognames; 
    adstring_array workdirs;
    ivector itmp=get_slave_information(hostnames,prognames,exedirs,workdirs,
      multiplier);
    cout << hostnames << endl;
    cout << prognames << endl;
    cout << exedirs << endl;
    cout << workdirs << endl;
    if(allocated(hostmatch))
      hostmatch.deallocate();
    
    minspawn=itmp.indexmin();
    maxspawn=itmp.indexmax();
    hostmatch.allocate(minspawn,maxspawn);
    hostmatch=itmp;

    if(allocated(ad_stid)) 
      ad_stid.deallocate();

    ad_stid.allocate(minspawn,maxspawn);

    // Here we have one slave per machine including the machine running the m aster
    // we will have to think about this
    // also where will the slave programs run

  
    ivector process_count(0,ad_nhost);
    process_count.initialize();
    
    for (i=hostmatch.indexmin(); i<=hostmatch.indexmax(); i++)	                           /* spawn processes on */
    {				/* all physical machines */
      //if (hostmatch(i)>-1)
      {
        //set up arguments for i'th spawned process
        // "scum:/home/dave/src/mfcl/trunk/src/multi2/sub", 
        adstring tmpexedir=hostnames(i) + ":" + exedirs(i);
        cout << tmpexedir << endl;

        int hostindex=hostmatch(i);
        process_count(hostmatch(i)++);
       
        int nopt,nnopt;
        if ((nopt=option_match(ad_argc,ad_argv,"-ddd",nnopt))>0)
        {
          myargv[0]=new char[100];
          strcpy(myargv[0],"--args");
          myargv[1]=new char[100];
          strcpy(myargv[1],prognames(i));
          myargv[2]=new char[100];
          strcpy(myargv[2],full_datafile_path);
          myargv[3]=new char[100];
          strcpy(myargv[3],full_input_parfile_path);
          myargv[4]=new char[100];
          strcpy(myargv[4],full_output_parfile_path);
          myargv[5]=new char[100];
          strcpy(myargv[5],"-slave");
          myargv[6]=0;
        
          // run slave with debugger
          //strcpy(myargv[1],"mftimecatchslave.exe")
          //check = pvm_spawn("drunslave", myargv, 
          //  PvmTaskHost | PvmTaskDebug,mf_pvm->ad_hostp[hostindex].hi_name, 1, 
          char ts[100];
          strcpy(ts,(char*)(tmpexedir));
	  check = pvm_spawn(prognames(i),  /*"/usr/bin/mfclo64", */
                 myargv,  /*(char**)0, */
                    PvmTaskDebug,
                    ts, /* "scum:/home/dave/src/mfcl/trunk/src/multi2/sub", */
                    1, 
                &((ad_stid)(i)));
          /*
          check = pvm_spawn("/usr/bin/mfclo64", myargv, 
            PvmTaskHost | PvmTaskDebug,mf_pvm->ad_hostp[0].hi_name, 1, 
            &((ad_stid)(1)));
          */
        }
        else
        {
          //myargv[0]=new char[100];
          //strcpy(myargv[0],prognames(i));
          myargv[0]=new char[100];
          strcpy(myargv[0],full_datafile_path);
          myargv[1]=new char[100];
          strcpy(myargv[1],full_input_parfile_path);
          myargv[2]=new char[100];
          strcpy(myargv[2],full_output_parfile_path);
          myargv[3]=new char[100];
          strcpy(myargv[3],"-slave");
          myargv[4]=0;
        
          char ts[100];
          strcpy(ts,(char*)(tmpexedir));
	  check = pvm_spawn(prognames(i),  /*"/usr/bin/mfclo64", */
                 myargv,  /*(char**)0, */
                    0, /*PvmTaskDebug,*/
                    ts, /* "scum:/home/dave/src/mfcl/trunk/src/multi2/sub", */
                    1, 
                &((ad_stid)(i)));
          /*
          check = pvm_spawn("/usr/bin/mfclo64", myargv, 
            PvmTaskHost | PvmTaskDebug,mf_pvm->ad_hostp[0].hi_name, 1, 
            &((ad_stid)(1)));
          */
        }

        if (!check) 
        {
          printf("Couldn't start process %d on %s\n",
            process_count(hostindex),mf_pvm->ad_hostp[hostindex].hi_name); 
          printf("arguments are\n");
          int icount=0;
          do
          {
            if (myargv[icount]==0) break;
            printf("%s\n",myargv[icount++]);
          }
          while(1);
          
        }

      }

      cout << hostmatch << endl;
        // so we can find out if slaves crash (I hope)        
      int ierr=pvm_notify(PvmTaskExit,301,hostmatch.size(),
            &((ad_stid)(i)));
          //  &((ad_stid)(hostmatch.indexmin())));
      if (ierr<0)
      {
        cerr << "Error in pvm_notify -- return value = " << ierr << endl;
        //ad_exit(1);
      }

    }
  }
}
ivector match_slave_names(int nhost,pvmhostinfo * hinf,adstring_array& hnames)
{
  int mmin=hnames.indexmin();
  int mmax=hnames.indexmax();
  ivector hostmatch(mmin,mmax);
  for (int j=mmin;j<=mmax;j++)
  {
    hostmatch(j)=-1;
    for (int i=0;i<nhost;i++)
    {
      if (adstring(hinf[i].hi_name)==hnames(j)) hostmatch(j)=i;
    }
  }
  return hostmatch;
}

ivector get_slave_information(adstring_array & hostnames,
  adstring_array& prognames,adstring_array & exedirs,
  adstring_array & workdirs,ivector& multiplier)
{
  const int maxnslaves=100;
  ivector nspawn(1,maxnslaves);
  ivector tmpmultiplier(1,maxnslaves);
  const int maxnprocesses=1000;
  ivector tmphostmatch(1,maxnprocesses);
  int nlist=0;
  int nprocesses=0;
  adstring tmpname;
  adstring tmpprogname;
  adstring tmpedir;
  adstring tmpwdir;
  MY_DOUBLE_TYPE tmpmult;
  adstring  hostfilename;
  int nopt,nnopt;
  if ((nopt=option_match(ad_argc,ad_argv,"-slinfo",nnopt))>0)
    hostfilename=ad_argv[nopt+1];
  else
    hostfilename="mfhostfile";
  cifstream cif((char*)hostfilename);

  if (!cif)
  {
    cerr << "Error trying to open slaveinfo file with name "
      << ad_argv[nopt+1] << endl;
    ad_exit(1);
  }

  do
  {
    cif >> tmpname;
    if (cif.eof()) break;
    if (tmpname.size()==0) break;
    cif >> tmpprogname;
    if (cif.eof()) break;
    cif >> tmpedir;
    if (cif.eof()) break;
    cif >> tmpwdir;
    if (cif.eof()) break;
    cif >> tmpmult;
    if (cif.eof()) break;
    // can we match the name with pvm hostname
    int hm=-1;
    for (int i=0;i<mf_pvm->ad_nhost;i++)
    {
      if (adstring(mf_pvm->ad_hostp[i].hi_name)==tmpname) 
      {
        if (++nprocesses>maxnprocesses)
        {
          cerr << "need to increase maxnprocesses" << endl;
          ad_exit(1);
        }
        tmphostmatch(nprocesses)=i;
        tmpmultiplier(nprocesses)=tmpmult;
        prognames.append_distinct(tmpprogname);
        hostnames.append_distinct(tmpname);
        exedirs.append_distinct(tmpedir);
        workdirs.append_distinct(tmpwdir);
        break;
      }
    }
  }
  while(1);
  //return nspawn(1,nlist);
  ivector hostmatch(1,nprocesses);
  hostmatch=tmphostmatch(1,nprocesses);
  if (allocated(multiplier)) multiplier.deallocate();
  multiplier.allocate(1,nprocesses);
  multiplier=tmpmultiplier(1,nprocesses);
  return hostmatch;
}


int check_number_on_master(void)
{
  int nspawn=1;
  int nopt,nnopt;
  if ((nopt=option_match(ad_argc,ad_argv,"-file",nnopt))>0)
  {
    if (nnopt==1)
      nspawn=atoi(ad_argv[nopt+1]);
    else 
      nspawn=1;
  }
  return nspawn;
}

void set_tag_group_assignments(dvar_len_fish_stock_history& fsh)
{
  int mmin=mf_pvm->minspawn;
  int mmax=mf_pvm->maxspawn;
  // number of slaves plus master
  int nprocs=mmax-mmin+2;
  if (mf_pvm->pvm_switch==1)
  {
    int i;
    ivector bs(mmin,mmax+1);
    ivector min_group(mmin,mmax+1);
    ivector max_group(mmin,mmax+1);
    dvector weights(mmin,mmax+1);
    weights(mmin,mmax)=mf_pvm->multiplier;
    //weights(mmax+1)=1000.0;
    weights(mmax+1)=0.0;
    weights/=sum(weights);

    bs=ivector(fsh.num_tag_releases*weights);
    int g=sum(bs);
    int r=fsh.num_tag_releases-g;
    for (i=mf_pvm->minspawn;i<=mf_pvm->minspawn+r-1;i++) bs(i)+=1;

    int offset=0;
    for (i=mf_pvm->minspawn; i<=mf_pvm->maxspawn+1; i++)
    {
      min_group(i)=1+offset;
      max_group(i)=bs(i)+offset;
      offset+=bs(i);
    }
    if (max_group(mf_pvm->maxspawn+1) !=fsh.num_tag_releases)
    {
      cerr << "error in sum of group assignments" << endl;
      ad_exit(1);
    }
    send_taggroup_assignments(min_group,max_group,fsh);
    fsh.min_tag_group=min_group(mf_pvm->maxspawn+1);
    fsh.max_tag_group=max_group(mf_pvm->maxspawn+1);
    //fsh.min_tag_group=0;
    //fsh.max_tag_group=-1;
  }
  else
  {
    get_taggroup_assignments(fsh);
  }
}


void mf_timer::initialize(void) { ftime(&t0); }

int mf_timer::report(void) 
{ 
  ftime(&t1);
  return  1000.*(t1.time-t0.time) + t1.millitm - t0.millitm;
}
void check_slave_crash()
{
  int bufid = pvm_nrecv(-1,301);

  int tid=0;
  if (bufid>0)
  {
    do
    {
      int info=pvm_upkint(&tid,1,1);
      if (info<0) break;
      cerr << "slave " << tid << " crashed" << endl;
    }
    while(1);
    ad_exit(1);
  }
}

#endif //#if defined(USE_ADPVM)
