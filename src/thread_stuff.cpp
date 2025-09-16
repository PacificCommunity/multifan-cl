/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#include <cstddef>
#include "adthread.h"
#if !defined(OPT_LIB)
#  if !defined(CHK_ID_STRING)
#    define CHK_ID_STRING
#  endif
#endif

/* prototype for thread routine */
void admb_thread(void * ptr);

pthread_mutex_t mutex_print;

//const int NSLAVES=1;

typedef char * pchar;
typedef pchar *  ppchar;
typedef ofstream * pofstream;
__ADMBTHREAD__ int adpthread_manager::slave_number;

int adpthread_manager::old_buffer_flag=0;


adpthread_manager::adpthread_manager(int ns,int bs) : buffer_size(1,ns),
  mflag(1,ns), sflag(1,ns), num_in_group(1,1),ngroups(1),
  initial_buffer_size(bs)
{
  num_in_group=ns;
  gmin.allocate(1,ngroups);
  gmax.allocate(1,ngroups);
  gmin(1)=1;
  gmax(1)=gmin(1)+num_in_group(1)-1;
  logflag=0;
  mflag.initialize();
  sflag.initialize();
  int i;
  buffer_size=bs;
  nslaves=ns;
  transfer_buffer=new pchar[ns];
  current_bptr=new pchar[ns];
  current_bptr--;
  buffend=new pchar[ns];
  transfer_buffer--;
  buffend--;
  // for timing 0 is for master 1-n for slaves
  if (logflag)
  {
    adt=new adtimer[ns+1];
    logfiles=new pofstream[ns+1]; 
    adstring fname="log_master";
    logfiles[0]=new ofstream(fname);
    for (i=1;i<=ns;i++)
    {
      adstring fname="log_slave_" + str(i);
      logfiles[i]=new ofstream(fname);
    }
  }
  ssflag.allocate(0,ns-1);
  smflag.allocate(0,ns-1);
  for (i=0;i<ns;i++)
  {
    ssflag(i).allocate(i+1,ns);
    smflag(i).allocate(i+1,ns);
  }
  smflag.initialize();
  ssflag.initialize();
  
  ssmutex=new ppthread_mutex_t[ns];
  for (i=0;i<ns;i++)
  {
    ssmutex[i]=new pthread_mutex_t[ns-i];
    ssmutex[i]-=i+1;
  }

  sbuffer_size.allocate(0,ns-1);
  for (i=0;i<ns;i++)
  {
    sbuffer_size(i).allocate(i+1,ns);
  }
  sbuffer_size.initialize();

 stransfer_buffer=new ppchar[ns];

  for (i=0;i<ns;i++)
  {
    stransfer_buffer[i]=new pchar[ns-i];
    stransfer_buffer[i]-=i+1;
  }
  scurrent_bptr=new ppchar[ns];
  for (i=0;i<ns;i++)
  {
    scurrent_bptr[i]=new pchar[ns-i];
    scurrent_bptr[i]-=i+1;
  }

 sbuffend=new ppchar[ns];
  for (i=0;i<ns;i++)
  {
    sbuffend[i]=new pchar[ns-i];
    sbuffend[i]-=i+1;
  }

  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      stransfer_buffer[i][j]=0;
      scurrent_bptr[i][j]=0;
      sbuffend[i][j]=0;
    }
  }

  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      pthread_mutex_init(ssmutex[i]+j,NULL);
    }
  }
  sscondition=new ppthread_cond_t[ns];
  for (i=0;i<ns;i++)
  {
    sscondition[i]=new pthread_cond_t[ns-i];
    sscondition[i]-=i+1;
  }
  
  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      pthread_cond_init(sscondition[i]+j,NULL);
    }
  }

  smcondition=new ppthread_cond_t[ns];
  for (i=0;i<ns;i++)
  {
    smcondition[i]=new pthread_cond_t[ns-i];
    smcondition[i]-=i+1;
  }
  
  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      pthread_cond_init(smcondition[i]+j,NULL);
    }
  }
}
adpthread_manager::adpthread_manager(int _ngroups,ivector& _num_in_group,int bs)
  : num_in_group(_num_in_group),ngroups(_ngroups),initial_buffer_size(bs)
{
  int ns=sum(num_in_group);
  gmin.allocate(1,ngroups);
  gmax.allocate(1,ngroups);
  gmin(1)=1;
  gmax(1)=gmin(1)+num_in_group(1)-1;
  for (int i=2;i<=ngroups;i++)
  {
    gmin(i)=gmax(i-1)+1;
    gmax(i)=gmin(i)+num_in_group(i)-1;
  }
    
  nslaves=ns;
  buffer_size.allocate(1,ns),
  mflag.allocate(1,ns);
  sflag.allocate(1,ns); 
  logflag=1;
  mflag.initialize();
  sflag.initialize();
  int i;
  buffer_size=bs;
  transfer_buffer=new pchar[ns];
  current_bptr=new pchar[ns];
  current_bptr--;
  buffend=new pchar[ns];
  transfer_buffer--;
  buffend--;
  // for timing 0 is for master 1-n for slaves
  if (logflag)
  {
    adt=new adtimer[ns+1];
    logfiles=new pofstream[ns+1]; 
    adstring fname="log_master";
    logfiles[0]=new ofstream(fname);
    for (i=1;i<=ns;i++)
    {
      adstring fname="log_slave_" + str(i);
      logfiles[i]=new ofstream(fname);
    }
  }
  
  for (i=1;i<=ns;i++)
  {
    transfer_buffer[i]=new char[bs];
    current_bptr[i]=transfer_buffer[i];
    buffend[i]=transfer_buffer[i]+bs-1;
  }
  smutex=new pthread_mutex_t[ns];
  smutex--;
  pthread_mutex_init(&copy_mutex,NULL);
  pthread_mutex_init(&start_mutex,NULL);
  for (i=1;i<=ns;i++)
  {
    pthread_mutex_init(smutex+i,NULL);
  }
  scondition=new pthread_cond_t[ns];
  scondition--;
  mcondition=new pthread_cond_t[ns];
  mcondition--;
  for (i=1;i<=ns;i++)
  {
    pthread_cond_init(scondition+i,0);
    pthread_cond_init(mcondition+i,0);
  }
  thread1=new pthread_t[ns];
  thread1--;
  ppf = new pthreadfun[ngroups];
  ppf--;
  for (i=1;i<=ngroups;i++)
  {
    ppf[i]= (pthreadfun)(0);
  }
  num_code=0;
  ssflag.allocate(0,ns-1);
  smflag.allocate(0,ns-1);
  for (i=0;i<ns;i++)
  {
    ssflag(i).allocate(i+1,ns);
    smflag(i).allocate(i+1,ns);
  }
  smflag.initialize();
  ssflag.initialize();

  ssmutex=new ppthread_mutex_t[ns];
  for (i=0;i<ns;i++)
  {
    ssmutex[i]=new pthread_mutex_t[ns-i];
    ssmutex[i]-=i+1;
  }
  sbuffer_size.allocate(0,ns-1);
  for (i=0;i<ns;i++)
  {
    sbuffer_size(i).allocate(i+1,ns);
  }
  sbuffer_size.initialize();

  stransfer_buffer=new ppchar[ns];

  for (i=0;i<ns;i++)
  {
    stransfer_buffer[i]=new pchar[ns-i];
    stransfer_buffer[i]-=i+1;
  }
  scurrent_bptr=new ppchar[ns];
  for (i=0;i<ns;i++)
  {
    scurrent_bptr[i]=new pchar[ns-i];
    scurrent_bptr[i]-=i+1;
  }

  sbuffend=new ppchar[ns];
  for (i=0;i<ns;i++)
  {
    sbuffend[i]=new pchar[ns-i];
    sbuffend[i]-=i+1;
  }

  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      stransfer_buffer[i][j]=0;
      scurrent_bptr[i][j]=0;
      sbuffend[i][j]=0;
    }
  }
  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      pthread_mutex_init(ssmutex[i]+j,NULL);
    }
  }
  sscondition=new ppthread_cond_t[ns];
  for (i=0;i<ns;i++)
  {
    sscondition[i]=new pthread_cond_t[ns-i];
    sscondition[i]-=i+1;
  }
  
  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      pthread_cond_init(sscondition[i]+j,NULL);
    }
  }

  smcondition=new ppthread_cond_t[ns];
  for (i=0;i<ns;i++)
  {
    smcondition[i]=new pthread_cond_t[ns-i];
    smcondition[i]-=i+1;
  }
  
  for (i=0;i<ns;i++)
  {
    for (int j=i+1;j<=ns;j++)
    {
      pthread_cond_init(smcondition[i]+j,NULL);
    }
  }
}

void adpthread_manager::attach_code(pthreadfun pf)
{
  ppf[++num_code]=pf;
}

void adpthread_manager::check_buffer_size_read(int nbytes,int sno)
{
  
  // if the buffer is too small on read there is an error
  if (current_bptr[sno]+nbytes>buffend[sno])
  {
    std::ptrdiff_t pd=current_bptr[sno]-transfer_buffer[sno];
    std::ptrdiff_t pd1=buffend[sno]-transfer_buffer[sno];
    cout << "current offset is " << pd << " bytes " << endl;
    cout << "trying to read  " << nbytes << " bytes " << endl;
    cout << "buffend is at     " << pd1 << " bytes " << endl;
    cerr << "error reading from master slave buffer for slave "
         << sno << endl;
    ad_exit(1);
  }
}

void adpthread_manager::check_buffer_size_read(int nbytes,int s1,int s2)
{
  if (stransfer_buffer[s1][s2]== 0)
  {
    cerr << "Need to fix this" << endl;
    //transfer_buffer[s1][s2]=new char[bs];
    //scurrent_bptr[s1][s2]=transfer_buffer[s1][s2];
    //sbuffend[s1][s2]=stransfer_buffer[s1][s2]+bs-1;
    cout << "Initialized transfer buffer for pair " 
         << s1 << "  " << s2 << endl;
  }
  
  // if the buffer is too small on read theree is an error
  if (scurrent_bptr[s1][s2]+nbytes>sbuffend[s1][s2])
  {
    std::ptrdiff_t pd=scurrent_bptr[s1][s2]-stransfer_buffer[s1][s2];
    std::ptrdiff_t pd1=sbuffend[s1][s2]-stransfer_buffer[s1][s2];
    cout << "current offset is " << pd << " bytes " << endl;
    cout << "trying to read  " << nbytes << " bytes " << endl;
    cout << "sbuffend is at     " << pd1 << " bytes " << endl;
    cerr << "error reading from  buffer for pair "
         << s1 << "  " << s2 << endl;
    ad_exit(1);
  }
}
void adpthread_manager::check_buffer_size(int nbytes,int s1,int s2)
{
  // if the buffer is too small make it bigger and copy old
  // buffer contents
  if (scurrent_bptr[s1][s2]+nbytes>sbuffend[s1][s2])
  {
    std::ptrdiff_t pd=scurrent_bptr[s1][s2]-stransfer_buffer[s1][s2];
    std::ptrdiff_t pd1=sbuffend[s1][s2]-stransfer_buffer[s1][s2];
    cout << "scurrent offset is " << pd << " bytes " << endl;
    cout << "trying to write " << nbytes << " bytes " << endl;
    cout << "scurrent sbuffend is at " << pd1 << " bytes " << endl;
    cerr << "master increasing master slave buffer for pair "
         << s1  << "  "  << s2 << endl;
    if (sbuffer_size(s1,s2)==0)
    {
      sbuffer_size(s1,s2)=initial_buffer_size;
      stransfer_buffer[s1][s2]=new char[sbuffer_size(s1,s2)];
      scurrent_bptr[s1][s2]=stransfer_buffer[s1][s2];
    }
    else
    {
      char * tmp=stransfer_buffer[s1][s2];
      sbuffer_size(s1,s2)*=2;
      stransfer_buffer[s1][s2]=new char[sbuffer_size(s1,s2)];
      ::memcpy(stransfer_buffer[s1][s2],tmp,pd+1);
      scurrent_bptr[s1][s2]=stransfer_buffer[s1][s2]+pd;
      delete [] tmp;
      tmp=0;
    }
    sbuffend[s1][s2]=stransfer_buffer[s1][s2]+sbuffer_size(s1,s2)-1;
    std::ptrdiff_t pd2=sbuffend[s1][s2]-stransfer_buffer[s1][s2];
    cout << " new sbuffend is at " << pd2 << " bytes " << endl;
  }
}


/*
void adpthread_manager::check_buffer_size(int nbytes,int sno)
{
  // if the buffer is too small make it bigger and copy old
  // buffer contents
  if (current_bptr[sno]+nbytes>buffend[sno])
  {
    std::ptrdiff_t pd=current_bptr[sno]-transfer_buffer[sno];
    std::ptrdiff_t pd1=buffend[sno]-transfer_buffer[sno];
    cout << "current offset is " << pd << " bytes " << endl;
    cout << "trying to write " << nbytes << " bytes " << endl;
    cout << "current buffend is at " << pd1 << " bytes " << endl;
    cerr << "master increasing master slave buffer for slave "
         << sno << endl;
    char * tmp=transfer_buffer[sno];
    buffer_size[sno]*=2;
    transfer_buffer[sno]=new char[buffer_size[sno]];
    ::memcpy(transfer_buffer[sno],tmp,pd+1);
    current_bptr[sno]=transfer_buffer[sno]+pd;
    delete [] tmp;
    tmp=0;
    buffend[sno]=transfer_buffer[sno]+buffer_size[sno]-1;
    std::ptrdiff_t pd2=buffend[sno]-transfer_buffer[sno];
    cout << " new buffend is at " << pd2 << " bytes " << endl;
  }
}
*/

void adjoint_write_lock_buffer_master(void)
{
  verify_identifier_string("CT");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_write_lock_buffer_master();
}

void adpthread_manager::adjoint_write_lock_buffer_master(void)
{
  //unlock read buffer slave
  // cout  << "adjoint_write_lock_buffer_master(void)" << endl;
  verify_identifier_string("OD");
  int sno=restore_int_value();
  current_bptr[sno]=transfer_buffer[sno];
  verify_identifier_string("YX");
  pthread_mutex_unlock(smutex+sno) ;
  // !!!!!!!!!!!!!!!!!!!!
  //pthread_cond_signal(mcondition+sno);
  pthread_cond_signal(scondition+sno);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[0].get_elapsed_time_and_reset();
    *(logfiles[0]) << tm << " adjoint_write_lock_buffer_master"   << endl;
  }
}
void adpthread_manager::write_lock_buffer_master(int sno)
{
  if (logflag)
  {
    *(logfiles[0]) << " write_lock_buffer_master"   << endl;
    adt[0].get_elapsed_time_and_reset();
  }
  // we should be using the mutexes from the class
  // to get rid of global variables
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  while (mflag[sno] == 1 || sflag[sno] ==1 ) 
    pthread_cond_wait(mcondition+sno,smutex+sno);
//    save_identifier_string("YX");
  const char * str1;
  str1="YX";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  save_int_value(sno);
//    save_identifier_string("OD");
  const char * str2;
  str2="OD";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  save_pointer_value(this);
//    save_identifier_string("CT");
  const char * str3;
  str3="CT";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(::adjoint_write_lock_buffer_master);
}   

void adjoint_read_lock_buffer_master(void)
{
  verify_identifier_string("ZE");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_read_lock_buffer_master();
}
void adpthread_manager::adjoint_read_lock_buffer_master(void)
{
  // cout  << "adjoint_read_lock_buffer_master(void)" << endl;
  verify_identifier_string("XE");
  int sno=restore_int_value();
  verify_identifier_string("WE");
  //write_unlock_buffer_master(sno);
  current_bptr[sno]=transfer_buffer[sno];
  mflag[sno] = 1;
  pthread_mutex_unlock(smutex+sno);
  pthread_cond_signal(scondition+sno);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[0].get_elapsed_time_and_reset();
    *(logfiles[0]) << tm << " adjoint_read_lock_buffer_master" << endl;
  }
}


void adpthread_manager::read_lock_buffer_master(int sno)
{
  // we should be using the mutexes from the class
  // to get rid of global variables
  if (logflag)
  {
    *(logfiles[0]) << "starting read_lock_buffer_master " << sno   << endl;
    adt[0].get_elapsed_time_and_reset();
  }

  pthread_mutex_lock(smutex+sno);
  // only read if buffer full and data is for you
  int sf=sflag[sno]; 
  while (sflag[sno] == 0) 
    pthread_cond_wait(mcondition+sno,smutex+sno);
  // cout  << "A sflag=0" << endl;
  sflag[sno] = 0;
//    save_identifier_string("WE");
  const char * str4;
  str4="WE";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  save_int_value(sno);
//    save_identifier_string("XE");
  const char * str5;
  str5="XE";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  save_pointer_value(this);
//    save_identifier_string("ZE");
  const char * str6;
  str6="ZE";
  char* strx6=const_cast <char*> (str6);
  save_identifier_string(strx6);
  gradient_structure::GRAD_STACK1->
     set_gradient_stack(::adjoint_read_lock_buffer_master);
}   
void adjoint_write_unlock_buffer_master(void)
{
  verify_identifier_string("ZD");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_write_unlock_buffer_master();
}
void adpthread_manager::adjoint_write_unlock_buffer_master(void)
{
  if (logflag)
  {
    adt[0].get_elapsed_time_and_reset();
  }
  // cout  << "adjoint_write_unlock_buffer_master(void)" << endl;
  verify_identifier_string("XD");
  int sno=restore_int_value();
  verify_identifier_string("TD");
  //read_lock_buffer_master(sno);
  pthread_mutex_lock(smutex+sno);
  // only read if buffer full and data is for you
  // cout  << endl << "XX sflag[sno] = " << sflag[sno] << endl;
  int sf=sflag[sno];
  while (sflag[sno] == 0) 
    pthread_cond_wait(mcondition+sno,smutex+sno);
  // cout  << endl << "VV sflag[sno] = " << sflag[sno] << endl;
  // cout  << "B sflag=0" << endl;
  sflag[sno] = 0;
}

void adpthread_manager::write_unlock_buffer_master(int sno)
{
  // we should be using the mutexes from the class
  // to get rid of global variables
  current_bptr[sno]=transfer_buffer[sno];
  mflag[sno] = 1;
  pthread_mutex_unlock(smutex+sno);
  pthread_cond_signal(scondition+sno);
//    save_identifier_string("TD");
  const char * str7;
  str7="TD";
  char* strx7=const_cast <char*> (str7);
  save_identifier_string(strx7);
  save_int_value(sno);
//    save_identifier_string("XD");
  const char * str8;
  str8="XD";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
  save_pointer_value(this);
//    save_identifier_string("ZD");
  const char * str9;
  str9="ZD";
  char* strx9=const_cast <char*> (str9);
  save_identifier_string(strx9);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_write_unlock_buffer_master);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[0].get_elapsed_time_and_reset();
    *(logfiles[0]) << tm << " write_unlock_buffer_master" << endl;
  }
}   
void adjoint_write_unlock_buffer_slave(void)
{
  verify_identifier_string("CZ");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_write_unlock_buffer_slave();
}

void adpthread_manager::adjoint_write_unlock_buffer_slave(void)
{
  // cout  << "adjoint_write_unlock_buffer_slave(void)" << endl;
  verify_identifier_string("CD");
  int sno=restore_int_value();
  if (logflag)
  {
    adt[sno].get_elapsed_time_and_reset();
  }
  verify_identifier_string("UD");
  // we should be using the mutexes from the class
  // to get rid of global variables
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  while (mflag[sno] == 0) 
    pthread_cond_wait(scondition+sno,smutex+sno);
  mflag[sno]=0;
}
void adpthread_manager::write_unlock_buffer_slave(int sno)
{
  current_bptr[sno]=transfer_buffer[sno];
  sflag[sno] = 1;
  pthread_mutex_unlock(smutex+sno);
  pthread_cond_signal(mcondition+sno);
//    save_identifier_string("UD");
  const char * str10;
  str10="UD";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
  save_int_value(sno);
//    save_identifier_string("CD");
  const char * str11;
  str11="CD";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
  save_pointer_value(this);
//    save_identifier_string("CZ");
  const char * str12;
  str12="CZ";
  char* strx12=const_cast <char*> (str12);
  save_identifier_string(strx12);
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(::adjoint_write_unlock_buffer_slave);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[sno].get_elapsed_time_and_reset();
    *(logfiles[sno]) << tm  <<  " write_unlock_buffer_slave" << endl;
  }

}   

void adjoint_write_lock_buffer_slave(void)
{
  verify_identifier_string("CY");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_write_lock_buffer_slave();
}

void adpthread_manager::adjoint_write_lock_buffer_slave(void)
{
  //unlock read buffer slave
  // cout  << "adjoint_write_lock_buffer_slave(void)" << endl;
  verify_identifier_string("ED");
  int sno=restore_int_value();
  current_bptr[sno]=transfer_buffer[sno];
  verify_identifier_string("FN");
  pthread_mutex_unlock(smutex+sno) ;
  pthread_cond_signal(mcondition+sno);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[sno].get_elapsed_time_and_reset();
    *(logfiles[sno]) << tm << " adjoint_write_lock_buffer_slave"   << endl;
  }
}
void adpthread_manager::write_lock_buffer_slave(int sno)
{
  if (logflag)
  {
    adt[sno].get_elapsed_time_and_reset();
    *(logfiles[sno]) << "starting write_lock_buffer_slave " << sno  << endl;
  }
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  while (mflag[sno] == 1 || sflag[sno] ==1 ) 
    pthread_cond_wait(scondition+sno,smutex+sno);
//    save_identifier_string("FN");
  const char * str13;
  str13="FN";
  char* strx13=const_cast <char*> (str13);
  save_identifier_string(strx13);
  save_int_value(sno);
//    save_identifier_string("ED");
  const char * str14;
  str14="ED";
  char* strx14=const_cast <char*> (str14);
  save_identifier_string(strx14);
  save_pointer_value(this);
//    save_identifier_string("CY");
  const char * str15;
  str15="CY";
  char* strx15=const_cast <char*> (str15);
  save_identifier_string(strx15);
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(::adjoint_write_lock_buffer_slave);
}

void adjoint_read_lock_buffer_slave(void)
{
  verify_identifier_string("D1");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_read_lock_buffer_slave();
}
void adpthread_manager::adjoint_read_lock_buffer_slave(void)
{
  // cout  << "adjoint_read_lock_buffer_slave(void)" << endl;
  //write_unlock_buffer_slave(sno);
  verify_identifier_string("E1");
  int sno=restore_int_value();
  verify_identifier_string("H1");
  // we should be using the mutexes from the class
  // to get rid of global variables
  current_bptr[sno]=transfer_buffer[sno];
  // cout  << "adpthread_manager::adjoint_read_lock_buffer_slave(void)=1" << endl;
  sflag[sno] = 1;
  pthread_mutex_unlock(smutex+sno);
  pthread_cond_signal(mcondition+sno);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[sno].get_elapsed_time_and_reset();
    *(logfiles[sno]) << tm   << " adjoint_read_lock_buffer_slave" << endl;
  }
}

void adpthread_manager::read_lock_buffer_slave(int sno)
{
  // we should be using the mutexes from the class
  // to get rid of global variables
  if (logflag)
  {
    *(logfiles[sno]) << "starting read_lock_buffer_slave "  << sno  << endl;
    adt[sno].get_elapsed_time_and_reset();
  }
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  while (mflag[sno] == 0) 
    pthread_cond_wait(scondition+sno,smutex+sno);
  mflag[sno]=0;
//    save_identifier_string("H1");
  const char * str16;
  str16="H1";
  char* strx16=const_cast <char*> (str16);
  save_identifier_string(strx16);
  save_int_value(sno);
//    save_identifier_string("E1");
  const char * str17;
  str17="E1";
  char* strx17=const_cast <char*> (str17);
  save_identifier_string(strx17);
  save_pointer_value(this);
//    save_identifier_string("D1");
  const char * str18;
  str18="D1";
  char* strx18=const_cast <char*> (str18);
  save_identifier_string(strx18);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_read_lock_buffer_slave);
}
void adjoint_master_write_lock(void)
{
  verify_identifier_string("G5");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  verify_identifier_string("G7");
  int sno=restore_int_value();
  verify_identifier_string("G9");
  ptr->adjoint_master_write_lock(sno);
}

void  adpthread_manager::adjoint_master_write_lock(int sno)
{
  // cout  << "adpthread_manager::adjoint_master_write_lock(int sno)" << endl;
  if (logflag)
  {
    adt[0].get_elapsed_time_and_reset();
    *(logfiles[0]) << "starting adjoint_master_write_lock " << sno   << endl;
  }
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  current_bptr[sno]=transfer_buffer[sno];
  while (mflag[sno] == 1 || sflag[sno] ==1 )
    pthread_cond_wait(mcondition+sno,smutex+sno);
}

void adjoint_read_unlock_buffer_master(void)
{
  verify_identifier_string("G5");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_read_unlock_buffer_master();
}

void  adpthread_manager::adjoint_read_unlock_buffer_master(void)
{
  // cout  << "adpthread_manager::adjoint_read_unlock_buffer_master(void);" << endl;
  verify_identifier_string("G7");
  int sno=restore_int_value();
  verify_identifier_string("G9");
  if (logflag)
  {
    *(logfiles[0]) << "starting adjoint_read_unlock_buffer_master " << sno 
                   << endl;
    adt[0].get_elapsed_time_and_reset();
  }
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  while (mflag[sno] == 1 || sflag[sno] ==1 )
    pthread_cond_wait(mcondition+sno,smutex+sno);
}

void adpthread_manager::read_unlock_buffer_master(int sno)
{
  // we should be using the mutexes from the class
  // to get rid of global variables
  current_bptr[sno]=transfer_buffer[sno];
  pthread_mutex_unlock(smutex+sno) ;
  pthread_cond_signal(scondition+sno);
  
//    save_identifier_string("G9");
  const char * str19;
  str19="G9";
  char* strx19=const_cast <char*> (str19);
  save_identifier_string(strx19);
  save_int_value(sno);
//    save_identifier_string("G7");
  const char * str20;
  str20="G7";
  char* strx20=const_cast <char*> (str20);
  save_identifier_string(strx20);
  save_pointer_value(this);
//    save_identifier_string("G5");
  const char * str21;
  str21="G5";
  char* strx21=const_cast <char*> (str21);
  save_identifier_string(strx21);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_read_unlock_buffer_master);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[0].get_elapsed_time_and_reset();
    *(logfiles[0]) << tm << " read_unlock_buffer_master"  << endl;
  }
}
void adjoint_read_unlock_buffer_slave(void)
{
  verify_identifier_string("K5");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_read_unlock_buffer_slave();
}
void  adpthread_manager::adjoint_read_unlock_buffer_slave(void)
{ 
  // cout  << "adpthread_manager::adjoint_read_unlock_buffer_slave(void)" << endl;
  // write_lock_slave_buffer
  verify_identifier_string("L7");
  int sno=restore_int_value();
  if (logflag)
  {
    *(logfiles[sno]) << "starting adjoint_read_unlock_buffer_slave " << sno  
                     << endl;
    adt[sno].get_elapsed_time_and_reset();
  }
  pthread_mutex_lock(smutex+sno);
  // only write if buffer empty
  int mf=mflag[sno];
  int sf=sflag[sno];
  while (mflag[sno] == 1 || sflag[sno] ==1 )
    pthread_cond_wait(scondition+sno,smutex+sno);
  verify_identifier_string("J5");
}


void adpthread_manager::read_unlock_buffer_slave(int sno)
{
  // we should be using the mutexes from the class
  // to get rid of global variables
//    save_identifier_string("J5");
  const char * str22;
  str22="J5";
  char* strx22=const_cast <char*> (str22);
  save_identifier_string(strx22);
  current_bptr[sno]=transfer_buffer[sno];
  pthread_mutex_unlock(smutex+sno) ;
  pthread_cond_signal(mcondition+sno);
  save_int_value(sno);
//    save_identifier_string("L7");
  const char * str23;
  str23="L7";
  char* strx23=const_cast <char*> (str23);
  save_identifier_string(strx23);
  save_pointer_value(this);
//    save_identifier_string("K5");
  const char * str24;
  str24="K5";
  char* strx24=const_cast <char*> (str24);
  save_identifier_string(strx24);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_read_unlock_buffer_slave);
  if (logflag)
  {
    MY_DOUBLE_TYPE tm=adt[sno].get_elapsed_time_and_reset();
    *(logfiles[sno]) << tm << " read_unlock_buffer_slave" << endl;
  }
}


void adpthread_manager::create_all(pthreadfun pf ,new_thread_data * ptr)
{
  for (int i=1;i<=nslaves;i++)
  {
    cerr << "Need to fix this" << endl;
    //pthread_create(thread1+i,NULL,pf,ptr+i);
  }
}

void adpthread_manager::create_all(pthreadfun pf ,thread_data * ptr)
{
  for (int i=1;i<=nslaves;i++)
  {
    cerr << "Need to fix this" << endl;
    //pthread_create(thread1+i,NULL,pf,ptr+i);
  }
}

void adpthread_manager::create_all(void * ptr)
{
  new_thread_data * dptr = (new_thread_data *)ptr; 
  int ii=0;
  for (int i=1;i<=ngroups;i++)
  {
    if (i>ngroups)
    {
      cerr << " index too high in "
        "adpthread_manager::create_all(void * dptr)" << endl;
      ad_exit(1);
    }
    for (int j=1;j<=num_in_group(i);j++)
    {
      cout << ppf[i] << endl;
      ++ii;
      cerr << "Need to fix this" << endl;
      //pthread_create(thread1+ii,NULL,ppf[i],dptr+ii);
    }
  }
}

void adpthread_manager::pthread_join_all(void)
{
  for (int i=1;i<=nslaves;i++)
  {
    pthread_join(thread1[i], NULL);
  }
}

void adpthread_manager::writebuffer(const void *x,int nbytes,int sno)
{
 /*
  if (old_buffer_flag==1)
  {
    check_buffer_size(nbytes,sno);
    ::memcpy(current_bptr[sno],x,nbytes);
    current_bptr[sno]+=nbytes;
  }
  else
 */
  {
    int s1,s2;
    int tn2=ad_comm::pthread_manager->get_slave_number();
    if (tn2==sno)
    {
      cerr << "This can't happen" << endl;
      ad_exit(1);
    }
    if (tn2<sno)
    {
      s1=tn2;
      s2=sno;
    }
    else
    {
      s1=sno;
      s2=tn2;
    }
    check_buffer_size(nbytes,s1,s2);
    ::memcpy(scurrent_bptr[s1][s2],x,nbytes);
    scurrent_bptr[s1][s2]+=nbytes;
  }
}

long int adpthread_manager::get_offset(int sno)
{
  int s1,s2;
  int tn2=ad_comm::pthread_manager->get_slave_number();
  if (tn2==sno)
  {
    cerr << "This can't happen" << endl;
    ad_exit(1);
  }
  if (tn2<sno)
  {
    s1=tn2;
    s2=sno;
  }
  else
  {
    s1=sno;
    s2=tn2;
  }
  std::ptrdiff_t pd=scurrent_bptr[s1][s2]-stransfer_buffer[s1][s2];
  return pd;
}

void adpthread_manager::readbuffer(const void *_x,int nbytes,int sno)
{
  void * x= (void *)(_x); 
 /*
  if (old_buffer_flag==1)
  {
    check_buffer_size_read(nbytes,sno);
    ::memcpy(x,current_bptr[sno],nbytes);
    current_bptr[sno]+=nbytes;
  }
  else
 */
  {
    int s1,s2;
    int tn2=ad_comm::pthread_manager->get_slave_number();
    if (tn2==sno)
    {
      cerr << "This can't happen" << endl;
      ad_exit(1);
    }
    if (tn2<sno)
    {
      s1=tn2;
      s2=sno;
    }
    else
    {
      s1=sno;
      s2=tn2;
    }
    check_buffer_size_read(nbytes,s1,s2);
    ::memcpy(x,scurrent_bptr[s1][s2],nbytes);
    scurrent_bptr[s1][s2]+=nbytes;
  }
}


void adpthread_manager::verify_id_string_from_master(const char * s,int sno)
{
#if defined(CHK_ID_STRING)
  char s1[10];
  int sz=strlen(s);
  s1[sz]='\0';
  readbuffer(s1,sz,sno);
  if (strcmp(s,s1))
  {
    cerr << "Error verifying master string " << s << endl;
    ad_exit(1);
  }
#endif
}

void adpthread_manager::verify_id_string_from_slave(const char * s,int sno)
{
#if defined(CHK_ID_STRING)
  char s1[10];
  int sz=strlen(s);
  s1[sz]='\0';
  readbuffer(s1,sz,sno);
  if (strcmp(s,s1))
  {
    cerr << "Error verifying slave string " << s <<  "  got " 
         << s1  << " instead" << endl;
    ad_exit(1);
  }
#endif
}

void adpthread_manager::send_double(const MY_DOUBLE_TYPE &x,int sno)
{
  send_id_string_to_slave("TY",sno);
  writebuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
}

void adpthread_manager::send_double_to_slave(const MY_DOUBLE_TYPE &x,int sno)
{
  send_id_string_to_slave("TY",sno);
  writebuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
}
void adpthread_manager::send_int_to_slave(int x,int sno)
{
  send_id_string_to_slave("RY",sno);
  writebuffer(&x,sizeof(int),sno);
}

void adpthread_manager::send_int(int x,int sno)
{
  send_id_string_to_slave("ES",sno);
  writebuffer(&x,sizeof(int),sno);
}

void adpthread_manager::send_int_to_master(int x,int sno)
{
  send_id_string_to_slave("YY",sno);
  writebuffer(&x,sizeof(int),sno);
}

int adpthread_manager::get_int_from_master(int sno)
{
  verify_id_string_from_master("RY",sno);
  int x;
  readbuffer(&x,sizeof(int),sno);
  return x;
}

int adpthread_manager::get_int(int sno)
{
  verify_id_string_from_master("ES",sno);
  int x;
  readbuffer(&x,sizeof(int),sno);
  return x;
}

int adpthread_manager::get_int_from_slave(int sno)
{
  verify_id_string_from_master("YY",sno);
  int x;
  readbuffer(&x,sizeof(int),sno);
  return x;
}

MY_DOUBLE_TYPE adpthread_manager::get_double_from_master(int sno)
{
  verify_id_string_from_master("TY",sno);
  MY_DOUBLE_TYPE x;
  readbuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
  return x;
}
MY_DOUBLE_TYPE adpthread_manager::get_double(int sno)
{
  verify_id_string_from_master("TY",sno);
  MY_DOUBLE_TYPE x;
  readbuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
  return x;
}

void adpthread_manager::send_double_to_master(const MY_DOUBLE_TYPE &x,int sno)
{
  send_id_string_to_master("EY",sno);
  writebuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
}

MY_DOUBLE_TYPE adpthread_manager::get_double_from_slave(int sno)
{
  verify_id_string_from_master("EY",sno);
  MY_DOUBLE_TYPE x;
  readbuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
  return x;
}

void adpthread_manager::send_id_string_to_slave(const char * x,int sno)
{
#if defined(CHK_ID_STRING)
  int sz=strlen(x);
  writebuffer(x,sz,sno);
#endif
}

void adpthread_manager::send_id_string_to_master(const char * x,int sno)
{
#if defined(CHK_ID_STRING)
  int sz=strlen(x);
  writebuffer(x,sz,sno);
#endif
}


void adjoint_get_dvar_vector_from_master(void)
{
  verify_identifier_string("D4");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_get_dvar_vector_from_master();
}

void adpthread_manager::adjoint_get_dvar_vector_from_master(void)
{
  verify_identifier_string("K6");
  int sno=restore_int_value();
  verify_identifier_string("Y7");
  dvar_vector_position dvpos=restore_dvar_vector_position();
  verify_identifier_string("C");
  send_id_string_to_master("SUX",sno);
  dvector dv=restore_dvar_vector_derivatives(dvpos);
  int mmin=dv.indexmin();
  int mmax=dv.indexmax();
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  int sz=mmax-mmin+1;
  writebuffer(&(dv(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
}
dvar_vector adpthread_manager::get_dvar_vector_from_master(int sno)
{
  int mmin;
  int mmax;
  verify_id_string_from_master("DCY",sno);
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  // cout  << "In dvar_vector get_dvar_vector_from_master " << endl;
  // cout  << " mmin = " << mmin   << " mmax = " << mmax  << endl;
  dvar_vector x(mmin,mmax);
  int sz=mmax-mmin+1;
  readbuffer(&(value(x(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("C");
  const char * str25;
  str25="C";
  char* strx25=const_cast <char*> (str25);
  save_identifier_string(strx25);
  x.save_dvar_vector_position();
//    save_identifier_string("Y7");
  const char * str26;
  str26="Y7";
  char* strx26=const_cast <char*> (str26);
  save_identifier_string(strx26);
  save_int_value(sno);
//    save_identifier_string("K6");
  const char * str27;
  str27="K6";
  char* strx27=const_cast <char*> (str27);
  save_identifier_string(strx27);
  save_pointer_value(this);
//    save_identifier_string("D4");
  const char * str28;
  str28="D4";
  char* strx28=const_cast <char*> (str28);
  save_identifier_string(strx28);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar_vector_from_master);
  return x;
}
dvar_vector adpthread_manager::get_dvar_vector(int sno)
{
  int mmin;
  int mmax;
  verify_id_string_from_master("DCY",sno);
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  // cout  << "In dvar_vector get_dvar_vector_from_master " << endl;
  // cout  << " mmin = " << mmin   << " mmax = " << mmax  << endl;
  dvar_vector x(mmin,mmax);
  int sz=mmax-mmin+1;
  readbuffer(&(value(x(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("C");
  const char * str29;
  str29="C";
  char* strx29=const_cast <char*> (str29);
  save_identifier_string(strx29);
  x.save_dvar_vector_position();
//    save_identifier_string("Y7");
  const char * str30;
  str30="Y7";
  char* strx30=const_cast <char*> (str30);
  save_identifier_string(strx30);
  save_int_value(sno);
//    save_identifier_string("K6");
  const char * str31;
  str31="K6";
  char* strx31=const_cast <char*> (str31);
  save_identifier_string(strx31);
  save_pointer_value(this);
//    save_identifier_string("D4");
  const char * str32;
  str32="D4";
  char* strx32=const_cast <char*> (str32);
  save_identifier_string(strx32);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar_vector_from_master);
  return x;
}

/*
dvariable adpthread_manager::get_dvariable_from_master(int sno)
{
  return ::get_dvariable_from_master(current_bptr[sno],sno);
}
*/

void adjoint_get_dvariable_from_slave(void)
{
  verify_identifier_string("G2");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_get_dvariable_from_slave();
}

void adpthread_manager::adjoint_get_dvariable_from_slave(void)
{
  verify_identifier_string("F7");
  int sno=restore_int_value();
  verify_identifier_string("D");
  prevariable_position dvpos=restore_prevariable_position();
  MY_DOUBLE_TYPE dv=restore_prevariable_derivative(dvpos);
  verify_identifier_string("C");
  send_id_string_to_slave("RUX",sno);
  writebuffer(&dv,sizeof(MY_DOUBLE_TYPE),sno);
}

dvariable adpthread_manager::get_dvariable_from_slave(int sno)
{
  dvariable x;
  MY_DOUBLE_TYPE cx;
  verify_id_string_from_slave("VP",sno);
  readbuffer(&cx,sizeof(MY_DOUBLE_TYPE),sno);
  value(x)=cx;
//    save_identifier_string("C");
  const char * str33;
  str33="C";
  char* strx33=const_cast <char*> (str33);
  save_identifier_string(strx33);
  x.save_prevariable_position();
//    save_identifier_string("D");
  const char * str34;
  str34="D";
  char* strx34=const_cast <char*> (str34);
  save_identifier_string(strx34);
  save_int_value(sno);
//    save_identifier_string("F7");
  const char * str35;
  str35="F7";
  char* strx35=const_cast <char*> (str35);
  save_identifier_string(strx35);
  save_pointer_value(this);
//    save_identifier_string("G2");
  const char * str36;
  str36="G2";
  char* strx36=const_cast <char*> (str36);
  save_identifier_string(strx36);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvariable_from_slave);
  return x;
}
dvariable adpthread_manager::get_dvariable(int sno)
{
  dvariable x;
  MY_DOUBLE_TYPE cx;
  verify_id_string_from_slave("VP",sno);
  readbuffer(&cx,sizeof(MY_DOUBLE_TYPE),sno);
  value(x)=cx;
//    save_identifier_string("C");
  const char * str37;
  str37="C";
  char* strx37=const_cast <char*> (str37);
  save_identifier_string(strx37);
  x.save_prevariable_position();
//    save_identifier_string("D");
  const char * str38;
  str38="D";
  char* strx38=const_cast <char*> (str38);
  save_identifier_string(strx38);
  save_int_value(sno);
//    save_identifier_string("F7");
  const char * str39;
  str39="F7";
  char* strx39=const_cast <char*> (str39);
  save_identifier_string(strx39);
  save_pointer_value(this);
//    save_identifier_string("G2");
  const char * str40;
  str40="G2";
  char* strx40=const_cast <char*> (str40);
  save_identifier_string(strx40);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvariable_from_slave);
  return x;
}
void adjoint_send_dvariable_to_master(void)
{
  verify_identifier_string("FH");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_send_dvariable_to_master();
}

void adpthread_manager::adjoint_send_dvariable_to_master(void)
{
  verify_identifier_string("KK");
  int sno=restore_int_value();
  verify_id_string_from_master("RUX",sno);
  verify_identifier_string("C");
  //double dv=restore_prevariable_derivative(dvpos);
  MY_DOUBLE_TYPE dv;
  readbuffer(&dv,sizeof(MY_DOUBLE_TYPE),sno);
  prevariable_position dvpos=restore_prevariable_position();
  save_double_derivative(dv,dvpos);
  verify_identifier_string("S5");
}

void adpthread_manager::send_dvariable_to_master(const prevariable& _x,int sno)
{
  ADUNCONST(prevariable,x)
//    save_identifier_string("S5");
  const char * str41;
  str41="S5";
  char* strx41=const_cast <char*> (str41);
  save_identifier_string(strx41);
  x.save_prevariable_position();
  send_id_string_to_master("VP",sno);
  MY_DOUBLE_TYPE cx=value(x);
  writebuffer(&cx,sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("C");
  const char * str42;
  str42="C";
  char* strx42=const_cast <char*> (str42);
  save_identifier_string(strx42);
  save_int_value(sno);
//    save_identifier_string("KK");
  const char * str43;
  str43="KK";
  char* strx43=const_cast <char*> (str43);
  save_identifier_string(strx43);
  save_pointer_value(this);
//    save_identifier_string("FH");
  const char * str44;
  str44="FH";
  char* strx44=const_cast <char*> (str44);
  save_identifier_string(strx44);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_send_dvariable_to_master);
}
void adpthread_manager::send_dvariable(const prevariable& _x,int sno)
{
  ADUNCONST(prevariable,x)
//    save_identifier_string("S5");
  const char * str45;
  str45="S5";
  char* strx45=const_cast <char*> (str45);
  save_identifier_string(strx45);
  x.save_prevariable_position();
  send_id_string_to_master("VP",sno);
  MY_DOUBLE_TYPE cx=value(x);
  writebuffer(&cx,sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("C");
  const char * str46;
  str46="C";
  char* strx46=const_cast <char*> (str46);
  save_identifier_string(strx46);
  save_int_value(sno);
//    save_identifier_string("KK");
  const char * str47;
  str47="KK";
  char* strx47=const_cast <char*> (str47);
  save_identifier_string(strx47);
  save_pointer_value(this);
//    save_identifier_string("FH");
  const char * str48;
  str48="FH";
  char* strx48=const_cast <char*> (str48);
  save_identifier_string(strx48);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_send_dvariable_to_master);
}

// ***********************************************
// ***********************************************
// ***********************************************
void adjoint_send_dvariable_to_slave(void)
{
  verify_identifier_string("RP");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_send_dvariable_to_slave();
}

void adpthread_manager::adjoint_send_dvariable_to_slave(void)
{
  verify_identifier_string("GH");
  int sno=restore_int_value();
  verify_identifier_string("UP");
  verify_id_string_from_slave("SUA",sno);
  MY_DOUBLE_TYPE x;
  readbuffer(&x,sizeof(MY_DOUBLE_TYPE),sno);
  prevariable_position dxpos=restore_prevariable_position();
  verify_identifier_string("B");
  save_double_derivative(x,dxpos);
}

void adpthread_manager::send_dvariable_to_slave(const prevariable &x,int sno)
{
  writebuffer(&(value(x)),sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("B");
  const char * str49;
  str49="B";
  char* strx49=const_cast <char*> (str49);
  save_identifier_string(strx49);
  x.save_prevariable_position();
//    save_identifier_string("UP");
  const char * str50;
  str50="UP";
  char* strx50=const_cast <char*> (str50);
  save_identifier_string(strx50);
  save_int_value(sno);
//    save_identifier_string("GH");
  const char * str51;
  str51="GH";
  char* strx51=const_cast <char*> (str51);
  save_identifier_string(strx51);
  save_pointer_value(this);
//    save_identifier_string("RP");
  const char * str52;
  str52="RP";
  char* strx52=const_cast <char*> (str52);
  save_identifier_string(strx52);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvariable_to_slave);
}

// ***********************************************
// ***********************************************

void adjoint_send_dvar_vector_to_slave(void)
{
  verify_identifier_string("UP");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_send_dvar_vector_to_slave();
}

void adpthread_manager::adjoint_send_dvar_vector_to_slave(void)
{
  verify_identifier_string("HH");
  int sno=restore_int_value();
  verify_identifier_string("OP");
  verify_id_string_from_slave("SUX",sno);
  int mmin;
  int mmax;
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  int sz=mmax-mmin+1;
  dvector x(mmin,mmax);
  readbuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  dvar_vector_position dvpos=restore_dvar_vector_position();
  verify_identifier_string("T");
  x.save_dvector_derivatives(dvpos);
}

void adpthread_manager::send_dvar_vector_to_slave(const dvar_vector &x,int sno)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int sz=mmax-mmin+1;
  send_id_string_to_slave("DCY",sno);
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  writebuffer(&(value(x(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("T");
  const char * str53;
  str53="T";
  char* strx53=const_cast <char*> (str53);
  save_identifier_string(strx53);
  x.save_dvar_vector_position();
//    save_identifier_string("OP");
  const char * str54;
  str54="OP";
  char* strx54=const_cast <char*> (str54);
  save_identifier_string(strx54);
  save_int_value(sno);
//    save_identifier_string("HH");
  const char * str55;
  str55="HH";
  char* strx55=const_cast <char*> (str55);
  save_identifier_string(strx55);
  save_pointer_value(this);
//    save_identifier_string("UP");
  const char * str56;
  str56="UP";
  char* strx56=const_cast <char*> (str56);
  save_identifier_string(strx56);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar_vector_to_slave);
    //set_gradient_stack(pthread_master_unpack_vector_derivatives);
}
void adpthread_manager::send_dvar_vector(const dvar_vector &x,int sno)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int sz=mmax-mmin+1;
  send_id_string_to_slave("DCY",sno);
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  writebuffer(&(value(x(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
//    save_identifier_string("T");
  const char * str57;
  str57="T";
  char* strx57=const_cast <char*> (str57);
  save_identifier_string(strx57);
  x.save_dvar_vector_position();
//    save_identifier_string("OP");
  const char * str58;
  str58="OP";
  char* strx58=const_cast <char*> (str58);
  save_identifier_string(strx58);
  save_int_value(sno);
//    save_identifier_string("HH");
  const char * str59;
  str59="HH";
  char* strx59=const_cast <char*> (str59);
  save_identifier_string(strx59);
  save_pointer_value(this);
//    save_identifier_string("UP");
  const char * str60;
  str60="UP";
  char* strx60=const_cast <char*> (str60);
  save_identifier_string(strx60);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar_vector_to_slave);
    //set_gradient_stack(pthread_master_unpack_vector_derivatives);
}

void adpthread_manager::send_dvector_to_slave(const dvector &x,int sno)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int sz=mmax-mmin+1;
  send_id_string_to_slave("CCX",sno);
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  writebuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  //::send_dvector_to_slave(x,current_bptr[sno],sno);
}

void adpthread_manager::send_dvector(const dvector &x,int sno)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int sz=mmax-mmin+1;
  send_id_string_to_slave("CCX",sno);
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  writebuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  //::send_dvector_to_slave(x,current_bptr[sno],sno);
}

void adpthread_manager::send_dvector_to_master(const dvector &x,int sno)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int sz=mmax-mmin+1;
  send_id_string_to_master("TYU",sno);
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  writebuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
}
dvector adpthread_manager::get_dvector_from_slave(int sno)
{
  int mmin;
  int mmax;
  verify_id_string_from_slave("TYU",sno);
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  dvector x(mmin,mmax);
  int sz=mmax-mmin+1;
  readbuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  return x;
}

dvector adpthread_manager::get_dvector_from_master(int sno)
{
  int mmin;
  int mmax;
  verify_id_string_from_master("CCX",sno);
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  dvector x(mmin,mmax);
  int sz=mmax-mmin+1;
  readbuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  return x;
}
dvector adpthread_manager::get_dvector(int sno)
{
  int mmin;
  int mmax;
  verify_id_string_from_master("CCX",sno);
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  dvector x(mmin,mmax);
  int sz=mmax-mmin+1;
  readbuffer(&(x(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  return x;
}
adpthread_manager * test_thread_manager=0;

void adjoint_send_dvar_matrix_to_slave(void)
{
  verify_identifier_string("UN");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_send_dvar_matrix_to_slave();
}

void adpthread_manager::adjoint_send_dvar_matrix_to_slave(void)
{
  verify_identifier_string("HH");
  int sno=restore_int_value();
  verify_identifier_string("OP");
  verify_id_string_from_slave("WVX",sno);
  int rmin;
  int rmax;
  readbuffer(&rmin,sizeof(int),sno);
  readbuffer(&rmax,sizeof(int),sno);
  dmatrix M(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin,mmax;
    readbuffer(&mmin,sizeof(int),sno);
    readbuffer(&mmax,sizeof(int),sno);
    M(i).allocate(mmin,mmax);
    int sz=mmax-mmin+1;
    readbuffer(&(M(i)(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
  dvar_matrix_position dmpos=restore_dvar_matrix_position();
  verify_identifier_string("Y");
  M.save_dmatrix_derivatives(dmpos);
}

void adpthread_manager::send_dvar_matrix_to_slave(const dvar_matrix &x,int sno)
{
  int rmin=x.indexmin();
  int rmax=x.indexmax();
  send_id_string_to_slave("WTZ",sno);
  writebuffer(&rmin,sizeof(int),sno);
  writebuffer(&rmax,sizeof(int),sno);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin=x(i).indexmin();
    int mmax=x(i).indexmax();
    int sz=mmax-mmin+1;
    writebuffer(&mmin,sizeof(int),sno);
    writebuffer(&mmax,sizeof(int),sno);
    writebuffer(&(value(x(i)(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
//    save_identifier_string("Y");
  const char * str61;
  str61="Y";
  char* strx61=const_cast <char*> (str61);
  save_identifier_string(strx61);
  // !!! should we optimize this ?
  x.save_dvar_matrix_position();
//    save_identifier_string("OP");
  const char * str62;
  str62="OP";
  char* strx62=const_cast <char*> (str62);
  save_identifier_string(strx62);
  save_int_value(sno);
//    save_identifier_string("HH");
  const char * str63;
  str63="HH";
  char* strx63=const_cast <char*> (str63);
  save_identifier_string(strx63);
  save_pointer_value(this);
//    save_identifier_string("UN");
  const char * str64;
  str64="UN";
  char* strx64=const_cast <char*> (str64);
  save_identifier_string(strx64);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar_matrix_to_slave);
}

void adpthread_manager::send_dvar_matrix(const dvar_matrix &x,int sno)
{
  int rmin=x.indexmin();
  int rmax=x.indexmax();
  send_id_string_to_slave("RAZ",sno);
  writebuffer(&rmin,sizeof(int),sno);
  writebuffer(&rmax,sizeof(int),sno);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin=x(i).indexmin();
    int mmax=x(i).indexmax();
    int sz=mmax-mmin+1;
    writebuffer(&mmin,sizeof(int),sno);
    writebuffer(&mmax,sizeof(int),sno);
    writebuffer(&(value(x(i)(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
//    save_identifier_string("Y");
  const char * str65;
  str65="Y";
  char* strx65=const_cast <char*> (str65);
  save_identifier_string(strx65);
  // !!! should we optimize this ?
  x.save_dvar_matrix_position();
//    save_identifier_string("OP");
  const char * str66;
  str66="OP";
  char* strx66=const_cast <char*> (str66);
  save_identifier_string(strx66);
  save_int_value(sno);
//    save_identifier_string("HH");
  const char * str67;
  str67="HH";
  char* strx67=const_cast <char*> (str67);
  save_identifier_string(strx67);
  save_pointer_value(this);
//    save_identifier_string("UN");
  const char * str68;
  str68="UN";
  char* strx68=const_cast <char*> (str68);
  save_identifier_string(strx68);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar_matrix_to_slave);
}
void adpthread_manager::send_dmatrix(const dmatrix &x,int sno)
{
  int rmin=x.indexmin();
  int rmax=x.indexmax();
  send_id_string_to_slave("SAY",sno);
  writebuffer(&rmin,sizeof(int),sno);
  writebuffer(&rmax,sizeof(int),sno);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin=x(i).indexmin();
    int mmax=x(i).indexmax();
    int sz=mmax-mmin+1;
    writebuffer(&mmin,sizeof(int),sno);
    writebuffer(&mmax,sizeof(int),sno);
    writebuffer(&(x(i)(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
}
void adjoint_get_dvar_matrix_from_master(void)
{
  verify_identifier_string("T4");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_get_dvar_matrix_from_master();
}

void adpthread_manager::adjoint_get_dvar_matrix_from_master(void)
{
  verify_identifier_string("K6");
  int sno=restore_int_value();
  verify_identifier_string("Y7");
  dvar_matrix_position dvpos=restore_dvar_matrix_position();
  verify_identifier_string("C2");
  dmatrix dv=restore_dvar_matrix_derivatives(dvpos);
  send_id_string_to_master("WVX",sno);
  int rmin=dv.indexmin();
  int rmax=dv.indexmax();
  writebuffer(&rmin,sizeof(int),sno);
  writebuffer(&rmax,sizeof(int),sno);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin=dv(i).indexmin();
    int mmax=dv(i).indexmax();
    int sz=mmax-mmin+1;
    writebuffer(&mmin,sizeof(int),sno);
    writebuffer(&mmax,sizeof(int),sno);
    writebuffer(&(dv(i)(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
}
dvar_matrix adpthread_manager::get_dvar_matrix_from_master(int sno)
{
  verify_id_string_from_master("WTZ",sno);
  int rmin;
  int rmax;
  readbuffer(&rmin,sizeof(int),sno);
  readbuffer(&rmax,sizeof(int),sno);
  // cout  << "In dvar_vector get_dvar_vector_from_master " << endl;
  // cout  << " mmin = " << mmin   << " mmax = " << mmax  << endl;
  dvar_matrix x(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin;
    int mmax;
    readbuffer(&mmin,sizeof(int),sno);
    readbuffer(&mmax,sizeof(int),sno);
    int sz=mmax-mmin+1;
    x(i).allocate(mmin,mmax);
    readbuffer(&(value(x(i)(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
//    save_identifier_string("C2");
  const char * str69;
  str69="C2";
  char* strx69=const_cast <char*> (str69);
  save_identifier_string(strx69);
  x.save_dvar_matrix_position();
//    save_identifier_string("Y7");
  const char * str70;
  str70="Y7";
  char* strx70=const_cast <char*> (str70);
  save_identifier_string(strx70);
  save_int_value(sno);
//    save_identifier_string("K6");
  const char * str71;
  str71="K6";
  char* strx71=const_cast <char*> (str71);
  save_identifier_string(strx71);
  save_pointer_value(this);
//    save_identifier_string("T4");
  const char * str72;
  str72="T4";
  char* strx72=const_cast <char*> (str72);
  save_identifier_string(strx72);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar_matrix_from_master);
  return x;
}
dvar_matrix adpthread_manager::get_dvar_matrix(int sno)
{
  verify_id_string_from_master("RAZ",sno);
  int rmin;
  int rmax;
  readbuffer(&rmin,sizeof(int),sno);
  readbuffer(&rmax,sizeof(int),sno);
  // cout  << "In dvar_vector get_dvar_vector_from_master " << endl;
  // cout  << " mmin = " << mmin   << " mmax = " << mmax  << endl;
  dvar_matrix x(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin;
    int mmax;
    readbuffer(&mmin,sizeof(int),sno);
    readbuffer(&mmax,sizeof(int),sno);
    int sz=mmax-mmin+1;
    x(i).allocate(mmin,mmax);
    readbuffer(&(value(x(i)(mmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
//    save_identifier_string("C2");
  const char * str73;
  str73="C2";
  char* strx73=const_cast <char*> (str73);
  save_identifier_string(strx73);
  x.save_dvar_matrix_position();
//    save_identifier_string("Y7");
  const char * str74;
  str74="Y7";
  char* strx74=const_cast <char*> (str74);
  save_identifier_string(strx74);
  save_int_value(sno);
//    save_identifier_string("K6");
  const char * str75;
  str75="K6";
  char* strx75=const_cast <char*> (str75);
  save_identifier_string(strx75);
  save_pointer_value(this);
//    save_identifier_string("T4");
  const char * str76;
  str76="T4";
  char* strx76=const_cast <char*> (str76);
  save_identifier_string(strx76);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar_matrix_from_master);
  return x;
}

dmatrix adpthread_manager::get_dmatrix(int sno)
{
  verify_id_string_from_master("SAY",sno);
  int rmin;
  int rmax;
  readbuffer(&rmin,sizeof(int),sno);
  readbuffer(&rmax,sizeof(int),sno);
  dmatrix x(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin;
    int mmax;
    readbuffer(&mmin,sizeof(int),sno);
    readbuffer(&mmax,sizeof(int),sno);
    int sz=mmax-mmin+1;
    x(i).allocate(mmin,mmax);
    readbuffer(&(x(i)(mmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
  }
  return x;
}

void add_slave_suffix(const adstring& _tmpstring)
{
  ADUNCONST(adstring,tmpstring)
  if (ad_comm::pthread_manager)
  {
    pthread_mutex_lock(&ad_comm::pthread_manager->copy_mutex);
    if (ad_comm::pthread_manager->is_slave())
    {
      tmpstring += "_slave_";
      tmpstring += str(ad_comm::pthread_manager->get_slave_number());
       //cout << "In slave " << tmpstring << endl;
    }
    else
    {
      tmpstring += "_master";
       //cout << "In master " << tmpstring << endl;
    }
    pthread_mutex_unlock(&ad_comm::pthread_manager->copy_mutex);
  }
}
new_thread_data::new_thread_data(void)
{
 // id=0;
  thread_no=0;
  pfsh=0;
  m=0;
}

imatrix adpthread_manager::get_imatrix(int sno)
{
  verify_id_string_from_master("FKY",sno);
  int rmin;
  int rmax;
  readbuffer(&rmin,sizeof(int),sno);
  readbuffer(&rmax,sizeof(int),sno);
  imatrix x(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin;
    int mmax;
    readbuffer(&mmin,sizeof(int),sno);
    readbuffer(&mmax,sizeof(int),sno);
    int sz=mmax-mmin+1;
    x(i).allocate(mmin,mmax);
    readbuffer(&(x(i)(mmin)),sz*sizeof(int),sno);
  }
  return x;
}

void adpthread_manager::send_imatrix(const imatrix &x,int sno)
{
  int rmin=x.indexmin();
  int rmax=x.indexmax();
  send_id_string_to_slave("FKY",sno);
  writebuffer(&rmin,sizeof(int),sno);
  writebuffer(&rmax,sizeof(int),sno);
  for (int i=rmin;i<=rmax;i++)
  {
    int mmin=x(i).indexmin();
    int mmax=x(i).indexmax();
    int sz=mmax-mmin+1;
    writebuffer(&mmin,sizeof(int),sno);
    writebuffer(&mmax,sizeof(int),sno);
    writebuffer(&(x(i)(mmin)),sz*sizeof(int),sno);
  }
}
ivector adpthread_manager::get_ivector(int sno)
{
  int mmin;
  int mmax;
  verify_id_string_from_master("EAC",sno);
  readbuffer(&mmin,sizeof(int),sno);
  readbuffer(&mmax,sizeof(int),sno);
  ivector x(mmin,mmax);
  int sz=mmax-mmin+1;
  readbuffer(&(x(mmin)),sz*sizeof(int),sno);
  return x;
}
void adpthread_manager::send_ivector(const ivector &x,int sno)
{
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int sz=mmax-mmin+1;
  send_id_string_to_slave("EAC",sno);
  writebuffer(&mmin,sizeof(int),sno);
  writebuffer(&mmax,sizeof(int),sno);
  writebuffer(&(x(mmin)),sz*sizeof(int),sno);
}
