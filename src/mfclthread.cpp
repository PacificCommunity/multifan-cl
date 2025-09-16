/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#ifdef USE_PTHREADS
#include <stdio.h>
#include "all.hpp"
#include "admsthread.h"

//#define NUM_THREADS     25
static const int NUM_THREADS=1;


class threadinfo
{
public:
  int threadid;
  int nvar;
  dvector &x;
  dvar_len_fish_stock_history& lfsh;
  static MY_DOUBLE_TYPE tsum;
  pthread_mutex_t * pmutex;
  threadinfo(dvar_len_fish_stock_history& lfsh,pthread_mutex_t *,
    int nvar,dvector& x);
};

MY_DOUBLE_TYPE threadinfo::tsum=0;

threadinfo::threadinfo(dvar_len_fish_stock_history&  _lfsh,
  pthread_mutex_t * _pmutex,int _nvar,dvector& _x) : lfsh(_lfsh), 
  pmutex(_pmutex),nvar(_nvar),x(_x) 
{
  threadid=0;
}

void * fit_catchability_devs(void * t)
{
  threadinfo * tinfo = (threadinfo*) (t);
  int threadid=tinfo->threadid;
  dvar_len_fish_stock_history& lfsh=tinfo->lfsh;
  dvector& x=tinfo->x;
  int nvar=tinfo->nvar;
  //lfsh.grouped_implicit_catchability_deviations_calc_part1(nvar,x);
   printf("\n%d: Hello World!\n", threadid);
   pthread_mutex_lock (tinfo->pmutex+threadid);
   tinfo->tsum=tinfo->tsum+threadid;
   pthread_mutex_unlock (tinfo->pmutex+threadid);
   pthread_exit(NULL);
   return (void*)(0);
}
//ad_master_slave_synchro psynch;

void * start_catchdev_part1(void * t)
{
  adtimer mytimer;
  int * p =new int;
  *p=2;
  int rc = pthread_setspecific(admb_pthread_key, (void*)(p));

  do
  {
    MY_DOUBLE_TYPE tt=mytimer.get_elapsed_time();
    //psynch.slave_get_signal(1);  
    dvar_fish_stock_history * pfsh = (dvar_fish_stock_history*)(t);
    // *******************************************8
    // *******************************************8
    // need to call part1 here
    pfsh->grouped_implicit_catchability_deviations_calc_part1_thread();
  
    // *******************************************8
    // *******************************************8
    //psynch.slave_send_signal();  
    MY_DOUBLE_TYPE tt1=mytimer.get_elapsed_time();
    cout << "slave before = " << tt << "  after = " << tt1 << endl;
  }
  while(1);

  delete p;
  pthread_exit(NULL);
}

#if defined(USE_ADPVM)
void start_catchdev_part1(dvar_fish_stock_history & fsh )
{
  adtimer mytimer;
  //int * p =new int;
  //*p=2;
  //int rc = pthread_setspecific(admb_pthread_key, (void*)(p));

  do
  {
    MY_DOUBLE_TYPE tt=mytimer.get_elapsed_time();
    //psynch.slave_get_signal(1);  
    // *******************************************8
    // *******************************************8
    // need to call part1 here
    fsh.grouped_implicit_catchability_deviations_calc_part1_pvm();
  
    // *******************************************8
    // *******************************************8
    //psynch.slave_send_signal();  
    MY_DOUBLE_TYPE tt1=mytimer.get_elapsed_time();
    cout << "slave before = " << tt << "  after = " << tt1 << endl;
  }
  while(1);

  //delete p;
  //pthread_exit(NULL);
}
#endif  //#if defined(USE_ADPVM)

void pthreads_master_send_signal_to_slave(void)
{
  //psynch.master_send_signal();  
}


void dvar_len_fish_stock_history::create_catchability_deviation_thread(void)
{
 pthread_t threads[NUM_THREADS];
 int rc;
  
  
  printf("Creating thread \n");
  rc = pthread_create(&threads[0], NULL, start_catchdev_part1, 
     (void *)this);
  if (rc){
     printf("ERROR; return code from pthread_create() is %d\n", rc);
     exit(-1);
  }
}

  // void dvar_len_fish_stock_history::create_catchability_deviation_thread(void)
  // {
  //    int rc, status;
  //    threadinfo * pt;
  //    pthread_mutex_t pmutex;
  //    pthread_mutex_init(&pmutex, NULL);
  // 
  // /* Initialize and set thread detached attribute */
  //    pthread_t threads;
  //    pthread_attr_t attr;
  //    pthread_attr_init(&attr);
  //    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  // 
  //    int n=sum(num_real_grouped_fish_times);
  //    dvector x(1,n);
  //    pt = new threadinfo(*this,pmutex,n,x);
  //    rc = pthread_create(&threads,&attr, fit_catchability_devs, (void *)(pt));
  //    _sleep(1);
  //    if (rc)
  //    {
  //      printf("ERROR; return code from pthread_create() is %d\n", rc);
  //      //exit(-1);
  //    }
  // /* Free attribute and wait for the other threads */
  //   
  //    pthread_attr_destroy(&attr);
  //   /*  don't want this
  //    rc = pthread_join(threads, (void **)&status);
  //    if (rc)
  //    {
  //       printf("ERROR; return code from pthread_join() is %d\n", rc);
  //       exit(-1);
  //    }
  //    printf("Completed join with thread  status= %d\n", status);
  //   */
  //   
  //    
  //    cout <<"mfclthread.cpp " << threadinfo::tsum << endl;
  //    pthread_exit(NULL);
  // }


dvariable dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part2_thread_with_signal(void)
{
  //psynch.master_get_signal();  
  dvariable vf=
    grouped_implicit_catchability_deviations_calc_part2_thread();
  return vf;
}

#endif //USE_PTHREADS
