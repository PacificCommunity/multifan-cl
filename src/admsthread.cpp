/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE



#include <stdio.h>
#include <admodel.h>

#include <pthread.h>
#include "admsthread.h"

void ad_master_slave_synchro::master_lock(void)  
{
 /*
  if (flag !=0)
  {
    cout << "this can't happen" << endl;
    ad_exit(1);
  }
 */
  pthread_mutex_lock(&mutex);
}

void ad_master_slave_synchro::slave_lock(void)  
{
  if (flag !=1)
  {
    cout << "this can't happen" << endl;
    ad_exit(1);
  }
  pthread_mutex_lock(&mutex);
}

void ad_master_slave_synchro::slave_get_signal(int i)  
{
  pthread_mutex_lock(&mutex);
  if (flag==0)
  {
    cout << "master not ready for slave " << i << endl;
    pthread_cond_wait(&cond1,&mutex);
  }
  pthread_mutex_unlock(&mutex);
}

void ad_master_slave_synchro::master_get_signal(void)  
{
  pthread_mutex_lock(&mutex);
  if (flag==1)
  {
    cout << "slave not ready"  << endl;
    pthread_cond_wait(&cond2,&mutex);
  }
  pthread_mutex_unlock(&mutex);
}

void ad_master_slave_synchro::master_send_signal(void)  
{
  flag=1;
  pthread_cond_signal(&cond1);
  pthread_mutex_unlock(&mutex);
}

void ad_master_slave_synchro::slave_send_signal(void)  
{
  flag=0;
  pthread_cond_signal(&cond2);
  pthread_mutex_unlock(&mutex);
}

ad_master_slave_synchro::ad_master_slave_synchro(void)  
{
  threadid=-1;  // not set yet
  flag=0;
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init (&cond1, NULL);
  pthread_cond_init (&cond2, NULL);
}
