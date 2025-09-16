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


void adjoint_send_dvar3_array_to_master(void)
{
  verify_identifier_string("UN");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_send_dvar3_array_to_master();
}

void adpthread_manager::adjoint_send_dvar3_array_to_master(void)
{
  verify_identifier_string("HH");
  int sno=restore_int_value();
  verify_id_string_from_master("JKH",sno);
  int smin;
  int smax;
  int i;
  readbuffer(&smin,sizeof(int),sno);
  readbuffer(&smax,sizeof(int),sno);
  d3_array M(smin,smax);
  for (i=smin;i<=smax;i++)
  {
    int rmin,rmax;
    readbuffer(&rmin,sizeof(int),sno);
    readbuffer(&rmax,sizeof(int),sno);
    M(i).allocate(rmin,rmax);
    for (int j=rmin;j<=rmax;j++)
    { 
      int cmin,cmax;
      readbuffer(&cmin,sizeof(int),sno);
      readbuffer(&cmax,sizeof(int),sno);
      M(i,j).allocate(cmin,cmax);
      int sz=cmax-cmin+1;
      readbuffer(&(M(i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }
  verify_identifier_string("OP");
  for (int i=smax;i>=smin;i--)
  {
    dvar_matrix_position dmpos=restore_dvar_matrix_position();
    M(i).save_dmatrix_derivatives(dmpos);
  }
  verify_identifier_string("TB");
}

void adpthread_manager::send_dvar3_array_to_master(const dvar3_array &x,int sno)
{
  send_id_string_to_master("FTZ",sno);
  int smin=x.indexmin();
  int smax=x.indexmax();
  int i;
  writebuffer(&smin,sizeof(int),sno);
  writebuffer(&smax,sizeof(int),sno);
  for (i=smin;i<=smax;i++)
  {
    int rmin=x(i).indexmin();
    int rmax=x(i).indexmax();
    writebuffer(&rmin,sizeof(int),sno);
    writebuffer(&rmax,sizeof(int),sno);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin=x(i,j).indexmin();
      int cmax=x(i,j).indexmax();
      writebuffer(&cmin,sizeof(int),sno);
      writebuffer(&cmax,sizeof(int),sno);
      int sz=cmax-cmin+1;
      writebuffer(&(value(x(i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }

//    save_identifier_string("TB");
  const char * str1;
  str1="TB";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  // !!! should we optimize this ?
  for (i=smin;i<=smax;i++)
  {
    x(i).save_dvar_matrix_position();
  }
//    save_identifier_string("OP");
  const char * str2;
  str2="OP";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  save_int_value(sno);
//    save_identifier_string("HH");
  const char * str3;
  str3="HH";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  save_pointer_value(this);
//    save_identifier_string("UN");
  const char * str4;
  str4="UN";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar3_array_to_master);
}

void adpthread_manager::send_dvar3_array(const dvar3_array &x,int sno)
{
  send_id_string_to_master("FTZ",sno);
  int smin=x.indexmin();
  int smax=x.indexmax();
  int i;
  writebuffer(&smin,sizeof(int),sno);
  writebuffer(&smax,sizeof(int),sno);
  for (i=smin;i<=smax;i++)
  {
    int rmin=x(i).indexmin();
    int rmax=x(i).indexmax();
    writebuffer(&rmin,sizeof(int),sno);
    writebuffer(&rmax,sizeof(int),sno);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin=x(i,j).indexmin();
      int cmax=x(i,j).indexmax();
      writebuffer(&cmin,sizeof(int),sno);
      writebuffer(&cmax,sizeof(int),sno);
      int sz=cmax-cmin+1;
      writebuffer(&(value(x(i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }

//    save_identifier_string("TB");
  const char * str5;
  str5="TB";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  // !!! should we optimize this ?
  for (i=smin;i<=smax;i++)
  {
    x(i).save_dvar_matrix_position();
  }
//    save_identifier_string("OP");
  const char * str6;
  str6="OP";
  char* strx6=const_cast <char*> (str6);
  save_identifier_string(strx6);
  save_int_value(sno);
//    save_identifier_string("HH");
  const char * str7;
  str7="HH";
  char* strx7=const_cast <char*> (str7);
  save_identifier_string(strx7);
  save_pointer_value(this);
//    save_identifier_string("UN");
  const char * str8;
  str8="UN";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar3_array_to_master);
}

void adpthread_manager::send_d3_array(const d3_array &x,int sno)
{
  send_id_string_to_master("GRH",sno);
  int smin=x.indexmin();
  int smax=x.indexmax();
  int i;
  writebuffer(&smin,sizeof(int),sno);
  writebuffer(&smax,sizeof(int),sno);
  for (i=smin;i<=smax;i++)
  {
    int rmin=x(i).indexmin();
    int rmax=x(i).indexmax();
    writebuffer(&rmin,sizeof(int),sno);
    writebuffer(&rmax,sizeof(int),sno);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin=x(i,j).indexmin();
      int cmax=x(i,j).indexmax();
      writebuffer(&cmin,sizeof(int),sno);
      writebuffer(&cmax,sizeof(int),sno);
      int sz=cmax-cmin+1;
      writebuffer(&(x(i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }
}

void adjoint_get_dvar3_array_from_slave(void)
{
  verify_identifier_string("C5");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_get_dvar3_array_from_slave();
}

void adpthread_manager::adjoint_get_dvar3_array_from_slave(void)
{
  verify_identifier_string("LG");
  int sno=restore_int_value();
  verify_identifier_string("D4");
  int smin=restore_int_value();
  int smax=restore_int_value();
  d3_array dv(smin,smax);
  for (int i=smax;i>=smin;i--)
  {
    dvar_matrix_position dvpos=restore_dvar_matrix_position();
    dv(i)=restore_dvar_matrix_derivatives(dvpos);
  }
  verify_identifier_string("A3");
  send_id_string_to_slave("JKH",sno);
  int i;
  writebuffer(&smin,sizeof(int),sno);
  writebuffer(&smax,sizeof(int),sno);
  for (i=smin;i<=smax;i++)
  {
    int rmin=dv(i).indexmin();
    int rmax=dv(i).indexmax();
    writebuffer(&rmin,sizeof(int),sno);
    writebuffer(&rmax,sizeof(int),sno);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin=dv(i,j).indexmin();
      int cmax=dv(i,j).indexmax();
      writebuffer(&cmin,sizeof(int),sno);
      writebuffer(&cmax,sizeof(int),sno);
      int sz=cmax-cmin+1;
      writebuffer(&(dv(i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }
}
dvar3_array adpthread_manager::get_dvar3_array_from_slave(int sno)
{
  verify_id_string_from_slave("FTZ",sno);
   //*************************************************
   //*************************************************

  int smin;
  int smax;
  int i;
  readbuffer(&smin,sizeof(int),sno);
  readbuffer(&smax,sizeof(int),sno);
  dvar3_array M(smin,smax);
  for (i=smin;i<=smax;i++)
  {
    int rmin,rmax;
    readbuffer(&rmin,sizeof(int),sno);
    readbuffer(&rmax,sizeof(int),sno);
    M(i).allocate(rmin,rmax);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin,cmax;
      readbuffer(&cmin,sizeof(int),sno);
      readbuffer(&cmax,sizeof(int),sno);
      M(i,j).allocate(cmin,cmax);
      int sz=cmax-cmin+1;
      readbuffer(&(value(M(i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }

//    save_identifier_string("A3");
  const char * str9;
  str9="A3";
  char* strx9=const_cast <char*> (str9);
  save_identifier_string(strx9);
  for (i=smin;i<=smax;i++)
  {
    M(i).save_dvar_matrix_position();
  }
  save_int_value(smax);
  save_int_value(smin);
//    save_identifier_string("D4");
  const char * str10;
  str10="D4";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
  save_int_value(sno);
//    save_identifier_string("LG");
  const char * str11;
  str11="LG";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
  save_pointer_value(this);
//    save_identifier_string("C5");
  const char * str12;
  str12="C5";
  char* strx12=const_cast <char*> (str12);
  save_identifier_string(strx12);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar3_array_from_slave);
  return M;
}

dvar3_array adpthread_manager::get_dvar3_array(int sno)
{
  verify_id_string_from_slave("FTZ",sno);
   //*************************************************
   //*************************************************

  int smin;
  int smax;
  int i;
  readbuffer(&smin,sizeof(int),sno);
  readbuffer(&smax,sizeof(int),sno);
  dvar3_array M(smin,smax);
  for (i=smin;i<=smax;i++)
  {
    int rmin,rmax;
    readbuffer(&rmin,sizeof(int),sno);
    readbuffer(&rmax,sizeof(int),sno);
    M(i).allocate(rmin,rmax);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin,cmax;
      readbuffer(&cmin,sizeof(int),sno);
      readbuffer(&cmax,sizeof(int),sno);
      M(i,j).allocate(cmin,cmax);
      int sz=cmax-cmin+1;
      readbuffer(&(value(M(i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }

//    save_identifier_string("A3");
  const char * str13;
  str13="A3";
  char* strx13=const_cast <char*> (str13);
  save_identifier_string(strx13);
  for (i=smin;i<=smax;i++)
  {
    M(i).save_dvar_matrix_position();
  }
  save_int_value(smax);
  save_int_value(smin);
//    save_identifier_string("D4");
  const char * str14;
  str14="D4";
  char* strx14=const_cast <char*> (str14);
  save_identifier_string(strx14);
  save_int_value(sno);
//    save_identifier_string("LG");
  const char * str15;
  str15="LG";
  char* strx15=const_cast <char*> (str15);
  save_identifier_string(strx15);
  save_pointer_value(this);
//    save_identifier_string("C5");
  const char * str16;
  str16="C5";
  char* strx16=const_cast <char*> (str16);
  save_identifier_string(strx16);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar3_array_from_slave);
  return M;
}

d3_array adpthread_manager::get_d3_array(int sno)
{
  verify_id_string_from_slave("GRH",sno);

  int smin;
  int smax;
  int i;
  readbuffer(&smin,sizeof(int),sno);
  readbuffer(&smax,sizeof(int),sno);
  d3_array M(smin,smax);
  for (i=smin;i<=smax;i++)
  {
    int rmin,rmax;
    readbuffer(&rmin,sizeof(int),sno);
    readbuffer(&rmax,sizeof(int),sno);
    M(i).allocate(rmin,rmax);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin,cmax;
      readbuffer(&cmin,sizeof(int),sno);
      readbuffer(&cmax,sizeof(int),sno);
      M(i,j).allocate(cmin,cmax);
      int sz=cmax-cmin+1;
      readbuffer(&(M(i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
    }
  }
  return M;
}

i3_array adpthread_manager::get_i3_array(int sno)
{
  verify_id_string_from_slave("EWS",sno);

  int smin;
  int smax;
  int i;
  readbuffer(&smin,sizeof(int),sno);
  readbuffer(&smax,sizeof(int),sno);
  i3_array M(smin,smax);
  for (i=smin;i<=smax;i++)
  {
    int rmin,rmax;
    readbuffer(&rmin,sizeof(int),sno);
    readbuffer(&rmax,sizeof(int),sno);
    M(i).allocate(rmin,rmax);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin,cmax;
      readbuffer(&cmin,sizeof(int),sno);
      readbuffer(&cmax,sizeof(int),sno);
      M(i,j).allocate(cmin,cmax);
      int sz=cmax-cmin+1;
      readbuffer(&(M(i,j)(cmin)),sz*sizeof(int),sno);
    }
  }
  return M;
}

void adpthread_manager::send_i3_array(const i3_array &x,int sno)
{
  send_id_string_to_master("EWS",sno);
  int smin=x.indexmin();
  int smax=x.indexmax();
  int i;
  writebuffer(&smin,sizeof(int),sno);
  writebuffer(&smax,sizeof(int),sno);
  for (i=smin;i<=smax;i++)
  {
    int rmin=x(i).indexmin();
    int rmax=x(i).indexmax();
    writebuffer(&rmin,sizeof(int),sno);
    writebuffer(&rmax,sizeof(int),sno);
    for (int j=rmin;j<=rmax;j++)
    {
      int cmin=x(i,j).indexmin();
      int cmax=x(i,j).indexmax();
      writebuffer(&cmin,sizeof(int),sno);
      writebuffer(&cmax,sizeof(int),sno);
      int sz=cmax-cmin+1;
      writebuffer(&(x(i,j)(cmin)),sz*sizeof(int),sno);
    }
  }
}

