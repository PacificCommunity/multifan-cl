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


void adjoint_send_dvar4_array_to_master(void)
{
  verify_identifier_string("CV");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_send_dvar4_array_to_master();
}

void adpthread_manager::adjoint_send_dvar4_array_to_master(void)
{
  verify_identifier_string("HH");
  int sno=restore_int_value();
  verify_id_string_from_master("FDG",sno);
  int hsmin;
  int hsmax;
  readbuffer(&hsmin,sizeof(int),sno);
  readbuffer(&hsmax,sizeof(int),sno);
  d4_array M(hsmin,hsmax);
  int l;
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    readbuffer(&(smin(l)),sizeof(int),sno);
    readbuffer(&(smax(l)),sizeof(int),sno);
    M(l).allocate(smin(l),smax(l));
    for (int i=smin(l);i<=smax(l);i++)
    {
      int rmin,rmax;
      readbuffer(&rmin,sizeof(int),sno);
      readbuffer(&rmax,sizeof(int),sno);
      M(l,i).allocate(rmin,rmax);
      for (int j=rmin;j<=rmax;j++)
      { 
        int cmin,cmax;
        readbuffer(&cmin,sizeof(int),sno);
        readbuffer(&cmax,sizeof(int),sno);
        M(l,i,j).allocate(cmin,cmax);
        int sz=cmax-cmin+1;
        readbuffer(&(M(l,i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }
  verify_identifier_string("OP");
  for (int l=hsmax;l>=hsmin;l--)
  {
    for (int i=smax(l);i>=smin(l);i--)
    {
      dvar_matrix_position dmpos=restore_dvar_matrix_position();
      M(l,i).save_dmatrix_derivatives(dmpos);
    }
  }
  verify_identifier_string("VC");
}

void adpthread_manager::send_dvar4_array_to_master(const dvar4_array &x,int sno)
{
  send_id_string_to_master("HYD",sno);
  int hsmin=x.indexmin();
  int hsmax=x.indexmax();
  writebuffer(&hsmin,sizeof(int),sno);
  writebuffer(&hsmax,sizeof(int),sno);
  int l;
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    smin(l)=x(l).indexmin();
    smax(l)=x(l).indexmax();
    int i;
    writebuffer(&(smin(l)),sizeof(int),sno);
    writebuffer(&(smax(l)),sizeof(int),sno);
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin=x(l,i).indexmin();
      int rmax=x(l,i).indexmax();
      writebuffer(&rmin,sizeof(int),sno);
      writebuffer(&rmax,sizeof(int),sno);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin=x(l,i,j).indexmin();
        int cmax=x(l,i,j).indexmax();
        writebuffer(&cmin,sizeof(int),sno);
        writebuffer(&cmax,sizeof(int),sno);
        int sz=cmax-cmin+1;
        writebuffer(&(value(x(l,i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }

//    save_identifier_string("VC");
  const char * str1;
  str1="VC";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  // !!! should we optimize this ?
  for (l=hsmin;l<=hsmax;l++)
  {
    for (int i=smin(l);i<=smax(l);i++)
    {
      x(l,i).save_dvar_matrix_position();
    }
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
//    save_identifier_string("CV");
  const char * str4;
  str4="CV";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar4_array_to_master);
}
void adpthread_manager::send_dvar4_array(const dvar4_array &x,int sno)
{
  send_id_string_to_master("HYD",sno);
  int hsmin=x.indexmin();
  int hsmax=x.indexmax();
  writebuffer(&hsmin,sizeof(int),sno);
  writebuffer(&hsmax,sizeof(int),sno);
  int l;
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    smin(l)=x(l).indexmin();
    smax(l)=x(l).indexmax();
    int i;
    writebuffer(&(smin(l)),sizeof(int),sno);
    writebuffer(&(smax(l)),sizeof(int),sno);
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin=x(l,i).indexmin();
      int rmax=x(l,i).indexmax();
      writebuffer(&rmin,sizeof(int),sno);
      writebuffer(&rmax,sizeof(int),sno);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin=x(l,i,j).indexmin();
        int cmax=x(l,i,j).indexmax();
        writebuffer(&cmin,sizeof(int),sno);
        writebuffer(&cmax,sizeof(int),sno);
        int sz=cmax-cmin+1;
        writebuffer(&(value(x(l,i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }

//    save_identifier_string("VC");
  const char * str5;
  str5="VC";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  // !!! should we optimize this ?
  for (l=hsmin;l<=hsmax;l++)
  {
    for (int i=smin(l);i<=smax(l);i++)
    {
      x(l,i).save_dvar_matrix_position();
    }
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
//    save_identifier_string("CV");
  const char * str8;
  str8="CV";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(::adjoint_send_dvar4_array_to_master);
}

void adpthread_manager::send_d4_array(const d4_array &x,int sno)
{
  send_id_string_to_master("KYS",sno);
  int hsmin=x.indexmin();
  int hsmax=x.indexmax();
  writebuffer(&hsmin,sizeof(int),sno);
  writebuffer(&hsmax,sizeof(int),sno);
  int l;
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    smin(l)=x(l).indexmin();
    smax(l)=x(l).indexmax();
    int i;
    writebuffer(&(smin(l)),sizeof(int),sno);
    writebuffer(&(smax(l)),sizeof(int),sno);
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin=x(l,i).indexmin();
      int rmax=x(l,i).indexmax();
      writebuffer(&rmin,sizeof(int),sno);
      writebuffer(&rmax,sizeof(int),sno);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin=x(l,i,j).indexmin();
        int cmax=x(l,i,j).indexmax();
        writebuffer(&cmin,sizeof(int),sno);
        writebuffer(&cmax,sizeof(int),sno);
        int sz=cmax-cmin+1;
        writebuffer(&(x(l,i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }
}

void adjoint_get_dvar4_array_from_slave(void)
{
  verify_identifier_string("I5");
  adpthread_manager * ptr=(adpthread_manager*)(restore_pointer_value());
  ptr->adjoint_get_dvar4_array_from_slave();
}

void adpthread_manager::adjoint_get_dvar4_array_from_slave(void)
{
  verify_identifier_string("LG");
  int sno=restore_int_value();
  verify_identifier_string("W8");
  int hsmin=restore_int_value();
  int hsmax=restore_int_value();
  verify_identifier_string("H3");
  int l;
  d4_array dv(hsmin,hsmax);
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmax;l>=hsmin;l--)
  {
    smin(l)=restore_int_value();
    smax(l)=restore_int_value();
    dv(l).allocate(smin(l),smax(l));
    for (int i=smax(l);i>=smin(l);i--)
    {
      dvar_matrix_position dvpos=restore_dvar_matrix_position();
      dv(l,i)=restore_dvar_matrix_derivatives(dvpos);
    }
  }
  verify_identifier_string("G9");
  send_id_string_to_slave("FDG",sno);
  writebuffer(&hsmin,sizeof(int),sno);
  writebuffer(&hsmax,sizeof(int),sno);
  for (l=hsmin;l<=hsmax;l++)
  {
    writebuffer(&(smin(l)),sizeof(int),sno);
    writebuffer(&(smax(l)),sizeof(int),sno);
    int i;
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin=dv(l,i).indexmin();
      int rmax=dv(l,i).indexmax();
      writebuffer(&rmin,sizeof(int),sno);
      writebuffer(&rmax,sizeof(int),sno);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin=dv(l,i,j).indexmin();
        int cmax=dv(l,i,j).indexmax();
        writebuffer(&cmin,sizeof(int),sno);
        writebuffer(&cmax,sizeof(int),sno);
        int sz=cmax-cmin+1;
        writebuffer(&(dv(l,i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }
}
dvar4_array adpthread_manager::get_dvar4_array_from_slave(int sno)
{
  verify_id_string_from_slave("HYD",sno);

  int hsmin;
  int hsmax;
  int l;
  readbuffer(&hsmin,sizeof(int),sno);
  readbuffer(&hsmax,sizeof(int),sno);
  dvar4_array M(hsmin,hsmax);
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    int i;
    readbuffer(&(smin(l)),sizeof(int),sno);
    readbuffer(&(smax(l)),sizeof(int),sno);
    M(l).allocate(smin(l),smax(l));
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin,rmax;
      readbuffer(&rmin,sizeof(int),sno);
      readbuffer(&rmax,sizeof(int),sno);
      M(l,i).allocate(rmin,rmax);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin,cmax;
        readbuffer(&cmin,sizeof(int),sno);
        readbuffer(&cmax,sizeof(int),sno);
        M(l,i,j).allocate(cmin,cmax);
        int sz=cmax-cmin+1;
        readbuffer(&(value(M(l,i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }

//    save_identifier_string("G9");
  const char * str9;
  str9="G9";
  char* strx9=const_cast <char*> (str9);
  save_identifier_string(strx9);
  for (l=hsmin;l<=hsmax;l++)
  {
    for (int i=smin(l);i<=smax(l);i++)
    {
      M(l,i).save_dvar_matrix_position();
    }
    save_int_value(smax(l));
    save_int_value(smin(l));
  }
//    save_identifier_string("H3");
  const char * str10;
  str10="H3";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
  save_int_value(hsmax);
  save_int_value(hsmin);
//    save_identifier_string("W8");
  const char * str11;
  str11="W8";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
  save_int_value(sno);
//    save_identifier_string("LG");
  const char * str12;
  str12="LG";
  char* strx12=const_cast <char*> (str12);
  save_identifier_string(strx12);
  save_pointer_value(this);
//    save_identifier_string("I5");
  const char * str13;
  str13="I5";
  char* strx13=const_cast <char*> (str13);
  save_identifier_string(strx13);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar4_array_from_slave);
  return M;
}

dvar4_array adpthread_manager::get_dvar4_array(int sno)
{
  verify_id_string_from_slave("HYD",sno);

  int hsmin;
  int hsmax;
  int l;
  readbuffer(&hsmin,sizeof(int),sno);
  readbuffer(&hsmax,sizeof(int),sno);
  dvar4_array M(hsmin,hsmax);
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    int i;
    readbuffer(&(smin(l)),sizeof(int),sno);
    readbuffer(&(smax(l)),sizeof(int),sno);
    M(l).allocate(smin(l),smax(l));
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin,rmax;
      readbuffer(&rmin,sizeof(int),sno);
      readbuffer(&rmax,sizeof(int),sno);
      M(l,i).allocate(rmin,rmax);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin,cmax;
        readbuffer(&cmin,sizeof(int),sno);
        readbuffer(&cmax,sizeof(int),sno);
        M(l,i,j).allocate(cmin,cmax);
        int sz=cmax-cmin+1;
        readbuffer(&(value(M(l,i,j)(cmin))),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }

//    save_identifier_string("G9");
  const char * str14;
  str14="G9";
  char* strx14=const_cast <char*> (str14);
  save_identifier_string(strx14);
  for (l=hsmin;l<=hsmax;l++)
  {
    for (int i=smin(l);i<=smax(l);i++)
    {
      M(l,i).save_dvar_matrix_position();
    }
    save_int_value(smax(l));
    save_int_value(smin(l));
  }
//    save_identifier_string("H3");
  const char * str15;
  str15="H3";
  char* strx15=const_cast <char*> (str15);
  save_identifier_string(strx15);
  save_int_value(hsmax);
  save_int_value(hsmin);
//    save_identifier_string("W8");
  const char * str16;
  str16="W8";
  char* strx16=const_cast <char*> (str16);
  save_identifier_string(strx16);
  save_int_value(sno);
//    save_identifier_string("LG");
  const char * str17;
  str17="LG";
  char* strx17=const_cast <char*> (str17);
  save_identifier_string(strx17);
  save_pointer_value(this);
//    save_identifier_string("I5");
  const char * str18;
  str18="I5";
  char* strx18=const_cast <char*> (str18);
  save_identifier_string(strx18);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(::adjoint_get_dvar4_array_from_slave);
  return M;
}

d4_array adpthread_manager::get_d4_array(int sno)
{
  verify_id_string_from_slave("KYS",sno);

  int hsmin;
  int hsmax;
  int l;
  readbuffer(&hsmin,sizeof(int),sno);
  readbuffer(&hsmax,sizeof(int),sno);
  d4_array M(hsmin,hsmax);
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    int i;
    readbuffer(&(smin(l)),sizeof(int),sno);
    readbuffer(&(smax(l)),sizeof(int),sno);
    M(l).allocate(smin(l),smax(l));
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin,rmax;
      readbuffer(&rmin,sizeof(int),sno);
      readbuffer(&rmax,sizeof(int),sno);
      M(l,i).allocate(rmin,rmax);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin,cmax;
        readbuffer(&cmin,sizeof(int),sno);
        readbuffer(&cmax,sizeof(int),sno);
        M(l,i,j).allocate(cmin,cmax);
        int sz=cmax-cmin+1;
        readbuffer(&(M(l,i,j)(cmin)),sz*sizeof(MY_DOUBLE_TYPE),sno);
      }
    }
  }
  return M;
}

i4_array adpthread_manager::get_i4_array(int sno)
{
  verify_id_string_from_slave("PDN",sno);

  int hsmin;
  int hsmax;
  int l;
  readbuffer(&hsmin,sizeof(int),sno);
  readbuffer(&hsmax,sizeof(int),sno);
  i4_array M(hsmin,hsmax);
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    int i;
    readbuffer(&(smin(l)),sizeof(int),sno);
    readbuffer(&(smax(l)),sizeof(int),sno);
    M(l).allocate(smin(l),smax(l));
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin,rmax;
      readbuffer(&rmin,sizeof(int),sno);
      readbuffer(&rmax,sizeof(int),sno);
      M(l,i).allocate(rmin,rmax);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin,cmax;
        readbuffer(&cmin,sizeof(int),sno);
        readbuffer(&cmax,sizeof(int),sno);
        M(l,i,j).allocate(cmin,cmax);
        int sz=cmax-cmin+1;
        readbuffer(&(M(l,i,j)(cmin)),sz*sizeof(int),sno);
      }
    }
  }
  return M;
}

void adpthread_manager::send_i4_array(const i4_array &x,int sno)
{
  send_id_string_to_master("PDN",sno);
  int hsmin=x.indexmin();
  int hsmax=x.indexmax();
  writebuffer(&hsmin,sizeof(int),sno);
  writebuffer(&hsmax,sizeof(int),sno);
  int l;
  ivector smin(hsmin,hsmax);
  ivector smax(hsmin,hsmax);
  for (l=hsmin;l<=hsmax;l++)
  {
    smin(l)=x(l).indexmin();
    smax(l)=x(l).indexmax();
    int i;
    writebuffer(&(smin(l)),sizeof(int),sno);
    writebuffer(&(smax(l)),sizeof(int),sno);
    for (i=smin(l);i<=smax(l);i++)
    {
      int rmin=x(l,i).indexmin();
      int rmax=x(l,i).indexmax();
      writebuffer(&rmin,sizeof(int),sno);
      writebuffer(&rmax,sizeof(int),sno);
      for (int j=rmin;j<=rmax;j++)
      {
        int cmin=x(l,i,j).indexmin();
        int cmax=x(l,i,j).indexmax();
        writebuffer(&cmin,sizeof(int),sno);
        writebuffer(&cmax,sizeof(int),sno);
        int sz=cmax-cmin+1;
        writebuffer(&(x(l,i,j)(cmin)),sz*sizeof(int),sno);
      }
    }
  }
}
