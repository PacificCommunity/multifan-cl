/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
void set_value_partial(const dvar_matrix& w,const dvar_vector& x,const int& ii,const ivector& range,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,const ivector& group,
  MY_DOUBLE_TYPE scale);

void set_value_inv_partial(const dvar_matrix& w,const dvector& x,const int& ii,const ivector& range,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale);

void set_value(const prevariable& x,const dvar_vector& v, const int& ii,MY_DOUBLE_TYPE s); 

void set_value(const dvar_vector& x,const dvar_vector& v, const int& ii,MY_DOUBLE_TYPE s);

void set_value(const dvar3_array& x,const dvar_vector& v, const int& ii,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, dvariable fpen,MY_DOUBLE_TYPE s);

void set_value_partial(const dvar_vector& x,const dvar_vector& v, const int& ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen, MY_DOUBLE_TYPE s);

void set_value_partial(const dvar_vector& x,const dvar_vector& v, const int& ii, int n,
  MY_DOUBLE_TYPE s);

void set_value_inv_partial(const dvar_vector& x,const dvector& v, const int& ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s);

void set_value_inv_partial(const dvar_vector& x,const dvector& v, const int& ii, int n,
  MY_DOUBLE_TYPE s);

void set_value(const dvar3_array& x,const dvar_vector& v, const int& ii,MY_DOUBLE_TYPE s);

void set_value_inv(const dvar3_array& x,const dvector& v, const int& ii, MY_DOUBLE_TYPE s);

MY_DOUBLE_TYPE boundp( MY_DOUBLE_TYPE xx, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const MY_DOUBLE_TYPE& fpen,
  MY_DOUBLE_TYPE s);

MY_DOUBLE_TYPE boundpin(MY_DOUBLE_TYPE x, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);


void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale);

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale);

//void set_value_inv_partial(const d3_array& w,const dvector& x,
//  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
//  const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)

void set_value_inv_partial(const d3_array& w,const dvector& x,const int& ii,
      const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,
      const ivector& group,MY_DOUBLE_TYPE scale);
int num_active_partial(const dvar3_array& w,const ivector& flags,
    const ivector& group,const ivector& range);
void set_value_partial(const dvar_matrix _x,const dvar_vector& v,
  const int& _ii, int n,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen, 
  MY_DOUBLE_TYPE s);

void set_value_inv_partial(const dvector& x,const dvector& _v, const int& _ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s);


void set_value_inv_partial(const dvar_vector& x,const dvector& _v, 
  const int& _ii, int nl, int nu,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s);

void set_value_inv_partial(const dvar_matrix& x,const dvector& _v,
  const int& _ii, int nl,int nu,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s);
void set_value_partial(const dvar_vector& _x,const dvar_vector& v,
  const int& _ii, int nl,int nu,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen, MY_DOUBLE_TYPE s);

void set_value_partial(const dvar_matrix _x,const dvar_vector& v,
  const int& _ii,int nl,int nu,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,
  const dvariable& fpen,MY_DOUBLE_TYPE s);
int size_count_partial(_CONST dvar_vector& x,int nl,int nu);

int size_count_partial(_CONST dvar_matrix& x, int nl,int nu);

void set_value_partial(const dvar_matrix& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& group,MY_DOUBLE_TYPE scale);
void set_value_inv_partial(const dvar_matrix& _w,const dvector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& group,MY_DOUBLE_TYPE scale);
int num_active_partial(const dvar_matrix& w,const ivector& group,
    const ivector& range);
