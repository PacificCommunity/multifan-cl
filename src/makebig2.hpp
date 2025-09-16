/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#if !defined(__MAKEBIG2__)
#define __MAKEBIG2__

#include <admodel.h>
typedef  cubic_spline_function * pcubic_spline_function; 
typedef  pcubic_spline_function * ppcubic_spline_function; 

//void add_slave_suffix(adstring const&) {;}

class cubic_spline_array
{
  int fsblen;
  int * ncopies;
  int rowmin;
  int rowmax;
  int colmin;
  int colmax;
  cubic_spline_function *** ptr;
  dvector x;
  dvar_matrix choleski_inv;
  dvar_matrix tcholeski_inv;
public:
  void set_z(const prevariable& z);
  void set_z(MY_DOUBLE_TYPE z);
  cubic_spline_array(int _maxindex,dvector& x, d3_array& y);
  void allocate(int _maxindex,dvector& x, d3_array& y);
  cubic_spline_array(void);
  cubic_spline_array(unsigned char fsb[]);
  ~cubic_spline_array();
  cubic_spline_array(const cubic_spline_array & csf);
  void allocate(const cubic_spline_array & csf);
  dmatrix operator () (MY_DOUBLE_TYPE& x);
  int get_rowmin() {return rowmin; }
  dvector get_x() { return x; }
  int get_rowmax() {return rowmax; }
  int get_colmin() {return colmin; }
  int get_colmax() {return colmax; }
  dvar_matrix get_choleski_inv() { return choleski_inv;}
  dvar_matrix get_tcholeski_inv() { return tcholeski_inv;}
  dvar_matrix get_tcholeski_inv(int n);
  banded_lower_triangular_dvar_matrix get_ltcholeski_inv(int n);
  dvariable get_ln_det_choleski_inv(int n);
  void myread(MY_DOUBLE_TYPE & n,unsigned char fsb[],int & offset);
  void myread(dvector& alpha,unsigned char fsb[],int &offset);
  void myread(dmatrix& M,unsigned char fsb[],int &offset);
  void myread(int & n,unsigned char fsb[],int & offset);
};

class ccubic_spline_array  // constant version
{
  int fsblen;
  int * ncopies;
  int rowmin;
  int rowmax;
  int colmin;
  int colmax;
  cubic_spline_function *** ptr;
  dvector x;
  dmatrix choleski_inv;
  dmatrix tcholeski_inv;
public:
  void set_z(MY_DOUBLE_TYPE z);
  ccubic_spline_array(int _maxindex,dvector& x, d3_array& y);
  void allocate(int _maxindex,dvector& x, d3_array& y);
  ccubic_spline_array(void);
  ccubic_spline_array(unsigned char fsb[]);
  ~ccubic_spline_array();
  ccubic_spline_array(const cubic_spline_array & csf);
  void allocate(const cubic_spline_array & csf);
  dmatrix operator () (MY_DOUBLE_TYPE& x);
  int get_rowmin() {return rowmin; }
  dvector get_x() { return x; }
  int get_rowmax() {return rowmax; }
  int get_colmin() {return colmin; }
  int get_colmax() {return colmax; }
  dmatrix get_choleski_inv() { return choleski_inv;}
  dmatrix get_tcholeski_inv() { return tcholeski_inv;}
  dmatrix get_tcholeski_inv(int n);
  banded_lower_triangular_dmatrix get_ltcholeski_inv(int n);
  MY_DOUBLE_TYPE get_ln_det_choleski_inv(int n);
  void myread(MY_DOUBLE_TYPE & n,unsigned char fsb[],int & offset);
  void myread(dvector& alpha,unsigned char fsb[],int &offset);
  void myread(dmatrix& M,unsigned char fsb[],int &offset);
  void myread(int & n,unsigned char fsb[],int & offset);
};


  
#endif 
