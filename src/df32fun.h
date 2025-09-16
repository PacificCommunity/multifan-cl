/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
/*
 * $Id$
 *
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 *
 * ADModelbuilder and associated libraries and documentations are
 * provided under the general terms of the "BSD" license.
 *
 * License:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2.  Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3.  Neither the name of the  University of California, Otter Research,
 * nor the ADMB Foundation nor the names of its contributors may be used
 * to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/**
 * \file
 * Description not yet available.
 */


#if !defined(__DF32FUN__)
#  define __DF32FUN__

/**
 * Description not yet available.
 * \param
 */
class df3_two_variable
{
   MY_DOUBLE_TYPE v[10];
 public:
   static int num_ind_var;
   MY_DOUBLE_TYPE *get_u(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[0]));
   }
   MY_DOUBLE_TYPE *get_u_x(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[1]));
   }
   MY_DOUBLE_TYPE *get_u_y(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[2]));
   }
   MY_DOUBLE_TYPE *get_u_xx(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[3]));
   }
   MY_DOUBLE_TYPE *get_u_xy(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[4]));
   }
   MY_DOUBLE_TYPE *get_u_yy(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[5]));
   }
   MY_DOUBLE_TYPE *get_u_xxx(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[6]));
   }
   MY_DOUBLE_TYPE *get_u_xxy(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[7]));
   }
   MY_DOUBLE_TYPE *get_u_xyy(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[8]));
   }
   MY_DOUBLE_TYPE *get_u_yyy(void) const
   {
      return (MY_DOUBLE_TYPE *)  (&(v[9]));
   }
   df3_two_variable & operator =(const df3_two_variable & v);
    df3_two_variable & operator =(MY_DOUBLE_TYPE v);
    df3_two_variable & operator +=(const df3_two_variable & v);
    df3_two_variable & operator *=(const df3_two_variable & v);
    df3_two_variable & operator *=(MY_DOUBLE_TYPE v);
    df3_two_variable & operator +=(MY_DOUBLE_TYPE v);
    df3_two_variable & operator -=(MY_DOUBLE_TYPE v);
    df3_two_variable & operator -=(const df3_two_variable & v);
    df3_two_variable & operator /=(const df3_two_variable & v);
    df3_two_variable(void);
    df3_two_variable(const df3_two_variable &);
    void set_independent_1(void);
    void set_independent_2(void);
};

/**
 * Description not yet available.
 * \param
 */
inline MY_DOUBLE_TYPE value(const df3_two_variable & x)
{
   return (MY_DOUBLE_TYPE) (*x.get_u());
}

/**
 * Description not yet available.
 * \param
 */
class init_df3_two_variable:public df3_two_variable
{
 public:
    init_df3_two_variable(MY_DOUBLE_TYPE);
};

/**
 * Description not yet available.
 * \param
 */
class df3_two_vector
{
   int index_min;
   int index_max;
   vector_shapex *shape;
   df3_two_variable *v;
 public:
   int indexmin(void) const
   {
      return int (index_min);
   }
   int indexmax(void) const
   {
      return int (index_max);
   }
   df3_two_vector(int min, int max);
    df3_two_vector(void);
   void allocate(void);
   void allocate(int min, int max);
    df3_two_variable & operator () (int i) const
   {
      return (df3_two_variable &) (*(v + i));
   }
   df3_two_variable & operator [] (int i) const
   {
      return (df3_two_variable &) (*(v + i));
   }
   void initialize(void);
   void deallocate(void);
   ~df3_two_vector();
    df3_two_vector(const df3_two_vector & m2);
};



dvector value(const df3_two_vector & v);

dvector first_derivatives(const df3_two_vector & v);

dvector second_derivatives(const df3_two_vector & v);

dvector third_derivatives(const df3_two_vector & v);

/**
 * Description not yet available.
 * \param
 */
class df3_two_matrix
{
   int index_min;
   int index_max;
   mat_shapex *shape;
   df3_two_vector *v;
 public:
   int indexmin(void) const
   {
      return int (index_min);
   }
   int indexmax(void) const
   {
      return int (index_max);
   }
   df3_two_matrix(int rmin, int rmax, int cmin, int cmax);
    df3_two_vector & operator () (int i) const
   {
      return (df3_two_vector &) * (v + i);
   }
   df3_two_vector & operator [] (int i) const
   {
      return (df3_two_vector &) * (v + i);
   }
   df3_two_variable & operator () (int i, int j) const
   {
      return (df3_two_variable &) (*(v + i)) (j);
   }
   void initialize(void);
   //df3_two_variable& operator () (int i,int j) const
   //  { return *((v+i)->(v+j)); }
   void deallocate(void);
   ~df3_two_matrix();
    df3_two_matrix(const df3_two_matrix & m2);
};

dmatrix value(const df3_two_matrix & v);

dmatrix first_derivatives(const df3_two_matrix & v);
dmatrix second_derivatives(const df3_two_matrix & v);
dmatrix third_derivatives(const df3_two_matrix & v);

/*
  df3_two_variable operator F(const df3_two_variable& x)
  {
    df3_two_variable z;

    *z.get_u() = ::F(*x.get_u());

    *z.get_udot() = ::D1F(*x.get_u())* *x.get_udot();

    *z.get_udot2() = ::D2F(*x.get_u())* square(*x.get_udot())
                   + ::D1F(*x.get_u())* *x.get_udot2();

    *z.get_udot3() = ::D3F(*x.get_u()) * cube(*x.get_udot())
                   + 3.0 * ::D2F(*x.get_u()) * *x.get_udot() * *x.get_udot2()
                   + ::D1F(*x.get_u()) * *x.get_udot3();
    return z;
  }

*/

df3_two_variable sin(const df3_two_variable & x);
df3_two_variable fabs(const df3_two_variable & x);
df3_two_variable sqrt(const df3_two_variable & x);
df3_two_variable atan(const df3_two_variable & x);
df3_two_variable cos(const df3_two_variable & x);
df3_two_variable tan(const df3_two_variable & x);
df3_two_variable log(const df3_two_variable & x);
df3_two_variable square(const df3_two_variable & x);
df3_two_variable cube(const df3_two_variable & x);
df3_two_variable pow(const df3_two_variable & x,
                     const df3_two_variable & y);
df3_two_variable sqrt(const df3_two_variable & x);
df3_two_variable exp(const df3_two_variable & x);
df3_two_variable inv(const df3_two_variable & x);
df3_two_variable operator *(const df3_two_variable & x,
                            const df3_two_variable & y);
df3_two_variable operator *(MY_DOUBLE_TYPE x, const df3_two_variable & y);
df3_two_variable operator *(const df3_two_variable & x, MY_DOUBLE_TYPE y);
df3_two_variable operator /(const df3_two_variable & x,
                            const df3_two_variable & y);
df3_two_variable operator /(const MY_DOUBLE_TYPE x, const df3_two_variable & y);

df3_two_variable operator /(const df3_two_variable & x, const MY_DOUBLE_TYPE y);

df3_two_variable operator +(const MY_DOUBLE_TYPE x, const df3_two_variable & y);

df3_two_variable operator +(const df3_two_variable & x, const MY_DOUBLE_TYPE y);

df3_two_variable operator +(const df3_two_variable & x,
                            const df3_two_variable & y);
df3_two_variable operator -(MY_DOUBLE_TYPE x, const df3_two_variable & y);
df3_two_variable operator -(const df3_two_variable & x, MY_DOUBLE_TYPE y);
df3_two_variable operator -(const df3_two_variable & x,
                            const df3_two_variable & y);
int operator <(const df3_two_variable & x, MY_DOUBLE_TYPE n);
int operator >(const df3_two_variable & x, MY_DOUBLE_TYPE n);
int operator >=(const df3_two_variable & x, MY_DOUBLE_TYPE n);
int operator ==(const df3_two_variable & x, const df3_two_variable & n);
int operator ==(const df3_two_variable & x, MY_DOUBLE_TYPE n);
int operator ==(MY_DOUBLE_TYPE x, const df3_two_variable & n);
int operator <(const df3_two_variable & x, const df3_two_variable & n);
int operator >(const df3_two_variable & x, const df3_two_variable & n);

df3_two_variable operator -(const df3_two_variable & v);
df3_two_matrix choleski_decomp(const df3_two_matrix & MM);

df3_two_variable cumd_gamma(const df3_two_variable & x,
                            const df3_two_variable & a);

df3_two_variable gammln(const df3_two_variable & xx);
#endif
