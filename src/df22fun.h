/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
/**
 * \file
 * Description not yet available.
 */


#if !defined(__DF22FUN__)
#  define __DF22FUN__

/**
 * Description not yet available.
 * \param
 */
class df2_two_variable
{
   MY_DOUBLE_TYPE v[6];
 public:
   static int num_ind_var;
   MY_DOUBLE_TYPE &get_u(void) 
   {
      return v[0];
   }
   MY_DOUBLE_TYPE &get_u_x(void) 
   {
      return v[1];
   }

   MY_DOUBLE_TYPE &get_u_y(void) 
   {
      return v[2];
   }

   MY_DOUBLE_TYPE &get_u_xx(void) 
   {
      return v[3];
   }

   MY_DOUBLE_TYPE &get_u_xy(void) 
   {
      return v[4];
   }

   MY_DOUBLE_TYPE &get_u_yy(void) 
   {
      return v[5];
   }

   df2_two_variable & operator =(const df2_two_variable & v);
    df2_two_variable & operator =(MY_DOUBLE_TYPE v);
    df2_two_variable & operator +=(const df2_two_variable & v);
    df2_two_variable & operator *=(const df2_two_variable & v);
    df2_two_variable & operator *=(MY_DOUBLE_TYPE v);
    df2_two_variable & operator +=(MY_DOUBLE_TYPE v);
    df2_two_variable & operator -=(MY_DOUBLE_TYPE v);
    df2_two_variable & operator -=(const df2_two_variable & v);
    df2_two_variable & operator /=(const df2_two_variable & v);
    df2_two_variable(void);
    df2_two_variable(const df2_two_variable &);
    void set_independent_1(void);
    void set_independent_2(void);
};

void mp(df2_two_variable & x);
/**
 * Description not yet available.
 * \param
 */
inline MY_DOUBLE_TYPE & value(const df2_two_variable & _x)
{
   ADUNCONST(df2_two_variable,x)
   return (MY_DOUBLE_TYPE&)(x.get_u());
}

/**
 * Description not yet available.
 * \param
 */
class init_df2_two_variable:public df2_two_variable
{
 public:
    init_df2_two_variable(MY_DOUBLE_TYPE);
};

/**
 * Description not yet available.
 * \param
 */
class df2_two_vector
{
   int index_min;
   int index_max;
   vector_shapex *shape;
   df2_two_variable *v;
 public:
   int indexmin(void) const
   {
      return int (index_min);
   }
   int indexmax(void) const
   {
      return int (index_max);
   }
   df2_two_vector(int min, int max);
    df2_two_vector(void);
   void allocate(void);
   void allocate(int min, int max);
    df2_two_variable & operator () (int i) const
   {
      return (df2_two_variable &) (*(v + i));
   }
   df2_two_variable & operator [] (int i) const
   {
      return (df2_two_variable &) (*(v + i));
   }
   void initialize(void);
   void deallocate(void);
   ~df2_two_vector();
    df2_two_vector(const df2_two_vector & m2);
};



dvector value(const df2_two_vector & v);

dvector first_derivatives(const df2_two_vector & v);

dvector second_derivatives(const df2_two_vector & v);

dvector third_derivatives(const df2_two_vector & v);

/**
 * Description not yet available.
 * \param
 */
class df2_two_matrix
{
   int index_min;
   int index_max;
   mat_shapex *shape;
   df2_two_vector *v;
 public:
   int indexmin(void) const
   {
      return int (index_min);
   }
   int indexmax(void) const
   {
      return int (index_max);
   }
   df2_two_matrix(int rmin, int rmax, int cmin, int cmax);
    df2_two_vector & operator () (int i) const
   {
      return (df2_two_vector &) * (v + i);
   }
   df2_two_vector & operator [] (int i) const
   {
      return (df2_two_vector &) * (v + i);
   }
   df2_two_variable & operator () (int i, int j) const
   {
      return (df2_two_variable &) (*(v + i)) (j);
   }
   void initialize(void);
   //df2_two_variable& operator () (int i,int j) const
   //  { return *((v+i)->(v+j)); }
   void deallocate(void);
   ~df2_two_matrix();
    df2_two_matrix(const df2_two_matrix & m2);
};

dmatrix value(const df2_two_matrix & v);

dmatrix first_derivatives(const df2_two_matrix & v);
dmatrix second_derivatives(const df2_two_matrix & v);
dmatrix third_derivatives(const df2_two_matrix & v);

/*
  df2_two_variable operator F(const df2_two_variable& x)
  {
    df2_two_variable z;

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

df2_two_variable sin(const df2_two_variable & x);
df2_two_variable fabs(const df2_two_variable & x);
df2_two_variable sqrt(const df2_two_variable & x);
df2_two_variable atan(const df2_two_variable & x);
df2_two_variable cos(const df2_two_variable & x);
df2_two_variable tan(const df2_two_variable & x);
df2_two_variable log(const df2_two_variable & x);
df2_two_variable square(const df2_two_variable & x);
df2_two_variable cube(const df2_two_variable & x);
df2_two_variable pow(const df2_two_variable & x,
                     const df2_two_variable & y);
df2_two_variable sqrt(const df2_two_variable & x);
df2_two_variable exp(const df2_two_variable & x);
df2_two_variable inv(const df2_two_variable & x);
df2_two_variable operator *(const df2_two_variable & x,
                            const df2_two_variable & y);
df2_two_variable operator *(MY_DOUBLE_TYPE x, const df2_two_variable & y);
df2_two_variable operator *(const df2_two_variable & x, MY_DOUBLE_TYPE y);
df2_two_variable operator /(const df2_two_variable & x,
                            const df2_two_variable & y);
df2_two_variable operator /(const MY_DOUBLE_TYPE x, const df2_two_variable & y);

df2_two_variable operator /(const df2_two_variable & x, const MY_DOUBLE_TYPE y);

df2_two_variable operator +(const MY_DOUBLE_TYPE x, const df2_two_variable & y);

df2_two_variable operator +(const df2_two_variable & x, const MY_DOUBLE_TYPE y);

df2_two_variable operator +(const df2_two_variable & x,
                            const df2_two_variable & y);
df2_two_variable operator -(MY_DOUBLE_TYPE x, const df2_two_variable & y);
df2_two_variable operator -(const df2_two_variable & x, MY_DOUBLE_TYPE y);
df2_two_variable operator -(const df2_two_variable & x,
                            const df2_two_variable & y);
int operator <(const df2_two_variable & x, MY_DOUBLE_TYPE n);
int operator >(const df2_two_variable & x, MY_DOUBLE_TYPE n);
int operator >=(const df2_two_variable & x, MY_DOUBLE_TYPE n);
int operator ==(const df2_two_variable & x, const df2_two_variable & n);
int operator ==(const df2_two_variable & x, MY_DOUBLE_TYPE n);
int operator ==(MY_DOUBLE_TYPE x, const df2_two_variable & n);
int operator <(const df2_two_variable & x, const df2_two_variable & n);
int operator >(const df2_two_variable & x, const df2_two_variable & n);

df2_two_variable operator -(const df2_two_variable & v);
df2_two_matrix choleski_decomp(const df2_two_matrix & MM);

df2_two_variable cumd_gamma(const df2_two_variable & x,
                            const df2_two_variable & a);

df2_two_variable gammln(const df2_two_variable & xx);
#endif
