/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#if !defined(__DF3FUN__)
#  define __DF3FUN__

  class df3_one_variable
  {
    MY_DOUBLE_TYPE v[4];
  public:
    static df1b2variable * ind_var;
    MY_DOUBLE_TYPE * get_u(void) const {return (MY_DOUBLE_TYPE*)  (&(v[0]));}
    MY_DOUBLE_TYPE * get_udot(void)const {return (MY_DOUBLE_TYPE*)  (&(v[1]));}
    MY_DOUBLE_TYPE * get_udot2(void) const {return (MY_DOUBLE_TYPE*)  (&(v[2]));}
    MY_DOUBLE_TYPE * get_udot3(void) const {return (MY_DOUBLE_TYPE*)  (&(v[3]));}
    df3_one_variable& operator = (const df3_one_variable& v);
    df3_one_variable& operator = (MY_DOUBLE_TYPE v);
    df3_one_variable& operator += (const df3_one_variable& v);
    df3_one_variable& operator -= (const df3_one_variable& v);
    df3_one_variable& operator -= (MY_DOUBLE_TYPE);
    df3_one_variable& operator *= (const df3_one_variable& v);
    df3_one_variable& operator /= (const df3_one_variable& v);
    df3_one_variable(void);
    df3_one_variable(const df3_one_variable& );
  };

  inline  MY_DOUBLE_TYPE value(const df3_one_variable& x) { return(MY_DOUBLE_TYPE)( *x.get_u()); }

  class init_df3_one_variable : public df3_one_variable
  {
  public:
    init_df3_one_variable(const df1b2variable& );
    init_df3_one_variable(MY_DOUBLE_TYPE );
  };



  class df3_one_vector
  {
    int index_min;
    int index_max;
    vector_shapex * shape;
    df3_one_variable * v;
  public:
    int indexmin(void) const { return int(index_min); }
    int indexmax(void) const { return int(index_max); }
    df3_one_vector(int min,int max);
    df3_one_vector(void);
    void allocate(void);
    void allocate(int min,int max);
    df3_one_variable& operator () (int i) const 
    { 
      return (df3_one_variable&) (*(v+i));
    }
    df3_one_variable& operator [] (int i) const 
    { 
      return (df3_one_variable&) (*(v+i));
    }
    void initialize(void);
    void deallocate(void);
    ~df3_one_vector();
    df3_one_vector(const df3_one_vector& m2);
  };



 dvector value(const df3_one_vector& v);
   
 dvector first_derivatives(const df3_one_vector& v);

 dvector second_derivatives(const df3_one_vector& v);
   
 dvector third_derivatives(const df3_one_vector& v);
    
  class df3_one_matrix
  {
    int index_min;
    int index_max;
    mat_shapex * shape;
    df3_one_vector * v;
  public:
    int indexmin(void) const { return int(index_min); }
    int indexmax(void) const { return int(index_max); }
    df3_one_matrix(int rmin,int rmax,int cmin,int cmax);
    df3_one_vector& operator () (int i) const 
    { 
      return (df3_one_vector&) *(v+i); 
    }
    df3_one_vector& operator [] (int i) const 
    { 
      return (df3_one_vector&) *(v+i); 
    }
    df3_one_variable& operator () (int i,int j) const 
    { 
      return (df3_one_variable&) (*(v+i))(j); 
    }
    void initialize(void);
    //df3_one_variable& operator () (int i,int j) const { return *((v+i)->(v+j)); }
    void deallocate(void);
    ~df3_one_matrix();
    df3_one_matrix(const df3_one_matrix& m2);
  };

 dmatrix value(const df3_one_matrix& v);
   
 dmatrix first_derivatives(const df3_one_matrix& v);
 dmatrix second_derivatives(const df3_one_matrix& v);
 dmatrix third_derivatives(const df3_one_matrix& v);
    
/*
  df3_one_variable operator F(const df3_one_variable& x)
  {
    df3_one_variable z;

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

  df3_one_variable sqrt(const df3_one_variable& x);
  df3_one_variable exp(const df3_one_variable& x);
  df3_one_variable log(const df3_one_variable& x);
  df3_one_variable inv(const df3_one_variable& x);
  df3_one_variable operator * (const df3_one_variable& x,
    const df3_one_variable& y);
  df3_one_variable operator * (MY_DOUBLE_TYPE x,
    const df3_one_variable& y);
  df3_one_variable operator * (const df3_one_variable& x,
    MY_DOUBLE_TYPE y);
  df3_one_variable operator / (const df3_one_variable& x,
    const df3_one_variable& y);

  df3_one_variable operator / (const MY_DOUBLE_TYPE x,const df3_one_variable& y);

  df3_one_variable operator / (const df3_one_variable& x,
    const MY_DOUBLE_TYPE y);

  df3_one_variable operator + (const MY_DOUBLE_TYPE x,
    const df3_one_variable& y);

  df3_one_variable operator + (const df3_one_variable& x,
    const MY_DOUBLE_TYPE y);

  df3_one_variable operator + (const df3_one_variable& x,
    const df3_one_variable& y);

  df3_one_variable operator - (const df3_one_variable& x,
    MY_DOUBLE_TYPE);
  df3_one_variable operator - (double,
    const df3_one_variable& y);
  df3_one_variable operator - (const df3_one_variable& x,
    const df3_one_variable& y);
  df3_one_variable operator - (const df3_one_variable& x,
    const df3_one_variable& y);
  df3_one_variable operator / (const df3_one_variable& x,
    const df3_one_variable& y);
  df3_one_variable operator * (const df3_one_variable& x,
    const df3_one_variable& y);


  df3_one_variable operator - (const df3_one_variable& v);
  df3_one_matrix choleski_decomp(const df3_one_matrix& MM);
 
#endif  // __DF3FUN__

