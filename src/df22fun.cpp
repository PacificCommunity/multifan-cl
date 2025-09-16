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
 */
/**
 * \file
 * Description not yet available.
 */

/**
 * Description not yet available.
 * \param
 */
  #include <fvar.hpp>
  #include "df22fun.h"
  df2_two_variable::df2_two_variable(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    v[0]=x.v[0];
    v[1]=x.v[1];
    v[2]=x.v[2];
    v[3]=x.v[3];
    v[4]=x.v[4];
    v[5]=x.v[5];
  }

  void df2_two_variable::set_independent_1()
  {
    v[1]=1.0;
    v[2]=0.0;
    v[3]=0.0;
    v[4]=0.0;
    v[5]=0.0;
  }
  void df2_two_variable::set_independent_2()
  {
    v[1]=0.0;
    v[2]=1.0;
    v[3]=0.0;
    v[4]=0.0;
    v[5]=0.0;
  }

/**
 * Description not yet available.
 * \param
 */
 df2_two_vector::df2_two_vector(const df2_two_vector& m2)
 {
   index_min=m2.index_min;
   index_max=m2.index_max;
   shape=m2.shape;
   if (shape)
   {
     (shape->ncopies)++;
   }
   v = m2.v;
 }

/**
 * Description not yet available.
 * \param
 */
 df2_two_vector::~df2_two_vector()
 {
   deallocate();
 }

/**
 * Description not yet available.
 * \param
 */
 void df2_two_vector::deallocate(void)
 {
   if(shape)
   {
     if (shape->ncopies)
     {
       (shape->ncopies)--;
     }
     else
     {
       v = (df2_two_variable*) (shape->trueptr);
       delete [] v;
       v = NULL;
       delete shape;
       shape=0;
     }
   }
 }

/**
 * Description not yet available.
 * \param
 */
 dvector value(const df2_two_vector& v)
 {
   int mmin=v.indexmin();
   int mmax=v.indexmax();
   dvector cv(mmin,mmax);
   for (int i=mmin;i<=mmax;i++)
   {
     cv(i)=value(v(i));
   }
   return cv;
 }

/**
 * Description not yet available.
 * \param
 */
  void df2_two_vector::initialize(void)
  {
    int mmin=indexmin();
    int mmax=indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      (*this)(i)=0.0;
    }
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_vector::df2_two_vector(void)
  {
    allocate();
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_vector::df2_two_vector(int min,int max)
  {
    allocate(min,max);
  }

/**
 * Description not yet available.
 * \param
 */
  void df2_two_vector::allocate(int min,int max)
  {
    index_min=min;
    index_max=max;
    v=new df2_two_variable[max-min+1];
    if (v==0)
    {
      cerr << "error allocating memory in df2_two_vector" << endl;
      ad_exit(1);
    }
    if ( (shape=new vector_shapex(min,max,v)) == NULL)
    {
      cerr << "Error trying to allocate memory for df2_two_vector"
           << endl;;
      ad_exit(1);
    }
    v-=min;
  }

/**
 * Description not yet available.
 * \param
 */
  void df2_two_vector::allocate(void)
  {
    index_min=0;
    index_max=-1;
    v=0;
    shape=0;
  }

/**
 * Description not yet available.
 * \param
 */
 dmatrix value(const df2_two_matrix& v)
 {
   int rmin=v.indexmin();
   int rmax=v.indexmax();
   dmatrix cm(rmin,rmax);
   for (int i=rmin;i<=rmax;i++)
   {
     int cmin=v(i).indexmin();
     int cmax=v(i).indexmax();
     cm(i).allocate(cmin,cmax);
     for (int j=cmin;j<=cmax;j++)
     {
       cm(i,j)=value(v(i,j));
     }
   }
   return cm;
 }

/**
 * Description not yet available.
 * \param
 */
 df2_two_matrix::df2_two_matrix(const df2_two_matrix& m2)
 {
   index_min=m2.index_min;
   index_max=m2.index_max;
   shape=m2.shape;
   if (shape)
   {
     (shape->ncopies)++;
   }
   v = m2.v;
 }

/**
 * Description not yet available.
 * \param
 */
 df2_two_matrix::~df2_two_matrix()
 {
   deallocate();
 }

/**
 * Description not yet available.
 * \param
 */
 void df2_two_matrix::deallocate(void)
 {
   if (shape)
   {
     if (shape->ncopies)
     {
       (shape->ncopies)--;
     }
     else
     {
       v = (df2_two_vector*) (shape->get_pointer());
       delete [] v;
       v=0;
       delete shape;
       shape=0;
     }
   }
 }

/**
 * Description not yet available.
 * \param
 */
  void df2_two_matrix::initialize(void)
  {
    int mmin=indexmin();
    int mmax=indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      (*this)(i).initialize();
    }
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_matrix::df2_two_matrix(int rmin,int rmax,int cmin,int cmax)
  {
    index_min=rmin;
    index_max=rmax;
    v=new df2_two_vector[rmax-rmin+1];
    if (v==0)
    {
      cerr << "error allocating memory in df2_two_matrix" << endl;
      ad_exit(1);
    }
    if ( (shape=new mat_shapex(v)) == NULL)
    {
      cerr << "Error trying to allocate memory for df2_two_vector"
           << endl;;
    }
    v-=rmin;

    for (int i=rmin;i<=rmax;i++)
    {
      v[i].allocate(cmin,cmax);
    }
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator -= (const df2_two_variable& _v)
  {
    ADUNCONST(df2_two_variable,v)
    get_u() -= v.get_u();
    get_u_x() -= v.get_u_x();
    get_u_y() -= v.get_u_y();
    get_u_xx() -= v.get_u_xx();
    get_u_xy() -= v.get_u_xy();
    get_u_yy() -= v.get_u_yy();
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
df2_two_variable operator-(const df2_two_variable& _v)
{
  ADUNCONST(df2_two_variable,v)
  df2_two_variable z;

  z.get_u() = - v.get_u();
  z.get_u_x() = -v.get_u_x();
  z.get_u_y() = -v.get_u_y();
  z.get_u_xx() = -v.get_u_xx();
  z.get_u_xy() = -v.get_u_xy();
  z.get_u_yy() = -v.get_u_yy();

  return z;
}

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator += (const df2_two_variable& _v)
  {
    ADUNCONST(df2_two_variable,v)
    get_u() += v.get_u();
    get_u_x() += v.get_u_x();
    get_u_y() += v.get_u_y();
    get_u_xx() += v.get_u_xx();
    get_u_xy() += v.get_u_xy();
    get_u_yy() += v.get_u_yy();

    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator += (MY_DOUBLE_TYPE v)
  {
    get_u() += v;
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator -= (MY_DOUBLE_TYPE v)
  {
    get_u() -= v;
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator *= (const df2_two_variable& _v)
  {
    ADUNCONST(df2_two_variable,v)
    df2_two_variable x=*this * v;
    *this=x;
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator *= (MY_DOUBLE_TYPE v)
  {
    get_u() *= v;
    get_u_x() *= v;
    get_u_y() *= v;
    get_u_xx() *= v;
    get_u_xy() *= v;
    get_u_yy() *= v;
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator /= (const df2_two_variable& y)
  {
    df2_two_variable x=*this / y;
    *this=x;
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
int operator <(const df2_two_variable & x, MY_DOUBLE_TYPE n)
{
   return value(x) < n;
}

/**
 * Description not yet available.
 * \param
 */
int operator >(const df2_two_variable & x, MY_DOUBLE_TYPE n)
{
   return value(x) > n;
}

/**
 * Description not yet available.
 * \param
 */
int operator >=(const df2_two_variable & x, MY_DOUBLE_TYPE n)
{
   return value(x) >= n;
}

/**
 * Description not yet available.
 * \param
 */
int operator ==(const df2_two_variable & x, const df2_two_variable & n)
{
   return value(x) == value(n);
}

/**
 * Description not yet available.
 * \param
 */
int operator ==(const df2_two_variable & x, MY_DOUBLE_TYPE n)
{
   return value(x) == n;
}

/**
 * Description not yet available.
 * \param
 */
int operator ==(MY_DOUBLE_TYPE x, const df2_two_variable & n)
{
   return x == value(n);
}

/**
 * Description not yet available.
 * \param
 */
int operator <(const df2_two_variable & x, const df2_two_variable & n)
{
   return value(x) < value(n);
}

/**
 * Description not yet available.
 * \param
 */
int operator >(const df2_two_variable & x, const df2_two_variable & n)
{
   return value(x) > value(n);
}

/**
 * Description not yet available.
 * \param
 */
void set_derivatives( df2_two_variable& z,const df2_two_variable& _x,MY_DOUBLE_TYPE u,
  MY_DOUBLE_TYPE zp,MY_DOUBLE_TYPE zp2)
{
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;

    z.get_u_x() = zp* x.get_u_x();

    z.get_u_y() = zp* x.get_u_y();

    z.get_u_xx() = zp2 * square(x.get_u_x())
                   + zp * x.get_u_xx();

    z.get_u_xy() = zp2 * x.get_u_x() * x.get_u_y()
                   + zp * x.get_u_xy();

    z.get_u_yy() = zp2 * square(x.get_u_y())
                   + zp * x.get_u_yy();

}

/**
 * Description not yet available.
 * \param
 */
void set_derivatives( df2_two_variable& z, const df2_two_variable& _x,
  const df2_two_variable& _y, MY_DOUBLE_TYPE u,
  MY_DOUBLE_TYPE f_u,MY_DOUBLE_TYPE f_v,MY_DOUBLE_TYPE f_uu,MY_DOUBLE_TYPE f_uv,MY_DOUBLE_TYPE f_vv)
{
    ADUNCONST(df2_two_variable,x)
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    //*z.get_u() = u;
    z.get_u() = u;

    z.get_u_x() = f_u* x.get_u_x()
                 + f_v* y.get_u_x();

    z.get_u_y() = f_u* x.get_u_y()
                 + f_v* y.get_u_y();

    z.get_u_xx() = f_uu * square(x.get_u_x())
                  + f_u  * x.get_u_xx()
                  + f_vv * square(y.get_u_x())
                  + f_v  * y.get_u_xx()
            + 2.0 * f_uv * x.get_u_x() * y.get_u_x();

    z.get_u_xy() = f_uu * x.get_u_x() * x.get_u_y()
                  + f_u  * x.get_u_xy()
                  + f_vv * y.get_u_x() * y.get_u_y()
                  + f_v  * y.get_u_xy()
                  + f_uv * (x.get_u_x() * y.get_u_y()
                         +  x.get_u_y() * y.get_u_x());

    z.get_u_yy() = f_uu * square(x.get_u_y())
                  + f_u  * x.get_u_yy()
                  + f_vv * square(y.get_u_y())
                  + f_v  * y.get_u_yy()
            + 2.0 * f_uv * x.get_u_y() * y.get_u_y();
}

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable sqrt(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u=sqrt(x.get_u());
    z.get_u()=u;
    MY_DOUBLE_TYPE xinv=1.0/(x.get_u());
    MY_DOUBLE_TYPE zp=0.5/u;
    MY_DOUBLE_TYPE zp2=-0.5*zp*xinv;


    set_derivatives(z,x,u,zp,zp2);

    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable atan(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE cx=value(x);
    MY_DOUBLE_TYPE d=1.0/(1+square(cx));
    MY_DOUBLE_TYPE d2=square(d);
    MY_DOUBLE_TYPE u=atan(cx);
    z.get_u()=u;
    MY_DOUBLE_TYPE zp=d;
    MY_DOUBLE_TYPE zp2=-2.0*cx*d2;

    set_derivatives(z,x,u,zp,zp2);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable square(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u=value(x);
    z.get_u()=u*u;
    MY_DOUBLE_TYPE zp=2.0*u;
    MY_DOUBLE_TYPE zp2=2.0;

    set_derivatives(z,x,u,zp,zp2);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable tan(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u=tan(x.get_u());
    z.get_u()=u;
    MY_DOUBLE_TYPE v=1.0/cos(x.get_u());
    MY_DOUBLE_TYPE w=sin(x.get_u());
    MY_DOUBLE_TYPE v2=v*v;
    MY_DOUBLE_TYPE zp=v2;
    MY_DOUBLE_TYPE zp2=2.0*w*v2*v;

    set_derivatives(z,x,u,zp,zp2);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable sin(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u=sin(x.get_u());
    z.get_u()=u;
    MY_DOUBLE_TYPE zp=cos(x.get_u());
    MY_DOUBLE_TYPE zp2=-u;

    set_derivatives(z,x,u,zp,zp2);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable fabs(const df2_two_variable& _v)
  {
    ADUNCONST(df2_two_variable,v)
    //*z.get_u() = u;
    df2_two_variable z;
    if (value(v)>=0)
    {
      z.get_u() = v.get_u();
      z.get_u_x() = v.get_u_x();
      z.get_u_y() = v.get_u_y();
      z.get_u_xx() = v.get_u_xx();
      z.get_u_xy() = v.get_u_xy();
      z.get_u_yy() = v.get_u_yy();
    }
    else
    {
      z.get_u() = -v.get_u();
      z.get_u_x() = -v.get_u_x();
      z.get_u_y() = -v.get_u_y();
      z.get_u_xx() = -v.get_u_xx();
      z.get_u_xy() = -v.get_u_xy();
      z.get_u_yy() = -v.get_u_yy();
    }

    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable log(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u=log(x.get_u());
    z.get_u()=u;
    MY_DOUBLE_TYPE zp=1/(x.get_u());
    MY_DOUBLE_TYPE zp2=-zp*zp;

    set_derivatives(z,x,u,zp,zp2);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable exp(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u=exp(x.get_u());
    z.get_u()=u;
    MY_DOUBLE_TYPE zp=u;
    MY_DOUBLE_TYPE zp2=u;

    set_derivatives(z,x,u,zp,zp2);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable pow(const df2_two_variable& _x,
                       const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,x)
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    //*z.get_u() = u;
    df2_two_variable z;
    z=exp(y*log(x));
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable inv(const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE xinv=1.0/(x.get_u());
    z.get_u()=xinv;
    MY_DOUBLE_TYPE zp=-xinv*xinv;
    MY_DOUBLE_TYPE zp2=-2.0*zp*xinv;

    set_derivatives(z,x,xinv,zp,zp2);

    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator = (const df2_two_variable& _x)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    get_u() = x.get_u();
    get_u_x() = x.get_u_x();
    get_u_y() = x.get_u_y();
    get_u_xx() = x.get_u_xx();
    get_u_xy() = x.get_u_xy();
    get_u_yy() = x.get_u_yy();
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable& df2_two_variable::operator = (MY_DOUBLE_TYPE x)
  {
    get_u() = x;
    get_u_x() =0.0;
    get_u_y() =0.0;
    get_u_xx() =0.0;
    get_u_xy() =0.0;
    get_u_yy() =0.0;
    return *this;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator * (const df2_two_variable& _x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,x)
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    //*z.get_u() = u;
    df2_two_variable z;
    MY_DOUBLE_TYPE u= x.get_u() * y.get_u();
    z.get_u() = u;
    MY_DOUBLE_TYPE f_u=y.get_u();
    MY_DOUBLE_TYPE f_v=x.get_u();
    MY_DOUBLE_TYPE f_uu=0.0;
    MY_DOUBLE_TYPE f_uv=1.0;
    MY_DOUBLE_TYPE f_vv=0.0;
    set_derivatives(z,x,y,u,
      f_u, f_v,
      f_uu, f_uv, f_vv);
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator * (MY_DOUBLE_TYPE x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() = x *  y.get_u();
    z.get_u_x() = x * y.get_u_x();
    z.get_u_y() = x * y.get_u_y();
    z.get_u_xx() = x * y.get_u_xx();
    z.get_u_xy() = x * y.get_u_xy();
    z.get_u_yy() = x * y.get_u_yy();

    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator * (const df2_two_variable& _y,
    MY_DOUBLE_TYPE x)
  {
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() = x *  y.get_u();
    z.get_u_x() = x * y.get_u_x();
    z.get_u_y() = x * y.get_u_y();
    z.get_u_xx() = x * y.get_u_xx();
    z.get_u_xy() = x * y.get_u_xy();
    z.get_u_yy() = x * y.get_u_yy();

    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator / (const df2_two_variable& _x,
    MY_DOUBLE_TYPE y)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    MY_DOUBLE_TYPE u=1/y;
    return x*u;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator / (const df2_two_variable& _x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,x)
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    //*z.get_u() = u;
    df2_two_variable u=inv(y);
    return x*u;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator / (const MY_DOUBLE_TYPE x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    df2_two_variable u=inv(y);
    return x*u;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator + (const MY_DOUBLE_TYPE x,const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() =  x + y.get_u();
    z.get_u_x() = y.get_u_x();
    z.get_u_y() = y.get_u_y();
    z.get_u_xx() = y.get_u_xx();
    z.get_u_xy() = y.get_u_xy();
    z.get_u_yy() = y.get_u_yy();
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator + (const df2_two_variable& _x,const MY_DOUBLE_TYPE y)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() =  x.get_u() + y;
    z.get_u_x() = x.get_u_x();
    z.get_u_y() = x.get_u_y();
    z.get_u_xx() = x.get_u_xx();
    z.get_u_xy() = x.get_u_xy();
    z.get_u_yy() = x.get_u_yy();
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator + (const df2_two_variable& _x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,x)
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() = x.get_u() + y.get_u();
    z.get_u_x() = x.get_u_x() + y.get_u_x();
    z.get_u_y() = x.get_u_y()+y.get_u_y();
    z.get_u_xx() = x.get_u_xx()+ y.get_u_xx();
    z.get_u_xy() = x.get_u_xy()+ y.get_u_xy();
    z.get_u_yy() = x.get_u_yy()+ y.get_u_yy();
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator - (const df2_two_variable& _x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,x)
    ADUNCONST(df2_two_variable,y)
    df2_two_variable z;
    z.get_u() = x.get_u() - y.get_u();
    z.get_u_x() = x.get_u_x()  - y.get_u_x();
    z.get_u_y() = x.get_u_y() - y.get_u_y();
    z.get_u_xx() = x.get_u_xx() - y.get_u_xx();
    z.get_u_xy() = x.get_u_xy() - y.get_u_xy();
    z.get_u_yy() = x.get_u_yy() - y.get_u_yy();
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator - (MY_DOUBLE_TYPE x,
    const df2_two_variable& _y)
  {
    ADUNCONST(df2_two_variable,y)
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() = x - y.get_u();
    z.get_u_x() = - y.get_u_x();
    z.get_u_y() = - y.get_u_y();
    z.get_u_xx() = - y.get_u_xx();
    z.get_u_xy() = - y.get_u_xy();
    z.get_u_yy() = - y.get_u_yy();
    return z;
  }

/**
 * Description not yet available.
 * \param
 */
  df2_two_variable operator - (const df2_two_variable& _x,
    MY_DOUBLE_TYPE y)
  {
    ADUNCONST(df2_two_variable,x)
    //*z.get_u() = u;
    //*z.get_u() = u;
    df2_two_variable z;
    z.get_u() = x.get_u()-y;
    z.get_u_x() = x.get_u_x();
    z.get_u_y() = x.get_u_y();
    z.get_u_xx() = x.get_u_xx();
    z.get_u_xy() = x.get_u_xy();
    z.get_u_yy() = x.get_u_yy();
    return z;
  }

/**
 * Description not yet available.
 * \param
 */

/**
 * Description not yet available.
 * \param
 */
  init_df2_two_variable::init_df2_two_variable(MY_DOUBLE_TYPE v)
  {
    get_u() =  v;
    get_u_x() = 0.0;
    get_u_y() = 0.0;
    get_u_xx() = 0.0;
    get_u_xy() = 0.0;
    get_u_yy() = 0.0;
  }

  df2_two_variable::df2_two_variable(void)
  {
  }

/**
 * Description not yet available.
 * \param
 */
df2_two_matrix choleski_decomp(const df2_two_matrix& MM)
{
  // kludge to deal with constantness
  df2_two_matrix & M= (df2_two_matrix &) MM;
  int rmin=M.indexmin();
  int cmin=M(rmin).indexmin();
  int rmax=M.indexmax();
  int cmax=M(rmin).indexmax();
  if (rmin !=1 || cmin !=1)
  {
    cerr << "minimum row and column inidices must equal 1 in "
      "df1b2matrix choleski_decomp(const df2_two_atrix& MM)"
         << endl;
    ad_exit(1);
  }
  if (rmax !=cmax)
  {
    cerr << "Error in df1b2matrix choleski_decomp(const df2_two_matrix& MM)"
      " Matrix not square" << endl;
    ad_exit(1);
  }

  int n=rmax-rmin+1;
  df2_two_matrix L(1,n,1,n);
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  int i,j,k;
  df2_two_variable tmp;

    if (value(M(1,1))<=0)
    {
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;
      ad_exit(1);
    }

  L(1,1)=sqrt(M(1,1));
  for (i=2;i<=n;i++)
  {
    L(i,1)=M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=2;j<=i-1;j++)
    {
      tmp=M(i,j);
      for (k=1;k<=j-1;k++)
      {
        tmp-=L(i,k)*L(j,k);
      }
      L(i,j)=tmp/L(j,j);
    }
    tmp=M(i,i);
    for (k=1;k<=i-1;k++)
    {
      tmp-=L(i,k)*L(i,k);
    }

    if (value(tmp)<=0)
    {
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;
      ad_exit(1);
    }

    L(i,i)=sqrt(tmp);
  }

  return L;
}

void mp(df2_two_variable& x) 
{ 
  cout << x.get_u() << " " << x.get_u_x() << " " << x.get_u_y() << "  "
   << x.get_u_xx() << " " << x.get_u_xy() << " " << x.get_u_yy() <<  endl;
}
/**
 * Description not yet available.
 * \param
 */

/**
 * Description not yet available.
 * \param
 */
