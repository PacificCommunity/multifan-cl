/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#include <admodel.h>
void gram_schmidt_qr(dmatrix& M,dmatrix& Q,dmatrix& R);
void gram_schmidt_qr(dvar_matrix& M,dvar_matrix& Q,dvar_matrix& R);

dvar_vector max(double x,const dvar_vector& v);
  

class dfils_manager
{
public:
  int get_nobs(){return nobs;}
  int get_niter(){return niter;}
  dmatrix& get_M(){return M;}
  dvar_vector& get_obs(){return obs;}
  dvar_vector& get_pred_obs(){return pred_obs;}
  dvar_vector& get_weights(){return weights;}
  dmatrix& get_Q(){return Q;}
  dmatrix& get_R(){return R;}
  dvar_matrix& get_Q1(){return Q1;}
  dvar_matrix& get_R1(){return R1;}
  dvar_matrix& get_theta_hat(){return theta_hat;}
  //dvariable& get_vhat(){return vhat;}
  dfils_manager(dmatrix& _M,dvar_vector& _obs,int _niter);
  void fit_data(void);
private:
  int nobs;
  int niter;
  dmatrix M;
  dvector minval;
  dvar3_array M1;
  dvar_vector obs;
  dvar_vector obs1;
  dvar_vector pred_obs;
  dvar_vector weights;
  dmatrix Q;
  dmatrix R;
  dvar_matrix Q1;
  dvar_matrix R1;
  dvar_matrix theta_hat;
  dmatrix silly_pred;
  //dvariable vhat;
};

class ils_manager
{
public:
  int get_nobs(){return nobs;}
  int get_niter(){return niter;}
  dmatrix& get_M(){return M;}
  dvector& get_obs(){return obs;}
  dvector& get_pred_obs(){return pred_obs;}
  dvector& get_weights(){return weights;}
  dmatrix& get_Q(){return Q;}
  dmatrix& get_R(){return R;}
  dvector& get_theta_hat(){return theta_hat;}
  double& get_vhat(){return vhat;}
  ils_manager(dmatrix& _M,dvector& _obs,int _niter);
  double fit_data(void);

private:

  int nobs;
  int niter;
  dmatrix M;
  d3_array M1;
  dvector obs;
  dvector obs1;
  dvector pred_obs;
  dvector weights;
  dmatrix Q;
  dmatrix R;
  dvector theta_hat;
  double vhat;
};

