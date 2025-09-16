/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#if !defined(__PHI_STUFF__)
#define __PHI_STUFF__
#include <admodel.h>
MY_DOUBLE_TYPE psiprime(MY_DOUBLE_TYPE p);
MY_DOUBLE_TYPE psiprime(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE e);
MY_DOUBLE_TYPE psi(MY_DOUBLE_TYPE p);
MY_DOUBLE_TYPE psiprime(MY_DOUBLE_TYPE p);
MY_DOUBLE_TYPE psi(MY_DOUBLE_TYPE p);
MY_DOUBLE_TYPE d2phi_deta2_1(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE d2phi_deta2_2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE d2phi_deta2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE dphi_deta_2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE dphi_deta_1(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE dphi_deta(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE dphi_dp(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
MY_DOUBLE_TYPE phi(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta);
#endif
