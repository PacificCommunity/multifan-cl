/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#if defined(NO_MY_DOUBLE_TYPE)
#  if defined(MY_DOUBLE_TYPE)
#    undef MY_DOUBLE_TYPE
#  endif
#  define MY_DOUBLE_TYPE double
#  if defined(MY_REAL_DOUBLE)
#    undef MY_REAL_DOUBLE_double
#  endif
#  define MY_REAL_DOUBLE_double
#endif
