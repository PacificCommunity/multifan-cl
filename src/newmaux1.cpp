/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"

    int operator > (fishery_header_record& fr1,
                                            fishery_header_record& fr)
    {
      if (fr1.year<fr.year)
      {
        return 0;
      }
      if (fr1.year>fr.year)
      {
        return 1;
      }
      if (fr1.month<fr.month)
      {
        return 0;
      }
      if (fr1.month>fr.month)
      {
        return 1;
      }
      if (fr1.week<fr.week)
      {
        return 0;
      }
      if (fr1.week>fr.week)
      {
        return 1;
      }
    // for now don't want fishery in the ordering
    /*
      if (fr1.fishery<fr.fishery)
      {
        return 0;
      }
      if (fr1.fishery>fr.fishery)
      {
        return 1;
      }
    */
      return 0;
    }

    int operator < (fishery_header_record& fr1,
                                            fishery_header_record& fr)
    {
      if (fr1.year>fr.year)
      {
        return 0;
      }
      if (fr1.year<fr.year)
      {
        return 1;
      }
      if (fr1.month>fr.month)
      {
        return 0;
      }
      if (fr1.month<fr.month)
      {
        return 1;
      }
      if (fr1.week>fr.week)
      {
        return 0;
      }
      if (fr1.week<fr.week)
      {
        return 1;
      }
    // for now don't want fishery in the ordering
    /*
      if (fr1.fishery>fr.fishery)
      {
        return 0;
      }
      if (fr1.fishery<fr.fishery)
      {
        return 1;
      }
    */
      return 0;
    }

#undef HOME_VERSION

