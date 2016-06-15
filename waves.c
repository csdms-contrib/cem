/*
 * =====================================================================
 *
 *       Filename:  waves.c
 *
 *    Description:  Generate wave angles
 *
 *        Version:  1.0
 *        Created:  04/16/2010 02:42:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eric Hutton (eh), huttone@colorado.edu
 *        Company:  Community Surface Dynamics Modeling System
 *
 * =====================================================================
 */

#include <math.h>
#include <glib.h>

#include "waves.h"


#define WAVE_ANGLE_SIGN (1.0)


double
waves_next_angle (GRand * rand, double asymmetry, double highness)
{
  double angle;

  g_assert (highness < 1. && highness >= 0);
  g_assert (asymmetry < 1. && asymmetry >= 0);

  {
    double f = g_rand_double (rand);

    /*
     * Variable Asym will determine fractional distribution of waves coming
     * from the positive direction (positive direction coming from left)
     * -i.e. fractional wave asymmetry
     */
    f = g_rand_double (rand);
    if (f > highness)
      angle = (f - highness) / (1. - highness) * M_PI * .25;
    else
      angle = ((f / highness) + 1.) * M_PI * .25;

    if (g_rand_double (rand) > asymmetry)
      angle *= -1.;
  }

  return WAVE_ANGLE_SIGN * angle;
}
