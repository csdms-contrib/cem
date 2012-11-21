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

#define WAVE_ANGLE_SIGN (0.7)
#define DEFAULT_ASYMMETRY (0.5)
#define DEFAULT_HIGHNESS (0.5)
#define DEFAULT_HEIGHT (2.)
#define DEFAULT_PERIOD (7.)
#define DEFAULT_SEED (1945)

double
waves_next_angle (GRand * rand, double asymmetry, double highness)
{
  double angle;

  g_assert (highness < 1. && highness >= 0);
  g_assert (asymmetry < 1. && asymmetry >= 0);

  {
    double f = g_rand_double (rand);

    double sign;

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

int
_waves_initialize (State * s)
{
  if (s)
    return TRUE;
  else
    return FALSE;
}

int
_waves_run_until (State * s, int until)
{
  int status = FALSE;

  g_assert (s);
  if (s)
  {
    const gint len = until - s->now;

    if (len > 0)
    {
      gint i;

      if (s->angles)
        s->angles = g_renew (double, s->angles, len);

      else
        s->angles = g_new (double, len);

      s->len = len;
      for (i = 0; i < len; i++)
        s->angles[i] = waves_next_angle (s->rand, s->asymmetry, s->highness);
      s->now = until;
      status = TRUE;
    }
    else if (len == 0 && s->angles == NULL)
    {
      s->angles = g_new (double, 1);

      s->len = 1;
      s->angles[0] = waves_next_angle (s->rand, s->asymmetry, s->highness);
      s->now = until;
      status = TRUE;
    }
    else
      status = FALSE;
  }
  else
    status = FALSE;

  return status;
}

int
_waves_finalize (State * s)
{
  return TRUE;
}

void
waves_init_state (State * s)
{
  g_assert (s);

  if (s)
  {
    s->asymmetry = DEFAULT_ASYMMETRY;
    s->highness = DEFAULT_HIGHNESS;
    s->height = DEFAULT_HEIGHT;
    s->period = DEFAULT_PERIOD;

    s->now = 0;
    s->end = 0;
    s->time_step = 1.;

    s->seed = DEFAULT_SEED;
    s->angles = NULL;
    s->len = 0;
    s->rand = g_rand_new_with_seed (s->seed);
  }

  return;
}

void
waves_free_state (State * s)
{
  if (s)
  {
    g_rand_free (s->rand);
    s->rand = NULL;
    g_free (s->angles);
    s->angles = NULL;
    s->len = 0;
  }
  return;
}
