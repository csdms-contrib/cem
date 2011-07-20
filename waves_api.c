#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "waves.h"
#include "waves_api.h"

#ifdef SWIG
% include waves_api.h
#endif
  Waves_state * waves_new (void)
{
  State *s = malloc (sizeof (State));

  waves_init_state (s);

  return (Waves_state *) s;
}

Waves_state *
waves_destroy (Waves_state * self)
{
  if (self)
  {
    waves_free_state ((State *) self);
    free (self);
  }
  return NULL;
}

Waves_state *
waves_init (Waves_state * self)
{
  if (!self)
    self = waves_new ();

  _waves_initialize ((State *) self);

  return self;
}

int
waves_run_until (Waves_state * self, double time_in_days)
{
  State *p = (State *) self;

  int until_time_step = time_in_days / p->time_step;

  return _waves_run_until (p, until_time_step);
}

Waves_state *
waves_finalize (Waves_state * self, int free)
{
  _waves_finalize ((State *) self);

  if (free)
    self = waves_destroy (self);

  return self;
}

Waves_state *
waves_set_angle_asymmetry (Waves_state * self, double asymmetry)
{
  State *p = (State *) self;

  p->asymmetry = asymmetry;
  return self;
}

Waves_state *
waves_set_angle_highness (Waves_state * self, double highness)
{
  State *p = (State *) self;

  p->highness = highness;
  return self;
}

Waves_state *
waves_set_height (Waves_state * self, double height_in_m)
{
  State *p = (State *) self;

  p->height = height_in_m;
  return self;
}

Waves_state *
waves_set_period (Waves_state * self, double period_in_s)
{
  State *p = (State *) self;

  p->period = period_in_s;
  return self;
}

double
waves_get_angle_asymmetry (Waves_state * self)
{
  State *p = (State *) self;

  return p->asymmetry;
}

double
waves_get_angle_highness (Waves_state * self)
{
  State *p = (State *) self;

  return p->highness;
}

double
waves_get_height (Waves_state * self)
{
  State *p = (State *) self;

  return p->height;
}

double
waves_get_period (Waves_state * self)
{
  State *p = (State *) self;

  return p->period;
}

double
waves_get_wave_angle (Waves_state * self)
{
  double angle = 0.;

  g_assert (self);
  {
    State *p = (State *) self;

    if (p->len > 0)
    {
      g_assert (p->len > 0 && p->angles != NULL);
      angle = p->angles[p->len - 1];
    }
    else
    {
      waves_run_until (self, 0.);
      angle = waves_get_wave_angle (self);
    }
  }
  return angle;
}

double
waves_get_wave_angle_max (Waves_state * self)
{
  double max = 0.;

  g_assert (self);
  {
    State *p = (State *) self;

    const gint len = p->len;

    if (len > 0)
    {
      gint i;

      max = 0;
      for (i = 0; i < len; i++)
        if (p->angles[i] > max)
          max = p->angles[i];
    }
    else
    {
      waves_run_until (self, 0.);
      max = waves_get_wave_angle_max (self);
    }
  }

  return max;
}

double
waves_get_wave_angle_min (Waves_state * self)
{
  double min = DBL_MAX;

  g_assert (self);
  {
    State *p = (State *) self;

    const gint len = p->len;

    if (len > 0)
    {
      gint i;

      for (i = 0; i < len; i++)
        if (p->angles[i] < min)
          min = p->angles[i];
    }
    else
    {
      waves_run_until (self, 0.);
      min = waves_get_wave_angle_min (self);
    }
  }

  return min;
}

double
waves_get_wave_angle_mean (Waves_state * self)
{
  double mean = 0.;

  g_assert (self);
  {
    State *p = (State *) self;

    const gint len = p->len;

    if (len > 0)
    {
      gint i;

      mean = 0.;
      for (i = 0; i < len; i++)
        mean += p->angles[i];
      mean /= len;
    }
    else
    {
      waves_run_until (self, 0.);
      mean = waves_get_wave_angle_mean (self);
    }
  }

  return mean;
}

const char *_waves_exchange_items[] = {
  "INCOMING_ANGLE",
  "INCOMING_ANGLE_MAX",
  "INCOMING_ANGLE_MIN",
  "INCOMING_ANGLE_MEAN",
  NULL
};

const char **
waves_get_exchange_items (void)
{
  return _waves_exchange_items;
}

double
waves_get_value (Waves_state * self, const char *value)
{
  double val = 0;

  if (strcasecmp (value, "INCOMING_ANGLE") == 0 ||
      strcasecmp (value, "sea_surface_wave_from_direction") == 0)
    val = waves_get_wave_angle (self);
  else if (strcasecmp (value, "INCOMING_ANGLE_MEAN") == 0)
    val = waves_get_wave_angle_mean (self);
  else if (strcasecmp (value, "INCOMING_ANGLE_MAX") == 0)
    val = waves_get_wave_angle_max (self);
  else if (strcasecmp (value, "INCOMING_ANGLE_MIN") == 0)
    val = waves_get_wave_angle_min (self);
  else if (strcasecmp (value, "sea_surface_wave_height") == 0)
    val = waves_get_height (self);
  else if (strcasecmp (value, "sea_surface_wave_period") == 0)
    val = waves_get_period (self);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return val;
}

double
waves_get_current_time (Waves_state * self)
{
  State *p = (State *) self;

  return p->now;
}

double
waves_get_start_time (Waves_state * self)
{
  return 0;
}

double
waves_get_end_time (Waves_state * self)
{
  return DBL_MAX;
}
