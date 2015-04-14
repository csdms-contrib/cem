#include "waves.h"
#include "waves_model.h"


#define DEFAULT_ASYMMETRY (0.5)
#define DEFAULT_HIGHNESS (0.5)
#define DEFAULT_HEIGHT (2.)
#define DEFAULT_PERIOD (7.)
#define DEFAULT_SEED (1945)


static int
_waves_run_until (WavesModel * s, int until)
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


WavesModel*
waves_new (void)
{
  WavesModel *s = (WavesModel*)g_malloc (sizeof (WavesModel));

  waves_init_state (s);

  return (WavesModel *) s;
}


void
waves_init_state (WavesModel * s)
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
waves_free_state (WavesModel * s)
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


WavesModel *
waves_destroy (WavesModel * self)
{
  if (self)
  {
    waves_free_state (self);
    g_free (self);
  }
  return NULL;
}


WavesModel *
waves_init (WavesModel * self)
{
  if (!self)
    self = waves_new ();

  return self;
}


int
waves_run_until (WavesModel * self, double time_in_days)
{
  int until_time_step = time_in_days / self->time_step;

  return _waves_run_until (self, until_time_step);
}

WavesModel *
waves_finalize (WavesModel * self, int free)
{
  if (free)
    self = waves_destroy (self);

  return self;
}


WavesModel *
waves_set_angle_asymmetry (WavesModel * self, double asymmetry)
{
  self->asymmetry = asymmetry;
  return self;
}


WavesModel *
waves_set_angle_highness (WavesModel * self, double highness)
{
  self->highness = highness;
  return self;
}


WavesModel *
waves_set_height (WavesModel * self, double height_in_m)
{
  self->height = height_in_m;
  return self;
}


WavesModel *
waves_set_period (WavesModel * self, double period_in_s)
{
  self->period = period_in_s;
  return self;
}


double
waves_get_angle_asymmetry (const WavesModel * self)
{
  return self->asymmetry;
}


double
waves_get_angle_highness (const WavesModel * self)
{
  return self->highness;
}


double
waves_get_height (const WavesModel * self)
{
  return self->height;
}


double
waves_get_period (const WavesModel * self)
{
  return self->period;
}


double
waves_get_wave_angle (WavesModel * self)
{
  double angle = 0.;

  g_assert (self);
  {
    if (self->len > 0)
    {
      g_assert (self->len > 0 && self->angles != NULL);
      angle = self->angles[self->len - 1];
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
waves_get_wave_angle_max (WavesModel * self)
{
  double max = 0.;

  g_assert (self);
  {
    const gint len = self->len;

    if (len > 0)
    {
      gint i;

      max = 0;
      for (i = 0; i < len; i++)
        if (self->angles[i] > max)
          max = self->angles[i];
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
waves_get_wave_angle_min (WavesModel * self)
{
  double min = DBL_MAX;

  g_assert (self);
  {
    const gint len = self->len;

    if (len > 0)
    {
      gint i;

      for (i = 0; i < len; i++)
        if (self->angles[i] < min)
          min = self->angles[i];
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
waves_get_wave_angle_mean (WavesModel * self)
{
  double mean = 0.;

  g_assert (self);
  {
    const gint len = self->len;

    if (len > 0)
    {
      gint i;

      mean = 0.;
      for (i = 0; i < len; i++)
        mean += self->angles[i];
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
