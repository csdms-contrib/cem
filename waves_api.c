#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "waves.h"
#include "waves_api.h"

#ifdef SWIG
% include waves_api.h
#endif

//typedef struct Waves_state _BMI_Model;

#define BMI_FAILURE_BAD_FILE_ERROR (2)
#define BMI_FAILURE_SCAN_ERROR (3)

int
BMI_Initialize (const char *file, BMI_Model **handle)
{
  int rtn = BMI_FAILURE;

  if (handle) {
    Waves_state * self = waves_new ();
    double end_time = 20.;
    double wave_height = 2.;
    double wave_period = 7.;
    double angle_highness = 0.2;
    double angle_asymmetry = 0.5;

    _waves_initialize ((State *) self);

    if (file) {
      FILE *fp = fopen (file, "r");
      if (fp) {
        int n_assigned;
        n_assigned = fscanf (fp, "%lf, %lf, %lf, %lf, %lf", &end_time, &wave_height, &wave_period, &angle_highness, &angle_asymmetry);
        if (n_assigned!=5)
          return BMI_FAILURE_SCAN_ERROR;
      }
      else
        return BMI_FAILURE_BAD_FILE_ERROR;
    }

    BMI_Set_double (self, "sea_water_surface_wave__height", &wave_height);
    BMI_Set_double (self, "sea_water_surface_wave__period", &wave_period);
    BMI_Set_double (self, "sea_water_surface_wave__model_from_direction_highness_constant", &angle_highness);
    BMI_Set_double (self, "sea_water_surface_wave__model_from_direction_asymmetry_constant", &angle_asymmetry);

    {
      State *p = (State *) self;
      p->end = end_time / p->time_step;
      fprintf (stderr, "Setting end time to %d\n", p->end);
      fflush (stderr);
    }

    *handle = self;

    rtn = BMI_SUCCESS;
  }

  return rtn;
}

int
BMI_Update (BMI_Model * self)
{
  State *p = (State *) self;
  double now;

  BMI_Get_current_time (self, &now);
  _waves_run_until (p, now+1.);

  return BMI_SUCCESS;
}

int
BMI_Update_until (BMI_Model * self, double time_in_days)
{
  State *p = (State *) self;
  int until_time_step = time_in_days / p->time_step;

  _waves_run_until (p, until_time_step);

  return BMI_SUCCESS;
}

int
BMI_Finalize (BMI_Model *self)
{
  int rtn = BMI_FAILURE;

  if (self) {
    _waves_finalize ((State *) self);
    waves_destroy (self);
    rtn = BMI_SUCCESS;
  }

  return rtn;
}

int
BMI_Get_component_name (BMI_Model * self, char *name)
{
  if (name) {
    strcpy (name, "Waves");
    return BMI_SUCCESS;
  }
  return BMI_FAILURE;
}

#define BMI_INPUT_VAR_NAME_COUNT (4)
const char *input_var_names[BMI_INPUT_VAR_NAME_COUNT] = {
  "sea_water_surface_wave__model_from_direction_asymmetry_constant",
  "sea_water_surface_wave__model_from_direction_highness_constant",
  "sea_water_surface_wave__height",
  "sea_water_surface_wave__period",
};

int
BMI_Get_input_var_name_count (BMI_Model * self, int *count)
{
  if (count) {
    *count = BMI_INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
  }
  return BMI_FAILURE;
}

int
BMI_Get_input_var_names (const BMI_Model * self, char ** names)
{
  int i;
  for (i=0; i<BMI_INPUT_VAR_NAME_COUNT; i++) {
    strncpy (names[i], input_var_names[i], BMI_VAR_NAME_MAX);
  }
  return BMI_SUCCESS;
}

#define BMI_OUTPUT_VAR_NAME_COUNT (6)
const char *output_var_names[BMI_OUTPUT_VAR_NAME_COUNT] = {
  "sea_water_surface_wave__from_direction",
  "max_over_increment_of_sea_water_surface_wave__from_direction",
  "min_over_increment_of_sea_water_surface_wave__from_direction",
  "mean_over_increment_of_sea_water_surface_wave__from_direction",
  "sea_water_surface_wave__height",
  "sea_water_surface_wave__period"
};

int
BMI_Get_output_var_name_count (BMI_Model * self, int *count)
{
  if (count) {
    *count = BMI_OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
  }
  return BMI_FAILURE;
}

int
BMI_Get_output_var_names (const BMI_Model * self, char ** names)
{
  int i;
  for (i=0; i<BMI_OUTPUT_VAR_NAME_COUNT; i++) {
    strncpy (names[i], output_var_names[i], BMI_VAR_NAME_MAX);
  }
  return BMI_SUCCESS;
}

int
BMI_Get_var_rank (const BMI_Model *self, const char * name, int *rank)
{
  if (rank) {
    *rank = 1;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_var_type (const BMI_Model *self, const char * name, BMI_Var_type *type)
{
  if (type) {
    *type = BMI_VAR_TYPE_DOUBLE;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_var_point_count (const BMI_Model *self, const char * name, int *count)
{
  if (count) {
    *count = 1;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_grid_shape (const BMI_Model *self, const char * name, int *shape)
{
  if (shape) {
    shape[0] = 1;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_grid_type (const BMI_Model *self, const char * name, BMI_Grid_type *type)
{
  if (type) {
    *type = BMI_GRID_TYPE_UNIFORM;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_double (BMI_Model * self, const char *value, double *dest)
{
  int rtn = BMI_FAILURE;

  if (dest) {
    double val = 0;

    if (strcasecmp (value, "sea_water_surface_wave__from_direction") == 0)
      val = waves_get_wave_angle (self);
    else if (strcasecmp (value, "mean_over_increment_of_sea_water_surface_wave__from_direction") == 0)
      val = waves_get_wave_angle_mean (self);
    else if (strcasecmp (value, "max_over_increment_of_sea_water_surface_wave__from_direction") == 0)
      val = waves_get_wave_angle_max (self);
    else if (strcasecmp (value, "min_over_increment_of_sea_water_surface_wave__from_direction") == 0)
      val = waves_get_wave_angle_min (self);
    else if (strcasecmp (value, "sea_water_surface_wave__height") == 0)
      val = waves_get_height (self);
    else if (strcasecmp (value, "sea_water_surface_wave__period") == 0)
      val = waves_get_period (self);
    else {
      return BMI_FAILURE;
    }

    dest[0] = val;
    rtn = BMI_SUCCESS;
  }

  return rtn;
}

int
BMI_Get_double_ptr (BMI_Model * self, const char *value, double **dest)
{
  return -BMI_FAILURE;
}

int
BMI_Set_double (BMI_Model *self, const char *value, double *src)
{
  int rtn = BMI_FAILURE;

  if (src) {
    double val = src[0];

    if (strcasecmp (value, "sea_water_surface_wave__model_from_direction_asymmetry_constant") == 0)
      waves_set_angle_asymmetry (self, val);
    else if (strcasecmp (value, "sea_water_surface_wave__model_from_direction_highness_constant") == 0)
      waves_set_angle_highness (self, val);
    else if (strcasecmp (value, "sea_water_surface_wave__height") == 0)
      waves_set_height (self, val);
    else if (strcasecmp (value, "sea_water_surface_wave__period") == 0)
      waves_set_period (self, val);
    else {
      return BMI_FAILURE;
    }

    rtn = BMI_SUCCESS;
  }

  return rtn;
}

int
BMI_Get_current_time (const BMI_Model * self, double *time)
{
  if (time) {
    State *p = (State *) self;
    *time = p->now;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_start_time (const BMI_Model * self, double *time)
{
  if (time) {
    *time = 0.;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_end_time (const BMI_Model * self, double *time)
{
  if (time) {
    State *p = (State *) self;
    *time = p->end / p->time_step;
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

int
BMI_Get_time_units (const BMI_Model * self, char *units)
{
  if (self && units) {
    strcpy (units, "d");
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}

Waves_state * waves_new (void)
{
  State *s = (State*)malloc (sizeof (State));

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
waves_get_angle_asymmetry (const Waves_state * self)
{
  State *p = (State *) self;

  return p->asymmetry;
}

double
waves_get_angle_highness (const Waves_state * self)
{
  State *p = (State *) self;

  return p->highness;
}

double
waves_get_height (const Waves_state * self)
{
  State *p = (State *) self;

  return p->height;
}

double
waves_get_period (const Waves_state * self)
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
  "sea_water_surface_wave__from_direction",
  "max_over_increment_of_sea_water_surface_wave__from_direction",
  "min_over_increment_of_sea_water_surface_wave__from_direction",
  "mean_over_increment_of_sea_water_surface_wave__from_direction",
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

  if (strcasecmp (value, "sea_water_surface_wave__from_direction") == 0)
    val = waves_get_wave_angle (self);
  else if (strcasecmp (value, "mean_over_increment_of_sea_water_surface_wave__from_direction") == 0)
    val = waves_get_wave_angle_mean (self);
  else if (strcasecmp (value, "max_over_increment_of_sea_water_surface_wave__from_direction") == 0)
    val = waves_get_wave_angle_max (self);
  else if (strcasecmp (value, "min_over_increment_of_sea_water_surface_wave__from_direction") == 0)
    val = waves_get_wave_angle_min (self);
  else if (strcasecmp (value, "sea_water_surface_wave__height") == 0)
    val = waves_get_height (self);
  else if (strcasecmp (value, "sea_water_surface_wave__period") == 0)
    val = waves_get_period (self);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return val;
}

double
waves_get_current_time (const Waves_state * self)
{
  State *p = (State *) self;

  return p->now;
}

double
waves_get_start_time (const Waves_state * self)
{
  return 0;
}

double
waves_get_end_time (const Waves_state * self)
{
  return DBL_MAX;
}
