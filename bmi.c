#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "bmi.h"
#include "globals.h"


int cem_initialize (void);
int cem_update (void);
int cem_update_until (int);
int cem_finalize (void);


struct _BMI_CEM_Model {
  double time;
};


int
BMI_CEM_Initialize (const char *config_file, BMI_CEM_Model ** handle)
{
  BMI_CEM_Model *self = NULL;
  int status;

  if (!handle)
    return BMI_FAILURE;

  self = malloc (sizeof (BMI_CEM_Model));

  status = cem_initialize ();

  *handle = self;

  if (status == 0)
    return BMI_SUCCESS;
  else
    return BMI_FAILURE;
}


int
BMI_CEM_Update (BMI_CEM_Model *self)
{
  if (cem_update () == 0)
    return BMI_SUCCESS;
  else
    return BMI_FAILURE;
}


int
BMI_CEM_Update_frac (BMI_CEM_Model *self, double f)
{
  if (f>0) {
    double dt;

    BMI_CEM_Get_time_step (self, &dt);

    TimeStep = f * dt;

    BMI_CEM_Update (self);

    TimeStep = dt;
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Update_until (BMI_CEM_Model *self, double t)
{
  {
    double dt;
    double now;

    BMI_CEM_Get_time_step (self, &dt);
    BMI_CEM_Get_current_time (self, &now);

    {
      int n;
      const double n_steps = (t - now) / dt;
      const int n_full_steps = (int)n_steps;
      for (n=0; n<n_full_steps; n++) {
        BMI_CEM_Update (self);
      }

      // BMI_CEM_Update_frac (self, n_steps - n_full_steps);
    }
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Finalize (BMI_CEM_Model *self)
{
  if (self) {
    cem_finalize ();
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_var_type (BMI_CEM_Model *self, const char *long_var_name, BMI_Var_type * type)
{
  if (strcasecmp (long_var_name, "surface_elevation") == 0) {
    *type = BMI_VAR_TYPE_DOUBLE;
    return BMI_SUCCESS;
  }
  else {
    *type = BMI_VAR_TYPE_UNKNOWN;
    return BMI_FAILURE;
  }
}


int
BMI_CEM_Get_var_units (BMI_CEM_Model *self, const char *long_var_name, char * units)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    strncpy (units, "meter", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
  }
  else {
    units[0] = '\0';
    return BMI_FAILURE;
  }
}


int
BMI_CEM_Get_var_rank (BMI_CEM_Model *self, const char *long_var_name, int * rank)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    *rank = 2;
    return BMI_SUCCESS;
  }
  else {
    *rank = -1;
    return BMI_FAILURE;
  }
}


int
BMI_CEM_Get_grid_shape (BMI_CEM_Model *self, const char *long_var_name, int * shape)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    shape[0] = Xmax;
    shape[1] = Ymax;
  }

  return BMI_SUCCESS;
}



int
BMI_CEM_Get_grid_spacing (BMI_CEM_Model *self, const char *long_var_name, double * spacing)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    spacing[0] = CellWidth;
    spacing[1] = CellWidth;
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_grid_origin (BMI_CEM_Model *self, const char *long_var_name, double * origin)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    origin[0] = 0.;
    origin[1] = 0.;
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_grid_type (BMI_CEM_Model *self, const char *long_var_name, BMI_Grid_type * type)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    *type = BMI_GRID_TYPE_UNIFORM;
    return BMI_SUCCESS;
  }
  else {
    *type = BMI_GRID_TYPE_UNKNOWN;
    return BMI_FAILURE;
  }
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_value (BMI_CEM_Model *self, const char *long_var_name, void *dest)
{
  void *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //src = (void*) self->z[0];
  }

  //memcpy (dest, src, sizeof (double) * self->n_x * self->n_y);

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_value_ptr (BMI_CEM_Model *self, const char *long_var_name, void **dest)
{
  void *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //src = (void*) self->z[0];
  }

  *dest = src;

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_value_at_indices (BMI_CEM_Model *self, const char *long_var_name, void *dest, int * inds, int len)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //src = self->z[0];
  }

  { /* Copy the data */
    int i;
    double * to = (double*) dest;
    for (i=0; i<len; i++) {
      to[i] = src[inds[i]];
    }
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_double (BMI_CEM_Model *self, const char *long_var_name, double *dest)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //src = self->z[0];
  }

  //memcpy (dest, src, sizeof (double) * self->n_x * self->n_y);

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_double_ptr (BMI_CEM_Model *self, const char *long_var_name, double **dest)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //src = self->z[0];
  }

  *dest = src;

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_double_at_indices (BMI_CEM_Model *self, const char *long_var_name, double *dest, int * inds, int len)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //src = self->z[0];
  }

  { /* Copy the data */
    int i;
    for (i=0; i<len; i++) {
      dest[i] = src[inds[i]];
    }
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Set_value (BMI_CEM_Model *self, const char *long_var_name, void *array)
{
  if (strcmp (long_var_name, "surface_elevation")==0) {
    //memcpy (self->z[0], array, sizeof (double) * self->n_x * self->n_y);
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Set_value_at_indices (BMI_CEM_Model *self, const char *long_var_name, int * inds, int len, void *src)
{
  double * to;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //to = self->z[0];
  }

  { /* Copy the data */
    int i;
    double * from = src;
    for (i=0; i<len; i++) {
      to[inds[i]] = from[i];
    }
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Set_double (BMI_CEM_Model *self, const char *long_var_name, double *array)
{
  if (strcmp (long_var_name, "surface_elevation")==0) {
    //memcpy (self->z[0], array, sizeof (double) * self->n_x * self->n_y);
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Set_double_at_indices (BMI_CEM_Model *self, const char *long_var_name, int * inds, int len, double *src)
{
  double * dest;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    //dest = self->z[0];
  }

  { /* Copy the data */
    int i;
    for (i=0; i<len; i++) {
      dest[inds[i]] = src[i];
    }
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_component_name (BMI_CEM_Model *self, char * name)
{
  strncpy (name, "Coastline Evolution Model", BMI_MAX_COMPONENT_NAME);
  return BMI_SUCCESS;
}


const char *input_var_names[BMI_CEM_INPUT_VAR_NAME_COUNT] = {
  "surface_elevation"
};

int
BMI_CEM_Get_input_var_names (BMI_CEM_Model *self, char ** names)
{
  int i;
  for (i=0; i<BMI_CEM_INPUT_VAR_NAME_COUNT; i++) {
    strncpy (names[i], input_var_names[i], BMI_MAX_VAR_NAME);
  }
  return BMI_SUCCESS;
}


const char *output_var_names[BMI_CEM_OUTPUT_VAR_NAME_COUNT] = {
  "surface_elevation"
};

int
BMI_CEM_Get_output_var_names (BMI_CEM_Model *self, char ** names)
{
  int i;
  for (i=0; i<BMI_CEM_OUTPUT_VAR_NAME_COUNT; i++) {
    strncpy (names[i], output_var_names[i], BMI_MAX_VAR_NAME);
  }
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_start_time (BMI_CEM_Model *self, double * time)
{
  if (time) {
    *time = 0.;
    return BMI_SUCCESS;
  }
  else {
    return BMI_FAILURE;
  }
}


int
BMI_CEM_Get_end_time (BMI_CEM_Model *self, double * time)
{
  *time = StopAfter;
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_current_time (BMI_CEM_Model *self, double * time)
{
  *time = CurrentTimeStep * TimeStep;
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_time_step (BMI_CEM_Model *self, double * dt)
{
  *dt = TimeStep;
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_time_units (BMI_CEM_Model *self, char * units)
{
  strncpy (units, "d", BMI_MAX_UNITS_NAME);
  return BMI_SUCCESS;
}
