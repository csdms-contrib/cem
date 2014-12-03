#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "bmi.h"


struct _BMI_CEM_Model {
  double dt;
  double t;
  double t_end;

  int n_x;
  int n_y;

  double dx;
  double dy;
  double **z;

  double **temp_z;
};


int
BMI_CEM_Initialize (const char *config_file, BMI_CEM_Model ** handle)
{
  BMI_CEM_Model *self = NULL;

  if (!handle)
    return BMI_FAILURE;

  self = malloc (sizeof (BMI_CEM_Model));

  if (config_file)
  { /* Read input file */
    FILE *fp = NULL;

    double dt = 0.;
    double t_end = 1.;
    int n_x = 0;
    int n_y = 0;

    fp = fopen (config_file, "r");
    if (!fp)
      return BMI_FAILURE;

    fscanf (fp, "%lf, %lf, %d, %d", &dt, &t_end, &n_x, &n_y);

    self->dt = dt;
    self->t_end = t_end;
    self->n_x = n_x;
    self->n_y = n_y;
  }
  else
  { /* Set to default values */
    self->dt = 1.;
    self->t_end = 10.;
    self->n_x = 10;
    self->n_y = 20;
  }

  self->dx = 1.;
  self->dy = 1.;

  { /* Initialize data */
    int i;
    const int len = self->n_x * self->n_y;
    double top_x = self->n_x - 1;

    /* Allocate memory */
    self->z = (double **)malloc (sizeof (double*) * self->n_y);
    self->temp_z = (double **)malloc (sizeof (double*) * self->n_y);

    if (!self->temp_z || !self->z)
      return BMI_FAILURE;

    self->z[0] = (double *)malloc (sizeof (double) * self->n_x * self->n_y);
    self->temp_z[0] = (double *)malloc (sizeof (double) * self->n_x * self->n_y);

    if (!self->temp_z[0] || !self->z[0])
      return BMI_FAILURE;

    for (i=1; i<self->n_y; i++) {
      self->z[i] = self->z[i-1] + self->n_x;
      self->temp_z[i] = self->temp_z[i-1] + self->n_x;
    }

    self->t = 0;
    for (i = 0; i < len; i++)
      self->z[0][i] = random ()*1./RAND_MAX * top_x * top_x * .5 - top_x * top_x * .25;
    for (i = 0; i < self->n_y; i++) {
      self->z[i][0] = 0.;
      self->z[i][self->n_x-1] = 0.;
    }
    for (i = 0; i < self->n_x; i++) {
      self->z[0][i] = 0.;
      self->z[self->n_y-1][i] = top_x*top_x*.25 - (i-top_x*.5) * (i-top_x*.5);
    }
    
    memcpy (self->temp_z[0], self->z[0], sizeof (double)*self->n_x*self->n_y);
  }

  *handle = self;

  return BMI_SUCCESS;
}


int
BMI_CEM_Update (BMI_CEM_Model *self)
{
  {
    int i, j;
    const double rho = 0.;
    const double dx2 = self->dx * self->dx;
    const double dy2 = self->dy * self->dy;
    const double dx2_dy2_rho = dx2 * dy2 * rho;
    const double denom = self->dt/(2 * (dx2 + dy2));
    double **z = self->z;

    for (i=1; i<self->n_y-1; i++)
      for (j=1; j<self->n_x-1; j++)
        self->temp_z[i][j] = denom * (dx2 * (z[i-1][j] + z[i+1][j]) +
                                      dy2 * (z[i][j-1] + z[i][j+1]) -
                                      dx2_dy2_rho);
    self->t += self->dt;
  }

  memcpy (self->z[0], self->temp_z[0], sizeof (double) * self->n_y * self->n_x);

  return BMI_SUCCESS;
}


int
BMI_CEM_Update_frac (BMI_CEM_Model *self, double f)
{
  if (f>0) {
    double dt;

    BMI_CEM_Get_time_step (self, &dt);

    self->dt = f * dt;

    BMI_CEM_Update (self);

    self->dt = dt;
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
      for (n=0; n<(int)n_steps; n++) {
        BMI_CEM_Update (self);
      }

      BMI_CEM_Update_frac (self, n_steps - (int)n_steps);
    }
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Finalize (BMI_CEM_Model *self)
{
  if (self)
  {
    free (self->temp_z[0]);
    free (self->temp_z);
    free (self->z[0]);
    free (self->z);
    free (self);
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
    shape[0] = self->n_y;
    shape[1] = self->n_x;
  }

  return BMI_SUCCESS;
}



int
BMI_CEM_Get_grid_spacing (BMI_CEM_Model *self, const char *long_var_name, double * spacing)
{
  if (strcmp (long_var_name, "surface_elevation") == 0) {
    spacing[0] = self->dy;
    spacing[1] = self->dx;
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
    src = (void*) self->z[0];
  }

  memcpy (dest, src, sizeof (double) * self->n_x * self->n_y);

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_value_ptr (BMI_CEM_Model *self, const char *long_var_name, void **dest)
{
  void *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    src = (void*) self->z[0];
  }

  *dest = src;

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_value_at_indices (BMI_CEM_Model *self, const char *long_var_name, void *dest, int * inds, int len)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    src = self->z[0];
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
    src = self->z[0];
  }

  memcpy (dest, src, sizeof (double) * self->n_x * self->n_y);

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_double_ptr (BMI_CEM_Model *self, const char *long_var_name, double **dest)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    src = self->z[0];
  }

  *dest = src;

  return BMI_SUCCESS;
}


int
BMI_CEM_Get_double_at_indices (BMI_CEM_Model *self, const char *long_var_name, double *dest, int * inds, int len)
{
  double *src = NULL;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    src = self->z[0];
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
    memcpy (self->z[0], array, sizeof (double) * self->n_x * self->n_y);
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Set_value_at_indices (BMI_CEM_Model *self, const char *long_var_name, int * inds, int len, void *src)
{
  double * to;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    to = self->z[0];
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
    memcpy (self->z[0], array, sizeof (double) * self->n_x * self->n_y);
  }

  return BMI_SUCCESS;
}


int
BMI_CEM_Set_double_at_indices (BMI_CEM_Model *self, const char *long_var_name, int * inds, int len, double *src)
{
  double * dest;

  if (strcmp (long_var_name, "surface_elevation")==0) {
    dest = self->z[0];
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
  strncpy (name, "Example C model", BMI_MAX_COMPONENT_NAME);
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
  *time = self->t_end;
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_current_time (BMI_CEM_Model *self, double * time)
{
  *time = self->t;
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_time_step (BMI_CEM_Model *self, double * dt)
{
  *dt = self->dt;
  return BMI_SUCCESS;
}


int
BMI_CEM_Get_time_units (BMI_CEM_Model *self, char * units)
{
  strncpy (units, "-", BMI_MAX_UNITS_NAME);
  return BMI_SUCCESS;
}
