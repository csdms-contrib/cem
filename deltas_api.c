#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "deltas.h"
#include "deltas_api.h"

#ifdef SWIG
% include deltas_api.h
#endif
/** \file
\brief The Caperiffic API
*/

char *
scan_next_non_comment_line (char *line, int len, FILE *fp)
{
  char * rtn = NULL;

  {
    char * buffer = (char*) malloc (sizeof (char) * len);
    int found = 0;

    while (fgets (buffer, len, fp) != NULL) {
      if (buffer[0] != '#')
        found = 1;
        break;
    }

    if (found) {
      strncpy (line, buffer, len);
      rtn = line;
    }

    free (buffer);
  }

  return rtn;
}

int
BMI_CEM_Initialize (const char *config_file, BMI_Model **handle)
{
  int rtn = BMI_FAILURE;

  if (handle) {
    Deltas_state *s = NULL;
    int n_rows, n_cols;
    int sed_flux_flag = 0;
    int shape[2];
    double dx;
    double shoreface_slope;
    double shoreface_depth;
    double shelf_slope;

    if (config_file) {
      FILE *fp = fopen (config_file, "r");
      char buffer[2048];

      fprintf (stderr, "CEM: trying to open file: %s\n", config_file);

      if (!fp)
        return BMI_FAILURE + 3;

      scan_next_non_comment_line (buffer, 2048, fp);
      fprintf (stderr, "CEM: line: %s\n", buffer);
      sscanf (buffer, "%d, %d, %lf, %d\n", &n_rows, &n_cols, &dx, &sed_flux_flag); 

      scan_next_non_comment_line (buffer, 2048, fp);
      sscanf (buffer, "%lf, %lf, %lf", &shoreface_slope, &shoreface_depth, &shelf_slope);

      fprintf (stderr, "CEM: number of rows, columns: %d, %d\n", n_rows, n_cols);

      //fscanf (fp, "%d, %d, %lf, %d\n", &n_rows, &n_cols, &dx, &sed_flux_flag);
      //fscanf (fp, "%lf, %lf, %lf", &shoreface_slope, &shoreface_depth, &shelf_slope);
    }
    else {
      n_rows = 50;
      n_cols = 100;
      dx = 100.;
      sed_flux_flag = 1;

      shoreface_slope = 0.01;
      shoreface_depth = 10.0;
      shelf_slope = 0.001;
    }

    shape[0] = n_rows;
    shape[1] = n_cols;

    s = deltas_new ();
    if (!s)
      return BMI_FAILURE + 4;

    deltas_init_grid_shape (s, shape);
    deltas_init_cell_width (s, dx);

    deltas_init (s);

    if (sed_flux_flag)
      deltas_use_sed_flux (s);

    deltas_set_shoreface_slope (s, shoreface_slope);
    deltas_set_shoreface_depth (s, shoreface_depth);
    deltas_set_shelf_slope (s, shelf_slope);

    deltas_set_save_file (s, "test.txt"); 

    *handle = s;
    rtn = BMI_SUCCESS;
  }
  else
    rtn = BMI_FAILURE + 1;

  return rtn;
}

Deltas_state * deltas_new (void)
{
  State *s = (State*) malloc (sizeof (State));

  deltas_init_state (s);

  return (Deltas_state *) s;
}

Deltas_state *
deltas_destroy (Deltas_state * s)
{
  if (s)
  {

    deltas_free_state ((State *) s);
    free (s);
  }
  return NULL;
}

Deltas_state *
deltas_init (Deltas_state * s)
{
  if (!s)
  {
    s = deltas_new ();
  }

  _cem_initialize ((State *) s);

  return s;
}

int
BMI_CEM_Get_component_name (BMI_Model *s, char *name)
{
  if (*name) {
    strncpy (name, "CEM", BMI_CEM_COMPONENT_NAME_MAX);
    return BMI_SUCCESS;
  }
}

int
BMI_CEM_Update (BMI_Model * s)
{
  State *p = (State *) s;
  double now;

  BMI_CEM_Get_current_time (s, &now);

  //fprintf (stderr, "Update until %f\n", now+1);
  _cem_run_until (p, (now+1) / p->time_step);

  BMI_CEM_Get_current_time (s, &now);
  //fprintf (stderr, "Current time is %f\n", now);

  return BMI_SUCCESS;
}

int
BMI_CEM_Update_until (Deltas_state * s, double time_in_days)
{
  State *p = (State *) s;
  int until_time_step = time_in_days / p->time_step;

  _cem_run_until (p, until_time_step);

  return BMI_SUCCESS;
}

int
deltas_run_until (Deltas_state * s, double time_in_days)
{
  State *p = (State *) s;

  int until_time_step = time_in_days / p->time_step;

  return _cem_run_until (p, until_time_step);
}

int
BMI_CEM_Finalize (Deltas_state * s)
{
  if (s) {
    State *p = (State *) s;
    _cem_finalize (p);
    deltas_destroy (s);
  }

  return BMI_SUCCESS;
}

Deltas_state *
deltas_finalize (Deltas_state * s, int free)
{
  State *p = (State *) s;

  _cem_finalize (p);

  if (free)
    s = deltas_destroy (s);

  return s;
}

Deltas_state *
deltas_set_save_file (Deltas_state * s, const char *name)
{
  State *p = (State *) s;

  p->savefilename = strdup (name);
  return s;
}

Deltas_state *
deltas_set_read_file (Deltas_state * s, char *name)
{
  State *p = (State *) s;

  p->readfilename = strdup (name);
  return s;
}

Deltas_state *
deltas_init_grid_shape (Deltas_state * s, int dimen[2])
{
  State *p = (State *) s;

  fprintf (stderr, "*** Grid size is (%d,%d)\n", p->nx, p->ny);
  fprintf (stderr, "*** Requested size is (%d,%d)\n", dimen[0], dimen[1]);
  if (s && (p->nx == 0 && p->ny == 0))
  {
    int i;

    //const int len = dimen[0] * (dimen[1] * 2);
    const int len = dimen[0] * dimen[1];

    //const int stride = dimen[1] * 2;
    const int stride = dimen[1];

    p->nx = dimen[0];
    p->ny = dimen[1]/2;
    p->max_beach_len = len;

    p->AllBeach = (char **)malloc (sizeof (char *) * p->nx);
    p->PercentFull = (double **)malloc (sizeof (double *) * p->nx);
    p->Age = (int **)malloc (sizeof (int *) * p->nx);
    p->CellDepth = (double **)malloc (sizeof (double *) * p->nx);
    p->InitDepth = (double **)malloc (sizeof (double *) * p->nx);

    p->AllBeach[0] = (char *)malloc (sizeof (char) * len);
    p->PercentFull[0] = (double *)malloc (sizeof (double) * len);
    p->Age[0] = (int *)malloc (sizeof (int) * len);
    p->CellDepth[0] = (double *)malloc (sizeof (double) * len);
    p->InitDepth[0] = (double *)malloc (sizeof (double) * len);

    for (i = 1; i < p->nx; i++)
    {
      p->AllBeach[i] = p->AllBeach[i - 1] + stride;
      p->PercentFull[i] = p->PercentFull[i - 1] + stride;
      p->Age[i] = p->Age[i - 1] + stride;
      p->CellDepth[i] = p->CellDepth[i - 1] + stride;
      p->InitDepth[i] = p->InitDepth[i - 1] + stride;
    }

    p->river_flux = (double *)malloc (sizeof (double) * len);
    p->river_x_ind = (int *)malloc (sizeof (int) * len);
    p->river_y_ind = (int *)malloc (sizeof (int) * len);
    p->river_x = (double *)malloc (sizeof (double) * len);
    p->river_y = (double *)malloc (sizeof (double) * len);
    p->n_rivers = 1;

    p->X = (int *)malloc (sizeof (int) * len);
    p->Y = (int *)malloc (sizeof (int) * len);
    p->InShadow = (char *)malloc (sizeof (char) * len);
    p->ShorelineAngle = (double *)malloc (sizeof (double) * len);
    p->SurroundingAngle = (double *)malloc (sizeof (double) * len);
    p->UpWind = (char *)malloc (sizeof (char) * len);
    p->VolumeIn = (double *)malloc (sizeof (double) * len);
    p->VolumeOut = (double *)malloc (sizeof (double) * len);
  }
  fprintf (stderr, "*** New grid size is (%d,%d)\n",
           deltas_get_nx (s), deltas_get_ny (s));

  return s;
}

Deltas_state *
deltas_init_cell_width (Deltas_state * s, double dx)
{
  State *p = (State *) s;

  if (s)
  {
    p->cell_width = dx;
  }
  return s;
}

Deltas_state *
deltas_init_grid (Deltas_state * s, double *z)
{
  State *p = (State *) s;

  if (s && (p->nx > 0 && p->ny > 0))
  {
    int i;

    const int len = p->nx * p->ny * 2;

    for (i = 0; i < len; i++)
      p->InitDepth[0][i] = z[i];
  }

  return s;
}

Deltas_state *
deltas_destroy_grid (Deltas_state * s)
{
  State *p = (State *) s;

  if (p)
  {
    p->nx = 0;
    p->ny = 0;

    free (p->AllBeach[0]);
    free (p->AllBeach);

    free (p->Age[0]);
    free (p->Age);

    free (p->PercentFull[0]);
    free (p->PercentFull);

    free (p->CellDepth[0]);
    free (p->CellDepth);

    free (p->InitDepth[0]);
    free (p->InitDepth);
  }

  return s;
}

double
deltas_get_sed_rate (Deltas_state * s)
{
  State *p = (State *) s;

  return p->SedRate;
}

Deltas_state *
deltas_find_river_mouth (Deltas_state * s, int n)
{
  State *p = (State *) s;

  int x,
    y;

  x = 0;
  //y = p->stream_spot;
  y = deltas_get_ny (s) / 2;

  while (p->AllBeach[x][y] == 'y')
  {
    x += 1;
  }

  p->river_x_ind[n] = x;
  p->river_y_ind[n] = y;
  p->river_x[n] = x*deltas_get_dx (s);
  p->river_y[n] = y*deltas_get_dy (s);

//fprintf (stderr, "x = %d\n", x);
//fprintf (stderr, "y = %d\n", y);

  return s;
}

void
deltas_avulsion (Deltas_state * s, double *qs, double river_flux)
{
  State *p = (State *) s;

  int x,
    y;

  int i;

  int qs_x,
    qs_y,
    qs_i;

  int len;

  deltas_find_river_mouth (s, 0);

  x = p->river_x_ind[0];
  y = p->river_y_ind[0];

  len = deltas_get_nx (s) * deltas_get_ny (s) / 2;
  for (i = 0; i < len; i++)
    qs[i] = 0;

  qs_x = x;
  //qs_y = y % deltas_get_ny (s) - deltas_get_ny (s) / 4;
  qs_y = y;
  qs_i = qs_x * deltas_get_ny (s) / 2 + qs_y;

  qs[qs_i] = river_flux;

  return;
}

Deltas_state *
deltas_set_sed_rate (Deltas_state * s, double rate)
{
  State *p = (State *) s;

  p->SedRate = rate;
  return s;
}

/** Set sediment flux in kg/s
*/
Deltas_state *
deltas_set_sed_flux (Deltas_state * s, double flux)
{
  State *p = (State *) s;

  p->SedFlux = flux;
  return s;
}

/** Set sediment flux in kg/s
*/
Deltas_state *
deltas_set_river_sed_flux (Deltas_state * s, double flux, int n)
{
  State *p = (State *) s;

  p->river_flux[n] = flux;
  return s;
}

/** Set sediment flux in kg/s
*/
Deltas_state *
deltas_set_river_position (Deltas_state * s, int x, int y, int n)
{
  State *p = (State *) s;

  p->river_x_ind[n] = x;
  p->river_y_ind[n] = y;
  p->river_x[n] = x*deltas_get_dx (s);
  p->river_y[n] = y*deltas_get_dy (s);
  return s;
}

Deltas_state *
deltas_set_sediment_flux_grid_old (Deltas_state * s, double *qs)
{
  State *p = (State *) s;

  int i;

  int n;

  int len;
  const int lower[2] = { deltas_get_ny (s) / 4, 0 };
  const int stride[2] = { 1, deltas_get_ny (s) / 2 };
  int dimen[3];
/*
  fprintf (stderr, "Set flux grid\n");
  fprintf (stderr, "Start.\n");
  fflush (stderr);
*/
  deltas_get_value_dimen (s, NULL, dimen);

  len = dimen[0] * dimen[1] * dimen[2];

  for (i = 0, n = 0; i < len; i++)
  {
    if (qs[i] > 0)
    {
      p->river_flux[n] = qs[i];
      p->river_x_ind[n] = i / stride[1];
      p->river_y_ind[n] = i % stride[1] + lower[0];
      p->river_x[n] = (i / stride[1])*deltas_get_dx (s);
      p->river_y[n] = (i % stride[1] + lower[0])*deltas_get_dy (s);
      /*
      fprintf (stderr, "  river position = %d, %d\n", p->river_x[n],
               p->river_y[n]);

      fprintf (stderr, "n = %d\n", n);
      fprintf (stderr, "i = %d\n", i);
      fprintf (stderr, "river_x = %d\n", p->river_x[n]);
      fprintf (stderr, "river_y = %d\n", p->river_y[n]);
      fprintf (stderr, "qs[%d] = %f\n", i, qs[i]);
      */
      n++;
    }
  }
  p->n_rivers = n;
  /*
  fprintf (stderr, "Number of rivers = %d\n", n);
  fprintf (stderr, "Done.\n");
  fflush (stderr);
  */

  return s;
}

Deltas_state *
deltas_set_sediment_flux_grid (Deltas_state * s, double *qs)
{
  State *p = (State *) s;

  int i;

  int n;

  int len;
  const int lower[2] = { 0, 0 };
  const int stride[2] = { 1, deltas_get_ny (s) };
  int dimen[3];
  const int qs_stride[2] = {deltas_get_ny (s)/2, 1};
  const int qs_lower[2] = {deltas_get_ny (s)/4, 0};

  deltas_get_value_dimen (s, NULL, dimen);

  len = dimen[0] * dimen[1] * dimen[2] / 2;

  fprintf (stderr, "Setting sediment flux grid\n");

  for (i = 0, n = 0; i < len; i++)
  {
    if (qs[i] > 0)
    {
      fprintf (stderr, "Found non-zero flux at %d\n", i);

      p->river_flux[n] = qs[i];

      p->river_x_ind[n] = i / qs_stride[0];
      p->river_y_ind[n] = i % qs_stride[0] + qs_lower[0];
      fprintf (stderr, "Found non-zero flux at x_ind = %d\n", p->river_x_ind[n]);
      fprintf (stderr, "Found non-zero flux at y_ind = %d\n", p->river_y_ind[n]);

      //p->river_x_ind[n] = i / stride[1];
      //p->river_y_ind[n] = i % stride[1] + lower[0];
      p->river_x[n] = p->river_x_ind[n]*deltas_get_dx (s);
      p->river_y[n] = p->river_y_ind[n]*deltas_get_dy (s);
      n++;
    }
  }
  p->n_rivers = n;

  return s;
}

Deltas_state *
deltas_set_rivers (Deltas_state * s, const double * x, const double * y,
                   double * qb, const int len)
{
  State *p = (State *) s;
  int i;
  const double dx = deltas_get_dx (s);
  const double dy = deltas_get_dy (s);
//fprintf (stderr, "DEBUG: n_rivers=%d\n", len);
  p->n_rivers = len;
  for (i=0; i<len; i++)
  {
//fprintf (stderr, "DEBUG: x[%d]=%d\n", i, (int)(x[i]/dx));
//fprintf (stderr, "DEBUG: y[%d]=%d\n", i, (int)(y[i]/dy));
//fprintf (stderr, "DEBUG: qb[%d]=%f\n", i, qb[i]);

    p->river_flux[i] = qb[i];
    p->river_x_ind[i] = x[i]/dx;
    p->river_y_ind[i] = y[i]/dy;
    p->river_x[i] = x[i];
    p->river_y[i] = y[i];
  }

  return s;
}

Deltas_state *
deltas_set_angle_asymmetry (Deltas_state * s, double angle_asymmetry)
{
  State *p = (State *) s;

  p->angle_asymmetry = angle_asymmetry;
  return s;
}

Deltas_state *
deltas_set_angle_highness (Deltas_state * s, double angle_highness)
{
  State *p = (State *) s;

  p->angle_highness = angle_highness;
  return s;
}

Deltas_state *
deltas_set_wave_angle (Deltas_state * s, double wave_angle)
{
  State *p = (State *) s;

  p->WaveAngle = wave_angle;
  return s;
}

Deltas_state *
deltas_set_depth (Deltas_state * s, double *depth)
{
  State *p = (State *) s;

  memcpy (p->CellDepth[0], depth, sizeof (p->CellDepth));
  return s;
}

Deltas_state *
deltas_set_wave_height (Deltas_state * s, double height_in_m)
{
  State *p = (State *) s;

  p->wave_height = height_in_m;
  return s;
}

Deltas_state *
deltas_set_wave_period (Deltas_state * s, double period_in_s)
{
  State *p = (State *) s;

  p->wave_period = period_in_s;
  return s;
}

Deltas_state *
deltas_set_shoreface_slope (Deltas_state * s, double shoreface_slope)
{
  State *p = (State *) s;

  p->shoreface_slope = shoreface_slope;
  return s;
}

Deltas_state *
deltas_set_shelf_slope (Deltas_state * s, double shelf_slope)
{
  State *p = (State *) s;

  p->shelf_slope = shelf_slope;
  return s;
}

Deltas_state *
deltas_set_shoreface_depth (Deltas_state * s, double shoreface_depth)
{
  State *p = (State *) s;

  p->shoreface_depth = shoreface_depth;
  return s;
}

Deltas_state *
deltas_set_cell_width (Deltas_state * s, double cell_width)
{
  State *p = (State *) s;

  p->cell_width = cell_width;
  return s;
}

const char *_deltas_exchange_items[] = {
  "sea_water__depth",
  "sea_water_to_sediment__depth_ratio",
  NULL
};

const char **
deltas_get_exchange_items (void)
{
  return _deltas_exchange_items;
}

const double *
deltas_get_value_grid (Deltas_state * s, const char *value)
{
  if (strcasecmp (value, "sea_water__depth") == 0)
    return deltas_get_depth (s);
  else if (strcasecmp (value, "sea_water_to_sediment__depth_ratio") == 0)
    return deltas_get_percent (s);
  else if (strcasecmp (value, "surface__elevation") == 0)
    return deltas_get_elevation_dup (s);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return NULL;
}

double *
deltas_get_value_grid_dup (Deltas_state * s, const char *value)
{
  if (strcasecmp (value, "sea_water__depth") == 0)
    return deltas_get_depth_dup (s);
  else if (strcasecmp (value, "sea_water_to_sediment__depth_ratio") == 0)
    return deltas_get_percent_dup (s);
  else if (strcasecmp (value, "surface__elevation") == 0)
    return deltas_get_elevation_dup (s);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return NULL;
}

const double *
deltas_get_value_data (Deltas_state * s, const char *value, int lower[2],
                       int upper[2], int stride[2])
{
  double *data = NULL;

  if (strncasecmp (value, "river_mouth", 11)==0)
  {
    lower[0] = 0;
    upper[0] = deltas_get_n_rivers (s) - 1;
    stride[0] = 1;

    lower[1] = 0;
    upper[1] = 0;
    stride[1] = 1;
  }
  else
  {
    lower[0] = 0;
    upper[0] = deltas_get_ny (s) - 1;
    stride[0] = 1;
    lower[1] = 0;
    upper[1] = deltas_get_nx (s) - 1;
    stride[1] = deltas_get_ny (s);
  }

  if (strcasecmp (value, "sea_water__depth") == 0)
    data = (double*)deltas_get_depth (s);
  else if (strcasecmp (value, "sea_water_to_sediment__depth_ratio") == 0)
    data = (double*)deltas_get_percent (s);
//  else if (strcasecmp (value, "ELEVATION") == 0)
//    data = (double*)deltas_get_elevation (s);
  else if (strcasecmp (value, "river_mouth__location_model_x_component") == 0)
  {
    data = (double*)deltas_get_river_x_position (s);
    fprintf (stderr, "cem_api: x[0]=%f\n", data[0]);
  }
  else if (strcasecmp (value, "river_mouth__location_model_y_component") == 0)
  {
    data = (double*)deltas_get_river_y_position (s);
    fprintf (stderr, "cem_api: y[0]=%f\n", data[0]);
  }
  else if (strcasecmp (value, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0)
    data = (double*)deltas_get_river_flux (s);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return (const double*)data;
}

int
BMI_CEM_Get_var_stride (Deltas_state *s, const char *name, int *stride)
{
  int rtn = BMI_FAILURE;

  if (stride) {
    State *p = (State *) s;

    fprintf (stderr, "getting stride for %s\n", name);
    if (strcmp (name, "surface_bed_load_sediment__mass_flow_rate") == 0 ||
        strcmp (name, "sea_water_to_sediment__depth_ratio") == 0 ||
        strcmp (name, "sea_water__depth") == 0) {
      stride[0] = p->ny * 2;
      stride[1] = 1;
    fprintf (stderr, "stride is %d, %d\n", stride[0], stride[1]);

    }
    else if (strcmp (name, "river_mouth__location_model_x_component") == 0 ||
             strcmp (name, "river_mouth__location_model_y_component") == 0 ||
             strcmp (name, "river_mouth_bed_load_sediment__mass_flow_rate") == 0) {
      stride[0] = 1;
    }

    rtn = BMI_SUCCESS;
  }

  return rtn;
}

#define BMI_OUTPUT_VAR_NAME_COUNT (6)

const char *output_var_names[BMI_OUTPUT_VAR_NAME_COUNT] = {
  "surface__elevation",
  "sea_water__depth",
  "sea_water_to_sediment__depth_ratio",
  "river_mouth__location_model_x_component",
  "river_mouth__location_model_y_component",
  "river_mouth_bed_load_sediment__mass_flow_rate",
  //"surface_bed_load_sediment__mass_flow_rate",
};

int
BMI_CEM_Get_output_var_names (Deltas_state *s, char **names)
{
  int i;
  for (i=0; i<BMI_OUTPUT_VAR_NAME_COUNT; i++)
    strncpy (names[i], output_var_names[i], BMI_CEM_VAR_NAME_MAX);
  return BMI_SUCCESS;
}

int
BMI_CEM_Get_output_var_name_count (Deltas_state *s, int *count)
{
  if (count) {
    *count = BMI_OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
  }
  return BMI_FAILURE;
}

#define BMI_INPUT_VAR_NAME_COUNT (6)

const char *input_var_names[BMI_INPUT_VAR_NAME_COUNT] = {
  "surface_bed_load_sediment__mass_flow_rate",
  "channel_outflow_end_bed_load_sediment__mass_flow_rate",
  "channel_outflow_end_suspended_load__mass_flow_rate",
  "sea_water_surface_wave__from_direction",
  "sea_water_surface_wave__height",
  "sea_water_surface_wave__period"
};

int
BMI_CEM_Get_input_var_names (Deltas_state *s, char **names)
{
  int i;
  for (i=0; i<BMI_INPUT_VAR_NAME_COUNT; i++) {
    strncpy (names[i], input_var_names[i], BMI_CEM_VAR_NAME_MAX);
  }
  return BMI_SUCCESS;
}

int
BMI_CEM_Get_input_var_name_count (Deltas_state *s, int *count)
{
  if (count) {
    *count = BMI_INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
  }
  return BMI_FAILURE;
}

int
BMI_CEM_Get_double (Deltas_state *s, const char *value, double *dest)
{
  int rtn = BMI_FAILURE;

  if (dest) {
    double * src = NULL;

    if (strcmp (value, "sea_water__depth") == 0 ||
        strcmp (value, "surface__elevation") == 0) {
      src = (double*)deltas_get_depth (s);
    }
    else if (strcmp (value, "sea_water_to_sediment__depth_ratio") == 0) {
      src = (double*)deltas_get_percent (s);
//  else if (strcasecmp (value, "ELEVATION") == 0)
//    data = (double*)deltas_get_elevation (s);
    }
    else if (strcmp (value, "river_mouth__location_model_x_component") == 0) {
      src = (double*)deltas_get_river_x_position (s);
    }
    else if (strcmp (value, "river_mouth__location_model_y_component") == 0) {
      src = (double*)deltas_get_river_y_position (s);
    }
    //else if (strcmp (value, "surface_bed_load_sediment__mass_flow_rate") == 0) {
    else if (strcmp (value, "river_mouth_bed_load_sediment__mass_flow_rate") == 0) {
      src = (double*)deltas_get_river_flux (s);
    }

    if (src) { /* Copy the subgrid to the destination array */
      if (strcmp (value, "river_mouth_bed_load_sediment__mass_flow_rate") == 0 ||
          strcmp (value, "river_mouth_bed_load_sediment__mass_flow_rate") == 0 ||
          strcmp (value, "river_mouth_bed_load_sediment__mass_flow_rate") == 0) {
        int i;
        const int len = deltas_get_n_rivers (s);
        for (i=0; i<len; i++)
          dest[i] = src[i];
      }
      else {
        State *p = (State *) s;
        int i, j;
        const int n_rows = p->nx;
        const int n_cols = p->ny;
        const int stride = p->ny * 2;
        double * src_row = src + p->ny / 2;
        double * dest_row = dest;

        if (strcmp (value, "surface__elevation") == 0) {
          const double scale = -1.;
          for (i=0; i<n_rows; i++) {
            for (j=0; j<n_cols; j++)
              dest_row[j] = scale * src_row[j];
            dest_row += n_cols;
            src_row += stride;
          }
        }
        else {
          for (i=0; i<n_rows; i++) {
            for (j=0; j<n_cols; j++)
              dest_row[j] = src_row[j];
            dest_row += n_cols;
            src_row += stride;
          }
        }
      }

      rtn = BMI_SUCCESS;
    }
  }

  return rtn;
}

int
BMI_CEM_Get_double_ptr (Deltas_state *s, const char *value, double **dest)
{
  int rtn = BMI_FAILURE;

  if (dest) {
    double * src = NULL;

    if (strcasecmp (value, "sea_water__depth") == 0) {
      src = (double*)deltas_get_depth (s);
    }
    else if (strcasecmp (value, "sea_water_to_sediment__depth_ratio") == 0) {
      src = (double*)deltas_get_percent (s);
    }
    else if (strcasecmp (value, "river_mouth_bed_load_sediment__mass_flow_rate") == 0) {
      src = (double*)deltas_get_river_flux (s);
    }
    else if (strcasecmp (value, "surface__elevation") == 0) {
      return -BMI_FAILURE;
    }

    if (src) { /* Get a pointer to the start of the data */
      State *p = (State *) s;
      *dest = src + p->ny / 2;
      rtn = BMI_SUCCESS;

      fprintf (stderr, "dest = %d\n", *dest);
      fprintf (stderr, "data = %f, %f, %f\n", (*dest)[0], (*dest)[1], (*dest)[2]);
    }
    else {
      if (strcasecmp (value, "river_mouth__location_model_x_component") == 0) {
        src = (double*)deltas_get_river_x_position (s);
      }
      else if (strcasecmp (value, "river_mouth__location_model_y_component") == 0) {
        src = (double*)deltas_get_river_y_position (s);
      }
      else
        return BMI_FAILURE;
      rtn = BMI_SUCCESS;
    }

  }

  return rtn;
}

int
BMI_CEM_Set_double (Deltas_state *s, const char * value, double *src)
{
  int rtn = BMI_FAILURE;

  if (src) {
    if (strcmp (value, "surface_bed_load_sediment__mass_flow_rate") == 0) {
      deltas_set_sediment_flux_grid (s, src);
    }
    else if (strcmp (value, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0 ||
        strcmp (value, "channel_outflow_end_suspended_load__mass_flow_rate") == 0) {
      deltas_set_sed_flux (s, src[0]);
    }
    else if (strcmp (value, "sea_water_surface_wave__from_direction") == 0) {
      deltas_set_wave_angle (s, src[0]);
    }
    else if (strcmp (value, "sea_water_surface_wave__height") == 0) {
      deltas_set_wave_height (s, src[0]);
    }
    else if (strcmp (value, "sea_water_surface_wave__period") == 0) {
      deltas_set_wave_period (s, src[0]);
    }
    else {
      return BMI_FAILURE;
    }

    rtn = BMI_SUCCESS;
  }

  return rtn;
}

double *
deltas_get_value_data_dup (Deltas_state * s, const char *value, int lower[2],
                       int upper[2], int stride[2])
{
  double *data = NULL;

/*
  lower[0] = deltas_get_ny (s)/4;
  lower[1] = 0;
  upper[0] = 3*deltas_get_ny (s)/4-1;
  upper[1] = deltas_get_nx (s)-1;
  stride[0] = 1;
  stride[1] = deltas_get_ny (s);
*/
  lower[0] = 0;
  lower[1] = 0;
  upper[0] = deltas_get_ny (s) / 2 - 1;
  upper[1] = deltas_get_nx (s) - 1;
  stride[0] = 1;
  stride[1] = deltas_get_ny (s) / 2;

  if (strcasecmp (value, "sea_water__depth") == 0)
    data = deltas_get_depth_dup (s);
  else if (strcasecmp (value, "sea_water_to_sediment__depth_ratio") == 0)
    data = deltas_get_percent_dup (s);
  else if (strcasecmp (value, "surface__elevation") == 0)
    data = deltas_get_elevation_dup (s);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return data;
}

int *
deltas_get_value_dimen_old (Deltas_state * s, const char *value, int shape[3])
{
  shape[0] = deltas_get_ny (s) / 2;
  //shape[0] = deltas_get_ny (s);
  shape[1] = deltas_get_nx (s);
  shape[2] = 1;

  return shape;
}

int *
deltas_get_value_dimen (Deltas_state * s, const char *value, int shape[3])
{
  if (value && strncasecmp (value, "river_mouth", 11)==0)
  {
    shape[0] = deltas_get_n_rivers (s);
    shape[1] = 1;
    shape[2] = 1;
  }
  else
  {
    shape[0] = deltas_get_ny (s);
    //shape[0] = deltas_get_ny (s);
    shape[1] = deltas_get_nx (s);
    shape[2] = 1;
  }

  return shape;
}

double *
deltas_get_value_res (Deltas_state * s, const char *value, double res[3])
{
  if (strncasecmp (value, "river_mouth", 11)==0)
  {
    res[0] = 1;
    res[1] = 1;
    res[2] = 1;
  }
  else
  {
    res[0] = deltas_get_dy (s);
    res[1] = deltas_get_dx (s);
    res[2] = 1;
  }

  return res;
}

const double *
deltas_get_depth (Deltas_state * s)
{
  State *p = (State *) s;

  return p->CellDepth[0];
}

double *
dup_subgrid (Deltas_state * s, double **src)
{
  double *dest = NULL;

  {
    int lower[2] = { deltas_get_ny (s) / 4, 0 };
    int upper[2] = { 3 * deltas_get_ny (s) / 4 - 1, deltas_get_nx (s) - 1 };
    int stride[2] = { 1, deltas_get_ny (s) };
    const int len = (upper[0] - lower[0] + 1) * (upper[1] - lower[1] + 1);

    dest = (double *)malloc (sizeof (double) * len);

    if (dest)
    {
      //int src_id;
      //int dst_id;
      int id;

      int i,
        j;

/*
      for (src_id=lower[0], dst_id=0, j=lower[1]; j<=upper[1];
           src_id+=(stride[1]-(upper[0]-lower[0])), j++)
        for (i=lower[0]; i<=upper[0]; src_id+=stride[0], dst_id++, i++)
          dest[dst_id] = src[src_id];
*/
      for (i = lower[1], id = 0; i <= upper[1]; i++)
        for (j = lower[0]; j <= upper[0]; j++, id++)
          dest[id] = src[i][j];
    }
  }

  return dest;
}

double *
deltas_get_depth_dup (Deltas_state * s)
{
  State *p = (State *) s;

  double *val = NULL;

  {
/*
    const int len = deltas_get_nx (s)*deltas_get_ny (s);
    const double* f_val = deltas_get_depth (s);

    val = (double*) malloc (sizeof(double)*len);
    if (val)
    {
      int i;
      for (i=0; i<len; i++)
        val[i] = f_val[i];
    }
*/
    val = dup_subgrid (s, p->CellDepth);
  }

  return val;
}

double *
deltas_get_elevation_dup (Deltas_state * s)
{
  State *p = (State *) s;

  double *val = NULL;

  {
    const int len = deltas_get_nx (s) * deltas_get_ny (s) / 2;

    int i;

    //val = deltas_get_depth_dup (s);
    val = dup_subgrid (s, p->CellDepth);
    for (i = 0; i < len; i++)
      val[i] *= -1;
  }

  return val;
}

const double *
deltas_get_percent (Deltas_state * s)
{
  State *p = (State *) s;

  return p->PercentFull[0];
}

double *
deltas_get_percent_dup (Deltas_state * s)
{
  State *p = (State *) s;

  double *val = NULL;

  {
/*
    const int len = deltas_get_nx (s)*deltas_get_ny (s);
    const double* f_val = deltas_get_percent (s);

    val = (double*) malloc (sizeof(double)*len);

    if (val)
    {
      int i;
      for (i=0; i<len; i++)
        val[i] = f_val[i];
    }
*/
    val = dup_subgrid (s, p->PercentFull);
  }

  return val;
}

const double*
deltas_get_river_x_position (Deltas_state * s)
{
  State *p = (State *) s;
  return p->river_x;
}

const double*
deltas_get_river_y_position (Deltas_state * s)
{
  State *p = (State *) s;
  return p->river_y;
}

const double*
deltas_get_river_flux (Deltas_state * s)
{
  State *p = (State *) s;
  return p->river_flux;
}

int
deltas_get_n_rivers (Deltas_state * s)
{
  State *p = (State *) s;
  return p->n_rivers;
}

double
deltas_get_angle_asymmetry (Deltas_state * s)
{
  State *p = (State *) s;

  return p->angle_asymmetry;
}

double
deltas_get_angle_highness (Deltas_state * s)
{
  State *p = (State *) s;

  return p->angle_highness;
}

double
deltas_get_wave_angle (Deltas_state * s)
{
  State *p = (State *) s;

  return p->WaveAngle;
}

double
deltas_get_start_time (Deltas_state * s)
{
  return 0.;
}

#include <float.h>

int
BMI_CEM_Get_end_time (Deltas_state * s, double * time)
{
  *time = DBL_MAX;
  return BMI_SUCCESS;
}

int
BMI_CEM_Get_current_time (Deltas_state * s, double * time)
{
  State *p = (State *) s;
  *time = p->CurrentTimeStep * p->time_step;
  return BMI_SUCCESS;
}

int
BMI_CEM_Get_start_time (Deltas_state * s, double * time)
{
  *time = 0.;
  return BMI_SUCCESS;
}

int
BMI_CEM_Get_time_step (Deltas_state * s, double * dt)
{
  State *p = (State *) s;
  *dt = p->time_step;
  return BMI_SUCCESS;
}

int
BMI_CEM_Get_time_units (Deltas_state * s, char *units)
{
  if (s && units) {
    strcpy (units, "d");
    return BMI_SUCCESS;
  }
  else
    return BMI_FAILURE;
}


double
deltas_get_end_time (Deltas_state * s)
{
  return DBL_MAX;
}

double
deltas_get_current_time (Deltas_state * s)
{
  State *p = (State *) s;

  return p->CurrentTimeStep * p->time_step;
}

int
BMI_CEM_Get_var_type (Deltas_state *s, const char *value, BMI_Var_type *type)
{
  int rtn = BMI_FAILURE;

  if (type) {
    *type = BMI_VAR_TYPE_DOUBLE;
    rtn = BMI_SUCCESS;
  }

  return rtn;
}

int
BMI_CEM_Get_var_rank (Deltas_state * s, const char *name, int *rank)
{
  int rtn = BMI_FAILURE;
  if (rank) {
    fprintf (stderr, "getting rank for %s\n", name);
    if (strcmp (name, "sea_water__depth") == 0 ||
        strcmp (name, "surface__elevation") == 0 ||
        strcmp (name, "sea_water_to_sediment__depth_ratio") == 0 ||
        strcmp (name, "surface_bed_load_sediment__mass_flow_rate") == 0)
      *rank = 2;
    else if (strcmp (name, "river_mouth__location_model_x_component") == 0 ||
             strcmp (name, "river_mouth__location_model_y_component") == 0 ||
             strcmp (name, "river_mouth_bed_load_sediment__mass_flow_rate") == 0)
             //strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0)
      *rank = 1;
    else if (strcmp (name, "channel_outflow_end_suspended_load__mass_flow_rate") == 0 ||
             strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0 ||
             strcmp (name, "sea_water_surface_wave__from_direction") == 0 ||
             strcmp (name, "sea_water_surface_wave__height") == 0 ||
             strcmp (name, "sea_water_surface_wave__period") == 0)
      *rank = 0;
    else
      return BMI_FAILURE;
    fprintf (stderr, "rank is %d\n", *rank);
    rtn = BMI_SUCCESS;
  }
  return rtn;
}

int
BMI_CEM_Get_var_point_count (Deltas_state * s, const char *name, int *count)
{
  if (count && s) {
    State *p = (State *) s;
    fprintf (stderr, "getting point count for %s\n", name);
    fprintf (stderr, "nx, ny = %d, %d\n", p->nx, p->ny);
    if (strcmp (name, "surface_bed_load_sediment__mass_flow_rate") == 0 ||
        strcmp (name, "sea_water__depth") == 0 ||
        strcmp (name, "surface__elevation") == 0 ||
        strcmp (name, "sea_water_to_sediment__depth_ratio") == 0) {
      *count = p->nx * p->ny;
    }
    else if (strcmp (name, "river_mouth__location_model_x_component") == 0 ||
             strcmp (name, "river_mouth__location_model_y_component") == 0 ||
             strcmp (name, "river_mouth_bed_load_sediment__mass_flow_rate") == 0)
             //strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0)
      *count = p->n_rivers;
    else if (strcmp (name, "channel_outflow_end_suspended_load__mass_flow_rate") == 0 ||
             strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0 ||
             strcmp (name, "sea_water_surface_wave__from_direction") == 0 ||
             strcmp (name, "sea_water_surface_wave__height") == 0 ||
             strcmp (name, "sea_water_surface_wave__period") == 0)
      *count = 1;
    else
      return BMI_FAILURE + 1;
    fprintf (stderr, "point count is %d\n", *count);
    return BMI_SUCCESS;
  }
  return BMI_FAILURE;
}

int
BMI_CEM_Get_grid_type (Deltas_state * s, const char *name, BMI_Grid_type *type)
{
  int rtn = BMI_FAILURE;
  if (type) {
    State *p = (State *) s;
    *type = BMI_GRID_TYPE_UNIFORM;
    rtn = BMI_SUCCESS;
  }
  return rtn;
}

int
BMI_CEM_Get_grid_shape (Deltas_state * s, const char *name, int *shape)
{
  int rtn = BMI_FAILURE;
  if (shape && s) {
    State *p = (State *) s;
    fprintf (stderr, "getting shape for %s\n", name);
    if (strcmp (name, "surface_bed_load_sediment__mass_flow_rate") == 0 ||
        strcmp (name, "sea_water__depth") == 0 ||
        strcmp (name, "surface__elevation") == 0 ||
        strcmp (name, "sea_water_to_sediment__depth_ratio") == 0) {
      shape[0] = p->nx;
      shape[1] = p->ny;
    fprintf (stderr, "shape is %d, %d\n", shape[0], shape[1]);
    }
    else if (strcmp (name, "river_mouth__location_model_x_component") == 0 ||
             strcmp (name, "river_mouth__location_model_y_component") == 0 ||
             strcmp (name, "river_mouth_bed_load_sediment__mass_flow_rate") == 0)
             //strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0)
      shape[0] = p->n_rivers;
    else if (strcmp (name, "channel_outflow_end_suspended_load__mass_flow_rate") == 0 ||
             strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0 ||
             strcmp (name, "sea_water_surface_wave__from_direction") == 0 ||
             strcmp (name, "sea_water_surface_wave__height") == 0 ||
             strcmp (name, "sea_water_surface_wave__period") == 0)
      shape[0] = 1;
    else
      return BMI_FAILURE;
    rtn = BMI_SUCCESS;
  }
  return rtn;
}

int
BMI_CEM_Get_grid_spacing (Deltas_state * s, const char *name, double *spacing)
{
  int rtn = BMI_FAILURE;
  if (spacing) {
    if (strcmp (name, "surface_bed_load_sediment__mass_flow_rate") == 0 ||
        strcmp (name, "sea_water__depth") == 0 ||
        strcmp (name, "surface__elevation") == 0 ||
        strcmp (name, "sea_water_to_sediment__depth_ratio") == 0) {
      spacing[0] = deltas_get_dy (s);
      spacing[1] = deltas_get_dx (s);
    }
    else if (strcmp (name, "river_mouth__location_model_x_component") == 0 ||
             strcmp (name, "river_mouth__location_model_y_component") == 0 ||
             strcmp (name, "river_mouth_bed_load_sediment__mass_flow_rate") == 0)
             //strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0)
      spacing[0] = 0.;
    else if (strcmp (name, "channel_outflow_end_suspended_load__mass_flow_rate") == 0 ||
             strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0 ||
             strcmp (name, "sea_water_surface_wave__from_direction") == 0 ||
             strcmp (name, "sea_water_surface_wave__height") == 0 ||
             strcmp (name, "sea_water_surface_wave__period") == 0)
      spacing[0] = 0.;
    else
      return BMI_FAILURE;
    rtn = BMI_SUCCESS;
  }
  return rtn;
}

int
BMI_CEM_Get_grid_origin (Deltas_state * s, const char *name, double *origin)
{
  int rtn = BMI_FAILURE;
  if (origin) {
    if (strcmp (name, "surface_bed_load_sediment__mass_flow_rate") == 0 ||
        strcmp (name, "sea_water__depth") == 0 ||
        strcmp (name, "surface__elevation") == 0 ||
        strcmp (name, "sea_water_to_sediment__depth_ratio") == 0) {
      origin[0] = 0.;
      origin[1] = 0.;
    }
    else if (strcmp (name, "river_mouth__location_model_x_component") == 0 ||
             strcmp (name, "river_mouth__location_model_y_component") == 0 ||
             strcmp (name, "river_mouth_bed_load_sediment__mass_flow_rate") == 0)
             //strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0)
      origin[0] = 0.;
    else if (strcmp (name, "channel_outflow_end_suspended_load__mass_flow_rate") == 0 ||
             strcmp (name, "channel_outflow_end_bed_load_sediment__mass_flow_rate") == 0 ||
             strcmp (name, "sea_water_surface_wave__from_direction") == 0 ||
             strcmp (name, "sea_water_surface_wave__height") == 0 ||
             strcmp (name, "sea_water_surface_wave__period") == 0)
      origin[0] = 0.;
    else
      return BMI_FAILURE;
    rtn = BMI_SUCCESS;
  }
  return rtn;
}

int
deltas_get_nx (Deltas_state * s)
{
  State *p = (State *) s;

  return p->nx;
}

int
deltas_get_ny (Deltas_state * s)
{
  State *p = (State *) s;

  return p->ny * 2;
  //return 2*Ymax;
}

int
deltas_get_stride_old (Deltas_state * s, int dimen)
{
  if (dimen == 0)
    return 1;
  else if (dimen == 1)
    return deltas_get_ny (s) / 2;
  else
    return 0;
}

int
deltas_get_stride (Deltas_state * s, int dimen)
{
  if (dimen == 0)
    return 1;
  else if (dimen == 1)
    return deltas_get_ny (s);
  else
    return 0;
}

int
deltas_get_len (Deltas_state * s, int dimen)
{
  if (dimen == 0)
    return deltas_get_nx (s);
  else if (dimen == 1)
    return deltas_get_ny (s);
  else
    return 0;
}

double
deltas_get_dx (Deltas_state * s)
{
  State *p = (State *) s;

  return p->cell_width;
  //return CellWidth;
}

double
deltas_get_dy (Deltas_state * s)
{
  State *p = (State *) s;

  return p->cell_width;
  //return CellWidth;
}

void
deltas_use_external_waves (Deltas_state * s)
{
  State *p = (State *) s;

  p->external_waves = TRUE;
  //fprintf (stderr, "*** Ignoring request for external waves\n");
  //p->external_waves = FALSE;
}

void
deltas_use_sed_flux (Deltas_state * s)
{
  State *p = (State *) s;

  p->use_sed_flux = TRUE;
  //fprintf (stderr, "*** Ignoring request for sediment flux\n");
  //p->use_sed_flux = FALSE;
}
