#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "cem_model.h"


static char *
_scan_next_non_comment_line (char *line, int len, FILE *fp)
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


static CemModel *
_deltas_new (void)
{
  CemModel *model = (CemModel*) malloc (sizeof (CemModel));

  deltas_init_state (model);

  return model;
}


static CemModel*
_deltas_init_grid_shape (CemModel * model, int dimen[2])
{
  fprintf (stderr, "*** Grid size is (%d,%d)\n", model->nx, model->ny);
  fprintf (stderr, "*** Requested size is (%d,%d)\n", dimen[0], dimen[1]);

  if (model && (model->nx == 0 && model->ny == 0))
  {
    int i;
    const int len = dimen[0] * dimen[1];
    const int stride = dimen[1];

    model->nx = dimen[0];
    model->ny = dimen[1] / 2;
    model->max_beach_len = len;

    model->AllBeach = (char **)malloc (sizeof (char *) * model->nx);
    model->PercentFull = (double **)malloc (sizeof (double *) * model->nx);
    model->Age = (int **)malloc (sizeof (int *) * model->nx);
    model->CellDepth = (double **)malloc (sizeof (double *) * model->nx);
    model->InitDepth = (double **)malloc (sizeof (double *) * model->nx);

    model->AllBeach[0] = (char *)malloc (sizeof (char) * len);
    model->PercentFull[0] = (double *)malloc (sizeof (double) * len);
    model->Age[0] = (int *)malloc (sizeof (int) * len);
    model->CellDepth[0] = (double *)malloc (sizeof (double) * len);
    model->InitDepth[0] = (double *)malloc (sizeof (double) * len);

    for (i = 1; i < model->nx; i++)
    {
      model->AllBeach[i] = model->AllBeach[i - 1] + stride;
      model->PercentFull[i] = model->PercentFull[i - 1] + stride;
      model->Age[i] = model->Age[i - 1] + stride;
      model->CellDepth[i] = model->CellDepth[i - 1] + stride;
      model->InitDepth[i] = model->InitDepth[i - 1] + stride;
    }

    model->river_flux = (double *)malloc (sizeof (double) * len);
    model->river_x_ind = (int *)malloc (sizeof (int) * len);
    model->river_y_ind = (int *)malloc (sizeof (int) * len);
    model->river_x = (double *)malloc (sizeof (double) * len);
    model->river_y = (double *)malloc (sizeof (double) * len);
    model->n_rivers = 1;

    model->X = (int *)malloc (sizeof (int) * len);
    model->Y = (int *)malloc (sizeof (int) * len);
    model->InShadow = (char *)malloc (sizeof (char) * len);
    model->ShorelineAngle = (double *)malloc (sizeof (double) * len);
    model->SurroundingAngle = (double *)malloc (sizeof (double) * len);
    model->UpWind = (char *)malloc (sizeof (char) * len);
    model->VolumeIn = (double *)malloc (sizeof (double) * len);
    model->VolumeOut = (double *)malloc (sizeof (double) * len);
  }
  fprintf (stderr, "*** New grid size is (%d,%d)\n",
           deltas_get_nx (model), deltas_get_ny (model));

  return model;
}


static CemModel *
_deltas_init_cell_width (CemModel * model, double dx)
{
  model->cell_width = dx;
  return model;
}


static CemModel *
_deltas_init (CemModel * model)
{
  if (!model)
  {
    model = _deltas_new ();
  }

  _cem_initialize (model);

  return model;
}


static void
_deltas_use_sed_flux (CemModel * model)
{
  model->use_sed_flux = 1;
}


int
cem_initialize (const char *config_file, CemModel **handle)
{
    CemModel *model = NULL;
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
        return 3;

      _scan_next_non_comment_line (buffer, 2048, fp);
      fprintf (stderr, "CEM: line: %s\n", buffer);
      sscanf (buffer, "%d, %d, %lf, %d\n", &n_rows, &n_cols, &dx,
          &sed_flux_flag); 

      _scan_next_non_comment_line (buffer, 2048, fp);
      sscanf (buffer, "%lf, %lf, %lf", &shoreface_slope, &shoreface_depth,
          &shelf_slope);

      fprintf (stderr, "CEM: number of rows, columns: %d, %d\n", n_rows,
          n_cols);
    }
    else {
      n_rows = 120;
      n_cols = 400;
      dx = 100.;
      sed_flux_flag = 1;

      shoreface_slope = 0.01;
      shoreface_depth = 10.0;
      shelf_slope = 0.001;
    }

    /* Double n_cols for the full grid with mirrors. */
    n_cols *= 2;

    shape[0] = n_rows;
    shape[1] = n_cols;

    model = _deltas_new ();
    if (!model)
      return 4;

    _deltas_init_grid_shape (model, shape);
    _deltas_init_cell_width (model, dx);

    _deltas_init (model);

    if (sed_flux_flag)
      _deltas_use_sed_flux (model);

    deltas_set_shoreface_slope (model, shoreface_slope);
    deltas_set_shoreface_depth (model, shoreface_depth);
    deltas_set_shelf_slope (model, shelf_slope);

    deltas_set_save_file (model, "test.txt"); 

    *handle = model;

    return 0;
}


int
cem_finalize (CemModel * model)
{
  if (model->savefilename)
    printf ("Run Complete.  Output file: %s\n", model->savefilename);

  if (model->savefilename)
    free(model->savefilename);

  if (model->readfilename)
    free(model->readfilename);

  free(model->river_flux);
  free(model->river_x);
  free(model->river_y);
  free(model->river_x_ind);
  free(model->river_y_ind);
  free (model->state);

  return 0;
}


int
deltas_get_n_rivers (CemModel * model)
{
  return model->n_rivers;
}


int
deltas_get_nx (CemModel * model)
{
  return model->nx;
}


int
deltas_get_ny (CemModel * model)
{
  return model->ny * 2;
  /* return 2*Ymax; */
}


double
deltas_get_dx (CemModel * model)
{
  return model->cell_width;
}


double
deltas_get_dy (CemModel * model)
{
  return model->cell_width;
}


double
deltas_get_current_time (CemModel * model)
{
  return model->CurrentTimeStep * model->time_step;
}


double
deltas_get_end_time (CemModel * model)
{
  return DBL_MAX;
}


double
deltas_get_time_step (CemModel * model)
{
  return model->time_step;
}


double*
deltas_get_river_x_position (CemModel * model)
{
  return model->river_x;
}


double*
deltas_get_river_y_position (CemModel * model)
{
  return model->river_y;
}


double*
deltas_get_river_flux (CemModel * model)
{
  return model->river_flux;
}


double *
deltas_get_depth (CemModel * model)
{
  return model->CellDepth[0];
}


double *
deltas_get_elevation_dup(CemModel * model, double *z)
{
  int i;
  int id;
  const int z_nrows = deltas_get_nx(model);
  const int z_ncols = deltas_get_ny(model) / 2;
  const int len = z_nrows * z_ncols;
  int row, col;
  double percent_full;

  for (i = 0; i < len; i++)
  {
    row = i / z_ncols;
    col = i % z_ncols;
    id = row * (2 * z_ncols) + col + z_ncols / 2;

    percent_full = model->PercentFull[0][id];
    if (percent_full < 1e-6)
      z[i] = - model->CellDepth[0][id];
    else if (percent_full > 1. - 1e-6)
      z[i] = - model->CellDepth[0][id];
    else
      z[i] = percent_full;
  }

  return z;
}


static double *
_dup_subgrid (CemModel * model, double **src, double *dest)
{
  {
    int lower[2] = { deltas_get_ny(model) / 4, 0 };
    int upper[2] = { 3 * deltas_get_ny(model) / 4 - 1, deltas_get_nx(model) - 1 };
    int stride[2] = { 1, deltas_get_ny(model) };
    const int len = (upper[0] - lower[0] + 1) * (upper[1] - lower[1] + 1);

    if (dest == NULL)
        dest = (double *)malloc(sizeof(double) * len);

    if (dest) {
      int id, i, j;

      for (i=lower[1], id=0; i <= upper[1]; i++)
        for (j=lower[0]; j <= upper[0]; j++, id++)
          dest[id] = src[i][j];
    }
  }

  return dest;
}


double *
deltas_get_depth_dup (CemModel * model, double *dest)
{
    return _dup_subgrid (model, model->CellDepth, dest);
}


CemModel *
deltas_set_shoreface_slope (CemModel * model, double shoreface_slope)
{
  model->shoreface_slope = shoreface_slope;
  return model;
}

CemModel *
deltas_set_shelf_slope (CemModel * model, double shelf_slope)
{
  model->shelf_slope = shelf_slope;
  return model;
}

CemModel *
deltas_set_shoreface_depth (CemModel * model, double shoreface_depth)
{
  model->shoreface_depth = shoreface_depth;
  return model;
}

CemModel *
deltas_set_save_file (CemModel * model, const char *name)
{
  model->savefilename = strdup (name);
  return model;
}


static void
_deltas_prograde(CemModel * model, int *row, int *col)
{
  const int max_row = deltas_get_nx(model) - 1;
  const int max_col = deltas_get_ny(model) - 1;

  if (*col < 0)
    *col = 0;
  if (*col > max_col)
    *col = max_col;
  if (*row < 0)
    *row = 0;
  if (*row >= max_row)
    *row = max_row;

  while (*row < max_row && model->AllBeach[*row][*col] == 'y') {
    *row += 1;
  }
}


static void
_deltas_set_river_mouth(CemModel * model, int river_id, int row, int col)
{
    model->river_x_ind[river_id] = row;
    model->river_y_ind[river_id] = col;
    model->river_x[river_id] = row * deltas_get_dx(model);
    model->river_y[river_id] = col * deltas_get_dy(model);
}


CemModel *
deltas_set_sediment_flux_grid (CemModel * model, double *qs)
{
  int i;
  int n_rivers;
  const int len = deltas_get_nx(model) * deltas_get_ny(model) / 2;
  const int qs_ncols = deltas_get_ny(model) / 2;

  for (i=0, n_rivers=0; i < len; i++) {
      if (qs[i] > 0) {
          int row, col;

          model->river_flux[n_rivers] = qs[i];

          row = i / qs_ncols;
          col = i % qs_ncols + deltas_get_ny(model) / 4;

          _deltas_prograde(model, &row, &col);
          _deltas_set_river_mouth(model, n_rivers, row, col);

          n_rivers ++;
      }
  }
  model->n_rivers = n_rivers;

  return model;
}


CemModel *
deltas_set_elevation_grid(CemModel * model, double * elevation)
{
  int i;
  int id;
  const int z_nrows = deltas_get_nx(model);
  const int z_ncols = deltas_get_ny(model) / 2;
  const int len = z_nrows * z_ncols;
  int row, col;
  double percent_full, cell_depth;
  char all_beach;

  for (i = 0; i < len; i++)
  {
    row = i / z_ncols;
    col = i % z_ncols + z_ncols / 2;
    id = row * z_ncols * 2 + col;

    if (elevation[i] > 0.) {
      if (elevation[i] < 1.) {
        percent_full = elevation[i];
        cell_depth = 0.;
      } else if (elevation[i] >= 1.) {
        percent_full = 1.;
        cell_depth = - elevation[i];
      }
      // cell_depth = 0.;
    } else {
      // percent_full = model->PercentFull[0][id];
      percent_full = 0.;
      cell_depth = - elevation[i];

      if (cell_depth < model->shoreface_depth)
        cell_depth = model->shoreface_depth;
    }

    all_beach = (percent_full < 1.) ? 'n' : 'y';

    model->CellDepth[0][id] = cell_depth;
    model->PercentFull[0][id] = percent_full;
    model->AllBeach[0][id] = all_beach;

    if (col < z_ncols)
      id += z_ncols;
    else
      id -= z_ncols;
    model->CellDepth[0][id] = cell_depth;
    model->PercentFull[0][id] = percent_full;
    model->AllBeach[0][id] = all_beach;
  }

  return model;
}


void
deltas_avulsion (CemModel * model, double *qs, double river_flux)
{
  const int len = deltas_get_nx(model) * deltas_get_ny(model) / 2;
  const int river_id = 0;
  int row = model->river_x_ind[river_id];
  int col = model->river_y_ind[river_id];

  _deltas_prograde(model, &row, &col);
  _deltas_set_river_mouth(model, river_id, row, col);

  memset(qs, 0, sizeof(double) * len);

  {
    const int qs_row = row;
    const int qs_col = col % deltas_get_ny(model) - deltas_get_ny(model) / 4;
    const int qs_i = qs_row * deltas_get_ny(model) / 2 + qs_col;

    qs[qs_i] = river_flux;
  }
}
