#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "deltas.h"
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
  return model->CurrentTimeStep * TimeStep;
}


double
deltas_get_end_time (CemModel * model)
{
  return DBL_MAX;
}


double
deltas_get_time_step (CemModel * model)
{
  return TimeStep;
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

CemModel * deltas_set_save_file (CemModel * model, const char *name)
{
  model->savefilename = strdup (name);
  return model;
}


static void
_deltas_find_river_mouth (CemModel * model, int n)
{
  int x, y;

  x = 0;
  //y = model->stream_spot;
  y = deltas_get_ny(model) / 2;

  while (model->AllBeach[x][y] == 'y') {
    x += 1;
  }

  model->river_x_ind[n] = x;
  model->river_y_ind[n] = y;
  model->river_x[n] = x * deltas_get_dx(model);
  model->river_y[n] = y * deltas_get_dy(model);

//fprintf (stderr, "x = %d\n", x);
//fprintf (stderr, "y = %d\n", y);
}


void
deltas_avulsion (CemModel * model, double *qs, double river_flux)
{
  int x, y;
  int i;
  int qs_x, qs_y, qs_i;
  int len;

  _deltas_find_river_mouth(model, 0);

  x = model->river_x_ind[0];
  y = model->river_y_ind[0];

  len = deltas_get_nx (model) * deltas_get_ny (model) / 2;
  for (i = 0; i < len; i++)
    qs[i] = 0;

  qs_x = x;
  //qs_y = y % deltas_get_ny (s) - deltas_get_ny (s) / 4;
  qs_y = y;
  qs_i = qs_x * deltas_get_ny(model) / 2 + qs_y;

  qs[qs_i] = river_flux;

  return;
}

