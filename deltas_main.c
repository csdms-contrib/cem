#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <deltas_cli.h>
#include "bmi_cem.h"

void print_matrix (double *x, int *shape);

int
main (int argc, char *argv[])
{
  BMI_Model *model = (BMI_Model*) malloc(sizeof(BMI_Model));
  int status;

  if (argc > 1) {
    if (strcmp (argv[1], "--version") == 0) {
      fprintf (stdout, "The Coastal Evolution Model version 0.1\n");
      exit (0);
    }
    else if (strcmp (argv[1], "--help") == 0) {
      fprintf (stdout, "Usage: run_deltas [--help] [--version] [FILE]\n");
      exit (0);
    }
  }

  register_bmi_cem(model);

  fprintf (stdout, "Initializing... ");
  if (argc>1)
    status = model->initialize(argv[1], &(model->self));
  else
    status = model->initialize(NULL, &(model->self));

  if (status == BMI_FAILURE) {
    fprintf (stdout, "FAIL.\n");
    fprintf (stderr, "Unable to initialize\n");
    return EXIT_FAILURE;
  }

  {
    int i;
    int grid;
    int rank;
    int *shape = NULL;
    int len;
    double *qs = NULL;
    double *z = NULL;
    double stop_time;
    const double river_flux = 250.;

    if (model->get_var_grid(model->self, "land_surface__elevation", &grid) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }
    if (model->get_grid_rank(model->self, grid, &rank) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }

    fprintf (stderr, "Grid rank: %d\n", rank);

    shape = (int*) malloc (sizeof (int)*rank);
    if (model->get_grid_shape(model->self, grid, shape) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }

    fprintf (stderr, "Grid shape: %d x %d\n", shape[0], shape[1]);

    if (model->get_grid_size(model->self, grid, &len) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }

    qs = (double *)malloc (sizeof (double) * len);
    z = (double *)malloc (sizeof (double) * len);
    fprintf (stderr, "len is %d\n", len);

    model->get_end_time(model->self, &stop_time);
    stop_time = 1000;

    for (i = 1; i <= stop_time; i++) {
      deltas_avulsion (model->self, qs, river_flux);

      if (model->set_value(model->self, "land_surface_water_sediment~bedload__mass_flow_rate", qs) == BMI_FAILURE) {
        fprintf(stderr, "unable to set qs\n");
        return EXIT_FAILURE;
      }

      if (model->update(model->self) == BMI_FAILURE) {
        fprintf(stderr, "unable to update\n");
        return EXIT_FAILURE;
      }

      if (model->get_value(model->self, "sea_water__depth", z) == BMI_FAILURE) {
        fprintf(stderr, "unable to get water depth\n");
        return EXIT_FAILURE;
      }

      if (i%100 == 0) {
        fprintf (stderr, "\n");
        fprintf (stderr, "Time: %d\n", i);
        fprintf (stderr, "Shape: %d x %d\n", shape[0], shape[1]);
        print_matrix (z, shape);
      }
    }

    free (qs);
    free (z);
    free (shape);

    {
      double time;
      int error;

      model->get_current_time(model->self, &time);
      if (error || fabs (time - stop_time) > 1e-6)
        return EXIT_FAILURE;
    }
  }

  model->finalize(model->self);

  return EXIT_SUCCESS;
}


void
print_matrix (double *x, int *shape)
{
  const int n_rows = shape[0];
  const int n_cols = shape[1];
  double *row = x;
  int i, j;

  for (i=0; i<n_rows; i++, row += n_cols) {
    for (j=0; j<n_cols; j++)
      fprintf (stderr, "%2.1f ", row[j]);
    fprintf (stderr, "\n");
  }

  return;
}
