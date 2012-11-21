#include <stdio.h>
#include <stdlib.h>
#include <deltas_api.h>
#include <deltas_cli.h>

void print_matrix (double *x, int *shape);

int
main (int argc, char *argv[])
{
  BMI_Model *self = NULL;
  int err;

  fprintf (stderr, "Initializing... ");
  if (argc>1)
    err = BMI_Initialize (argv[1], &self);
  else
    err = BMI_Initialize (NULL, &self);

  if (err) {
    fprintf (stderr, "Error: %d: Unable to initialize\n", err);
    return EXIT_FAILURE;
  }

  fprintf (stderr, "PASS\n");

  {
    int i;
    int rank;
    int *shape = NULL;
    int len;
    double *qs = NULL;
    double *z = NULL;
    double stop_time;
    const double river_flux = 250.;

    BMI_Get_var_rank (self, "SedimentFlux", &rank);
    fprintf (stderr, "Grid rank: %d\n", rank);

    shape = (int*) malloc (sizeof (int)*rank);
    BMI_Get_grid_shape (self, "SedimentFlux", shape);
    fprintf (stderr, "Grid shape: %d x %d\n", shape[0], shape[1]);

    BMI_Get_var_point_count (self, "SedimentFlux", &len);
    //for (i=0, len=1; i<rank; i++)
    //  len *= shape[i];

    qs = (double *)malloc (sizeof (double) * len);
    z = (double *)malloc (sizeof (double) * len);
    fprintf (stderr, "len is %d\n", len);

    //BMI_Get_double (self, "PERCENT_FILLED", z);
    //print_matrix (z, shape);
    //return 0;

    BMI_Get_end_time (self, &stop_time);
    stop_time = 1000;
    for (i = 1; i <= stop_time; i++) {
      deltas_avulsion (self, qs, river_flux);

      //fprintf (stderr, "Update until %d... ", i);
      err = BMI_Set_double (self, "SedimentFlux", qs);
      if (err)
        fprintf (stderr, "Error %d\n", err);

      if (BMI_Update (self) == BMI_SUCCESS) {
        //fprintf (stderr, "PASS\n");
      }

      BMI_Get_double (self, "PERCENT_FILLED", z);
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
  }

  fprintf (stderr, "Finalizing... ");
  BMI_Finalize (self);
  fprintf (stderr, "PASS");
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

#if 0
int
main (int argc, char *argv[])
{
  cem_args_st *args;

  args = parse_command_line (argc, argv);

  {
    //Deltas_state* s = deltas_init (NULL);
    Deltas_state *s = deltas_new ();
    int dimen[2] = { 200, 1000 };

    fprintf (stderr, "Set grid shape\n");
    deltas_init_grid_shape (s, dimen);

    fprintf (stderr, "Set cell width\n");
    deltas_init_cell_width (s, 100.);

    fprintf (stderr, "Initialize\n");
    deltas_init (s);

    fprintf (stderr, "Set save file\n");
    //deltas_set_save_file( s, "fileout_0" );
    deltas_set_save_file (s, args->out_prefix);

//    deltas_run_until( s, args->stop_time );
    {
      int i;

      const double dt = 5.2;

      const double river_flux = 250.;

      const int len = deltas_get_nx (s) * deltas_get_ny (s);

      double *qs = (double *)malloc (sizeof (double) * len);

      deltas_use_sed_flux (s);

      for (i = 0; i <= args->stop_time; i++)
      {
        fprintf (stderr, "Run (%d)\n", i);
        //deltas_set_river_sed_flux (s, river_flux, 0);
        //deltas_find_river_mouth (s, 0);

        deltas_avulsion (s, qs, river_flux);
        deltas_set_sediment_flux_grid (s, qs);

        deltas_run_until (s, i);
      }
    }

    deltas_finalize (s, TRUE);
  }

  free (args);

  return EXIT_SUCCESS;
}
#endif

