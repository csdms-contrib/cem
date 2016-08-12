#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <deltas_cli.h>
#include "cem_model.h"
#include "waves_model.h"
#include "bmi_cem.h"
#include "bmi_waves.h"

void print_model_info(BMI_Model *model);
void print_matrix (double *x, int *shape, FILE *fp);

int
main (int argc, char *argv[])
{
  BMI_Model *cem = (BMI_Model*) malloc(sizeof(BMI_Model));
  BMI_Model *waves = (BMI_Model*) malloc(sizeof(BMI_Model));
  FILE *output_file = NULL;
  int status;

  if (argc > 1) {
    if (strcmp (argv[1], "--version") == 0) {
      fprintf (stdout, "The Coastal Evolution Model version 0.1.2\n");
      exit (0);
    }
    else if (strcmp (argv[1], "--help") == 0) {
      fprintf (stdout, "Usage: run_cem [--help] [--version] WAVE_FILE DELTA_FILE [OUTPUT_FILE]\n");
      exit (0);
    }
  }

  if (argc != 3 && argc != 4) {
      fprintf (stdout, "Usage: run_cem [--help] [--version] WAVE_FILE DELTA_FILE [OUTPUT_FILE]\n");
    exit (1);
  }

  if (argc == 4)
    output_file = fopen(argv[3], "wb");
  else
    output_file = stdout;

  register_bmi_waves(waves);
  register_bmi_cem(cem);

  fprintf (stdout, "Initializing... ");

  status = waves->initialize(argv[1], &(waves->self));
  if (status == BMI_FAILURE) {
    fprintf (stdout, "FAIL.\n");
    fprintf (stderr, "Unable to initialize\n");
    return EXIT_FAILURE;
  }

  status = cem->initialize(argv[2], &(cem->self));
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
    double stop_time = 2500;
    double time;
    const double river_flux = 1750.;
    double angle, wave_height, wave_period;
    double time_step;


    cem->get_value(cem->self, "model__time_step", &time_step);
    fprintf(stderr, "Time step: %f\n", time_step);

    if (cem->get_var_grid(cem->self, "land_surface__elevation", &grid) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }
    if (cem->get_grid_rank(cem->self, grid, &rank) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }
    else
      fprintf (stderr, "Grid rank: %d\n", rank);

    shape = (int*) malloc (sizeof (int)*rank);
    if (cem->get_grid_shape(cem->self, grid, shape) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }
    else
      fprintf (stderr, "Grid shape: %d x %d\n", shape[0], shape[1]);

    if (cem->get_grid_size(cem->self, grid, &len) == BMI_FAILURE) {
      fprintf(stderr, "unable to get var grid\n");
      return EXIT_FAILURE;
    }

    qs = (double *)malloc (sizeof (double) * len);
    z = (double *)malloc (sizeof (double) * len);

    waves->get_current_time(waves->self, &time);
    waves->get_end_time(waves->self, &stop_time);
    i = 0;

    while (time < stop_time) {
    //for (i = 1; i <= stop_time; i++) {
      deltas_avulsion (cem->self, qs, river_flux);

      if (cem->set_value(cem->self, "land_surface_water_sediment~bedload__mass_flow_rate", qs) == BMI_FAILURE) {
        fprintf(stderr, "unable to set qs\n");
        return EXIT_FAILURE;
      }

      if (waves->update(waves->self) == BMI_FAILURE) {
        fprintf(stderr, "unable to update waves\n");
        return EXIT_FAILURE;
      }
      waves->get_current_time(waves->self, &time);

      status += waves->get_value(waves->self, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity", &angle);
      status += waves->get_value(waves->self, "sea_surface_water_wave__height", &wave_height);
      status += waves->get_value(waves->self, "sea_surface_water_wave__period", &wave_period);

      status += cem->set_value(cem->self, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity", &angle);
      status += cem->set_value(cem->self, "sea_surface_water_wave__height", &wave_height);
      status += cem->set_value(cem->self, "sea_surface_water_wave__period", &wave_period);

      if (cem->update_until(cem->self, time) == BMI_FAILURE) {
        fprintf(stderr, "unable to update cem\n");
        return EXIT_FAILURE;
      }

      if (cem->get_value(cem->self, "sea_water__depth", z) == BMI_FAILURE) {
        fprintf(stderr, "unable to get water depth\n");
        return EXIT_FAILURE;
      }

      i += 1;

      if (i%1000 == 0) {
        print_model_info(cem);
        print_matrix (z, shape, output_file);
      }

    }

    free (qs);
    free (z);
    free (shape);

    {
      double time;

      cem->get_current_time(cem->self, &time);
      if (fabs (time - stop_time) > 1e-6)
        return EXIT_FAILURE;
    }
  }

  cem->finalize(cem->self);
  waves->finalize(waves->self);

  fclose(output_file);

  return EXIT_SUCCESS;
}


void
print_model_info(BMI_Model *model)
{
    double angle = 0.;
    double wave_height = 0., wave_period = 0.;
    double time = 0;
    int status = 0;

    fprintf (stderr, "\n");

    status += model->get_value(model->self, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity", &angle);
    status += model->get_value(model->self, "sea_surface_water_wave__height", &wave_height);
    status += model->get_value(model->self, "sea_surface_water_wave__period", &wave_period);
    status += model->get_current_time(model->self, &time);

    if (status != 0) {
        fprintf(stderr, "Error getting model info.\n");
        return;
    }

    fprintf (stderr, "Time: %f\n", time);
    fprintf (stderr, "Angle: %f\n", angle);
    fprintf (stderr, "Wave height: %f\n", wave_height);
    fprintf (stderr, "Wave period: %f\n", wave_period);
}


void
print_matrix (double *x, int *shape, FILE *fp)
{
  const int n_rows = shape[0];
  const int n_cols = shape[1];
  double *row = x;
  int i, j;

  for (i=0; i<n_rows; i++, row += n_cols) {
    for (j=0; j<n_cols - 1; j++)
      fprintf(fp, "%f, ", row[j]);
    fprintf(fp, "%f\n", row[n_cols - 1]);
  }

  return;
}
