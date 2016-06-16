/*
 * =====================================================================
 *
 *       Filename:  waves_main.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  04/16/2010 04:04:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eric Hutton (eh), huttone@colorado.edu
 *        Company:  Community Surface Dynamics Modeling System
 *
 * =====================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>

#include "waves_cli.h"
#include "bmi_waves.h"

#define WAVES_PROGRAM_STR   "waves"
#define WAVES_MAJOR_VERSION (0)
#define WAVES_MINOR_VERSION (1)
#define WAVES_MICRO_VERSION (1)

double waves_next_angle (GRand * rand, double asymmetry, double highness);

int
main (int argc, char *argv[])
{
  waves_args_st *args = NULL;
  GError *error = NULL;

  args = parse_command_line (argc, argv, &error);
  g_assert (args);

  if (args->version)
  {
    fprintf (stdout, "%s %d.%d.%d\n", WAVES_PROGRAM_STR, WAVES_MAJOR_VERSION,
             WAVES_MINOR_VERSION, WAVES_MICRO_VERSION);
    return EXIT_SUCCESS;
  }

  {
    gint i;
    double a;
    GRand *rand = g_rand_new_with_seed (args->seed);
    const gint len = args->n_waves;
    const double asymmetry = args->asymmetry;
    const double highness = args->highness;
    const double scale = (args->radians) ? 1. : 180. / M_PI;

    BMI_Model *model = (BMI_Model*) malloc(sizeof(BMI_Model));

    if (args->verbose)
    {
      const gchar *units = (args->radians) ? "Radians" : "Degrees";

      fprintf (stderr, "Highness factor: %f\n", highness);
      fprintf (stderr, "Asymmetry factor: %f\n", asymmetry);
      fprintf (stderr, "Units: %s\n", units);
    }

    register_bmi_waves(model);

    fprintf (stdout, "Initializing... ");
    if (model->initialize(NULL, &(model->self)) == BMI_FAILURE) {
      fprintf (stdout, "FAIL.\n");
      fprintf (stderr, "Unable to initialize\n");
      return EXIT_FAILURE;
    }
    else
      fprintf (stdout, "PASS.\n");

    if (model->set_value(model->self, "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_asymmetry_parameter", &(args->asymmetry)) == BMI_FAILURE)
      return EXIT_FAILURE;
    if (model->set_value(model->self, "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_highness_parameter", &(args->highness)) == BMI_FAILURE)
      return EXIT_FAILURE;
    
    for (i = 0; i < len; i++) {
      fprintf (stdout, "Updating... ");
      if (model->update(model->self) == BMI_SUCCESS)
        fprintf (stdout, "PASS\n");
      else {
        fprintf (stdout, "FAIL.\n");
        fprintf (stderr, "Unable to update\n");
      }

      if (model->get_value(model->self, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity", &a) == BMI_SUCCESS)
        fprintf (stdout, "Wave angle: %f\n", a * scale);
      else
        fprintf (stderr, "Unable to get wave angle\n");
    }

    fprintf (stdout, "Finalizing... ");
    if (model->finalize(model->self) == BMI_SUCCESS)
      fprintf (stdout, "PASS\n");
    else {
      fprintf (stdout, "FAIL.\n");
      fprintf (stderr, "Unable to finalize\n");
    }
    g_rand_free (rand);
  }

  g_free (args);

  return EXIT_SUCCESS;
}
