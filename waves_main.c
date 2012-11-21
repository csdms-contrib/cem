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
#include "waves_api.h"

#define WAVES_PROGRAM_STR   "waves"
#define WAVES_MAJOR_VERSION (0)
#define WAVES_MINOR_VERSION (1)
#define WAVES_MICRO_VERSION (0)

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

    Waves_state *waves_state = NULL;
    BMI_Model *w = NULL;
    int err;

    if (args->verbose)
    {
      gchar *units = (args->radians) ? "Radians" : "Degrees";

      fprintf (stderr, "Highness factor: %f\n", highness);
      fprintf (stderr, "Asymmetry factor: %f\n", asymmetry);
      fprintf (stderr, "Units: %s\n", units);
    }

    fprintf (stdout, "Initializing... ");
    err = BMI_Initialize (NULL, &w);

    if (err) {
      fprintf (stdout, "FAIL.\n");
      fprintf (stderr, "Error: %d: Unable to initialize\n", err);
    }
    else
      fprintf (stdout, "PASS.\n");

    err = BMI_Set_double (w, "wave_asymmetry", &(args->asymmetry));
    err = BMI_Set_double (w, "wave_highness", &(args->highness));
    for (i = 0; i < len; i++) {
      fprintf (stdout, "Updating... ");
      if (BMI_Update (w) == BMI_SUCCESS)
        fprintf (stdout, "PASS\n");
      else {
        fprintf (stdout, "FAIL.\n");
        fprintf (stderr, "Error: %d: Unable to update\n", err);
      }

      if (BMI_Get_double (w, "sea_surface_wave_from_direction", &a) == BMI_SUCCESS)
        fprintf (stdout, "Wave angle: %f\n", a * scale);
      else
        fprintf (stderr, "Error: %d: Unable to get wave angle\n", err);
    }

    fprintf (stdout, "Finalizing... ");
    if (BMI_Finalize (w) == BMI_SUCCESS)
      fprintf (stdout, "PASS\n");
    else {
      fprintf (stdout, "FAIL.\n");
      fprintf (stderr, "Error: %d: Unable to finalize\n", err);
    }
/*
    waves_state = waves_init (NULL);
    waves_set_angle_asymmetry (waves_state, args->asymmetry);
    waves_set_angle_highness (waves_state, args->highness);

    for (i = 0; i < len; i++)
    {
      waves_run_until (waves_state, i);
      a = waves_get_wave_angle (waves_state) * scale;
      //a = waves_next_angle (rand, asymmetry, highness)*scale;
      fprintf (stdout, "%f\n", a);
    }

    waves_state = waves_finalize (waves_state, TRUE);
*/
    g_rand_free (rand);
  }

  g_free (args);

  return EXIT_SUCCESS;
}
