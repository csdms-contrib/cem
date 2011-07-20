#include <stdio.h>
#include <stdlib.h>
#include <deltas_api.h>
#include <deltas_cli.h>

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
