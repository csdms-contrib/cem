#include <stdio.h>
#include <stdlib.h>
#include <deltas_api.h>
#include <deltas_cli.h>

int
main( int argc, char *argv[] )
{
  cem_args_st* args;

  args = parse_command_line( argc, argv );

  {
    Deltas_state* s = deltas_init( NULL );

    deltas_set_save_file( s, "fileout_0" );

    {
      int i;
      const double dt = 5.2;

      for ( i=0; i<=100; i++ )
        deltas_run_until( s, i*dt );
    }

    deltas_finalize( s, TRUE );
  }

  free( args );

  return EXIT_SUCCESS;
}

