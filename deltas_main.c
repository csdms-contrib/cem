#include <stdio.h>
#include <stdlib.h>
#include <deltas_api.h>

int
main( void )
{
  Deltas_state* s = deltas_init( NULL );

  deltas_set_save_file( s, "fileout_0" );

  {
    int i;
    const double dt = 5.2;

    for ( i=0; i<=100; i++ )
      deltas_run_until( s, i*dt );
  }

  deltas_finalize( s );

  return EXIT_SUCCESS;
}

