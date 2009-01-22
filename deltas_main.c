#include <stdlib.h>
#include <deltas_api.h>

int
main( void )
{
   Deltas_state* s_0 = deltas_new();
   Deltas_state* s_1 = deltas_new();

   deltas_init( s_0 );
   deltas_init( s_1 );

   deltas_set_save_file( s_0, "fileout_0" );
   deltas_set_save_file( s_1, "fileout_1" );

   {
      int i;
      for ( i=0; i<=100; i++ )
      {
         deltas_run_until( s_0, i*26.F );
         deltas_run_until( s_1, i*52.F );
      }
   }

   deltas_finalize( s_0 );
   deltas_finalize( s_1 );

   deltas_destroy( s_0 );
   deltas_destroy( s_1 );

   return EXIT_SUCCESS;
}

