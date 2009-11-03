#include <string.h>
#include <stdlib.h>

#include "deltas.h"
#include "deltas_api.h"

#ifdef SWIG
%include deltas_api.h
#endif

/** \file
\brief The Caperiffic API
*/

Deltas_state*
deltas_new( void )
{
   State* s = malloc( sizeof(State) );

   deltas_init_state( s );

   return (Deltas_state*)s;
}

Deltas_state*
deltas_destroy( Deltas_state* s )
{
   if ( s )
   {
      deltas_free_state( (State*)s );
      free( s );
   }
   return NULL;
}

Deltas_state*
deltas_init( Deltas_state* s )
{
   if ( !s )
      s = deltas_new();

   _cem_initialize( (State*)s );

   return s;
}

int
deltas_run_until( Deltas_state* s, float time_in_days )
{
   State* p = (State*)s;
   int until_time_step = time_in_days / TimeStep;
   return _cem_run_until( p, until_time_step );
}

Deltas_state*
deltas_finalize( Deltas_state* s, int free )
{
   _cem_finalize( (State*)s );

   if ( free )
      s = deltas_destroy( s );

   return s;
}

Deltas_state*
deltas_set_save_file( Deltas_state* s, char* name )
{
   State* p = (State*)s;
   p->savefilename = strdup( name );
   return s;
}

Deltas_state*
deltas_set_read_file( Deltas_state* s, char* name )
{
   State* p = (State*)s;
   p->readfilename = strdup( name );
   return s;
}

float
deltas_get_sed_rate( Deltas_state* s )
{
   State* p = (State*)s;
   return p->SedRate;
}

Deltas_state*
deltas_set_sed_rate( Deltas_state* s, float rate )
{
   State* p = (State*)s;
   p->SedRate = rate;
   return s;
}

Deltas_state*
deltas_set_angle_asymmetry (Deltas_state* s, float angle_asymmetry)
{
   State* p = (State*)s;
   p->angle_asymmetry = angle_asymmetry;
   return s;
}

Deltas_state*
deltas_set_angle_highness (Deltas_state* s, float angle_highness)
{
   State* p = (State*)s;
   p->angle_highness = angle_highness;
   return s;
}

Deltas_state*
deltas_set_depth( Deltas_state* s, float* depth )
{
   State* p = (State*)s;
   memcpy( p->CellDepth[0], depth, sizeof(p->CellDepth) );
   return s;
}

float*
deltas_get_depth( Deltas_state* s )
{
   State* p = (State*)s;
   return p->CellDepth[0];
}

float*
deltas_get_percent (Deltas_state* s)
{
   State* p = (State*)s;
   return p->PercentFull[0];
}

double
deltas_get_angle_asymmetry (Deltas_state* s)
{
  State* p = (State*)s;
  return p->angle_asymmetry;
}

double
deltas_get_angle_highness (Deltas_state* s)
{
  State* p = (State*)s;
  return p->angle_highness;
}

int
deltas_get_nx (Deltas_state* s)
{
  return Xmax;
}

int
deltas_get_ny (Deltas_state* s)
{
  return 2*Ymax;
}

