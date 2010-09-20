#include <stdio.h>
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
   {
      s = deltas_new();
   }

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

/** Set sediment flux in kg/s
*/
Deltas_state*
deltas_set_sed_flux (Deltas_state* s, float flux)
{
   State* p = (State*)s;
   p->SedFlux = flux;
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
deltas_set_wave_angle (Deltas_state* s, float wave_angle)
{
  State* p = (State*)s;
  p->WaveAngle = wave_angle;
  return s;
}

Deltas_state*
deltas_set_depth( Deltas_state* s, float* depth )
{
   State* p = (State*)s;
   memcpy( p->CellDepth[0], depth, sizeof(p->CellDepth) );
   return s;
}

Deltas_state*
deltas_set_wave_height (Deltas_state* s, float height_in_m)
{
  State* p = (State*)s;
  p->wave_height = height_in_m;
  return s;
}

Deltas_state*
deltas_set_wave_period (Deltas_state* s, float period_in_s)
{
  State* p = (State*)s;
  p->wave_period = period_in_s;
  return s;
}

Deltas_state*
deltas_set_shoreface_slope (Deltas_state* s, float shoreface_slope)
{
  State* p = (State*)s;
  p->shoreface_slope = shoreface_slope;
  return s;
}

Deltas_state*
deltas_set_shelf_slope (Deltas_state* s, float shelf_slope)
{
  State* p = (State*)s;
  p->shelf_slope = shelf_slope;
  return s;
}

Deltas_state*
deltas_set_shoreface_depth (Deltas_state* s, float shoreface_depth)
{
  State* p = (State*)s;
  p->shoreface_depth = shoreface_depth;
  return s;
}

const char* _deltas_exchange_items[] =
{
  "DEPTH",
  "PERCENT_FILLED",
  NULL
};

const const char**
deltas_get_exchange_items (void)
{
  return _deltas_exchange_items;
}

double*
deltas_get_value_grid (Deltas_state* s, const char* value)
{
  if (strcasecmp (value, "DEPTH")==0)
    return deltas_get_depth_dup (s);
  else if (strcasecmp (value, "PERCENT_FILLED")==0)
    return deltas_get_percent_dup (s);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return NULL;
}

double*
deltas_get_value_data (Deltas_state* s, const char* value, int lower[2],
                       int upper[2], int stride[2])
{
  double* data = NULL;
  lower[0] = deltas_get_ny (s)/4;
  lower[1] = 0;
  upper[0] = 3*deltas_get_ny (s)/4-1;
  upper[1] = deltas_get_nx (s)-1;
  stride[0] = 1;
  stride[1] = deltas_get_ny (s);
/*
  lower[0] = deltas_get_ny (s)/4;
  lower[1] = 0;
  upper[0] = 3*deltas_get_ny (s)/4-1;
  upper[1] = deltas_get_nx (s)-1;
  stride[0] = 1;
  stride[1] = deltas_get_ny (s);
*/

  fprintf (stderr, "nx: %d\n", deltas_get_nx (s));
  fprintf (stderr, "ny: %d\n", deltas_get_ny (s));
  fflush (stderr);

  if (strcasecmp (value, "DEPTH")==0)
    data = deltas_get_depth_dup (s);
  else if (strcasecmp (value, "PERCENT_FILLED")==0)
    data = deltas_get_percent_dup (s);
  else
    fprintf (stderr, "ERROR: %s: Bad value string.", value);

  return data;
}

int*
deltas_get_value_dimen (Deltas_state* s, const char* value, int shape[3])
{
  shape[0] = deltas_get_ny (s)/2;
  //shape[0] = deltas_get_ny (s);
  shape[1] = deltas_get_nx (s);
  shape[2] = 1;

  return shape;
}

double*
deltas_get_value_res (Deltas_state* s, const char* value, double res[3])
{
  res[0] = deltas_get_dy (s);
  res[1] = deltas_get_dx (s);
  res[2] = 1;

  return res;
}

const float*
deltas_get_depth( Deltas_state* s )
{
   State* p = (State*)s;
   return p->CellDepth[0];
}

double*
deltas_get_depth_dup (Deltas_state* s)
{
  State* p = (State*)s;
  double* val = NULL;

  {
    const int len = deltas_get_nx (s)*deltas_get_ny (s);
    const float* f_val = deltas_get_depth (s);

    val = (double*) malloc (sizeof(double)*len);
    if (val)
    {
      int i;
      for (i=0; i<len; i++)
        val[i] = f_val[i];
    }
  }

  return val;
}

const float*
deltas_get_percent (Deltas_state* s)
{
   State* p = (State*)s;
   return p->PercentFull[0];
}

double*
deltas_get_percent_dup (Deltas_state* s)
{
  State* p = (State*)s;
  double* val = NULL;

  {
    const int len = deltas_get_nx (s)*deltas_get_ny (s);
    const float* f_val = deltas_get_percent (s);

    val = (double*) malloc (sizeof(double)*len);

    if (val)
    {
      int i;
      for (i=0; i<len; i++)
        val[i] = f_val[i];
    }
  }

  return val;
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

double
deltas_get_wave_angle (Deltas_state* s)
{
  State* p = (State*)s;
  return p->WaveAngle;
}

double
deltas_get_start_time (Deltas_state* s)
{
  return 0.;
}

#include <float.h>
double
deltas_get_end_time (Deltas_state* s)
{
  return DBL_MAX;
}

double
deltas_get_current_time (Deltas_state* s)
{
  State* p = (State*)s;
  return p->CurrentTimeStep*TimeStep;
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

int
deltas_get_stride (Deltas_state* s, int dimen)
{
  if (dimen==0)
    return 1;
  else if (dimen==1)
    return deltas_get_nx (s);
  else
    return 0;
}

int
deltas_get_len (Deltas_state* s, int dimen)
{
  if (dimen==0)
    return deltas_get_nx (s);
  else if (dimen==1)
    return deltas_get_ny (s);
  else
    return 0;
}

double
deltas_get_dx (Deltas_state* s)
{
  return CellWidth;
}

double
deltas_get_dy (Deltas_state* s)
{
  return CellWidth;
}

void
deltas_use_external_waves (Deltas_state* s)
{
   State* p = (State*)s;
   p->external_waves = TRUE;
   //fprintf (stderr, "*** Ignoring request for external waves\n");
   //p->external_waves = FALSE;
}

void
deltas_use_sed_flux (Deltas_state* s)
{
  State* p = (State*)s;
  p->use_sed_flux = TRUE;
  //fprintf (stderr, "*** Ignoring request for sediment flux\n");
  //p->use_sed_flux = FALSE;
}

