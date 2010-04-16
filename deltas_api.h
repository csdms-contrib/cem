#if !defined( DELTAS_API_H )
#define DELTAS_API_H 

#if !defined(TRUE)
#define TRUE (1)
#endif

#if !defined(FALSE)
#define FALSE (0)
#endif

/** An opaque data structure that holds the state of an instance of a 
deltas model.
*/
typedef struct
{
}
Deltas_state;

Deltas_state* deltas_new     ( void );
Deltas_state* deltas_destroy ( Deltas_state* );

Deltas_state* deltas_init     ( Deltas_state* );
int           deltas_run_until( Deltas_state*, float );
Deltas_state* deltas_finalize ( Deltas_state*, int );

Deltas_state* deltas_set_save_file( Deltas_state*, char* );
Deltas_state* deltas_set_read_file( Deltas_state*, char* );
Deltas_state* deltas_set_depth    ( Deltas_state*, float*);
Deltas_state* deltas_set_sed_rate ( Deltas_state*, float );
Deltas_state* deltas_set_angle_asymmetry (Deltas_state* s,
                                          float angle_asymmetry);
Deltas_state* deltas_set_angle_highness (Deltas_state* s, float angle_highness);
Deltas_state* deltas_set_wave_angle (Deltas_state* s, float wave_angle);

//char*  deltas_get_save_file( Deltas_state* );
//char*  deltas_get_read_file( Deltas_state* );
float* deltas_get_depth    ( Deltas_state* );
float* deltas_get_percent  ( Deltas_state* );
float  deltas_get_sed_rate ( Deltas_state* );
double deltas_get_angle_asymmetry (Deltas_state*);
double deltas_get_angle_highness (Deltas_state*);
double deltas_get_wave_angle (Deltas_state* s);
int deltas_get_nx (Deltas_state*);
int deltas_get_ny (Deltas_state*);

#endif

