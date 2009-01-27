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

char*  deltas_get_save_file( Deltas_state* );
char*  deltas_get_read_file( Deltas_state* );
float* deltas_get_depth    ( Deltas_state* );
float  deltas_get_sed_rate ( Deltas_state* );

#endif

