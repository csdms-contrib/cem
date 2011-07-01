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
Deltas_state* deltas_init_grid_shape (Deltas_state* s, int dimen[2]);
Deltas_state* deltas_init_cell_width (Deltas_state* s, double dx);
Deltas_state* deltas_init_grid (Deltas_state* s, double* z);
Deltas_state* deltas_set_grid (Deltas_state*, float*, int[2]);
Deltas_state* deltas_set_depth    ( Deltas_state*, float*);
Deltas_state* deltas_set_sed_rate ( Deltas_state*, float );
Deltas_state* deltas_find_river_mouth (Deltas_state* s, int n);
void deltas_avulsion (Deltas_state* s, double* qs, double river_flux);
Deltas_state* deltas_set_sed_flux ( Deltas_state*, float );
Deltas_state* deltas_set_river_sed_flux (Deltas_state* s, float flux, int n);
Deltas_state* deltas_set_river_position (Deltas_state* s, int x, int y, int n);
Deltas_state* deltas_set_angle_asymmetry (Deltas_state* s,
                                          float angle_asymmetry);
Deltas_state* deltas_set_angle_highness (Deltas_state* s, float angle_highness);
Deltas_state* deltas_set_wave_angle (Deltas_state* s, float wave_angle);
Deltas_state* deltas_set_wave_height (Deltas_state* s, float wave_height);
Deltas_state* deltas_set_wave_period (Deltas_state* s, float wave_period);
Deltas_state* deltas_set_shoreface_slope (Deltas_state* s, float shoreface_slope);
Deltas_state* deltas_set_shelf_slope (Deltas_state* s, float shelf_slope);
Deltas_state* deltas_set_shoreface_depth (Deltas_state* s, float shoreface_depth);

const char** deltas_get_exchange_items (void);
double* deltas_get_value_grid (Deltas_state* s, const char* value);
double* deltas_get_value_data (Deltas_state* s, const char* value, int lower[2],
                              int upper[2], int stride[2]);
int* deltas_get_value_dimen (Deltas_state* s, const char* value, int shape[3]);
double* deltas_get_value_res (Deltas_state* s, const char* value,
                              double res[3]);

//char*  deltas_get_save_file( Deltas_state* );
//char*  deltas_get_read_file( Deltas_state* );
const float* deltas_get_depth (Deltas_state*);
const float* deltas_get_percent (Deltas_state*);
double* deltas_get_depth_dup (Deltas_state*);
double* deltas_get_elevation_dup (Deltas_state* s);
double* deltas_get_percent_dup (Deltas_state*);
float  deltas_get_sed_rate ( Deltas_state* );
double deltas_get_angle_asymmetry (Deltas_state*);
double deltas_get_angle_highness (Deltas_state*);
double deltas_get_wave_angle (Deltas_state* s);

double deltas_get_end_time (Deltas_state* s);
double deltas_get_current_time (Deltas_state* s);
double deltas_get_start_time (Deltas_state* s);

int deltas_get_nx (Deltas_state*);
int deltas_get_ny (Deltas_state*);
int deltas_get_stride (Deltas_state* s, int dimen);
int deltas_get_len (Deltas_state* s, int dimen);
double deltas_get_dx (Deltas_state*);
double deltas_get_dy (Deltas_state*);

void deltas_use_external_waves (Deltas_state* s);
void deltas_use_sed_flux (Deltas_state* s);
#endif

