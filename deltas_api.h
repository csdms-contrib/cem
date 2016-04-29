#if !defined( DELTAS_API_H )
#define DELTAS_API_H

#ifdef __cplusplus
extern "C" {
#endif

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

typedef Deltas_state BMI_Model;

#define BMI_SUCCESS (0)
#define BMI_FAILURE (1)

#define BMI_CEM_COMPONENT_NAME_MAX (2048)
#define BMI_CEM_VAR_NAME_MAX (2048)
#define BMI_CEM_UNIT_NAME_MAX (2048)

typedef enum {
  BMI_VAR_TYPE_UNKNOWN = 0,
  BMI_VAR_TYPE_CHAR,
  BMI_VAR_TYPE_UNSIGNED_CHAR,
  BMI_VAR_TYPE_INT,
  BMI_VAR_TYPE_LONG,
  BMI_VAR_TYPE_UNSIGNED_INT,
  BMI_VAR_TYPE_UNSIGNED_LONG,
  BMI_VAR_TYPE_FLOAT,
  BMI_VAR_TYPE_DOUBLE,
  BMI_VAR_TYPE_COUNT
}
BMI_Var_type;

typedef enum {
  BMI_GRID_TYPE_UNKNOWN = 0,
  BMI_GRID_TYPE_UNIFORM,
  BMI_GRID_TYPE_RECTILINEAR,
  BMI_GRID_TYPE_STRUCTURED,
  BMI_GRID_TYPE_UNSTRUCTURED,
  BMI_GRID_TYPE_COUNT
}
BMI_Grid_type;

/* BMI Function definitions */
int BMI_CEM_Initialize (const char *config_file, BMI_Model **handle);
int BMI_CEM_Update (BMI_Model * s);
int BMI_CEM_Update_until (BMI_Model * s, double time_in_days);
int BMI_CEM_Finalize (BMI_Model * s);

int BMI_CEM_Get_component_name (BMI_Model *s, char *name);
int BMI_CEM_Get_output_var_names (BMI_Model *s, char **names);
int BMI_CEM_Get_output_var_name_count (Deltas_state *s, int *count);
int BMI_CEM_Get_input_var_names (BMI_Model *s, char **names);
int BMI_CEM_Get_input_var_name_count (Deltas_state *s, int *count);

int BMI_CEM_Get_double (BMI_Model *s, const char *value, double *dest);
int BMI_CEM_Get_double_ptr (BMI_Model *s, const char *value, double **dest);
int BMI_CEM_Set_double (BMI_Model *s, const char * value, double *src);

int BMI_CEM_Get_end_time (BMI_Model * s, double * time);
int BMI_CEM_Get_current_time (BMI_Model * s, double * time);
int BMI_CEM_Get_start_time (BMI_Model * s, double * time);
int BMI_CEM_Get_time_step (BMI_Model * s, double * dt);
int BMI_CEM_Get_time_units (Deltas_state * s, char *units);

int BMI_CEM_Get_var_type (BMI_Model *s, const char *value, BMI_Var_type *type);
int BMI_CEM_Get_var_rank (BMI_Model * s, const char *name, int *rank);
int BMI_CEM_Get_var_stride (BMI_Model *s, const char *value, int *stride);
int BMI_CEM_Get_var_point_count (BMI_Model * model, const char *name, int *count);

int BMI_CEM_Get_grid_type (BMI_Model * s, const char *name, BMI_Grid_type *type);
int BMI_CEM_Get_grid_shape (BMI_Model * s, const char *name, int *shape);
int BMI_CEM_Get_grid_spacing (BMI_Model * s, const char *name, double *spacing);
int BMI_CEM_Get_grid_origin (BMI_Model * s, char const *name, double *origin);

#define NO_BMI_CEM_GET_GRID_CONNECTIVITY
#define NO_BMI_CEM_GET_GRID_OFFSET
#define NO_BMI_CEM_GET_GRID_X
#define NO_BMI_CEM_GET_GRID_Y
#define NO_BMI_CEM_GET_GRID_Z

Deltas_state *deltas_new (void);

Deltas_state *deltas_destroy (Deltas_state *);

Deltas_state *deltas_init (Deltas_state *);

int deltas_run_until (Deltas_state *, double);

Deltas_state *deltas_finalize (Deltas_state *, int);

Deltas_state *deltas_set_save_file (Deltas_state *, const char *);

Deltas_state *deltas_set_read_file (Deltas_state *, char *);

Deltas_state *deltas_init_grid_shape (Deltas_state * s, int dimen[2]);

Deltas_state *deltas_init_cell_width (Deltas_state * s, double dx);

Deltas_state *deltas_init_grid (Deltas_state * s, double *z);

Deltas_state *deltas_set_grid (Deltas_state *, double *, int[2]);

Deltas_state *deltas_set_depth (Deltas_state *, double *);

Deltas_state *deltas_set_sed_rate (Deltas_state *, double);

Deltas_state *deltas_find_river_mouth (Deltas_state * s, int n);

void deltas_avulsion (Deltas_state * s, double *qs, double river_flux);

Deltas_state *deltas_set_sed_flux (Deltas_state *, double);

Deltas_state *deltas_set_river_sed_flux (Deltas_state * s, double flux, int n);

Deltas_state *deltas_set_river_position (Deltas_state * s, int x, int y, int n);

Deltas_state *deltas_set_sediment_flux_grid (Deltas_state * s, double *qs);

Deltas_state * deltas_set_rivers (Deltas_state * s, const double * x,
                                  const double * y, double * qb,
                                  const int len);

Deltas_state *deltas_set_angle_asymmetry (Deltas_state * s,
                                          double angle_asymmetry);
Deltas_state *deltas_set_angle_highness (Deltas_state * s,
                                         double angle_highness);
Deltas_state *deltas_set_wave_angle (Deltas_state * s, double wave_angle);

Deltas_state *deltas_set_wave_height (Deltas_state * s, double wave_height);

Deltas_state *deltas_set_wave_period (Deltas_state * s, double wave_period);

Deltas_state *deltas_set_shoreface_slope (Deltas_state * s,
                                          double shoreface_slope);
Deltas_state *deltas_set_shelf_slope (Deltas_state * s, double shelf_slope);

Deltas_state *deltas_set_shoreface_depth (Deltas_state * s,
                                          double shoreface_depth);

const char **deltas_get_exchange_items (void);

const double *deltas_get_value_grid (Deltas_state * s, const char *value);

double *deltas_get_value_grid_dup (Deltas_state * s, const char *value);

const double *deltas_get_value_data (Deltas_state * s, const char *value,
                               int lower[2], int upper[2], int stride[2]);
double *deltas_get_value_data_dup (Deltas_state * s, const char *value,
                                   int lower[2], int upper[2], int stride[2]);
int *deltas_get_value_dimen (Deltas_state * s, const char *value, int shape[3]);

double *deltas_get_value_res (Deltas_state * s, const char *value,
                              double res[3]);

//char*  deltas_get_save_file( Deltas_state* );
//char*  deltas_get_read_file( Deltas_state* );
const double *deltas_get_depth (Deltas_state *);

const double *deltas_get_percent (Deltas_state *);

double *deltas_get_depth_dup (Deltas_state *);

double *deltas_get_elevation_dup (Deltas_state * s);

double *deltas_get_percent_dup (Deltas_state *);

const double* deltas_get_river_x_position (Deltas_state * s);
const double* deltas_get_river_y_position (Deltas_state * s);
const double* deltas_get_river_flux (Deltas_state * s);
int deltas_get_n_rivers (Deltas_state * s);

double deltas_get_sed_rate (Deltas_state *);

double deltas_get_angle_asymmetry (Deltas_state *);

double deltas_get_angle_highness (Deltas_state *);

double deltas_get_wave_angle (Deltas_state * s);

double deltas_get_end_time (Deltas_state * s);

double deltas_get_current_time (Deltas_state * s);

double deltas_get_start_time (Deltas_state * s);

int deltas_get_nx (Deltas_state *);

int deltas_get_ny (Deltas_state *);

int deltas_get_stride (Deltas_state * s, int dimen);

int deltas_get_len (Deltas_state * s, int dimen);

double deltas_get_dx (Deltas_state *);

double deltas_get_dy (Deltas_state *);

void deltas_use_external_waves (Deltas_state * s);

void deltas_use_sed_flux (Deltas_state * s);

#ifdef __cplusplus
}
#endif

#endif
