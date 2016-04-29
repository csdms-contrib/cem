#if !defined( WAVES_API_H )
#define WAVES_API_H

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(TRUE)
#define TRUE (1)
#endif

#if !defined(FALSE)
#define FALSE (0)
#endif

//typedef struct _BMI_Model BMI_Model;

/** An opaque data structure that holds the state of an instance of a
waves model.
*/
typedef struct
{
}
Waves_state;

typedef Waves_state BMI_Model;


#define BMI_SUCCESS (0)
#define BMI_FAILURE (1)

#define BMI_WAVES_COMPONENT_NAME_MAX (2048)
#define BMI_WAVES_VAR_NAME_MAX (2048)
#define BMI_WAVES_UNIT_NAME_MAX (2048)

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

int BMI_WAVES_Initialize (const char*, BMI_Model **);
int BMI_WAVES_Update (BMI_Model*);
int BMI_WAVES_Update_until (BMI_Model*, double);
int BMI_WAVES_Finalize (BMI_Model*);

int BMI_WAVES_Get_component_name (BMI_Model * self, char *name);
int BMI_WAVES_Get_input_var_name_count (BMI_Model * self, int *count);
int BMI_WAVES_Get_input_var_names (const BMI_Model * self, char ** names);
int BMI_WAVES_Get_output_var_name_count (BMI_Model * self, int *count);
int BMI_WAVES_Get_output_var_names (const BMI_Model * self, char ** names);

int BMI_WAVES_Get_var_rank (const BMI_Model *self, const char * name, int *rank);
int BMI_WAVES_Get_var_stride (BMI_Model * self, const char *value, int *stride);
int BMI_WAVES_Get_var_type (const BMI_Model *self, const char * name, BMI_Var_type *type);
int BMI_WAVES_Get_var_point_count (const BMI_Model *self, const char * name, int *count);
int BMI_WAVES_Get_grid_shape (const BMI_Model *self, const char * name, int *shape);
int BMI_WAVES_Get_grid_spacing (const BMI_Model *self, const char * name, double *spacing);
int BMI_WAVES_Get_grid_origin (const BMI_Model *self, const char * name, double *origin);

int BMI_WAVES_Get_double (BMI_Model * self, const char *value, double *dest);
int BMI_WAVES_Get_double_ptr (BMI_Model * self, const char *value, double **dest);
int BMI_WAVES_Set_double (BMI_Model *self, const char *value, double *src);

int BMI_WAVES_Get_grid_type (const BMI_Model *self, const char * name, BMI_Grid_type *type);
int BMI_WAVES_Get_current_time (const BMI_Model * self, double *time);
int BMI_WAVES_Get_start_time (const BMI_Model * self, double *time);
int BMI_WAVES_Get_end_time (const BMI_Model * self, double *time);
int BMI_WAVES_Get_time_units (const BMI_Model * self, char *units);

#define NO_BMI_WAVES_GET_GRID_CONNECTIVITY
#define NO_BMI_WAVES_GET_GRID_OFFSET
#define NO_BMI_WAVES_GET_GRID_X
#define NO_BMI_WAVES_GET_GRID_Y
#define NO_BMI_WAVES_GET_GRID_Z
//#define NO_BMI_WAVES_GET_GRID_SPACING
//#define NO_BMI_WAVES_GET_GRID_ORIGIN

Waves_state *waves_new (void);

Waves_state *waves_destroy (Waves_state *);

Waves_state *waves_init (Waves_state *);

int waves_run_until (Waves_state *, double);

Waves_state *waves_finalize (Waves_state * self, int free);

Waves_state *waves_set_angle_asymmetry (Waves_state *, double asymmetry);

Waves_state *waves_set_angle_highness (Waves_state *, double highness);

Waves_state *waves_set_height (Waves_state *, double height_in_m);

Waves_state *waves_set_period (Waves_state *, double height_in_s);

double waves_get_angle_asymmetry (const Waves_state *);

double waves_get_angle_highness (const Waves_state *);

double waves_get_wave_angle (Waves_state *);

double waves_get_wave_angle_mean (Waves_state *);

double waves_get_wave_angle_max (Waves_state *);

double waves_get_wave_angle_min (Waves_state *);

double waves_get_height (const Waves_state *);

double waves_get_period (const Waves_state *);

const char **waves_get_exchange_items (void);

double waves_get_value (Waves_state * self, const char *value);

double waves_get_current_time (const Waves_state * self);

double waves_get_start_time (const Waves_state * self);

double waves_get_end_time (const Waves_state * self);

#ifdef __cplusplus
}
#endif

#endif
