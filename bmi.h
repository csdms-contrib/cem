#ifndef BMI_CEM_INCLUDED
#define BMI_CEM_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi_const.h"

#define BMI_CEM_INPUT_VAR_NAME_COUNT 1
#define BMI_CEM_OUTPUT_VAR_NAME_COUNT 1

typedef struct _BMI_CEM_Model BMI_CEM_Model;

/* Model Control functions */
int BMI_CEM_Initialize (const char *, BMI_CEM_Model**);
int BMI_CEM_Update (BMI_CEM_Model *);
int BMI_CEM_Update_until (BMI_CEM_Model *, double);
int BMI_CEM_Finalize (BMI_CEM_Model *);
int BMI_CEM_Run_model (BMI_CEM_Model *);

/* Model information functions */
int BMI_CEM_Get_component_name (BMI_CEM_Model *, char *);
int BMI_CEM_Get_input_var_name_count (BMI_CEM_Model, int *);
int BMI_CEM_Get_output_var_name_count (BMI_CEM_Model, int *);
int BMI_CEM_Get_input_var_names (BMI_CEM_Model *, char **);
int BMI_CEM_Get_output_var_names (BMI_CEM_Model *, char **);

/* Variable information functions */
int BMI_CEM_Get_var_type (BMI_CEM_Model *, const char *, BMI_Var_type *);
int BMI_CEM_Get_var_units (BMI_CEM_Model *, const char *, char *);
int BMI_CEM_Get_var_rank (BMI_CEM_Model *, const char *, int *);
int BMI_CEM_Get_current_time (BMI_CEM_Model *, double *);
int BMI_CEM_Get_start_time (BMI_CEM_Model *, double *);
int BMI_CEM_Get_end_time (BMI_CEM_Model *, double *);
int BMI_CEM_Get_time_units (BMI_CEM_Model *, char *);
int BMI_CEM_Get_time_step (BMI_CEM_Model *, double *);

/* Variable getter and setter functions */
int BMI_CEM_Get_double (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Get_double_ptr (BMI_CEM_Model *, const char *, double **);
int BMI_CEM_Get_double_at_indices (BMI_CEM_Model *, const char *, double *, int *, int);

int BMI_CEM_Set_double (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Set_double_ptr (BMI_CEM_Model *, const char *, double **);
int BMI_CEM_Set_double_at_indices (BMI_CEM_Model *, const char *, int *, int, double *);

/* Grid information functions */
int BMI_CEM_Get_grid_shape (BMI_CEM_Model *, const char *, int *);
int BMI_CEM_Get_grid_spacing (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Get_grid_origin (BMI_CEM_Model *, const char *, double *);

int BMI_CEM_Get_grid_x (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Get_grid_y (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Get_grid_z (BMI_CEM_Model *, const char *, double *);

int BMI_CEM_Get_grid_cell_count (BMI_CEM_Model *, const char *, int *);
int BMI_CEM_Get_grid_point_count (BMI_CEM_Model *, const char *, int *);
int BMI_CEM_Get_grid_vertex_count (BMI_CEM_Model *, const char *, int *);

int BMI_CEM_Get_grid_connectivity (BMI_CEM_Model *, const char *, int *);
int BMI_CEM_Get_grid_offset (BMI_CEM_Model *, const char *, int *);






// Assumes arrays start at 0, and have contiguous elements (unit stride).
double BMI_CEM_Get_0d_double (BMI_CEM_Model *, const char *);
double *BMI_CEM_Get_1d_double (BMI_CEM_Model *, const char *, int [1]);
double *BMI_CEM_Get_2d_double (BMI_CEM_Model *, const char *, int [2]);
double *BMI_CEM_Get_3d_double (BMI_CEM_Model *, const char *, int [3]);
double *BMI_CEM_Get_1d_double_at_indices (BMI_CEM_Model *, const char *, int *,
    int , double *);
double *BMI_CEM_Get_2d_double_at_indices (BMI_CEM_Model *, const char *, int *, int);
// A more general getter
//double *BMI_CEM_Get_double (BMI_CEM_Model *, const char *, int *, int **);

void BMI_CEM_Set_0d_double (BMI_CEM_Model *, const char *, double);
void BMI_CEM_Set_1d_double (BMI_CEM_Model *, const char *, const double *);
void BMI_CEM_Set_2d_double (BMI_CEM_Model *, const char *, const double *);
void BMI_CEM_Set_3d_double (BMI_CEM_Model *, const char *, const double *);
void BMI_CEM_Set_2d_double_at_indices (BMI_CEM_Model *, const char *, int *,
    const double *, int);

int BMI_CEM_Get_0d_int (BMI_CEM_Model *, const char *);
int *BMI_CEM_Get_1d_int (BMI_CEM_Model *, const char *, int [1]);
int *BMI_CEM_Get_2d_int (BMI_CEM_Model *, const char *, int [2]);
int *BMI_CEM_Get_3d_int (BMI_CEM_Model *, const char *, int [3]);
int *BMI_CEM_Get_2d_int_at_indices (BMI_CEM_Model *, const char *, int *, int);

void BMI_CEM_Set_0d_int (BMI_CEM_Model *, const char *, int);
void BMI_CEM_Set_1d_int (BMI_CEM_Model *, const char *, const int *);
void BMI_CEM_Set_2d_int (BMI_CEM_Model *, const char *, const int *);
void BMI_CEM_Set_3d_int (BMI_CEM_Model *, const char *, const int *);
void BMI_CEM_Set_2d_int_at_indices (BMI_CEM_Model *, const char *, int *,
    const int *, int);

/*
int *BMI_CEM_Get_grid_dimen (BMI_CEM_Model *, const char *, int *);
double *BMI_CEM_Get_grid_res (BMI_CEM_Model *, const char *, int *);
double *BMI_CEM_Get_grid_corner (BMI_CEM_Model *, const char *, int *);
 */

/*
  IElementSet get_Element_Set (BMI_CEM_Model *handle);
  IValueSet get_Value_Set (BMI_CEM_Model *handle, char *long_var_name, ITimeStamp);
*/

// Since these are just wrappers for other BMI functions, I don't
// think they should be included in the interface definition. They
// could be CMI functions.
int BMI_CEM_Is_scalar (BMI_CEM_Model *, const char *);
int BMI_CEM_Is_vector (BMI_CEM_Model *, const char *);
int BMI_CEM_Is_grid (BMI_CEM_Model *, const char *);
int BMI_CEM_Has_var (BMI_CEM_Model *, const char *);

// However, something that indicates if the grid is raster, or
// uniform rectilinear would be needed.
int BMI_CEM_Is_raster_grid (BMI_CEM_Model *, const char *);

#if defined(__cplusplus)
}
#endif

#endif
