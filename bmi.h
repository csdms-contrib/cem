#ifndef BMI_CEM_INCLUDED
#define BMI_CEM_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi_const.h"

#define BMI_CEM_INPUT_VAR_NAME_COUNT 2
#define BMI_CEM_OUTPUT_VAR_NAME_COUNT 2

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
int BMI_CEM_Get_value (BMI_CEM_Model *, const char *, void *);
int BMI_CEM_Get_value_ptr (BMI_CEM_Model *, const char *, void **);
int BMI_CEM_Get_value_at_indices (BMI_CEM_Model *self, const char *, void *,
    int *, int);

int BMI_CEM_Set_value (BMI_CEM_Model *, const char *, void *);
int BMI_CEM_Set_value_ptr (BMI_CEM_Model *, const char *, void **);
int BMI_CEM_Set_value_at_indices (BMI_CEM_Model *, const char *, int *, int,
    void *);

int BMI_CEM_Get_double (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Get_double_ptr (BMI_CEM_Model *, const char *, double **);
int BMI_CEM_Get_double_at_indices (BMI_CEM_Model *, const char *, double *, int *, int);

int BMI_CEM_Set_double (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Set_double_ptr (BMI_CEM_Model *, const char *, double **);
int BMI_CEM_Set_double_at_indices (BMI_CEM_Model *, const char *, int *, int, double *);

/* Grid information functions */
int BMI_CEM_Get_grid_type (BMI_CEM_Model *, const char *, BMI_Grid_type *);
int BMI_CEM_Get_grid_shape (BMI_CEM_Model *, const char *, int *);
int BMI_CEM_Get_grid_spacing (BMI_CEM_Model *, const char *, double *);
int BMI_CEM_Get_grid_origin (BMI_CEM_Model *, const char *, double *);

#if defined(__cplusplus)
}
#endif

#endif
