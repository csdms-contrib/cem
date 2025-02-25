#ifndef BMI_H_INCLUDED
#define BMI_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#define BMI_SUCCESS (0)
#define BMI_FAILURE (1)

#define BMI_MAX_UNITS_NAME (2048)
#define BMI_MAX_COMPONENT_NAME (2048)
#define BMI_MAX_VAR_NAME (2048)


typedef struct BMI_Model {
  void * data;

  int (*initialize)(struct BMI_Model* self, const char* config_file);
  int (* update)(struct BMI_Model*);
  int (* update_until)(struct BMI_Model *, double);
  int (* update_frac)(struct BMI_Model *, double);
  int (* finalize)(struct BMI_Model *);
  int (* run_model)(struct BMI_Model *);

  int (* get_component_name)(struct BMI_Model *, char *);
  int (* get_input_var_name_count)(struct BMI_Model *, int *);
  int (* get_output_var_name_count)(struct BMI_Model *, int *);
  int (* get_input_var_names)(struct BMI_Model *, char **);
  int (* get_output_var_names)(struct BMI_Model *, char **);

  int (* get_var_grid)(struct BMI_Model *, const char *, int *);
  int (* get_var_type)(struct BMI_Model *, const char *, char *);
  int (* get_var_units)(struct BMI_Model *, const char *, char *);
  int (* get_var_itemsize)(struct BMI_Model *, const char *, int *);
  int (* get_var_nbytes)(struct BMI_Model *, const char *, int *);
  int (* get_var_location)(struct BMI_Model *, const char *, char *);
  int (* get_current_time)(struct BMI_Model *, double *);
  int (* get_start_time)(struct BMI_Model *, double *);
  int (* get_end_time)(struct BMI_Model *, double *);
  int (* get_time_units)(struct BMI_Model *, char *);
  int (* get_time_step)(struct BMI_Model *, double *);

  /* Variable getter and setter functions */
  int (* get_value)(struct BMI_Model *, const char *, void *);
  int (* get_value_ptr)(struct BMI_Model *, const char *, void **);
  int (* get_value_at_indices)(struct BMI_Model *, const char *, void *, int *, int);

  int (* set_value)(struct BMI_Model *, const char *, void *);
  int (* set_value_ptr)(struct BMI_Model *, const char *, void **);
  int (* set_value_at_indices)(struct BMI_Model *, const char *, int *, int, void *);

  /* Grid information functions */
  int (* get_grid_rank)(struct BMI_Model *, int, int *);
  int (* get_grid_size)(struct BMI_Model *, int, int *);
  int (* get_grid_type)(struct BMI_Model *, int, char *);
  int (* get_grid_shape)(struct BMI_Model *, int, int *);
  int (* get_grid_spacing)(struct BMI_Model *, int, double *);
  int (* get_grid_origin)(struct BMI_Model *, int, double *);

  int (* get_grid_x)(struct BMI_Model *, int, double *);
  int (* get_grid_y)(struct BMI_Model *, int, double *);
  int (* get_grid_z)(struct BMI_Model *, int, double *);

  int (* get_grid_face_count)(struct BMI_Model *, int, int *);
  int (* get_grid_point_count)(struct BMI_Model *, int, int *);
  int (* get_grid_vertex_count)(struct BMI_Model *, int, int *);

  int (* get_grid_connectivity)(struct BMI_Model *, int, int *);
  int (* get_grid_offset)(struct BMI_Model *, int, int *);
} BMI_Model;


#if defined(__cplusplus)
}
#endif

#endif
