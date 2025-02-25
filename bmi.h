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


typedef struct Bmi {
  void * data;

  int (*initialize)(struct Bmi* self, const char* config_file);
  int (* update)(struct Bmi*);
  int (* update_until)(struct Bmi *, double);
  int (* update_frac)(struct Bmi *, double);
  int (* finalize)(struct Bmi *);
  int (* run_model)(struct Bmi *);

  int (* get_component_name)(struct Bmi *, char *);
  int (* get_input_var_name_count)(struct Bmi *, int *);
  int (* get_output_var_name_count)(struct Bmi *, int *);
  int (* get_input_var_names)(struct Bmi *, char **);
  int (* get_output_var_names)(struct Bmi *, char **);

  int (* get_var_grid)(struct Bmi *, const char *, int *);
  int (* get_var_type)(struct Bmi *, const char *, char *);
  int (* get_var_units)(struct Bmi *, const char *, char *);
  int (* get_var_itemsize)(struct Bmi *, const char *, int *);
  int (* get_var_nbytes)(struct Bmi *, const char *, int *);
  int (* get_var_location)(struct Bmi *, const char *, char *);
  int (* get_current_time)(struct Bmi *, double *);
  int (* get_start_time)(struct Bmi *, double *);
  int (* get_end_time)(struct Bmi *, double *);
  int (* get_time_units)(struct Bmi *, char *);
  int (* get_time_step)(struct Bmi *, double *);

  /* Variable getter and setter functions */
  int (* get_value)(struct Bmi *, const char *, void *);
  int (* get_value_ptr)(struct Bmi *, const char *, void **);
  int (* get_value_at_indices)(struct Bmi *, const char *, void *, int *, int);

  int (* set_value)(struct Bmi *, const char *, void *);
  int (* set_value_ptr)(struct Bmi *, const char *, void **);
  int (* set_value_at_indices)(struct Bmi *, const char *, int *, int, void *);

  /* Grid information functions */
  int (* get_grid_rank)(struct Bmi *, int, int *);
  int (* get_grid_size)(struct Bmi *, int, int *);
  int (* get_grid_type)(struct Bmi *, int, char *);
  int (* get_grid_shape)(struct Bmi *, int, int *);
  int (* get_grid_spacing)(struct Bmi *, int, double *);
  int (* get_grid_origin)(struct Bmi *, int, double *);

  int (* get_grid_x)(struct Bmi *, int, double *);
  int (* get_grid_y)(struct Bmi *, int, double *);
  int (* get_grid_z)(struct Bmi *, int, double *);

  int (* get_grid_face_count)(struct Bmi *, int, int *);
  int (* get_grid_point_count)(struct Bmi *, int, int *);
  int (* get_grid_vertex_count)(struct Bmi *, int, int *);

  int (* get_grid_connectivity)(struct Bmi *, int, int *);
  int (* get_grid_offset)(struct Bmi *, int, int *);
} Bmi;


#if defined(__cplusplus)
}
#endif

#endif
