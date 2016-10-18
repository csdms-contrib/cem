#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bmi.h"
#include "bmi_cem.h"
#include "globals.h"

#define return_on_error(stmt)                 \
  {                                           \
    const int status = (stmt);                \
    if (status != BMI_SUCCESS) return status; \
  }

int cem_initialize(void);
int cem_update(void);
int cem_update_until(int);
int cem_finalize(void);

static int get_component_name(void *self, char *name) {
  strncpy(name, "Coastline Evolution Model", BMI_MAX_COMPONENT_NAME);
  return BMI_SUCCESS;
}

#define INPUT_VAR_NAME_COUNT (2)
const char *input_var_names[INPUT_VAR_NAME_COUNT] = {
    "sea_surface_water_wave__significant_height",
    "sea_surface_water_wave__azimuth_angle_of_group_velocity"};

#define OUTPUT_VAR_NAME_COUNT (2)
const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {
    "sea_water__depth", "land_surface__elevation"};

static int get_input_var_name_count(void *self, int *count) {
  *count = INPUT_VAR_NAME_COUNT;
  return BMI_SUCCESS;
}

static int get_input_var_names(void *self, char **names) {
  int i;
  for (i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
    strncpy(names[i], input_var_names[i], BMI_MAX_VAR_NAME);
  }
  return BMI_SUCCESS;
}

static int get_output_var_name_count(void *self, int *count) {
  *count = OUTPUT_VAR_NAME_COUNT;
  return BMI_SUCCESS;
}

static int get_output_var_names(void *self, char **names) {
  int i;
  for (i = 0; i < OUTPUT_VAR_NAME_COUNT; i++) {
    strncpy(names[i], output_var_names[i], BMI_MAX_VAR_NAME);
  }
  return BMI_SUCCESS;
}

static int get_start_time(void *self, double *time) {
  if (time) {
    *time = 0.;
    return BMI_SUCCESS;
  } else {
    return BMI_FAILURE;
  }
}

static int get_end_time(void *self, double *time) {
  *time = StopAfter;
  return BMI_SUCCESS;
}

static int get_current_time(void *self, double *time) {
  *time = CurrentTimeStep * TimeStep;
  return BMI_SUCCESS;
}

static int get_time_step(void *self, double *dt) {
  *dt = TimeStep;
  return BMI_SUCCESS;
}

static int get_time_units(void *self, char *units) {
  strncpy(units, "d", BMI_MAX_UNITS_NAME);
  return BMI_SUCCESS;
}

static int initialize(const char *config_file, void **handle) {
  CemModel *self = NULL;
  int status;

  if (!handle) return BMI_FAILURE;

  self = malloc(sizeof(CemModel));

  status = cem_initialize();

  *handle = self;

  if (status == 0)
    return BMI_SUCCESS;
  else
    return BMI_FAILURE;
}

static int update(void *self) {
  if (cem_update() == 0)
    return BMI_SUCCESS;
  else
    return BMI_FAILURE;
}

static int update_frac(void *self, double f) {
  if (f > 0) {
    double dt;

    get_time_step(self, &dt);

    TimeStep = f * dt;

    update(self);

    TimeStep = dt;
  }

  return BMI_SUCCESS;
}

static int update_until(void *self, double t) {
  {
    double dt;
    double now;

    get_time_step(self, &dt);
    get_current_time(self, &now);

    {
      int n;
      const double n_steps = (t - now) / dt;
      const int n_full_steps = (int)n_steps;
      for (n = 0; n < n_full_steps; n++) {
        update(self);
      }

      update_frac(self, n_steps - n_full_steps);
    }
  }

  return BMI_SUCCESS;
}

static int finalize(void *self) {
  if (self) {
    cem_finalize();
  }

  return BMI_SUCCESS;
}

static int get_grid_rank(void *self, int grid_id, int *rank) {
  if (grid_id == 0) {
    *rank = 2;
    return BMI_SUCCESS;
  } else {
    *rank = -1;
    return BMI_FAILURE;
  }
}

static int get_grid_size(void *self, int grid_id, int *size) {
  if (grid_id == 0) {
    *size = Xmax * Ymax;
    return BMI_SUCCESS;
  } else {
    *size = -1;
    return BMI_FAILURE;
  }
}

static int get_grid_shape(void *self, int grid_id, int *shape) {
  if (grid_id == 0) {
    shape[0] = Xmax;
    shape[1] = Ymax;
    return BMI_SUCCESS;
  } else {
    return BMI_FAILURE;
  }
}

static int get_grid_spacing(void *self, int grid_id, double *spacing) {
  if (grid_id == 0) {
    spacing[0] = CellWidth;
    spacing[1] = CellWidth;
    return BMI_SUCCESS;
  } else {
    return BMI_FAILURE;
  }
}

static int get_grid_origin(void *self, int grid_id, double *origin) {
  if (grid_id == 0) {
    origin[0] = 0.;
    origin[1] = 0.;
    return BMI_SUCCESS;
  } else {
    return BMI_FAILURE;
  }
}

static int get_grid_type(void *self, int grid_id, char *type) {
  if (grid_id == 0) {
    strncpy(type, "uniform_rectilinear", 2048);
    return BMI_SUCCESS;
  } else {
    *type = '\0';
    return BMI_FAILURE;
  }
}

static int get_var_grid(void *self, const char *name, int *grid_id) {
  if (strcmp(name, "sea_water__depth") == 0 ||
      strcmp(name, "land_surface__elevation") == 0) {
    *grid_id = 0;
    return BMI_SUCCESS;
  } else {
    *grid_id = -1;
    return BMI_FAILURE;
  }
}

static int get_var_type(void *self, const char *long_var_name, char *type) {
  if (strcmp(long_var_name, "sea_water__depth") == 0 ||
      strcmp(long_var_name, "land_surface__elevation") == 0) {
    strncpy(type, "double", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
  } else {
    *type = '\0';
    return BMI_FAILURE;
  }
}

static int get_var_units(void *self, const char *long_var_name, char *units) {
  if (strcmp(long_var_name, "sea_water__depth") == 0 ||
      strcmp(long_var_name, "land_surface__elevation") == 0) {
    strncpy(units, "meter", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
  } else {
    units[0] = '\0';
    return BMI_FAILURE;
  }
}

static int get_var_itemsize(void *self, const char *name, int *itemsize) {
  if (strcmp(name, "sea_water__depth") == 0 ||
      strcmp(name, "land_surface__elevation") == 0) {
    *itemsize = sizeof(double);
    return BMI_SUCCESS;
  } else {
    *itemsize = 0;
    return BMI_FAILURE;
  }
}

static int get_var_nbytes(void *self, const char *name, int *nbytes) {
  int id, size, itemsize;

  return_on_error(get_var_grid(self, name, &id));
  return_on_error(get_grid_size(self, id, &size));
  return_on_error(get_var_itemsize(self, name, &itemsize));

  *nbytes = itemsize * size;

  return BMI_SUCCESS;
}

static int get_value_ptr(void *self, const char *long_var_name, void **dest) {
  void *data = NULL;

  if (strcmp(long_var_name, "sea_water__depth") == 0)
    data = (void *)ShelfDepth[0];
  else if (strcmp(long_var_name, "land_surface__elevation") == 0)
    data = (void *)Topography[0];
  else if (strcmp(long_var_name,
                  "sea_surface_water_wave__significant_height") == 0)
    data = (void *)Hsig[0];
  else if (strcmp(long_var_name,
                  "sea_surface_water_wave__azimuth_angle_of_group_velocity") ==
           0)
    data = (void *)Dir[0];

  *dest = data;

  return (data == NULL) ? BMI_FAILURE : BMI_SUCCESS;
}

static int get_value(void *self, const char *name, void *dest) {
  void *src = NULL;

  if (get_value_ptr(self, name, &src) == BMI_FAILURE) {
    return BMI_FAILURE;
  } else {
    int nbytes;
    return_on_error(get_var_nbytes(self, name, &nbytes));
    memcpy(dest, src, nbytes);
  }

  return BMI_SUCCESS;
}

static int set_value(void *self, const char *name, void *array) {
  void *dest = NULL;

  if (get_value_ptr(self, name, &dest) == BMI_FAILURE) {
    return BMI_FAILURE;
  } else {
    int nbytes = 0;
    return_on_error(get_var_nbytes(self, name, &nbytes));
    memcpy(dest, array, nbytes);
    return BMI_SUCCESS;
  }
}

BMI_Model *register_bmi_cem(BMI_Model *model) {
  model->self = NULL;

  model->initialize = initialize;
  model->update = update;
  model->update_until = update_until;
  model->update_frac = update_frac;
  model->finalize = finalize;
  model->run_model = NULL;

  model->get_component_name = get_component_name;
  model->get_input_var_name_count = get_input_var_name_count;
  model->get_output_var_name_count = get_output_var_name_count;
  model->get_input_var_names = get_input_var_names;
  model->get_output_var_names = get_output_var_names;

  model->get_var_grid = get_var_grid;
  model->get_var_type = get_var_type;
  model->get_var_units = get_var_units;
  model->get_var_nbytes = get_var_nbytes;
  model->get_current_time = get_current_time;
  model->get_start_time = get_start_time;
  model->get_end_time = get_end_time;
  model->get_time_units = get_time_units;
  model->get_time_step = get_time_step;

  model->get_value = get_value;
  model->get_value_ptr = get_value_ptr;
  model->get_value_at_indices = NULL;

  model->set_value = set_value;
  model->set_value_ptr = NULL;
  model->set_value_at_indices = NULL;

  model->get_grid_rank = get_grid_rank;
  model->get_grid_size = get_grid_size;
  model->get_grid_type = get_grid_type;
  model->get_grid_shape = get_grid_shape;
  model->get_grid_spacing = get_grid_spacing;
  model->get_grid_origin = get_grid_origin;

  model->get_grid_x = NULL;
  model->get_grid_y = NULL;
  model->get_grid_z = NULL;

  return model;
}
