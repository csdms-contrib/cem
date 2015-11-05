#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "bmi.h"

/* Implement this: Add model-specific includes */
#include "cem_model.h"


#define return_on_error(stmt) { \
  const int status = (stmt); \
  if (status != BMI_SUCCESS) \
    return status; \
}

static int
get_component_name (void *self, char * name)
{
    strncpy (name, "cem", BMI_MAX_COMPONENT_NAME);
    return BMI_SUCCESS;
}


#define INPUT_VAR_NAME_COUNT (7)
static const char *input_var_names[INPUT_VAR_NAME_COUNT] = {
    "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity",
    "basin_outlet_water_sediment~bedload__mass_flow_rate",
    "land_surface_water_sediment~bedload__mass_flow_rate",
    "sea_surface_water_wave__period",
    "basin_outlet_water_sediment~suspended__mass_flow_rate",
    "sea_surface_water_wave__height",
    "land_surface__elevation"
};


static int
get_input_var_name_count(void *self, int *count)
{
    *count = INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}


static int
get_input_var_names(void *self, char **names)
{
    int i;
    for (i=0; i<INPUT_VAR_NAME_COUNT; i++) {
        strncpy(names[i], input_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}


#define OUTPUT_VAR_NAME_COUNT (5)
static const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {
    "basin_outlet~coastal_center__x_coordinate",
    "basin_outlet~coastal_water_sediment~bedload__mass_flow_rate",
    "land_surface__elevation",
    "sea_water__depth",
    "basin_outlet~coastal_center__y_coordinate"
};


static int
get_output_var_name_count(void *self, int *count)
{
    *count = OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}


static int
get_output_var_names(void *self, char **names)
{
    int i;
    for (i=0; i<OUTPUT_VAR_NAME_COUNT; i++) {
        strncpy(names[i], output_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}


static int
get_start_time(void * self, double *time)
{
    *time = 0.0;
    return BMI_SUCCESS;
}


static int
get_end_time(void * self, double *time)
{ /* Implement this: Set end time */
    *time = deltas_get_end_time((CemModel*)self);
    return BMI_SUCCESS;
}


static int
get_current_time(void * self, double *time)
{ /* Implement this: Set current time */
    *time = deltas_get_current_time((CemModel*)self);
    return BMI_SUCCESS;
}


static int
get_time_step(void * self, double *dt)
{ /* Implement this: Set time step */
    *dt = deltas_get_time_step((CemModel*)self);
    return BMI_SUCCESS;
}


static int
get_time_units(void * self, char *units)
{
    strncpy(units, "d", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}


static int
initialize(const char * file, void **handle)
{ /* Implement this: Create and initialize a model handle */
    return cem_initialize(file, (CemModel**)handle);
}


static int
update_frac(void * self, double f)
{ /* Implement this: Update for a fraction of a time step */
    return BMI_FAILURE;
}


static int
update(void * self)
{
    cem_advance_one_time_step((CemModel*)self);
    return BMI_SUCCESS;
}


static int
update_until(void * self, double then)
{
    double dt;
    double now;

    return_on_error(get_time_step(self, &dt));
    return_on_error(get_current_time(self, &now));

    {
        int n;
        const double n_steps = (then - now) / dt;
        for (n=0; n<(int)n_steps; n++) {
          return_on_error(update(self));
        }
    }

    return BMI_SUCCESS;
}


static int
finalize(void * self)
{ /* Implement this: Clean up */
    cem_finalize((CemModel*)self);
    return BMI_SUCCESS;
}


static int
get_grid_type(void *self, int id, char *type)
{
    if (id == 0) {
        strncpy(type, "scalar", 2048);
    } else if (id == 1) {
        strncpy(type, "vector", 2048);
    } else if (id == 2) {
        strncpy(type, "uniform_rectilinear", 2048);
    } else {
        type[0] = '\0'; return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_grid_rank(void *self, int id, int *rank)
{
    if (id == 0) {
        *rank = 0;
    } else if (id == 1) {
        *rank = 1;
    } else if (id == 2) {
        *rank = 2;
    } else {
        *rank = -1; return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_grid_shape(void *self, int id, int *shape)
{ /* Implement this: set shape of structured grids */
    if (id == 2) {
        shape[0] = deltas_get_nx((CemModel*)self);
        shape[1] = deltas_get_ny((CemModel*)self) / 2;
    } else {
        return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_grid_spacing(void *self, int id, double *spacing)
{ /* Implement this: set spacing of uniform rectilinear grids */
    if (id == 2) {
        spacing[0] = deltas_get_dx((CemModel*)self);
        spacing[1] = deltas_get_dy((CemModel*)self);
    } else {
        return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_grid_origin(void *self, int id, double *origin)
{ /* Implement this: set origin of uniform rectilinear grids */
    if (id == 2) {
        origin[0] = 0.;
        origin[1] = 0.;
    } else {
        return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_grid_size(void *self, int id, int *size)
{
    if (id == 0)
        *size = 1;
    else if (id == 1)
        *size = deltas_get_n_rivers((CemModel*)self);
    else if (id == 2)
        *size = deltas_get_nx((CemModel*)self) * deltas_get_ny((CemModel*)self);
    else
        return BMI_FAILURE;

    return BMI_SUCCESS;
}


static int
get_var_grid(void *self, const char *name, int *grid)
{
    if (strcmp(name, "basin_outlet~coastal_center__x_coordinate") == 0) {
        *grid = 1;
    } else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *grid = 0;
    } else if (strcmp(name, "basin_outlet_water_sediment~bedload__mass_flow_rate") == 0) {
        *grid = 0;
    } else if (strcmp(name, "basin_outlet~coastal_water_sediment~bedload__mass_flow_rate") == 0) {
        *grid = 1;
    } else if (strcmp(name, "land_surface_water_sediment~bedload__mass_flow_rate") == 0) {
        *grid = 2;
    } else if (strcmp(name, "sea_surface_water_wave__period") == 0) {
        *grid = 0;
    } else if (strcmp(name, "land_surface__elevation") == 0) {
        *grid = 2;
    } else if (strcmp(name, "sea_water__depth") == 0) {
        *grid = 2;
    } else if (strcmp(name, "basin_outlet_water_sediment~suspended__mass_flow_rate") == 0) {
        *grid = 0;
    } else if (strcmp(name, "sea_surface_water_wave__height") == 0) {
        *grid = 0;
    } else if (strcmp(name, "basin_outlet~coastal_center__y_coordinate") == 0) {
        *grid = 1;
    } else {
      fprintf(stderr, "bad grid. returning %d", BMI_FAILURE);
        *grid = -1; return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_var_type(void *self, const char *name, char *type)
{
    if (strcmp(name, "basin_outlet~coastal_center__x_coordinate") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet_water_sediment~bedload__mass_flow_rate") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet~coastal_water_sediment~bedload__mass_flow_rate") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "land_surface_water_sediment~bedload__mass_flow_rate") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_surface_water_wave__period") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "land_surface__elevation") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_water__depth") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet_water_sediment~suspended__mass_flow_rate") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_surface_water_wave__height") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet~coastal_center__y_coordinate") == 0) {
        strncpy(type, "double", BMI_MAX_UNITS_NAME);
    } else {
        type[0] = '\0'; return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_var_units(void *self, const char *name, char *units)
{
    if (strcmp(name, "basin_outlet~coastal_center__x_coordinate") == 0) {
        strncpy(units, "meters", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        strncpy(units, "radians", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet_water_sediment~bedload__mass_flow_rate") == 0) {
        strncpy(units, "kg / s", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet~coastal_water_sediment~bedload__mass_flow_rate") == 0) {
        strncpy(units, "kg / s", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "land_surface_water_sediment~bedload__mass_flow_rate") == 0) {
        strncpy(units, "kg / s", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_surface_water_wave__period") == 0) {
        strncpy(units, "seconds", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "land_surface__elevation") == 0) {
        strncpy(units, "meters", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_water__depth") == 0) {
        strncpy(units, "meters", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet_water_sediment~suspended__mass_flow_rate") == 0) {
        strncpy(units, "kg / s", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "sea_surface_water_wave__height") == 0) {
        strncpy(units, "meters", BMI_MAX_UNITS_NAME);
    } else if (strcmp(name, "basin_outlet~coastal_center__y_coordinate") == 0) {
        strncpy(units, "meters", BMI_MAX_UNITS_NAME);
    } else {
        units[0] = '\0'; return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_var_itemsize(void *self, const char *name, int *itemsize)
{
    if (strcmp(name, "basin_outlet~coastal_center__x_coordinate") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "basin_outlet_water_sediment~bedload__mass_flow_rate") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "basin_outlet~coastal_water_sediment~bedload__mass_flow_rate") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "land_surface_water_sediment~bedload__mass_flow_rate") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "sea_surface_water_wave__period") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "land_surface__elevation") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "sea_water__depth") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "basin_outlet_water_sediment~suspended__mass_flow_rate") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "sea_surface_water_wave__height") == 0) {
        *itemsize = sizeof(double);
    } else if (strcmp(name, "basin_outlet~coastal_center__y_coordinate") == 0) {
        *itemsize = sizeof(double);
    } else {
        *itemsize = 0; return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}


static int
get_var_nbytes(void *self, const char *name, int *nbytes)
{
    int id, size, itemsize;

    return_on_error(get_var_grid(self, name, &id));
    return_on_error(get_grid_size(self, id, &size));
    return_on_error(get_var_itemsize(self, name, &itemsize));

    *nbytes = itemsize * size;

    return BMI_SUCCESS;
}


static int
get_var_ndim(void *self, const char *name, int *ndim)
{
    int id, rank;

    return_on_error(get_var_grid(self, name, &id));
    return_on_error(get_grid_rank(self, id, &rank));

    *ndim = rank;

    return BMI_SUCCESS;
}


static int
get_var_stride(void *self, const char *name, int *stride)
{
    int id;

    return_on_error(get_var_grid(self, name, &id));

    if (id == 1) {
        stride[0] = 1;
    }
    else if (id == 2) {
        stride[0] = deltas_get_ny((CemModel*)self);
        stride[1] = 1;
    }

    return BMI_SUCCESS;
}


static int
get_value_ptr(void *self, const char *name, void **dest)
{
    if (strcmp(name, "basin_outlet~coastal_center__x_coordinate") == 0) {
        *dest = (double*)deltas_get_river_x_position((CemModel*)self);
    } else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *dest = &(((CemModel*)self)->WaveAngle);
    } else if (strcmp(name, "basin_outlet_water_sediment~bedload__mass_flow_rate") == 0) {
        *dest = NULL;
    } else if (strcmp(name, "basin_outlet~coastal_water_sediment~bedload__mass_flow_rate") == 0) {
        *dest = (double*)deltas_get_river_flux((CemModel*)self);
    } else if (strcmp(name, "land_surface_water_sediment~bedload__mass_flow_rate") == 0) {
        *dest = NULL;
    } else if (strcmp(name, "sea_surface_water_wave__period") == 0) {
        *dest = &(((CemModel*)self)->wave_period);
    } else if (strcmp(name, "land_surface__elevation") == 0) {
        *dest = NULL;
    } else if (strcmp(name, "sea_water__depth") == 0) {
        *dest = (double*)deltas_get_depth((CemModel*)self) + deltas_get_ny((CemModel*)self) / 2;
    } else if (strcmp(name, "basin_outlet_water_sediment~suspended__mass_flow_rate") == 0) {
        *dest = NULL;
    } else if (strcmp(name, "sea_surface_water_wave__height") == 0) {
        *dest = &(((CemModel*)self)->wave_height);
    } else if (strcmp(name, "basin_outlet~coastal_center__y_coordinate") == 0) {
        *dest = (double*)deltas_get_river_y_position((CemModel*)self);
    } else {
        *dest = NULL; return BMI_FAILURE;
    }

    if (*dest)
        return BMI_SUCCESS;
    else
        return BMI_FAILURE;
}


int
get_value(void * self, const char * name, void *dest)
{
    if (strcmp(name, "sea_water__depth") == 0) {
        deltas_get_depth_dup((CemModel*)self, dest);
        return BMI_SUCCESS;
    }
    else if (strcmp(name, "land_surface__elevation") == 0) {
        deltas_get_elevation_dup((CemModel*)self, dest);
        return BMI_SUCCESS;
    }
    else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0)
        *(double*)dest = ((CemModel*)self)->WaveAngle;
    else if (strcmp(name, "sea_surface_water_wave__height") == 0)
        *(double*)dest = ((CemModel*)self)->wave_height;
    else if (strcmp(name, "sea_surface_water_wave__period") == 0)
        *(double*)dest = ((CemModel*)self)->wave_period;
    else {
        void *src = NULL;
        int nbytes = 0;

        return_on_error(get_value_ptr(self, name, &src));
        return_on_error(get_var_nbytes(self, name, &nbytes));

        memcpy(dest, src, nbytes);
    }
    return BMI_SUCCESS;
}


static int
get_value_at_indices (void *self, const char *name, void *dest,
    int * inds, int len)
{
    void *src = NULL;
    int itemsize = 0;

    return_on_error(get_value_ptr(self, name, &src));
    return_on_error(get_var_itemsize(self, name, &itemsize));

    { /* Copy the data */
        int i;
        int offset;
        char * ptr;
        for (i=0, ptr=(char*)dest; i<len; i++, ptr+=itemsize) {
            offset = inds[i] * itemsize;
            memcpy (ptr, (char*)src + offset, itemsize);
        }
    }

    return BMI_SUCCESS;
}


static int
set_value (void *self, const char *name, void *array)
{
    void * dest = NULL;
    int nbytes = 0;

    if (strcmp(name, "land_surface_water_sediment~bedload__mass_flow_rate") == 0)
        deltas_set_sediment_flux_grid ((CemModel*)self, (double*)array);
    if (strcmp(name, "land_surface__elevation") == 0)
        deltas_set_elevation_grid ((CemModel*)self, (double*)array);
    else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0)
      ((CemModel*)self)->WaveAngle = *((double*)array);
    else if (strcmp(name, "sea_surface_water_wave__height") == 0)
      ((CemModel*)self)->wave_height = *((double*)array);
    else if (strcmp(name, "sea_surface_water_wave__period") == 0)
      ((CemModel*)self)->wave_period = *((double*)array);
    else {
        return_on_error(get_value_ptr(self, name, &dest));
        return_on_error(get_var_nbytes(self, name, &nbytes));

        memcpy (dest, array, nbytes);
    }

    return BMI_SUCCESS;
}


static int
set_value_at_indices (void *self, const char *name, int * inds, int len,
    void *src)
{
    void * to = NULL;
    int itemsize = 0;

    return_on_error(get_value_ptr(self, name, &to));
    return_on_error(get_var_itemsize(self, name, &itemsize));

    { /* Copy the data */
        int i;
        int offset;
        char * ptr;
        for (i=0, ptr=(char*)src; i<len; i++, ptr+=itemsize) {
            offset = inds[i] * itemsize;
            memcpy ((char*)to + offset, ptr, itemsize);
        }
    }
    return BMI_SUCCESS;
}


BMI_Model*
register_bmi_cem(BMI_Model *model)
{
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
    model->get_value_at_indices = get_value_at_indices;

    model->set_value = set_value;
    model->set_value_ptr = NULL;
    model->set_value_at_indices = set_value_at_indices;

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
