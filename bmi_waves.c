#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "bmi.h"
#include "bmi_waves.h"

#include "waves_model.h"


typedef struct {
    const char *name;
    const char *units;
    const char *type;
    int itemsize;
} VarInfo;


static const VarInfo variables[] = {
    {
        "sea_surface_water_wave__min_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity",
        "radians", "double", sizeof(double)
    },
    {
        "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity",
        "radians", "double", sizeof(double)
    },
    {
        "sea_surface_water_wave__mean_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity",
        "radians", "double", sizeof(double)
    },
    {
        "sea_surface_water_wave__max_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity",
        "radians", "double", sizeof(double)
    },
    {
        "sea_surface_water_wave__height",
        "meters", "double", sizeof(double)
    },
    {
        "sea_surface_water_wave__period",
        "seconds", "double", sizeof(double)
    },
    {
        "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_highness_parameter",
        "", "double", sizeof(double)
    },
    {
        "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_asymmetry_parameter",
        "", "double", sizeof(double)
    }
};

#define VAR_COUNT (sizeof(variables) / sizeof(variables[0]))

static const VarInfo*
find_variable(const char *name) {
	size_t i;
    for (i = 0; i < VAR_COUNT; i++) {
        if (strcmp(name, variables[i].name) == 0) {
            return &variables[i];
        }
    }
    return NULL;
}


static int
get_component_name (Bmi *self, char * name)
{
    strncpy (name, "waves", BMI_MAX_COMPONENT_NAME);
    return BMI_SUCCESS;
}


#define INPUT_VAR_NAME_COUNT (4)
static const char *input_var_names[INPUT_VAR_NAME_COUNT] = {
    "sea_surface_water_wave__height",
    "sea_surface_water_wave__period",
    "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_highness_parameter",
    "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_asymmetry_parameter"
};


static int
get_input_item_count(Bmi *self, int *count)
{
    *count = INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}


static int
get_input_var_names(Bmi *self, char **names)
{
    int i;
    for (i=0; i<INPUT_VAR_NAME_COUNT; i++) {
        strncpy(names[i], input_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}


#define OUTPUT_VAR_NAME_COUNT (6)
static const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {
    "sea_surface_water_wave__min_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity",
    "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity",
    "sea_surface_water_wave__mean_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity",
    "sea_surface_water_wave__max_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity",
    "sea_surface_water_wave__height",
    "sea_surface_water_wave__period"
};


static int
get_output_item_count(Bmi *self, int *count)
{
    *count = OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}


static int
get_output_var_names(Bmi *self, char **names)
{
    int i;
    for (i=0; i<OUTPUT_VAR_NAME_COUNT; i++) {
        strncpy(names[i], output_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}


static int
get_start_time(Bmi * self, double *time)
{
    *time = 0.0;
    return BMI_SUCCESS;
}


static int
get_end_time(Bmi * self, double *time)
{
    WavesModel *model = (WavesModel*)self->data;
    *time = model->end * model->time_step;
    return BMI_SUCCESS;
}


static int
get_current_time(Bmi * self, double *time)
{
    WavesModel *model = (WavesModel*)self->data;
    *time = model->now * model->time_step;
    return BMI_SUCCESS;
}


static int
get_time_step(Bmi * self, double *dt)
{
    WavesModel *model = (WavesModel*)self->data;
    *dt = model->time_step;
    return BMI_SUCCESS;
}


static int
get_time_units(Bmi * self, char *units)
{
    strncpy(units, "d", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}


static int
initialize(Bmi* handle, const char * file)
{
    WavesModel * self = (WavesModel*)handle->data;
    double end_time = 20.;
    double wave_height = 2.;
    double wave_period = 7.;
    double angle_highness = 0.2;
    double angle_asymmetry = 0.5;

    if (file) {
      FILE *fp = fopen(file, "r");
      if (fp) {
        int n_assigned;
        n_assigned = fscanf(fp, "%lf, %lf, %lf, %lf, %lf", &end_time,
            &wave_height, &wave_period, &angle_highness, &angle_asymmetry);
        if (n_assigned != 5)
          return BMI_FAILURE;
      }
      else
        return BMI_FAILURE;
    }

    waves_set_height(self, wave_height);
    waves_set_period(self, wave_period);
    waves_set_angle_highness(self, angle_highness);
    waves_set_angle_asymmetry(self, angle_asymmetry);

    {
      WavesModel *p = (WavesModel *) self;
      p->end = end_time / p->time_step;
    }

  return BMI_SUCCESS;
}


static int
update_frac(Bmi * self, double f)
{
    WavesModel *p = (WavesModel *) self->data;
    double now;
    //int until_time_step = p->time_step * f;

    get_current_time(self, &now);
    //until_time_step = now + p->time_step * f;

    //waves_run_until (p, until_time_step);
    waves_run_until (p, now + p->time_step * f);

    return BMI_SUCCESS;
}


static int
update(Bmi * self)
{
    return update_frac(self, 1.);
}


static int
update_until(Bmi * self, double then)
{
    double dt;
    double now;

    if (get_time_step(self, &dt) == BMI_FAILURE)
        return BMI_FAILURE;

    if (get_current_time(self, &now) == BMI_FAILURE)
        return BMI_FAILURE;

    {
        int n;
        const double n_steps = (then - now) / dt;
        for (n=0; n<(int)n_steps; n++) {
            if (update(self) == BMI_FAILURE)
                return BMI_FAILURE;
        }

        if (update_frac(self, n_steps - (int)n_steps) == BMI_FAILURE)
            return BMI_FAILURE;
    }

    return BMI_SUCCESS;
}


static int
finalize(Bmi * self)
{
    waves_destroy ((WavesModel*)self->data);

    return BMI_SUCCESS;
}


static int
get_var_grid(Bmi *self, const char *name, int *grid)
{
    *grid = -1;
    return BMI_FAILURE;
}


static int
get_var_type(Bmi *self, const char *name, char *type)
{
    const VarInfo *var = find_variable(name);
    if (!var) {
        type[0] = '\0';
        return BMI_FAILURE;
    }
    strncpy(type, var->type, BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}


static int
get_var_units(Bmi *self, const char *name, char *units)
{
    const VarInfo *var = find_variable(name);
    if (!var) {
        units[0] = '\0';
        return BMI_FAILURE;
    }
    strncpy(units, var->units, BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}


static int
get_var_itemsize(Bmi *self, const char *name, int *itemsize)
{
    const VarInfo *var = find_variable(name);
    if (!var) {
        *itemsize = 0;
        return BMI_FAILURE;
    }
    *itemsize = var->itemsize;
    return BMI_SUCCESS;
}


static int
get_var_nbytes(Bmi *self, const char *name, int *nbytes)
{
    int itemsize;
	const int size = 1;

    if (get_var_itemsize(self, name, &itemsize) == BMI_FAILURE)
        return BMI_FAILURE;

    *nbytes = itemsize * size;

    return BMI_SUCCESS;
}


static int
get_var_location(Bmi *self, const char *name, char *location)
{
    const VarInfo *var = find_variable(name);
    if (!var) {
        location[0] = '\0';
        return BMI_FAILURE;
    }
    strncpy(location, "none", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}


static int
get_value(Bmi *self, const char *name, void *dest)
{
    double *dptr = (double*)dest;

    if (strcmp(name, "sea_surface_water_wave__min_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *dptr = waves_get_wave_angle_min ((WavesModel*)self->data);
    } else if (strcmp(name, "sea_surface_water_wave__azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *dptr = waves_get_wave_angle ((WavesModel*)self->data);
    } else if (strcmp(name, "sea_surface_water_wave__mean_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *dptr =  waves_get_wave_angle_mean ((WavesModel*)self->data);
    } else if (strcmp(name, "sea_surface_water_wave__max_of_increment_of_azimuth_angle_of_opposite_of_phase_velocity") == 0) {
        *dptr = waves_get_wave_angle_max ((WavesModel*)self->data);
    } else if (strcmp(name, "sea_surface_water_wave__height") == 0) {
        *dptr = waves_get_height ((WavesModel*)self->data);
    } else if (strcmp(name, "sea_surface_water_wave__period") == 0) {
        *dptr = waves_get_period ((WavesModel*)self->data);
    } else {
        return BMI_FAILURE;
    }

    return BMI_SUCCESS;
}


static int
set_value (Bmi *self, const char *name, void *array)
{
    if (strcmp(name, "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_asymmetry_parameter") == 0)
      waves_set_angle_asymmetry((WavesModel*)self->data, *(double*)array);
    else if (strcmp(name, "sea_shoreline_wave~incoming~deepwater__ashton_et_al_approach_angle_highness_parameter") == 0)
      waves_set_angle_highness((WavesModel*)self->data, *(double*)array);
    else if (strcmp(name, "sea_surface_water_wave__height") == 0)
      waves_set_height((WavesModel*)self->data, *(double*)array);
    else if (strcmp(name, "sea_surface_water_wave__period") == 0)
      waves_set_period((WavesModel*)self->data, *(double*)array);

    return BMI_SUCCESS;
}


Bmi*
register_bmi_waves(Bmi *model)
{
    model->data = waves_new();

    model->initialize = initialize;
    model->update = update;
    model->update_until = update_until;
    model->finalize = finalize;

    model->get_component_name = get_component_name;
    model->get_input_item_count = get_input_item_count;
    model->get_output_item_count = get_output_item_count;
    model->get_input_var_names = get_input_var_names;
    model->get_output_var_names = get_output_var_names;

    model->get_var_grid = get_var_grid;
    model->get_var_type = get_var_type;
    model->get_var_units = get_var_units;
    model->get_var_nbytes = get_var_nbytes;
    model->get_var_location = get_var_location;
    model->get_var_itemsize = get_var_itemsize;
    model->get_current_time = get_current_time;
    model->get_start_time = get_start_time;
    model->get_end_time = get_end_time;
    model->get_time_units = get_time_units;
    model->get_time_step = get_time_step;

    model->get_value = get_value;
    model->get_value_ptr = NULL;
    model->get_value_at_indices = NULL;

    model->set_value = set_value;
    model->set_value_at_indices = NULL;

    model->get_grid_rank = NULL;
    model->get_grid_size = NULL;
    model->get_grid_type = NULL;
    model->get_grid_shape = NULL;
    model->get_grid_spacing = NULL;
    model->get_grid_origin = NULL;

    model->get_grid_x = NULL;
    model->get_grid_y = NULL;
    model->get_grid_z = NULL;

    return model;
}
