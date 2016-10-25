#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "bmi/bmi_cem.h"


static BMI_Model* setup_model(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  register_bmi_cem(model);
  return model;
}


static void test_bmi_name(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  int status;
  char name[2048];

  register_bmi_cem(model);

  status = model->get_component_name(model->self, name);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(name, ==, "Coastline Evolution Model");
}


static void test_bmi_in_var_count(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  int status;
  int n_names;

  register_bmi_cem(model);

  status = model->get_input_var_name_count(model->self, &n_names);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(n_names, ==, 2);
}


static void test_bmi_out_var_count(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  int status;
  int n_names;

  register_bmi_cem(model);

  status = model->get_output_var_name_count(model->self, &n_names);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(n_names, ==, 2);
}


static void test_bmi_in_vars(void) {
  BMI_Model *model = setup_model();
  int status;
  char **names;

  names = (char**)malloc(2 * sizeof(char*));
  names[0] = (char*)malloc(2 * 2048 * sizeof(char));
  names[1] = names[0] + 2048;

  status = model->get_input_var_names(model->self, names);
  g_assert_cmpint(status, ==, BMI_SUCCESS);

  g_assert_cmpstr(names[0], ==, "sea_surface_water_wave__significant_height");
  g_assert_cmpstr(names[1], ==, "sea_surface_water_wave__azimuth_angle_of_group_velocity");
}


static void test_bmi_out_vars(void) {
  BMI_Model *model = setup_model();
  int status;
  char **names;

  names = (char**)malloc(2 * sizeof(char*));
  names[0] = (char*)malloc(2 * 2048 * sizeof(char));
  names[1] = names[0] + 2048;

  status = model->get_output_var_names(model->self, names);
  g_assert_cmpint(status, ==, BMI_SUCCESS);

  g_assert_cmpstr(names[0], ==, "sea_water__depth");
  g_assert_cmpstr(names[1], ==, "land_surface__elevation");
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/bmi/meta/name", &test_bmi_name);
  g_test_add_func("/bmi/meta/in_var_count", &test_bmi_in_var_count);
  g_test_add_func("/bmi/meta/out_var_count", &test_bmi_out_var_count);
  g_test_add_func("/bmi/meta/in_vars", &test_bmi_in_vars);
  g_test_add_func("/bmi/meta/out_vars", &test_bmi_out_vars);

  return g_test_run();
}
