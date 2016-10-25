#include <stdlib.h>
#include <glib.h>

#include "bmi_cem.h"
/*
  model->get_var_nbytes = get_var_nbytes;
  model->get_current_time = get_current_time;
  model->get_start_time = get_start_time;
  model->get_end_time = get_end_time;
  model->get_time_units = get_time_units;
  model->get_time_step = get_time_step;
*/

static BMI_Model* setup_model(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  register_bmi_cem(model);
  return model;
}


static void test_bmi_var_grid(void) {
  BMI_Model *model = setup_model();
  int status;
  int grid;

  status = model->get_var_grid(model->self, "not-a-var-name", &grid);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_var_grid(model->self, "land_surface__elevation", &grid);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(grid, ==, 2);

  status = model->get_var_grid(model->self, "sea_water__depth", &grid);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(grid, ==, 2);
}


static void test_bmi_var_type(void) {
  BMI_Model *model = setup_model();
  int status;
  char type[2048];

  status = model->get_var_type(model->self, "not-a-var-name", type);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_var_type(model->self, "land_surface__elevation", type);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(type, ==, "double");

  status = model->get_var_type(model->self, "sea_water__depth", type);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(type, ==, "double");
}


static void test_bmi_var_units(void) {
  BMI_Model *model = setup_model();
  int status;
  char units[2048];

  status = model->get_var_units(model->self, "not-a-var-name", units);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_var_units(model->self, "land_surface__elevation", units);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(units, ==, "meter");

  status = model->get_var_units(model->self, "sea_water__depth", units);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(units, ==, "meter");
}


static void test_bmi_var_nbytes(void) {
  BMI_Model *model = setup_model();
  int status;
  int nbytes;

  status = model->get_var_nbytes(model->self, "not-a-var-name", &nbytes);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_var_nbytes(model->self, "land_surface__elevation", &nbytes);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(nbytes, ==, 80000);

  status = model->get_var_nbytes(model->self, "sea_water__depth", &nbytes);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(nbytes, ==, 80000);
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/bmi/var/grid", &test_bmi_var_grid);
  g_test_add_func("/bmi/var/type", &test_bmi_var_type);
  g_test_add_func("/bmi/var/units", &test_bmi_var_units);
  g_test_add_func("/bmi/var/nbytes", &test_bmi_var_nbytes);

  return g_test_run();
}
