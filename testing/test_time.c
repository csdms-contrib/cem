#include <stdlib.h>
#include <glib.h>

#include "bmi/bmi_cem.h"


static BMI_Model* setup_model(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  register_bmi_cem(model);
  return model;
}


static void test_bmi_current_time(void) {
  BMI_Model *model = setup_model();
  int status;
  double time;

  status = model->get_current_time(model->self, &time);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpfloat(time, ==, 0.);
}


static void test_bmi_start_time(void) {
  BMI_Model *model = setup_model();
  int status;
  double time;

  status = model->get_start_time(model->self, &time);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpfloat(time, ==, 0.);
}


static void test_bmi_time_step(void) {
  BMI_Model *model = setup_model();
  int status;
  double time_step;

  status = model->get_time_step(model->self, &time_step);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpfloat(time_step, ==, 1.);
}


static void test_bmi_time_units(void) {
  BMI_Model *model = setup_model();
  int status;
  char units[2048];

  status = model->get_time_units(model->self, units);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(units, ==, "d");
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/bmi/time/now", &test_bmi_current_time);
  g_test_add_func("/bmi/time/start", &test_bmi_start_time);
  g_test_add_func("/bmi/time/time_step", &test_bmi_time_step);
  g_test_add_func("/bmi/time/units", &test_bmi_time_units);

  return g_test_run();
}
