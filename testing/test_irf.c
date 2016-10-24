#include <stdlib.h>
#include <glib.h>

#include "bmi_cem.h"


static void test_bmi_initialize(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  int status;

  register_bmi_cem(model);

  status = model->initialize(NULL, &(model->self));

  g_assert_cmpint(status, ==, BMI_SUCCESS);
}


static void test_bmi_update(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  int status;

  register_bmi_cem(model);

  status = model->initialize(NULL, &(model->self));
  g_assert_cmpint(status, ==, BMI_SUCCESS);

  status = model->update(&(model->self));
  g_assert_cmpint(status, ==, BMI_SUCCESS);
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/bmi/initialize", &test_bmi_initialize);
  g_test_add_func("/bmi/update", &test_bmi_update);

  return g_test_run();
}
