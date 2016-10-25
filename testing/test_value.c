#include <stdlib.h>
#include <glib.h>

#include "bmi/bmi_cem.h"


static BMI_Model* setup_model(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  register_bmi_cem(model);
  return model;
}


static void test_bmi_get_value(void) {
  BMI_Model *model = setup_model();
  int status;
  int nbytes;
  double *buffer = NULL;

  status = model->initialize(NULL, &(model->self));
  g_assert_cmpint(status, ==, BMI_SUCCESS);

  status = model->get_value(model->self, "not-a-var-name", NULL);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_var_nbytes(model->self, "land_surface__elevation", &nbytes);
  buffer = (double*)malloc(nbytes);

  status = model->get_value(model->self, "land_surface__elevation", buffer);
  g_assert_cmpint(status, ==, BMI_SUCCESS);

  status = model->get_value(model->self, "sea_water__depth", buffer);
  g_assert_cmpint(status, ==, BMI_SUCCESS);

  free(buffer);
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/bmi/value/get", &test_bmi_get_value);

  return g_test_run();
}
