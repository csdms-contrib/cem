#include <stdlib.h>
#include <glib.h>

#include "bmi/bmi_cem.h"


static BMI_Model* setup_model(void) {
  BMI_Model *model = (BMI_Model *)malloc(sizeof(BMI_Model));
  register_bmi_cem(model);
  return model;
}


static void test_bmi_grid_rank(void) {
  BMI_Model *model = setup_model();
  int status;
  int rank;

  status = model->get_grid_rank(model->self, -1, &rank);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_grid_rank(model->self, 2, &rank);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(rank, ==, 2);
}


static void test_bmi_grid_size(void) {
  BMI_Model *model = setup_model();
  int status;
  int size;

  status = model->get_grid_size(model->self, -1, &size);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_grid_size(model->self, 2, &size);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(size, ==, 10000);
}


static void test_bmi_grid_type(void) {
  BMI_Model *model = setup_model();
  int status;
  char type[2048];

  status = model->get_grid_type(model->self, -1, type);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_grid_type(model->self, 2, type);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpstr(type, ==, "uniform_rectilinear");
}


static void test_bmi_grid_shape(void) {
  BMI_Model *model = setup_model();
  int status;
  int shape[2] = {-1, -1};

  status = model->get_grid_shape(model->self, -1, shape);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_grid_shape(model->self, 2, shape);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpint(shape[0], ==, 50);
  g_assert_cmpint(shape[1], ==, 200);
}


static void test_bmi_grid_spacing(void) {
  BMI_Model *model = setup_model();
  int status;
  double spacing[2] = {-1., -1.};

  status = model->get_grid_spacing(model->self, -1, spacing);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_grid_spacing(model->self, 2, spacing);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpfloat(spacing[0], ==, 100.);
  g_assert_cmpfloat(spacing[1], ==, 100.);
}


static void test_bmi_grid_origin(void) {
  BMI_Model *model = setup_model();
  int status;
  double origin[2] = {-1., -1.};

  status = model->get_grid_origin(model->self, -1, origin);
  g_assert_cmpint(status, ==, BMI_FAILURE);

  status = model->get_grid_origin(model->self, 2, origin);
  g_assert_cmpint(status, ==, BMI_SUCCESS);
  g_assert_cmpfloat(origin[0], ==, 0.);
  g_assert_cmpfloat(origin[1], ==, 0.);
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/bmi/grid/rank", &test_bmi_grid_rank);
  g_test_add_func("/bmi/grid/size", &test_bmi_grid_size);
  g_test_add_func("/bmi/grid/type", &test_bmi_grid_type);
  g_test_add_func("/bmi/grid/shape", &test_bmi_grid_shape);
  g_test_add_func("/bmi/grid/spacing", &test_bmi_grid_spacing);
  g_test_add_func("/bmi/grid/origin", &test_bmi_grid_origin);

  return g_test_run();
}
