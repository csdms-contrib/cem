#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "cem/utils.h"


static void test_set_periodic_multiple_of_4(void) {
  const int len = 12;
  int array[12] =    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  int expected[12] = {6, 7, 8, 3, 4, 5, 6, 7, 8, 3,  4,  5};

  apply_periodic_boundary(array, sizeof(int), len);

  g_assert_cmpmem(array, len, expected, len);
}


static void test_set_periodic_even_width(void) {
  const int len = 10;
  int array[10] =    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int expected[10] = {5, 6, 2, 3, 4, 5, 6, 2, 3, 4};

  apply_periodic_boundary(array, sizeof(int), len);

  g_assert_cmpmem(array, len, expected, len);
}


static void test_set_periodic_odd_width(void) {
  const int len = 11;
  int array[11] =    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int expected[11] = {5, 6, 7, 3, 4, 5, 6, 7, 3, 4,  5};

  apply_periodic_boundary(array, sizeof(int), len);

  g_assert_cmpmem(array, len, expected, len);
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/cem/utils/set_periodic_4", &test_set_periodic_multiple_of_4);
  g_test_add_func("/cem/utils/set_periodic_even", &test_set_periodic_even_width);
  g_test_add_func("/cem/utils/set_periodic_odd", &test_set_periodic_odd_width);

  return g_test_run();
}
