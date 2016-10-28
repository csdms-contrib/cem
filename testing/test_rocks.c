#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "cem/rocks.h"


static void test_set_rock_blocks(void) {
  const int shape[2] = {10, 20};
  double **z = (double**)malloc(sizeof(double*) * shape[0]);
  char **rock_type = (char**)malloc(sizeof(char*) * shape[0]);
  int i;
  char expected[20] = {
    'f', 'f', 's', 's', 's',
    'f', 'f', 's', 's', 's', 'f', 'f', 's', 's', 's',
    'f', 'f', 's', 's', 's'
  };

  z[0] = (double*)malloc(sizeof(double) * shape[0] * shape[1]);
  rock_type[0] = (char*)malloc(sizeof(char) * shape[0] * shape[1]);
  for (i=1; i < shape[0]; i++) {
    z[i] = z[i-1] + shape[1];
    rock_type[i] = rock_type[i-1] + shape[1];
  }

  set_rock_blocks(rock_type, z, shape[0], shape[1], 2);

  for (i=0; i < shape[0]; i++)
    g_assert_cmpmem(rock_type[i], shape[1], expected, shape[1]);

  free(rock_type[0]);
  free(rock_type);
  free(z[0]);
  free(z);
}


int main(int argc, char* argv[]) {
  g_test_init(&argc, &argv, NULL);
  g_test_add_func("/cem/rocks/set_blocks", &test_set_rock_blocks);

  return g_test_run();
}
