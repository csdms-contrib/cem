#include "utils.h"
#include "rocks.h"
#include "consts.h"

// LMV Assign fast and slow weathering portions
// Divide rows of rocks into fast and slow rock types.
void set_rock_blocks(char **rock_type, double **topography, const int n_rows,
    const int n_cols, const int n_blocks)
{
  const int chunk_length = (n_cols / 2) / n_blocks;
  const int n_fast = chunk_length / 2;
  const int n_slow = chunk_length - n_fast;
  int col;

  { // Set blocks of rock type.
    char *stripe = (char*)malloc(sizeof(char) * (n_fast + n_slow));
    for (col=0; col < n_fast; col++) stripe[col] = 'f';
    for (col=0; col < n_slow; col++) stripe[col + n_fast] = 's';
    stripe_cem_matrix((void**)rock_type, n_rows, n_cols, sizeof(char),
        stripe, n_fast + n_slow);
    free(stripe);
  }

  { // Set blocks of cliff height.
    double *stripe = (double*)malloc(sizeof(double) * (n_fast + n_slow));
    for (col=0; col < n_fast; col++) stripe[col] = kCliffHeightFast;
    for (col=0; col < n_slow; col++) stripe[col + n_fast] = kCliffHeightSlow;
    stripe_cem_matrix((void**)topography, n_rows, n_cols, sizeof(double),
        stripe, n_fast + n_slow);
    free(stripe);
  }

  return;
}
