#include <stdlib.h>
#include <string.h>

#include "utils.h"


void **malloc2d(size_t n_rows, size_t n_cols, size_t itemsize)
{
  size_t i;
  void **matrix = malloc(sizeof(void*) * n_rows);
  matrix[0] = malloc(itemsize * n_rows * n_cols);
  for (i = 1; i < n_rows ; i++)
    matrix[i] = matrix[i - 1] + n_cols * itemsize;

  return matrix;
}


void free2d(void **mem)
{
  free(mem[0]);
  free(mem);
}


void apply_periodic_boundary(void *array, const int itemsize, const int nitems)
{
  const int width_items = nitems / 2;
  const int left_items = (nitems - width_items) / 2;
  const int right_items = nitems - (width_items + left_items);
  const int nbytes = nitems * itemsize;
  const int width_nbytes = width_items * itemsize;
  const int left_nbytes = left_items * itemsize;
  const int right_nbytes = right_items * itemsize;

  memcpy(array, array + width_nbytes, left_nbytes);
  memcpy(array + nbytes - right_nbytes, array + left_nbytes, right_nbytes);
}


// Repeat a block within an array until the array is filled.
void repeat_mem(void *dst, const int len, void *block, int block_len) {
  const int n_blocks = len / block_len;
  const int last_block_len = len - n_blocks * block_len;
  int offset;

  for (offset=0; offset < len; offset += block_len) {
    memcpy(dst + offset, block, block_len);
  }
  memcpy(dst + len - last_block_len, block, last_block_len);
}


// Run vertical stripes down a CEM matrix, making certain to obey the
// periodic boundary conditions.
void stripe_cem_matrix(void **matrix, int n_rows, int n_cols, int itemsize,
    void *stripe, int nitems)
{
  const int offset_to_left = (n_cols / 4) * itemsize;
  const int width = (n_cols / 2) * itemsize;
  void *row = malloc(itemsize * n_cols);
  int i;

  repeat_mem(row + offset_to_left, width, stripe, nitems * itemsize);
  apply_periodic_boundary(row, itemsize, n_cols);
  for (i = 0; i < n_rows; i++)
    memcpy(matrix[i], row, itemsize * n_cols);

  free(row);
}
