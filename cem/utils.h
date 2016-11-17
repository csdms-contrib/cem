#ifndef UTILS_INCLUDED
#define UTILS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include <stdlib.h>


void **malloc2d(size_t n_rows, size_t n_cols, size_t itemsize, size_t pointersize);
void free2d(void **mem);
void apply_periodic_boundary(void *array, const int itemsize,
    const int nitems);
void repeat_mem(void *dst, const int len, void *block, int block_len);
void stripe_cem_matrix(void **matrix, int n_rows, int n_cols, int itemsize,
    void *stripe, int nitems);


#if defined(__cplusplus)
}
#endif

#endif
