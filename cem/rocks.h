#ifndef ROCKS_INCLUDED
#define ROCKS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif


void set_rock_blocks(char **rock_type, double **topography, const int n_rows,
    const int n_cols, const int n_blocks);


#if defined(__cplusplus)
}
#endif

#endif
