#ifndef CEM_GLOBALS_INCLUDED
#define CEM_GLOBALS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "consts.h"

extern double **topography;

// Special SWAN matrices.
extern double **shelf_depth; // SWAN bathymetry.
extern double **wave_h_sig; // SWAN wave heights.
extern double **wave_dir; // SWAN wave angles.
extern int current_time_step; // Time step of current calculation
extern double stop_after; // Stop after what number of time steps
extern double time_step; // days; reflects rate of sediment transport per time step

#if defined(__cplusplus)
}
#endif

#endif
