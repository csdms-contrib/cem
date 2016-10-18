#ifndef CEM_GLOBALS_INCLUDED
#define CEM_GLOBALS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "consts.h"

extern double **Topography;

// Special SWAN matrices.
extern double **ShelfDepth; // SWAN bathymetry.
extern double **Hsig; // SWAN wave heights.
extern double **Dir; // SWAN wave angles.
extern int CurrentTimeStep; // Time step of current calculation
extern double StopAfter; // Stop after what number of time steps
extern double TimeStep; // days; reflects rate of sediment transport per time step

#if defined(__cplusplus)
}
#endif

#endif
