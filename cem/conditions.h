#ifndef CONDITIONS_INCLUDED
#define CONDITIONS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

void InitConds(double** cell_depth, char** all_beach, char** all_rock, double** percent_full_sand, double** percent_full_rock, char** type_of_rock, double** topography);

#if defined(__cplusplus)
}
#endif

#endif
