#if !defined (WAVES_H)
#define WAVES_H

#include <glib.h>

typedef struct
{
  double asymmetry;
  double highness;
  double height;
  double period;

  gint now;
  double time_step;

  double* angles;
  gint len;
  GRand* rand;
  guint seed;
}
State;

int _waves_initialize (State* s);
int _waves_run_until (State* s, int until);
int _waves_finalize (State* s);

void waves_init_state (State* s);
void waves_free_state (State* s);

#endif

