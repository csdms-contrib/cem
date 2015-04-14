#if !defined (WAVES_MODEL_H)
#define WAVES_MODEL_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <glib.h>


typedef struct
{
  double asymmetry;
  double highness;
  double height;
  double period;

  gint now;
  gint end;
  double time_step;

  double *angles;
  gint len;
  GRand *rand;
  guint seed;
}
WavesModel;


int waves_run_until (WavesModel * s, double until);
WavesModel* waves_new (void);
void waves_init_state (WavesModel * s);
void waves_free_state (WavesModel * s);
WavesModel * waves_destroy (WavesModel * self);
WavesModel * waves_finalize (WavesModel * self, int free);

WavesModel * waves_set_angle_asymmetry (WavesModel * self, double asymmetry);
WavesModel * waves_set_angle_highness (WavesModel * self, double highness);
WavesModel * waves_set_height (WavesModel * self, double height_in_m);
WavesModel * waves_set_period (WavesModel * self, double period_in_s);
double waves_get_angle_asymmetry (const WavesModel * self);
double waves_get_angle_highness (const WavesModel * self);
double waves_get_height (const WavesModel * self);
double waves_get_period (const WavesModel * self);
double waves_get_wave_angle (WavesModel * self);
double waves_get_wave_angle_max (WavesModel * self);
double waves_get_wave_angle_min (WavesModel * self);
double waves_get_wave_angle_mean (WavesModel * self);

#if defined(__cplusplus)
}
#endif

#endif
