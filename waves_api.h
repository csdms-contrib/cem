#if !defined( WAVES_API_H )
#define WAVES_API_H 

#if !defined(TRUE)
#define TRUE (1)
#endif

#if !defined(FALSE)
#define FALSE (0)
#endif

/** An opaque data structure that holds the state of an instance of a 
waves model.
*/
typedef struct
{
}
Waves_state;

Waves_state* waves_new     (void);
Waves_state* waves_destroy (Waves_state*);

Waves_state* waves_init     (Waves_state*);
int          waves_run_until(Waves_state*, double);
Waves_state* waves_finalize (Waves_state * self, int free);

Waves_state* waves_set_angle_asymmetry (Waves_state*, double asymmetry);
Waves_state* waves_set_angle_highness (Waves_state*, double highness);
Waves_state* waves_set_height (Waves_state*, double height_in_m);
Waves_state* waves_set_period (Waves_state*, double height_in_s);

double waves_get_angle_asymmetry (Waves_state*);
double waves_get_angle_highness (Waves_state*);
double waves_get_wave_angle (Waves_state*);
double waves_get_wave_angle_mean (Waves_state*);
double waves_get_wave_angle_max (Waves_state*);
double waves_get_wave_angle_min (Waves_state*);
double waves_get_height (Waves_state*);
double waves_get_period (Waves_state*);

const char** waves_get_exchange_items (void);
double waves_get_value (Waves_state * self, const char* value);

double waves_get_current_time (Waves_state* self);
double waves_get_start_time (Waves_state* self);
double waves_get_end_time (Waves_state* self);

#endif

