#include "waves_cli.h"
#include <glib.h>

static gboolean verbose = FALSE;
static gboolean silent = TRUE;
static gboolean version = FALSE;
static gboolean radians = TRUE;
static gint n_waves = 1;
static double highness = 0.5;
static double asymmetry = 0.5;
static gint seed = 1945;
      
static GOptionEntry entries[] =
{
  {"verbose", 'V', 0, G_OPTION_ARG_NONE, &verbose, "Be verbose", NULL} ,
  {"silent", 'S', G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &silent, "Be silent", NULL} ,
  {"version", 'v', 0, G_OPTION_ARG_NONE, &version, "Version number", NULL} ,
  {"seed", 's', 0, G_OPTION_ARG_INT, &seed, "Seed for rng", "SEED"} ,
  {"n-waves", 'n', 0, G_OPTION_ARG_INT, &n_waves, "Number of waves", "N"} ,
  {"highness", 'h', 0, G_OPTION_ARG_DOUBLE, &highness,
   "Highness factor [0,1)", "VAL"} ,
  {"asymmetry", 'a', 0, G_OPTION_ARG_DOUBLE, &asymmetry,
   "Asymmetry [0,1)", "VAL"} ,
  {"radians", 'r', 0, G_OPTION_ARG_NONE, &radians, "Angles in radians", NULL} ,
  {"degrees", 'd', G_OPTION_FLAG_REVERSE, G_OPTION_ARG_NONE, &radians,
   "Angles in degrees", NULL} ,
  {NULL}
};

/** Parse command line options for waves.

Usage: waves [options]

Use: 'waves --help' for a full list of command line options.
*/    
waves_args_st*
parse_command_line (int argc, char* argv[], GError** error)
{
  waves_args_st* args = NULL;

  g_assert ((error==NULL || *error==NULL));
  g_assert (argv);

  {
    GError* tmp_err = NULL;
    GOptionContext* context = g_option_context_new (
                                "Generate random incoming wave angles" );

    g_option_context_add_main_entries (context, entries, NULL);
    if (!g_option_context_parse (context, &argc, &argv, &tmp_err))
    {
      g_propagate_error (error, tmp_err);
      args = NULL;
    }
    else
    {
      args = g_new (waves_args_st, 1);

      args->verbose = verbose;
      args->silent = silent;
      args->version = version;
      args->radians = radians;
      args->n_waves = n_waves;
      args->highness = highness;
      args->asymmetry = asymmetry;
    }
  }

  return args;
}

