/*
 * =====================================================================
 *
 *       Filename:  waves_cli.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/16/2010 04:50:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eric Hutton (eh), huttone@colorado.edu
 *        Company:  Community Surface Dynamics Modeling System
 *
 * =====================================================================
 */

#if !defined (WAVES_CLI_H)
#define WAVES_CLI_H

#include <glib.h>

typedef struct
{
  gboolean verbose;
  gboolean silent;
  gboolean version;
  gboolean radians;
  gint n_waves;
  double highness;
  double asymmetry;
  guint32 seed;
}
waves_args_st;

waves_args_st* parse_command_line (int argc, char* argv[], GError** error);

#endif

