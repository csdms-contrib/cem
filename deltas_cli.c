#include "deltas_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#define CEM_PROGRAM_STR   "cem"
#define CEM_MAJOR_VERSION (0)
#define CEM_MINOR_VERSION (1)
#define CEM_MICRO_VERSION (2)

#define CEM_DEFAULT_END_TIME (100)
#define CEM_DEFAULT_OUT_PREFIX "fileout"

static int verbose_flag = 0;

static int version_flag = 0;

static int help_flag = 0;

static struct option cem_long_opts[] = {
  {"verbose", no_argument, &verbose_flag, 1},
  {"brief", no_argument, &verbose_flag, 0},
  {"version", no_argument, &version_flag, 1},
  {"help", no_argument, &help_flag, 1},
  {"stop-time", required_argument, NULL, 's'},
  {"out-prefix", required_argument, NULL, 'o'},
  {NULL, 0, NULL, 0}
};

static const char *help_msg[] = {
  "Usage: deltas [options]",
  "Options:",
  "  -V or --verbose      Be verbose",
  "  --brief              Be terse",
  "  -v or --version      Print version number and exit",
  "  -?,-h or --help      Print this information and exit",
  "  --stop-time or -s    Model end time (in days)",
  "  --out-prefix or -o   Prefix for output files",
  NULL
};

/** Parse command line options for deltas.

Usage: deltas [options]

Use: 'deltas --help' for a full list of command line options.
*/
cem_args_st *
parse_command_line (int argc, char *argv[])
{
  cem_args_st *args = NULL;

  if (argv)
  {
    int ch;

    float end = CEM_DEFAULT_END_TIME;

    char *out_prefix = NULL;

    while ((ch = getopt_long (argc, argv, "vVh?", cem_long_opts, NULL)) != -1)
    {
      switch (ch)
      {
      case 's':
        sscanf (optarg, "%f", &end);
        break;
      case 'o':
        out_prefix = strdup (optarg);
        break;
      case 'v':
        version_flag = 1;
        break;
      case 'V':
        verbose_flag = 1;
        break;
      case '?':
      case 'h':
        help_flag = 1;
        break;
      case 0:
        break;
      }
    }

    if (help_flag)
    {
      const char **str;

      for (str = help_msg; *str; str++)
        fprintf (stdout, "%s\n", *str);
      exit (EXIT_SUCCESS);
    }

    if (version_flag)
    {
      fprintf (stdout, "%s version %d.%d.%d\n",
               CEM_PROGRAM_STR, CEM_MAJOR_VERSION,
               CEM_MINOR_VERSION, CEM_MICRO_VERSION);
      exit (EXIT_SUCCESS);
    }

    if (optind < argc)
    {
      fprintf (stderr, "Error: Extra command line arguments.\n");
      fprintf (stderr, "%s --help for help\n", CEM_PROGRAM_STR);
      exit (EXIT_FAILURE);
    }

    if (!out_prefix)
      out_prefix = strdup (CEM_DEFAULT_OUT_PREFIX);

    if (verbose_flag)
    {
      fprintf (stdout, "End time is %f\n", end);
      fprintf (stdout, "Output prefix is %s\n", out_prefix);
    }
    args = (cem_args_st*) malloc (sizeof (cem_args_st));
    args->verbose = verbose_flag;
    args->stop_time = end;
    args->out_prefix = out_prefix;
  }

  return args;
}
