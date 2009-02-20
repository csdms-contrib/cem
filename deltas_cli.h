#if !defined( CEM_CLI_H )
#define CEM_CLI_H

typedef struct
{
  int verbose;
  float stop_time;
  char* out_prefix;
}
cem_args_st;

cem_args_st *parse_command_line (int argc, char *argv[]);

#endif

