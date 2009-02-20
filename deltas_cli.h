#if !defined( CEM_CLI_H )
#define CEM_CLI_H

typedef struct
{
  int verbose;
}
cem_args_st;

cem_args_st *parse_command_line (int argc, char *argv[]);

#endif

