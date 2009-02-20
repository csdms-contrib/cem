ACLOCAL_AMFLAGS       = -I m4

EXTRA_DIST            = bootstrap

bin_PROGRAMS          = deltas

deltas_SOURCES        = deltas_main.c
deltas_DEPENDENCIES   = libdeltas.la

include_HEADERS       = deltas_api.h
noinst_HEADERS        = deltas.h deltas_cli.h

lib_LTLIBRARIES       = libdeltas.la
libdeltas_la_SOURCES  = ndelta4.c ndelta4_api.c deltas_cli.c

deltas_LDADD          = -ldeltas

