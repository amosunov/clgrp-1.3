AUTOMAKE_OPTIONS = foreign
CFLAGS=-Wall -std=gnu99 -O3
LDFLAGS=
bin_PROGRAMS = clgrp verify
lib_LIBRARIES = libclgrp.a
libclgrp_a_SOURCES = functions.c sieve.c clgrp.c verify.c
clgrp_SOURCES = functions.c sieve.c clgrp.c clgrp_main.c
verify_SOURCES = functions.c sieve.c verify.c verify_main.c
clgrpincludedir = $(includedir)/libclgrp
clgrpinclude_HEADERS = functions.h sieve.h clgrp.h verify.h
