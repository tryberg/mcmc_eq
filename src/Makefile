CC      = gcc

# -O3      : highest standard optimization level (the old "-O4" was silently clamped to -O3)
# -march=native : tune/vectorize for the build machine's CPU.
#           NOTE: ties the binary to this CPU family. On a heterogeneous cluster, either build
#           on the target node or drop -march=native (e.g. use -mtune=generic) for portability.
# -flto    : link-time optimization; lets the compiler inline find_in_cell/dst/traveltimet across
#           translation units.
CFLAGS  = -O3 -march=native -flto

all:	 mcmc_eq  analyse_eq fw fw_mod

clean:
	\rm -f *.o mcmc_eq  fw  core analyse_eq fw_mod

mod_grd.o:  mod_grd.c  mc.h
	$(CC) $(CFLAGS) -c mod_grd.c -o mod_grd.o

time_2d.o:  time_2d.c  
	$(CC) $(CFLAGS) -c time_2d.c -o time_2d.o

mcmc_eq:  mcmc_eq.c  forward.c  mc.h time_2d.o mod_grd.o
	$(CC) $(CFLAGS)  mcmc_eq.c time_2d.o mod_grd.o -o mcmc_eq -lm 

analyse_eq: analyse_eq.c mc.h mod_grd.o
	$(CC) $(CFLAGS) analyse_eq.c mod_grd.o -o analyse_eq -lm 

fw:  fw.c  mc.h time_2d.o mod_grd.o
	$(CC) $(CFLAGS) fw.c time_2d.o mod_grd.o -o fw -lm 

fw_mod:  fw_mod.c  mc.h time_2d.o mod_grd.o
	$(CC) $(CFLAGS) fw_mod.c time_2d.o mod_grd.o -o fw_mod -lm 
