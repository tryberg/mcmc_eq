CC      = gcc
CC      = /Library/Developer/CommandLineTools/usr/bin/gcc

all:	 mcmc_eq  analyse_eq fw fw_mod

clean:
	\rm -f *.o mcmc_eq  fw  core analyse_eq fw_mod

mod_grd.o:  mod_grd.c  mc.h
	$(CC) -O4 -c mod_grd.c -o mod_grd.o

time_2d.o:  time_2d.c  
	$(CC) -O4 -c time_2d.c -o time_2d.o

mcmc_eq:  mcmc_eq.c  mc.h time_2d.o mod_grd.o interpol.c misfit.c
	$(CC) -O4  mcmc_eq.c time_2d.o mod_grd.o -o mcmc_eq -lm 

analyse_eq: analyse_eq.c mc.h mod_grd.o
	$(CC) -O4 analyse_eq.c mod_grd.o -o analyse_eq -lm 

fw:  fw.c  mc.h time_2d.o mod_grd.o interpol.c misfit.c
	$(CC) -O4 fw.c time_2d.o mod_grd.o -o fw -lm 

fw_mod:  fw_mod.c  mc.h time_2d.o mod_grd.o interpol.c misfit.c
	$(CC) -O4 fw_mod.c time_2d.o mod_grd.o -o fw_mod -lm 
