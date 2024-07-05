// Program to calculate P&S traveltimes from gridded model

/* *******************************************
		     mcmc_eq
	
	Inversion of travel times to locate earthquakes
	and derive 1D velocity models using a statistical 
	Markov chain Monte Carlo method
	
	Trond Ryberg, Christian Haberland & Jeremy D. Pesicek
	
	Copyright (C) 2024
	- Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences 

	Version 2.0	 1. May  2024
	Version 3.0	 4. July 2024

   SPDX-FileCopyrightText: 2024 Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences
   SPDX-License-Identifier: GPL-3.0-only 	

   ******************************************* */
   
/*
 *
 * This file is part of the mcmc_eq package.
 *
 * mcmc_eq is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; version GPL-3.0-only.
 *
 * mcmc_eq is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with mcmc_eq; see the file COPYING.  If not, write to
 * the Free Software Foundation, 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 *
 */


// 190624 cleanup, misfit & interpol in common file

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>


#include <string.h>
#include "mc.h"

int TRIA, DR, aflag;

FILE *fpin;
int read_mcmcdata (FILE *f, struct DATA *d);
void write_mcmcdata (struct DATA *d, int ne);
float cal_fit_newx (struct Model *m, struct DATA *d, int ne, float ***tttp, float ***ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3,  int flag, int eikonal, int out); 
float traveltimet (float **ttt, int nx, int ny, int nz, float h, float dist, float z, float z0);
float dst (float x1, float x2, float y1, float y2);
void setup_table (struct Model *m, float **ttt, struct GRDHEAD gh, int ps);
void setup_table_new (struct Model *m, float ***ttt, struct GRDHEAD gh, int ps);
int output_ttt (float **ttt, int nx, int nz, float h);
void read_single_line(FILE *fp, char x[]);
void copy_model(struct Model *dest, struct Model *src);
int find_neighbor_cell(struct Model *model, int n);
int find_in_cell(struct Model *modx, float x);
void model_to_hsbuf(struct Model *modx, struct GRDHEAD gh, float *hsbuf);

#include "interpol.c"
#include "misfit.c"



// Method to allocate a 3D array of floats
float*** make_3d_array(int nx, int ny, int nz) {
    float*** arr;
    int i,j;

    arr = (float ***) malloc(nx*sizeof(float**));

    for (i = 0; i < nx; i++) {
        arr[i] = (float **) malloc(ny*sizeof(float*));

        for(j = 0; j < ny; j++) {
            arr[i][j] = (float *) malloc(nz * sizeof(float));
        }
    }

    return arr;
} 

// Method to allocate a 2D array of floats
float** make_2d_array(int nx, int ny) {
    float** arr;
    int i;

    arr = (float **) malloc(nx*sizeof(float*));

    for (i = 0; i < nx; i++) arr[i] = (float *) malloc(ny*sizeof(float));

    return arr;
} 

float rand_gauss (void)
// returns random float with zero mean, gaussian distributed, sdev=1 
{
  float v1,v2,s;
  do {
    v1 = 2.0 * ((float) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((float) rand()/RAND_MAX) - 1;
    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );
  if (s == 0.0)
    return 0.0;
  else
    return (v1*sqrt(-2.0 * log(s) / s));
}

float nexp (v)
float v;
{
    float h;
    if (v<log(FLT_MAX/1000.0)) {h=exp(v);} else {h=exp(log(FLT_MAX/1000.0));}
    return (h);
}

int out_of_bounds (v0, lower, upper)
float v0, lower, upper;
{
  if ((v0>=lower) && (v0<=upper)) return (0); else {return(1);}
}

float rand_gauss_bounded (v0, sdev, lower, upper)
float v0, sdev, lower, upper;
// returns gaussian distributed random number with zero mean & sdev, centered at v0 & limited by lower and upper
{
  float dv;
  do
  {
     dv=rand_gauss()*sdev;
  }
  while (!(((v0+dv)>lower) && ((v0+dv)<upper)));
  return (dv);
}


int rand_eq_int (int n) 
// returns random integer between 0 and n-1, eq distributed 
{
    return ((int) (((float) rand()/RAND_MAX*n)));
}

float rand_eq (void)
// returns random float between 0 and 1, eq distributed 
{
    return (((float) rand()/RAND_MAX));
}

float rand_eq_limited (min, max)
float min,max;
// returns random float between min and max, eq distributed 
{
    return (((float) min+(max-min)*rand()/RAND_MAX));
}


unsigned int get_seed (void)
{ 
	unsigned int seed;
	FILE *fp;
	fp = fopen("/dev/urandom", "r");
	if (!(fp = fopen("/dev/urandom", "r")))
	{
    		fprintf(stdout,"Could not open random device\n"); exit(0);
	}
	if (fread(&seed, sizeof seed, 1, fp)!=1)
	{
        	fprintf(stdout,"Could not read from random device\n"); exit(0);
	}
	fclose(fp);
	return(seed);
}


// --------------main-------------------------------------------------------------

int main(int argc, char *argv[])
{
 int i,j,start_cell_number,sdev_start_cell_number, topo_shift, max_dim;
 double old_ll;
 float sdevx, sdevy, sdevz, sdevvp, sdevn;
 float xmin, xmax, ymin, ymax, zmin, zmax, vpmin, vpmax;
 float start_vp, sdev_start_vp, start_noise;
 float dr_fac;
 float noise_min, noise_max;
 double old_misfit, old_rms;
 int deci,j_max_start, j_max_main,true_random,topo_flag;

 char		    line[1024], ds[1000];
 char   dstring_start[1000],dstring_main[1000], topo_file_name[1000];
 FILE 		*config_file;		/* config file	*/

// CH
 FILE 		*finp;			/* data file 			*/

 struct GRDHEAD  gh;			/* FD grid header structure 	*/
 struct DATA 	*s;		/*  data 			*/


 struct Model old_model;
 float mfp0, mfs0, mfp1, mfs1, mfp2, mfs2, mfp3, mfs3;

 float sdevxs, sdevys, sdevzs, sdevresidual;
 float residual_min, residual_max;

 int n_ppicks,n_spicks, n_ppicks0, n_spicks0, n_ppicks1, n_spicks1, n_ppicks2, n_spicks2, n_ppicks3, n_spicks3, ne;


 float 	***tttpr, ***tttsr;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */

 int nxmod,nos,eikonal;

 float start_delay, sdev_start_delay,r_start_eqh,r_start_eqv;
 float xz,vpm1, vpm2, vpmp, vsm1, vsm2, vsmp, df;
 int di,mod_from_file;
 float inv_control;

 float vpvsmin, vpvsmax, sdevvpvs, start_vpvs, sdev_start_vpvs;
 float zrmin, zrmax, xtest;


// read observations
//*  --- usage */
        if (argc<4)
        {
                fprintf(stderr, "usage: xxx config.dat output_file pick_file[s]\n");
                fprintf(stderr, "Ch. Haberland GFZ Potsdam July 2016, Lena Delta\n");
                fprintf(stderr, "V. %04.1f\n", VERSION);
                exit (0);
        }
// ini data s
	if (!(s=malloc(MAX_NOQ*sizeof(struct DATA)))) { fprintf(stderr, "malloc failed\n"); exit(0); }

/* ---  read config file */
	if (!(config_file = fopen (argv[1], "r"))) { fprintf(stderr, "could not open config file %s\n",argv[1]); exit (0);}
// read config file
        read_single_line(config_file, line); sscanf(line, "%f ", &gh.h); 
        read_single_line(config_file, line); sscanf(line, "%d ", &gh.nx); 
        read_single_line(config_file, line); sscanf(line, "%d ", &gh.ny); 
        read_single_line(config_file, line); sscanf(line, "%d ", &gh.nz); 
        read_single_line(config_file, line); sscanf(line, "%f ", &gh.x0);
        read_single_line(config_file, line); sscanf(line, "%f ", &gh.y0);
        read_single_line(config_file, line); sscanf(line, "%f ", &gh.z0); 
        read_single_line(config_file, line); sscanf(line, "%d ", &max_dim); 
	if (max_dim>MD) { fprintf(stderr, "max_dim > MD, recompile!\n"); exit (0);}
//	if (max_dim>gh.nz) { fprintf(stderr, "max_dim > NZ: -> max_dim=%d\n",gh.nz); max_dim=gh.nz;}
        read_single_line(config_file, line); sscanf(line, "%f ", &vpmin); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vpmax); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vpvsmin); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vpvsmax); 
        read_single_line(config_file, line); sscanf(line, "%f ", &noise_min); 
        read_single_line(config_file, line); sscanf(line, "%f ", &noise_max);
        read_single_line(config_file, line); sscanf(line, "%f ", &residual_min); 
        read_single_line(config_file, line); sscanf(line, "%f ", &residual_max); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevx); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevy); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevz);
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevvp); 
	read_single_line(config_file, line); sscanf(line, "%f ", &sdevvpvs); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevn);
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevxs); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevys); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevzs);
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevresidual);
        read_single_line(config_file, line); sscanf(line, "%f ", &inv_control);
        read_single_line(config_file, line); sscanf(line, "%f ", &dr_fac);
        read_single_line(config_file, line); sscanf(line, "%d ", &TRIA);
        read_single_line(config_file, line); sscanf(line, "%d %d ", &j_max_start, &j_max_main);
        read_single_line(config_file, line); sscanf(line, "%d ", &deci);
        read_single_line(config_file, line); sscanf(line, "%d %d", &true_random, &eikonal);
        read_single_line(config_file, line); sscanf(line, "%s %s", dstring_start, dstring_main);
        read_single_line(config_file, line); sscanf(line, "%d ", &aflag); mod_from_file=0; if (aflag==3) {mod_from_file=1; aflag=0;}

        read_single_line(config_file, line); sscanf(line, "%d %s %d ", &topo_flag, topo_file_name, &topo_shift);
        read_single_line(config_file, line); sscanf(line, "%f %f ", &start_vp, &sdev_start_vp); 
        read_single_line(config_file, line); sscanf(line, "%f %f ", &start_vpvs, &sdev_start_vpvs); 
        read_single_line(config_file, line); sscanf(line, "%d %d ", &start_cell_number, &sdev_start_cell_number);
        read_single_line(config_file, line); sscanf(line, "%f ", &start_noise);
        read_single_line(config_file, line); sscanf(line, "%f %f ", &start_delay, &sdev_start_delay);
        read_single_line(config_file, line); sscanf(line, "%f %f ", &r_start_eqh, &r_start_eqv);

	fclose (config_file);

 fpin = fopen(argv[2], "r");

// ini random number generator
 if (true_random<=0) {unsigned int seed; seed=get_seed(); srand(seed);  fprintf(stderr, "random seed %u\n", seed); }
 if (true_random>0)  {unsigned int seed; seed=true_random; srand(seed); fprintf(stderr, "fixed  seed %u\n", seed); }

// define boundaries
 xmin=gh.x0;
 xmax=gh.x0+(gh.nx-1)*gh.h;
 ymin=gh.y0;
 ymax=gh.y0+(gh.ny-1)*gh.h;
 zmin=gh.z0;
 zmax=gh.z0+(gh.nz-1)*gh.h;

 fprintf(stderr,"xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n",xmin, xmax, ymin, ymax, zmin, zmax);

// open data file
  if (!(finp=fopen(argv[3],"r"))) { fprintf (stderr, "could not open data file %s\n", argv[3]); exit (0);}       
  ne = read_mcmcdata(finp, s);
//  fprintf(stderr, "found %d events\n", ne);
  old_model.noq=ne;


// count # of picks

  n_ppicks0=0; for (i=0; i<old_model.noq; i++) n_ppicks0=n_ppicks0+s[i].nobs_p0;
  n_spicks0=0; for (i=0; i<old_model.noq; i++) n_spicks0=n_spicks0+s[i].nobs_s0;
  fprintf(stderr, "class0 %d P picks %d S picks\n",n_ppicks0,n_spicks0); 

  n_ppicks1=0; for (i=0; i<old_model.noq; i++) n_ppicks1=n_ppicks1+s[i].nobs_p1;
  n_spicks1=0; for (i=0; i<old_model.noq; i++) n_spicks1=n_spicks1+s[i].nobs_s1;
  fprintf(stderr, "class1 %d P picks %d S picks\n",n_ppicks1,n_spicks1); 

  n_ppicks2=0; for (i=0; i<old_model.noq; i++) n_ppicks2=n_ppicks2+s[i].nobs_p2;
  n_spicks2=0; for (i=0; i<old_model.noq; i++) n_spicks2=n_spicks2+s[i].nobs_s2;
  fprintf(stderr, "class2 %d P picks %d S picks\n",n_ppicks2,n_spicks2); 
  
  n_ppicks3=0; for (i=0; i<old_model.noq; i++) n_ppicks3=n_ppicks3+s[i].nobs_p3;
  n_spicks3=0; for (i=0; i<old_model.noq; i++) n_spicks3=n_spicks3+s[i].nobs_s3;
  fprintf(stderr, "class3 %d P picks %d S picks\n",n_ppicks3,n_spicks3); 
  
// sum of all pick classes
  n_ppicks=0; for (i=0; i<old_model.noq; i++) n_ppicks=n_ppicks+s[i].nobs_p;
  n_spicks=0; for (i=0; i<old_model.noq; i++) n_spicks=n_spicks+s[i].nobs_s;
  fprintf(stderr, "sum of all picks : %d P picks %d S picks\n",n_ppicks,n_spicks); 

//for (i=0; i<old_model.noq; i++) fprintf(stdout, "EQ = %d nobs_p=%d nobs_s=%d\n",i,s[i].nobs_p,s[i].nobs_s); 


// count # of stations
  nos=0;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].st_id>nos) nos=s[i].s_picks[j].st_id;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].st_id>nos) nos=s[i].p_picks[j].st_id;}
  old_model.nos=nos+1;

  fprintf(stderr, "Number of stations = %ld\n",old_model.nos);
  fprintf(stderr, "Number of quakes   = %ld\n",old_model.noq);

// check if stations inside search boundaries
  xtest=1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].x<xtest) xtest=s[i].p_picks[j].x;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].x<xtest) xtest=s[i].s_picks[j].x;}
  if (xtest<xmin) {fprintf(stderr, "station x position outside search boundaries %f %f\n",xtest,xmin); exit(0);}
  xtest=-1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].x>xtest) xtest=s[i].p_picks[j].x;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].x>xtest) xtest=s[i].s_picks[j].x;}
  if (xtest>xmax) {fprintf(stderr, "station x position outside search boundaries %f %f\n",xtest,xmax); exit(0);}
  xtest=1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].y<xtest) xtest=s[i].p_picks[j].y;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].y<xtest) xtest=s[i].s_picks[j].y;}
  if (xtest<ymin) {fprintf(stderr, "station y position outside search boundaries %f %f\n",xtest,ymin); exit(0);}
  xtest=-1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].y>xtest) xtest=s[i].p_picks[j].y;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].y>xtest) xtest=s[i].s_picks[j].y;}
  if (xtest>ymax) {fprintf(stderr, "station y position outside search boundaries %f %f\n",xtest,ymax); exit(0);}
  zrmin=1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].z<zrmin) zrmin=s[i].p_picks[j].z;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].z<zrmin) zrmin=s[i].s_picks[j].z;}
  if (zrmin<zmin) {fprintf(stderr, "station z position outside search boundaries %f %f\n",zrmin,zmin); exit(0);}
  zrmax=-1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].z>zrmax) zrmax=s[i].p_picks[j].z;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].z>zrmax) zrmax=s[i].s_picks[j].z;}
  if (zrmax>zmax) {fprintf(stderr, "station z position outside search boundaries %f %f\n",zrmax,zmax); exit(0);}

  fprintf(stderr, "found station elevations between %f and %f\n",zrmin,zrmax);

// calculate weights for receiver elevation interpolation
  for (i=0; i<old_model.noq; i++) 
  {
	for (j=0; j<s[i].nobs_p; j++) 
	{
		s[i].p_picks[j].layer=(int)((s[i].p_picks[j].z-gh.z0)/gh.h);
		s[i].p_picks[j].w2=-(s[i].p_picks[j].layer*gh.h+gh.z0-s[i].p_picks[j].z)/gh.h;
		s[i].p_picks[j].w1=1.0-s[i].p_picks[j].w2;
	}
	for (j=0; j<s[i].nobs_s; j++)
	{
		s[i].s_picks[j].layer=(int)((s[i].s_picks[j].z-gh.z0)/gh.h);
		s[i].s_picks[j].w2=-(s[i].s_picks[j].layer*gh.h+gh.z0-s[i].s_picks[j].z)/gh.h;
		s[i].s_picks[j].w1=1.0-s[i].s_picks[j].w2;
	}
  }

// expand nx (diagonal)
  nxmod=(int) sqrt(gh.nx*gh.nx+gh.ny*gh.ny);

/* ***  ALLOCATING TRAVELTIME TABLES  */

// tttpr ttts	
	tttpr = make_3d_array(gh.nz, gh.nz, nxmod);
	tttsr = make_3d_array(gh.nz, gh.nz, nxmod);

// sscanf(argv[3], "%d", &deci);


//-------------------------------------------------------------------//
// read results
  old_model.number=0;
  old_model.dimension=gh.nz;
  old_model.p_noise0=start_noise;
  old_model.s_noise0=start_noise;
// ini ALL station corrections
for (i=0; i<MAX_OBS; i++) old_model.pres[i]=-99999;
for (i=0; i<MAX_OBS; i++) old_model.sres[i]=-99999;

// Vp/Vs STAN   0.000   5.265   0.125   3.013   0.024   5.279   0.058   3.014   0.021   5.280   3.010
  for (i=0; i<gh.nz; i++)
  {
        read_single_line(fpin, line); sscanf(line, "%s %f %f %f %f %f %f %f %f %f %f %f ", ds, &xz, &vpm1, &df, &vsm1, &df, &vpm2, &df, &vsm2, &df, &vpmp, &vsmp); 
	old_model.z[i]=xz;
	old_model.vp[i]=vpm2;
	old_model.vpvs[i]=vsm2;
//fprintf(stderr,"model z=%f vp=%f vpvs=%f\n",old_model.z[i],old_model.vp[i],old_model.vp[i]/old_model.vpvs[i]);
  }


//  for (i=0; i<gh.nz; i++) fprintf(stderr,"model z=%f vp=%f vpvs=%f\n",old_model.z[i],old_model.vp[i],old_model.vpvs[i]);

// quakes EQ    0     7.000    26.162    41.680     0.230     0.246     1.633 1104637590.560  -0.067   0.547
  for (i=0; i<old_model.noq; i++)
  {
        read_single_line(fpin, line); sscanf(line, "%s %d %f %f %f %f %f %f %f %f %f ", ds, &di, &old_model.eq[i].x, &old_model.eq[i].y, &old_model.eq[i].z, &df, &df, &df, &df, &old_model.origin[i], &df); 
	if ((old_model.eq[i].x<xmin) || (old_model.eq[i].x>xmax)) {fprintf(stderr, "quake %d ouside model x %f\n",i,old_model.eq[i].x); exit(0);}
	if ((old_model.eq[i].y<ymin) || (old_model.eq[i].y>ymax)) {fprintf(stderr, "quake %d ouside model y %f\n",i,old_model.eq[i].y); exit(0);}
	if ((old_model.eq[i].z<zmin) || (old_model.eq[i].z>zmax)) {fprintf(stderr, "quake %d ouside model z %f\n",i,old_model.eq[i].z); exit(0);}

  }
//  for (i=0; i<old_model.noq; i++) fprintf(stderr,"EQ %d x=%f y=%f z=%f\n",i,old_model.eq[i].x,old_model.eq[i].y,old_model.eq[i].z);
// skip EZ
  for (i=0; i<old_model.noq; i++) read_single_line(fpin, line);
 
// res RES   29  -0.986   0.024   0.023   0.038 old_model.pres[i]
  for (i=0; i<old_model.nos; i++)
  {
        read_single_line(fpin, line); sscanf(line, "%s %d %f %f %f %f  ", ds, &di, &old_model.pres[i], &old_model.sres[i], &df, &df); 
//	fprintf(stdout,"statcorr=%d %f %f \n",i,old_model.pres[i],old_model.sres[i]);
  }
// noise level
	read_single_line(fpin, line); sscanf(line, "%f %f %f %f %f %f %f %f ", &old_model.p_noise0, &old_model.p_noise1, &old_model.p_noise2,  &old_model.p_noise3,
&old_model.s_noise0, &old_model.s_noise1, &old_model.s_noise2, &old_model.s_noise3); 

 fprintf(stderr,"Start misfit calc\n");


 old_misfit = cal_fit_newx(&old_model, s, ne, tttpr, tttsr, gh, 3, &mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,1);
// for (ix=0; ix<gh.nz; ix++) for (jx=0; jx<nxmod; jx++) {tttp_old[ix][jx]=tttp[ix][jx]; ttts_old[ix][jx]=ttts[ix][jx];}

// fprintf(stderr,"XXX %f %f %d\n",mfp0+mfp1+mfp2+mfp3, mfs0+mfs1+mfs2+mfs3,n_ppicks+n_spicks );
 old_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/(n_ppicks+n_spicks));
 old_ll=-old_misfit/2.0;






 fprintf(stderr,"Start model found with loglikelihood %f RMS=%f\n",old_ll,old_rms);



 fclose(fpin);
 exit(0);
}

/* ============================================= */
int read_mcmcdata (FILE *f, struct DATA *d)
{
	int i, j, ne, end, k;
	char buffer[4096];
	char dum[64], *c;
	int st_id;
	float x, y, z;
	double t;
	char ph[2];
	int np, ns;
	int cl;

	np=0; ns=0;
	i=0; end=1; j=0; k=0;
        while (!feof(f) && end)
        {		
		c = fgets (buffer, 4096, f);
//		fprintf(stderr, "%s\n", buffer);
	        if (!c) 
               		break;
			
		if (strchr(buffer, '#')==NULL)
		{

			memset (ph, '\0', 2);
			sscanf(buffer, "%s %d %s %f %f %f %lf %d\n", 
				dum, &st_id, ph, &x, &y, &z, &t, &cl);
			if (strchr(ph,'P')!=NULL)
			{
				d[i].p_picks[np].st_id = st_id;
				d[i].p_picks[np].x = x;
				d[i].p_picks[np].y = y;
				d[i].p_picks[np].z = z;
				d[i].p_picks[np].t = (float)t;
				d[i].p_picks[np].cl = cl;
				if (cl==0) d[i].nobs_p0=d[i].nobs_p0+1;
				if (cl==1) d[i].nobs_p1=d[i].nobs_p1+1;
				if (cl==2) d[i].nobs_p2=d[i].nobs_p2+1;
				if (cl==3) d[i].nobs_p3=d[i].nobs_p3+1;
				np++;
			}
			else
			{
				d[i].s_picks[ns].st_id = st_id;
				d[i].s_picks[ns].x = x;
				d[i].s_picks[ns].y = y;
				d[i].s_picks[ns].z = z;
				d[i].s_picks[ns].t = (float)t;
				d[i].s_picks[ns].cl = cl;
				if (cl==0) d[i].nobs_s0=d[i].nobs_s0+1;
				if (cl==1) d[i].nobs_s1=d[i].nobs_s1+1;
				if (cl==2) d[i].nobs_s2=d[i].nobs_s2+1;
				if (cl==3) d[i].nobs_s3=d[i].nobs_s3+1;
				ns++;
			}
			k++;
		}
		else
		{	
			if (j!=0) i++;
			d[i].nobs_p0=0; d[i].nobs_p1=0; d[i].nobs_p2=0; d[i].nobs_p3=0;
			d[i].nobs_s0=0; d[i].nobs_s1=0; d[i].nobs_s2=0; d[i].nobs_s3=0;
//			fprintf(stderr, "EVENT %d\n", i);
			sscanf(buffer, "%s %d %d %d %lf", dum, &d[i].eq_id,
					 &d[i].nobs_p, &d[i].nobs_s, &d[i].reftime);
			k=0;
			np=0; ns=0;
		}
		j++;
	}
	ne = i+1;
	
	return(ne);
}

/* -------------------------------------------------------------------- */
float dst (float x1, float x2, float y1, float y2)
{
	return (sqrt ( ( (x1-x2)*(x1-x2) )+ ( (y1-y2)*(y1-y2) ) ) );
}



