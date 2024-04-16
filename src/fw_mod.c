/*
     Program to calculate FD traveltimes for shot data (2D)
     
     model is defined on arbitrary points 
     performs delauney triangulation of these points (meshing)
     and generates FD grid from the mesh
     
     Written by Christian Haberland, GFZ Potsdam, July/August 2016
     
     uses Triangle routine by J.R Shewchuk and 
     time_2d routine by Podvin & Lecomte
     
     compile:
     cc -DTRILIBRARY -O -c triangle.c 
     gcc -O -c time_2d.c -o time_2d.o ; gcc -O -c mod_grd.c -o mod_grd.o ;  gcc -O mcmc_2d_tomo.c time_2d.o triangle.o mod_grd.o -o mcmc_2d_tomo -lm 


*/
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
float cal_fit_newx (struct Model *m, struct DATA *d, int ne, float ***tttp, float ***ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3,  int flag, int eikonal); 
float cal_fit_new (struct Model *m, struct DATA *d, int ne, float **tttp, float **ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3, int flag, int eikonal); 
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
//int gettimeofday(struct timeval *tv, struct timezone *tz);

// Method to allocate a 2D array of floats
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
    int i,j;

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
 int i,j,k,testcase,index,ideath,start_cell_number,sdev_start_cell_number, topo_shift, max_dim;
 double new_ll, new2_ll, old_ll;
 float sdevx, sdevy, sdevz, sdevvp, sdevvs, sdevn;
 float xmin, xmax, ymin, ymax, zmin, zmax, vpmin, vpmax, vsmin,vsmax;
 float dz, dy, dx, dvp, dvs, dns, dnp;
 float xx,yy,zz,start_vp, sdev_start_vp, start_noise;
 double log_fac;
 float alpha12,alpha13, alpha32, dr_fac, l2, q1;
 float newx, newz;
 float noise_min, noise_max, dn;
 double new_misfit, new2_misfit, old_misfit, new_rms, new2_rms, old_rms, best_rms;
 char *decision, ctestcase;
 int n_picks,no, acce;
 float x0,x1,x2,z0,z1,z2, dxnew, dynew, dznew, dvnew, sum;
 int deci,j_max_start, j_max_main,true_random,topo_flag;
 struct timeval tv;
// struct timezone tz;
 int nmod,nsa,nsr,npa,npr,nba,nbr,nda,ndr,nna,nnr,nra,nrr,nma,nmr,nqa,nqr;

 char		    line[1024], ds[1000];
 char   dstring_start[1000],dstring_main[1000], topo_file_name[1000];
 FILE 		*config_file;		/* config file	*/
 FILE 		*topo_file;		/* topo file	*/
// CH
 FILE 		*finp;			/* data file 			*/

 struct GRDHEAD  gh;			/* FD grid header structure 	*/
 struct DATA 	*s;		/*  data 			*/
 int 	nshot, err=0;
 int    *topo;
 int	fl_s_n, fl_s_v, fl_s_c;
struct Model old_model, new_model, best_model;
 float mfp0, mfs0, mfp1, mfs1, mfp2, mfs2, mfp3, mfs3;
 float xs, ys, zs, dt;
 float sdevxs, sdevys, sdevzs, sdevresidual;
 float residual_min, residual_max;

 int n_ppicks,n_spicks, n_ppicks0, n_spicks0, n_ppicks1, n_spicks1, n_ppicks2, n_spicks2, n_ppicks3, n_spicks3, iq, ne;


 float 	**tttp, **ttts;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */
 float 	**tttp_old, **ttts_old;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */
 float 	***tttpr, ***tttsr;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */
 float 	***tttpr_old, ***tttsr_old;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */

 int calct,nxmod,nos, model_not_valid, eikonal, jj, ix, jx;

 float start_delay, sdev_start_delay,r_start_eqh,r_start_eqv;
 float xz,vpm1, vpm2, vpmp, vsm1, vsm2, vsmp, df;
 int di,mod_from_file;
 float inv_control;

 float vpvsmin, vpvsmax, sdevvpvs, start_vpvs, sdev_start_vpvs;
 float zrmin, zrmax, xtest;
 long tx;
 char *buf;
 char *c,*c2;
 int l;
 
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
	tttpr_old = make_3d_array(gh.nz, gh.nz, nxmod);
	tttsr_old = make_3d_array(gh.nz, gh.nz, nxmod);

// sscanf(argv[3], "%d", &deci);
// read results  
  if (!(buf=malloc(5000*3*20*sizeof(char))))
  {
  	fprintf(stderr, "malloc failed for buf\n");
	exit(0);
  }
      
      
//  fprintf(stderr, "XXX open file %s\n", argv[2]);

  if (!(finp=fopen(argv[2],"r"))) { fprintf (stderr, "could not open model file %s\n", argv[2]); exit (0);}       
  if (fgets(buf, 5000*3*20*sizeof(char), finp)!=NULL) 
  {
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      old_model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      old_model.dimension = atoi(strtok (NULL, " "));/* read dimension */	
      old_rms = atof(strtok (NULL, " "));/* read rms */	
      old_model.p_noise0 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.p_noise1 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.p_noise2 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.p_noise3 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.s_noise0 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.s_noise1 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.s_noise2 = atof(strtok (NULL, " "));/* read model noise */ 
      old_model.s_noise3 = atof(strtok (NULL, " "));/* read model noise */ 
      l=0; 
      while (l<old_model.dimension)
      {
       old_model.z[l] = atof(strtok (NULL, " "));
       old_model.vp[l] = atof(strtok (NULL, " ")); 
       old_model.vpvs[l] = atof(strtok (NULL, " ")); 
//       fprintf(stderr, "XXX evaluate model %ld at z %f with vp %f vs %f\n", old_model.number,old_model.z[l],old_model.vp[l],old_model.vpvs[l]);
       l++; 
      }

  }
  
//-------------------------------------------------------------------//
// quakes
  for (i=0; i<old_model.noq; i++)
  {
        read_single_line(finp, line); sscanf(line, "%s %s %d %d %f %f %f %f %f %f ", ds, ds, &di, &di, &df, &old_model.eq[i].x, &old_model.eq[i].y, &old_model.eq[i].z, &old_model.origin[i], &df); 
	if ((old_model.eq[i].x<xmin) || (old_model.eq[i].x>xmax)) {fprintf(stderr, "quake %d ouside model x %f\n",i,old_model.eq[i].x); exit(0);}
	if ((old_model.eq[i].y<ymin) || (old_model.eq[i].y>ymax)) {fprintf(stderr, "quake %d ouside model y %f\n",i,old_model.eq[i].y); exit(0);}
	if ((old_model.eq[i].z<zmin) || (old_model.eq[i].z>zmax)) {fprintf(stderr, "quake %d ouside model z %f\n",i,old_model.eq[i].z); exit(0);}
//	fprintf(stderr, "XXX EQ %d at x %f y %f z %f t %f\n", i,old_model.eq[i].x, old_model.eq[i].y, old_model.eq[i].z,old_model.origin[i]);
  }

// station corrections  
  for (i=0; i<old_model.nos; i++)
  {
  	 read_single_line(finp, line); sscanf(line, "%s %s %d %d %f %f %f ", ds, ds, &di, &di, &df, &old_model.pres[i],&old_model.sres[i]); 
//	 fprintf(stderr, "XXX SC %d P %f S %f\n", i,old_model.pres[i], old_model.sres[i]);

  }
  
 fprintf(stderr,"Start misfit calc\n");

 old_misfit = cal_fit_newx(&old_model, s, ne, tttpr, tttsr, gh, 3, &mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal);
// old_misfit = cal_fit_new(&old_model, s, ne, tttp, ttts, gh, 3, &mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal);
// for (ix=0; ix<gh.nz; ix++) for (jx=0; jx<nxmod; jx++) {tttp_old[ix][jx]=tttp[ix][jx]; ttts_old[ix][jx]=ttts[ix][jx];}

// fprintf(stderr,"XXX %f %f %d\n",mfp0+mfp1+mfp2+mfp3, mfs0+mfs1+mfs2+mfs3,n_ppicks+n_spicks );
 old_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/(n_ppicks+n_spicks));
 old_ll=-old_misfit/2.0;






 fprintf(stderr,"Start model found with loglikelihood %f RMS=%f\n",old_ll,old_rms);



 fclose(fpin);
 exit(0);
}



//-----------------------------------------------------
// CCCCCCCCCCCCCCCCCCCCCCCCC
// HHHHHHHHHHHHHHHHHHHHHHHHH

/* ============================================= */
float cal_fit_newx (struct Model *m, struct DATA *d, int ne, float ***tttp, float ***ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3, int flag, int eikonal) 
{
	int i, j, k;
	float tp, ts, dist;

        float txp[MAX_OBS],txs[MAX_OBS],sum,pp[MAX_OBS],ss[MAX_OBS];
	float rmsp0, rmss0, rmsp1, rmss1, rmsp2, rmss2, rmsp3, rmss3; /* rms value (data fit) for this model */
	float station_correction;

	float xs, zs, zr, xr, w1, w2, tsoll;
	int i1, i2;

	rmsp0=0.0; rmss0=0.0;
	rmsp1=0.0; rmss1=0.0;
	rmsp2=0.0; rmss2=0.0;
	rmsp3=0.0; rmss3=0.0;


	if (flag==1) {*mfp0=rmsp0; *mfs0=rmss0; *mfp1=rmsp1; *mfs1=rmss1; *mfp2=rmsp2; *mfs2=rmss2; *mfp3=rmsp3; *mfs3=rmss3; return(1.0);}
	if ( eikonal == 1 )
	{
		if (calct==1)			/* calculate p table */
		{
			setup_table_new(m, tttp, gh, 1);
		}
		if (calct==2)			/* calculate s table */	
		{
			setup_table_new(m, ttts, gh, 2);
		}
		if (calct==3)			/* calculate p&s table */	
		{
			setup_table_new(m, tttp, gh, 1);
			setup_table_new(m, ttts, gh, 2);
		}
	}

	for (i=0; i<m->noq; i++)
	{

// P travel times
		for (j=0; j<d[i].nobs_p; j++)
		{
			dist = dst(d[i].p_picks[j].x, m->eq[i].x, d[i].p_picks[j].y, m->eq[i].y);
			if (eikonal==0) tp = sqrt(dist*dist+m->eq[i].z*m->eq[i].z)/m->vp[find_in_cell(m,0.0)];
 			if (eikonal==1) tp = traveltimet(tttp[d[i].p_picks[j].layer], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].p_picks[j].w1+traveltimet(tttp[d[i].p_picks[j].layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].p_picks[j].w2;
			station_correction=m->pres[d[i].p_picks[j].st_id];
//fprintf(stderr,"XXX P %d %f\n",d[i].p_picks[j].st_id,station_correction);
			if (station_correction<-1000) {fprintf(stderr, "ERROR points to invalid P station correction\n"); exit(0);}
			tp+=station_correction;	/* add residual static correction */
//			fprintf(stderr, "EQ=%d station %03d P dist: %f tpcalc: %f tpobs: %f vstat: %f\n", i, d[i].p_picks[j].st_id, dist, tp, d[i].p_picks[j].t,vpstat);
			txp[j]=tp;
			pp[j]=tp;

		}
// diff
//		for (k=0; k<d[i].nobs_p; k++) fprintf(stderr, "XXX P %f  %f %f\n", txp[k],d[i].p_picks[k].t,txp[k]-d[i].p_picks[k].t);
		for (k=0; k<d[i].nobs_p; k++) txp[k]=txp[k]-d[i].p_picks[k].t;

// sum
		sum=0;
		for (k=0; k<d[i].nobs_p; k++) sum=sum+txp[k];

// S travel times
		for (j=0; j<d[i].nobs_s; j++)
		{
			dist = dst(d[i].s_picks[j].x, m->eq[i].x, d[i].s_picks[j].y, m->eq[i].y);
			if (eikonal==0) ts=sqrt(dist*dist+m->eq[i].z*m->eq[i].z)/(m->vp[find_in_cell(m,0.0)]/m->vpvs[find_in_cell(m,0.0)]);
 			if (eikonal==1) ts = traveltimet (ttts[d[i].s_picks[j].layer], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].s_picks[j].w1+traveltimet(ttts[d[i].s_picks[j].layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].s_picks[j].w2;

			station_correction=m->sres[d[i].s_picks[j].st_id];
			if (station_correction<-1000) {fprintf(stderr, "ERROR points to invalid S station correction\n"); exit(0);}
			ts+=station_correction;	/* add residual static correction */
//			fprintf(stderr, "EQ=%d station %03d S dist: %f tscalc: %f tsobs: %f vstat: %f\n", i, d[i].s_picks[j].st_id, dist, ts, d[i].s_picks[j].t,vsstat);
			txs[j]=ts;
			ss[j]=ts;
		}
// diff
//		for (k=0; k<d[i].nobs_s; k++) fprintf(stderr, "XXX S %f  %f %f \n", txs[k],d[i].s_picks[k].t,txs[k]-d[i].s_picks[k].t);

		for (k=0; k<d[i].nobs_s; k++) txs[k]=txs[k]-d[i].s_picks[k].t;
		for (k=0; k<d[i].nobs_s; k++) sum=sum+txs[k];

		sum=sum/(d[i].nobs_s+d[i].nobs_p);

// origin time
		m->origin[i]=-sum;
		fprintf(stdout, "EVENT %d  %lf %f %f %f %f\n", i, d[i].reftime, m->eq[i].x, m->eq[i].y, m->eq[i].z, m->origin[i]);


//fprintf(stderr, "XXX picktime=%f elstat=%f vp=%f\n", d[i].p_picks[0].t,d[i].p_picks[0].z/vpstat,vpstat);

// de-mean
		for (k=0; k<d[i].nobs_p; k++) txp[k]=txp[k]-sum;
		for (k=0; k<d[i].nobs_s; k++) txs[k]=txs[k]-sum;

// output			
//		fprintf(stdout, "EQ %d\n", i);
		for (k=0; k<d[i].nobs_p; k++)
		{
			dist = dst(d[i].p_picks[k].x, m->eq[i].x, d[i].p_picks[k].y, m->eq[i].y);
			fprintf(stdout, "%f %f %f %f %f %f P\n", txp[k],dist,m->eq[i].z,m->origin[i],d[i].p_picks[k].t,pp[k]);
		}

		for (k=0; k<d[i].nobs_s; k++)
		{	
			dist = dst(d[i].s_picks[k].x, m->eq[i].x, d[i].s_picks[k].y, m->eq[i].y);
			fprintf(stdout, "%f %f %f %f %f %f S\n", txs[k],dist,m->eq[i].z,m->origin[i],d[i].s_picks[k].t,ss[k]);
		}
// rms
		
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==0) rmsp0=rmsp0+txp[k]*txp[k];
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==0) rmss0=rmss0+txs[k]*txs[k];
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==1) rmsp1=rmsp1+txp[k]*txp[k];
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==1) rmss1=rmss1+txs[k]*txs[k];
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==2) rmsp2=rmsp2+txp[k]*txp[k];
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==2) rmss2=rmss2+txs[k]*txs[k];
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==3) rmsp3=rmsp3+txp[k]*txp[k];
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==3) rmss3=rmss3+txs[k]*txs[k];
			
	}
	*mfp0=rmsp0; *mfs0=rmss0;
	*mfp1=rmsp1; *mfs1=rmss1;
	*mfp2=rmsp2; *mfs2=rmss2;
	*mfp3=rmsp3; *mfs3=rmss3; 
//fprintf(stderr, "XXX rmsp0 %f rmsp1 %f rmsp2 %f rmss0 %f rmss1 %f rmss2 %f\n",rmsp0, rmsp1, rmsp2, rmss0, rmss1, rmss2);

	return(1.0);
}

/* -------------------------------------------------------------------- */

void setup_table_new (struct Model *m, float ***ttt, struct GRDHEAD gh, int ps)
{
	int i, j, k;
	float z;
	int ix, iz;
	float *hsbuf, *tbuf;
	float v,xs,zs;
	int p,nxmod;

	struct Model temp_mod;
	float xx;
	float a,b;

	float *vp, *vs;

	nxmod=(int) sqrt(gh.nx*gh.nx+gh.ny*gh.ny);

	if (!(hsbuf=malloc(nxmod * gh.nz * sizeof(float))))
	{
		fprintf(stderr, "malloc failed\n");
		exit(0);
	}
	if (!(tbuf=malloc(nxmod * gh.nz * sizeof(float))))
	{
		fprintf(stderr, "malloc failed\n");
		exit(0);
	}

	if (!(vp=malloc(gh.nz * sizeof(float))))
	{
		fprintf(stderr, "malloc failed\n");
		exit(0);
	}
	if (!(vs=malloc(gh.nz * sizeof(float))))
	{
		fprintf(stderr, "malloc failed\n");
		exit(0);
	}

// ---- set-up v(z) table Voronoi
	if ( TRIA == 0 )
	{
		for (iz=0; iz<gh.nz; iz++)		/* depth */
		{
			z = gh.z0+(float)iz * gh.h;
			k=find_in_cell(m,z);   
			vp[iz]= m->vp[k];
//vp[iz]=6.0;
			vs[iz]= vp[iz]/m->vpvs[k];
		}
	}

// ---- set-up v(z) table tria
	if ( TRIA == 1 )
	{
// sort   
		copy_model(&temp_mod,m);

		j=1;	
		while ( j != 0 )
		{
			j=0;
			for (i=0; i<temp_mod.dimension-1; i++)
			{
				if (temp_mod.z[i]>temp_mod.z[i+1])
				{
					j=1;
					xx=temp_mod.z[i+1]; temp_mod.z[i+1]=temp_mod.z[i]; temp_mod.z[i]=xx;
					xx=temp_mod.vp[i+1]; temp_mod.vp[i+1]=temp_mod.vp[i]; temp_mod.vp[i]=xx;
					xx=temp_mod.vpvs[i+1]; temp_mod.vpvs[i+1]=temp_mod.vpvs[i]; temp_mod.vpvs[i]=xx;
				}
			}
   		}

// interpol
		for (iz=0; iz<gh.nz; iz++)		/* depth */
		{
			z = gh.z0+(float)iz * gh.h;

			for (i=0; i<m->dimension-1; i++) if ((z>=temp_mod.z[i]) && (z<temp_mod.z[i+1])) k=i;
//vp
			a=(temp_mod.vp[k+1]-temp_mod.vp[k])/(temp_mod.z[k+1]-temp_mod.z[k]);
			b=temp_mod.vp[k]-a*temp_mod.z[k];
			vp[iz]= a*z+b;
//vs
			a=(temp_mod.vp[k+1]/temp_mod.vpvs[k+1]-temp_mod.vp[k]/temp_mod.vpvs[k])/(temp_mod.z[k+1]-temp_mod.z[k]);
			b=temp_mod.vp[k]/temp_mod.vpvs[k]-a*temp_mod.z[k];
			vs[iz]= a*z+b;
		}
	}

// prep model
	j=0;
	for (ix=0; ix<nxmod; ix++)			/* epi-distance */
	{
		for (iz=0; iz<gh.nz; iz++)		/* depth */
		{
			if (ps==1) v= vp[iz];
			if (ps==2) v= vs[iz];
			hsbuf[j] = gh.h/v;
			j++;
		}
	}


/* ---- calculate travel times and set-up table */
//	for (i=0; i<nxmod*gh.nz; i++) fprintf(stderr, "%d %f\n", i, hsbuf[i]);
//fprintf(stderr, "start tt calc\n");
	for (iz=0; iz<gh.nz; iz++)
	{

		for (i=0; i<nxmod*gh.nz; i++) tbuf[i] = 0.;
		xs=0.0;
		zs=(float)iz;	
		z=zs*gh.h;

//		fprintf(stderr, "working on source at depth z= %f %f with nxmod=%d nz=%d\n", z, zs, nxmod, gh.nz);	
		time_2d(hsbuf, tbuf, nxmod, gh.nz, xs, zs, 0.001, 0);

// all receiver elevations
		for (j=0; j<gh.nz; j++)
		{
//fprintf(stderr, "sorting rec elev %d\n", j);	
			for (i=0; i<nxmod; i++)
			{
				p = i * gh.nz+j;  
				ttt[j][iz][i] = tbuf[p];
			}
		}
	}	
	free(hsbuf); free(tbuf);
	free(vp); free(vs);
	return;
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



/* -------------------------------------------------------------------- */
float traveltimet (float **ttt, int nx, int ny, int nz, float h, float dist, float z, float z0)
{
	int m1, iz1, m2, iz2;
	float x1, y1, x2, y2;
	float v, v1, v2, v3, v4;

	float x,y ;
	int nxmod;

	nxmod=(int) sqrt(nx*nx+ny*ny);



	m1  = (int)(dist / h);		
	iz1 = (int)((z-z0) / h);



	x = dist;
	y = z-z0;
	
	if (m1>=nxmod-1 || iz1 >= nz-1)
		return (1e30);
		
	m2 = m1 + 1;
	iz2 = iz1 + 1;
// bilinar interpolation wikipedia https://en.wikipedia.org/wiki/Bilinear_interpolation
	x1 = (float)m1 * h;
	y1 = (float)iz1 * h;
	x2 = (float)m2 * h;
	y2 = (float)iz2 * h;
	
	v1 = ttt[iz1][m1];
	v2 = ttt[iz1][m2];
	v3 = ttt[iz2][m1];
	v4 = ttt[iz2][m2];
//fprintf(stdout, "%f %f %d %d %d %d %f %f %f %f %f %f %f %f\n", x,y,m1,m2,iz1,iz2, v1, v2, v3, v4, x1, x2, y1, y2);
	v=1.0/(x2-x1)/(y2-y1)*(v1*(x2-x)*(y2-y)+v2*(x-x1)*(y2-y)+v3*(x2-x)*(y-y1)+v4*(x-x1)*(y-y1));

	return (v);
}



