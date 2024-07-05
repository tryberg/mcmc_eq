// Program to invert P & S picks for location & structure by MCMC

// don't use TRIA
// no priors, no weights
// fixed P, S & P/S station corrections
// 191023 modification for compiler
// 190624 cleanup, misfit & interpol in common file

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>

#include <string.h>
#include "mc.h"

int TRIA, DR, aflag;
FILE *fpout;

int read_mcmcdata (FILE *f, struct DATA *d);
void write_mcmcdata (struct DATA *d, int ne);

float cal_fit_newx (struct Model *m, struct DATA *d, int ne, float ***tttp, float ***ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3,  int flag, int eikonal, int out); 
float traveltimet (float **ttt, int nx, int ny, int nz, float h, float dist, float z, float z0);
float dst (float x1, float x2, float y1, float y2);

void setup_table_new (struct Model *m, float ***ttt, struct GRDHEAD gh, int ps);
void read_single_line(FILE *fp, char x[]);
void copy_model(struct Model *dest, struct Model *src);
int find_neighbor_cell(struct Model *modx, int n);
int find_in_cell(struct Model *modx, float x);
void model_to_hsbuf(struct Model *modx, struct GRDHEAD gh, float *hsbuf);

#include "interpol.c"
#include "misfit.c"


int check_string(char string[], char x)
{
 char * ptr;
 int i;
 ptr = strchr( string, x );
 i=0;
 if (ptr) i=1;
 return(i);
}

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

float nexp (float v)
{
    float h;
    if (v<log(FLT_MAX/1000.0)) {h=exp(v);} else {h=exp(log(FLT_MAX/1000.0));}
    return (h);
}

int out_of_bounds (float v0, float lower, float upper)
{
  if ((v0>=lower) && (v0<=upper)) return (0); else {return(1);}
}

float rand_gauss_bounded (float v0, float sdev, float lower, float upper)
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

float rand_eq_limited (float min, float max)
// returns random float between min and max, eq distributed 
{
    return (((float) min+(max-min)*rand()/RAND_MAX));
}

int model_valid (struct Model *model, float dz, float zmin, float zmax, float inv_control)
{
  int i,swap,iflag;
  float dzmin,help;
  float zz[MD],bd[MD],th[MD],vp[MD],vs[MD];
  if (model->dimension==1) return (0);
  for (i=0; i<model->dimension; i++) {zz[i]=model->z[i]; vp[i]=model->vp[i];vs[i]=vp[i]/model->vpvs[i]; }
  
// sort
  do
  {
  	swap=0;
  	for (i=1; i<model->dimension; i++)
  	{
		if (zz[i-1]>zz[i])
// swap
		{
			help=zz[i-1]; zz[i-1]=zz[i]; zz[i]=help; swap=1;
			help=vp[i-1]; vp[i-1]=vp[i]; vp[i]=help;
			help=vs[i-1]; vs[i-1]=vs[i]; vs[i]=help;	
		}
  	}
  }
  while (swap==1);

// layer boundaries
  for (i=0; i<model->dimension-1; i++) bd[i]=(zz[i]+zz[i+1])/2.0; 
  bd[model->dimension-1]=zmax;
// layer thickness
  for (i=1; i<model->dimension; i++) th[i]=bd[i]-bd[i-1]; 
  th[0]=bd[0]-zmin;
  dzmin=FLT_MAX;
// test for layer thickness
  for (i=0; i<model->dimension; i++) if (th[i]<dzmin) dzmin=th[i];

// ini model valid
  iflag=0;
// check for minimum layer thickness
  if (dzmin<(sqrt(inv_control*inv_control)*dz)) iflag=1;
  
  
// check for LVZ
  swap=0;
  for (i=0; i<model->dimension-1; i++) if (vp[i]>vp[i+1]) swap=swap+1;
  for (i=0; i<model->dimension-1; i++) if (vs[i]>vs[i+1]) swap=swap+1;
// check for LVZ 
  if (inv_control<0) if (swap>0) iflag=1;

  if (iflag==1) return (1); else {return(0);}
}



// print model raw
void print_model_raw(struct Model *model, struct DATA *d, float rms, char *text, char *deci)
{
 int i;
 fprintf(fpout,"%3s %2s %8ld %3ld %f %f %f %f %f %f %f %f %f",text,deci,model->number,model->dimension,rms,model->p_noise0,model->p_noise1, model->p_noise2, model->p_noise3, model->s_noise0, model->s_noise1, model->s_noise2, model->s_noise3);
 for (i=0; i<model->dimension; i++) fprintf(fpout," %f %f %f",model->z[i],model->vp[i],model->vpvs[i]);
 fprintf(fpout,"\n");


 for (i=0; i<model->noq; i++) fprintf(fpout,"EQ  %2s %8ld %d %f %f %f %f %lf %f\n",deci,model->number,i,rms,model->eq[i].x,model->eq[i].y,model->eq[i].z,d[i].reftime,model->origin[i]);
// for (i=0; i<model->noq; i++) fprintf(fpout,"EQ  %2s %8d %d %f %f %f %f %f\n",deci,model->number,i,model->eq[i].x,model->eq[i].y,model->eq[i].z,model->p_ori[i],model->s_ori[i]);
 for (i=0; i<model->nos; i++) fprintf(fpout,"RES %2s %8ld %d %f %f %f\n",deci,model->number,i,rms,model->pres[i],model->sres[i]);


 fflush(fpout);
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
 int i,j,k,testcase,index,ideath,start_cell_number,sdev_start_cell_number, max_dim;
 double new_ll, old_ll;
 float sdevx, sdevy, sdevz, sdevvp, sdevvpvs, sdevn;
 float xmin, xmax, ymin, ymax, zmin, zmax, vpmin, vpmax, vpvsmin,vpvsmax;
 float dz, dy, dx, dvp, dvs, dns, dnp;
 float start_vp, sdev_start_vp, start_vp_grad, start_noise;
 double log_fac;
 float alpha12;
 float newz;
 float noise_min, noise_max;
 float xrmin, xrmax, yrmin, yrmax, zrmin, zrmax;
 double new_misfit, old_misfit, new_rms, old_rms, best_rms;
 char *decision, ctestcase;
 int acce,reject;

 int deci,j_max_start, j_max_main,true_random,topo_shift,topo_flag;


 int nmod,nsa,nsr,npa,npr,nba,nbr,nda,ndr,nna,nnr,nra,nrr,nma,nmr,nqa,nqr;

 char		    line[1024];
 char   dstring_start[1000],dstring_main[1000],topo_file_name[1000];
 char	pstring_start[100000],pstring_main[100000];
 FILE 		*config_file;		/* config file	*/
// CH
 FILE 		*finp;			/* data file 			*/

 struct GRDHEAD  gh;			/* FD grid header structure 	*/
 struct DATA 	*s;		/*  data 			*/



struct Model old_model, new_model, best_model;

 float mfp0, mfs0, mfp1, mfs1, mfp2, mfs2, mfp3, mfs3;
 float sdevxs, sdevys, sdevzs, sdevresidual;
 float residual_min, residual_max;
 int n_ppicks, n_spicks,n_ppicks0, n_spicks0, n_ppicks1, n_spicks1, n_ppicks2, n_spicks2, n_ppicks3, n_spicks3, iq, ne;


 float 	***tttpr, ***tttsr;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */
 float 	***tttpr_old, ***tttsr_old;		/* travel time (in sec.) table. ttt[source depth km][distance km]  */

 int nxmod,nos, model_not_valid, eikonal, jj, ix, iz, jx;

 float start_delay, sdev_start_delay,r_start_eqh,r_start_eqv;
 float start_vpvs, sdev_start_vpvs,inv_control;
 int ff,reference_station;
 

 float fac, epi_search;

 int mod_from_file, scor_flag, xflag, lvz_flag, revert;
 FILE *fmodel;			/* model file			*/
 
 char sd[100], inp_model_switch[10];
 float df, f1, f2, f3, f4, f5, f6, f7, f8, value, ref_statcor_P, ref_statcor_S;
 int di, sum_of_picks;

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
        read_single_line(config_file, line); sscanf(line, "%f %f ", &sdevxs, &epi_search); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevys); 
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevzs);
        read_single_line(config_file, line); sscanf(line, "%f ", &sdevresidual);
        read_single_line(config_file, line); sscanf(line, "%f ", &inv_control); if (inv_control==0.0) { fprintf(stderr, "inv_control should be != 0!\n"); exit (0);}
	lvz_flag=0; if ( inv_control > 0 ) {lvz_flag=1; inv_control=-inv_control;}
        read_single_line(config_file, line); sscanf(line, "%d %d %f %f ", &reference_station, &scor_flag, &ref_statcor_P, &ref_statcor_S);
        read_single_line(config_file, line); sscanf(line, "%d ", &TRIA);
        read_single_line(config_file, line); sscanf(line, "%d %d ", &j_max_start, &j_max_main);
        read_single_line(config_file, line); sscanf(line, "%d ", &deci);
        read_single_line(config_file, line); sscanf(line, "%d %d", &true_random, &eikonal);
        read_single_line(config_file, line); sscanf(line, "%s %s", dstring_start, dstring_main);
        read_single_line(config_file, line); sscanf(line, "%d %s ", &aflag, inp_model_switch); mod_from_file=0; if (aflag==3) {mod_from_file=1; aflag=0;}
	read_single_line(config_file, line); sscanf(line, "%d %s %d ", &topo_flag, topo_file_name, &topo_shift);
	read_single_line(config_file, line); sscanf(line, "%f %f %f ", &start_vp, &sdev_start_vp, &start_vp_grad); 
        read_single_line(config_file, line); sscanf(line, "%f %f ", &start_vpvs, &sdev_start_vpvs); 
        read_single_line(config_file, line); sscanf(line, "%d %d ", &start_cell_number, &sdev_start_cell_number);
        read_single_line(config_file, line); sscanf(line, "%f ", &start_noise);
        read_single_line(config_file, line); sscanf(line, "%f %f ", &start_delay, &sdev_start_delay);
        read_single_line(config_file, line); sscanf(line, "%f %f ", &r_start_eqh, &r_start_eqv);

	fclose (config_file);
	
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

 fprintf(stderr,"xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f REF station %d\n",xmin, xmax, ymin, ymax, zmin, zmax, reference_station);

// open data file
  if (!(finp=fopen(argv[3],"r"))) { fprintf (stderr, "could not open data file %s\n", argv[3]); exit (0);}       
  ne = read_mcmcdata(finp, s);
//  fprintf(stderr, "found %d events\n", ne);
  old_model.noq=ne;


// count # of picks

  n_ppicks0=0; for (i=0; i<old_model.noq; i++) n_ppicks0=n_ppicks0+s[i].nobs_p0; 
  n_spicks0=0; for (i=0; i<old_model.noq; i++) n_spicks0=n_spicks0+s[i].nobs_s0; 
  fprintf(stderr, "class0 %7d P picks %7d S picks\n",n_ppicks0,n_spicks0); 


  n_ppicks1=0; for (i=0; i<old_model.noq; i++) n_ppicks1=n_ppicks1+s[i].nobs_p1;
  n_spicks1=0; for (i=0; i<old_model.noq; i++) n_spicks1=n_spicks1+s[i].nobs_s1; 
  fprintf(stderr, "class1 %7d P picks %7d S picks\n",n_ppicks1,n_spicks1); 


  n_ppicks2=0; for (i=0; i<old_model.noq; i++) n_ppicks2=n_ppicks2+s[i].nobs_p2;
  n_spicks2=0; for (i=0; i<old_model.noq; i++) n_spicks2=n_spicks2+s[i].nobs_s2;
  fprintf(stderr, "class2 %7d P picks %7d S picks\n",n_ppicks2,n_spicks2); 
  
  n_ppicks3=0; for (i=0; i<old_model.noq; i++) n_ppicks3=n_ppicks3+s[i].nobs_p3;
  n_spicks3=0; for (i=0; i<old_model.noq; i++) n_spicks3=n_spicks3+s[i].nobs_s3; 
  fprintf(stderr, "class3 %7d P picks %7d S picks\n",n_ppicks3,n_spicks3); 


  
  
// sum of all pick classes
  n_ppicks=0; for (i=0; i<old_model.noq; i++) n_ppicks=n_ppicks+s[i].nobs_p;
  n_spicks=0; for (i=0; i<old_model.noq; i++) n_spicks=n_spicks+s[i].nobs_s;
  sum_of_picks=n_ppicks+n_spicks;

  fprintf(stderr, "sum of all picks : %8d P picks %8d S picks\n",n_ppicks,n_spicks); 

//for (i=0; i<old_model.noq; i++) fprintf(stdout, "EQ = %d nobs_p=%d nobs_s=%d\n",i,s[i].nobs_p,s[i].nobs_s); 


// check for largest station number
  nos=0;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].st_id>nos) nos=s[i].s_picks[j].st_id;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].st_id>nos) nos=s[i].p_picks[j].st_id;}
  old_model.nos=nos+1;
  
// check for "missing" stations
  for (k=0; k<old_model.nos; k++)
  {
    xflag=0; 
    for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) if (s[i].s_picks[j].st_id==k) xflag=1;
    for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) if (s[i].p_picks[j].st_id==k) xflag=1;
    if (xflag==0) {fprintf(stdout, "WARNING station %d missing in pick file\n",k);}
  } 

// check for "missing" quakes
  for (i=0; i<old_model.noq; i++) if (s[i].eq_id!=i) {fprintf(stdout, "WARNING quakes not correctly sorted in pick file around EQ_ID %d \n",i);}



  fprintf(stderr, "Number of stations = %ld\n",old_model.nos);
  fprintf(stderr, "Number of quakes   = %ld\n",old_model.noq);
  
// check if reference station in data file
  if (reference_station>=old_model.nos) { fprintf (stderr, "reference station not in data file\n"); exit (0);}

// check if stations inside search boundaries
  xrmin=1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].x<xrmin) xrmin=s[i].p_picks[j].x;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].x<xrmin) xrmin=s[i].s_picks[j].x;}
  if (xrmin<xmin) {fprintf(stderr, "station x position outside search boundaries %f %f\n",xrmin,xmin); exit(0);}
  xrmax=-1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].x>xrmax) xrmax=s[i].p_picks[j].x;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].x>xrmax) xrmax=s[i].s_picks[j].x;}
  if (xrmax>xmax) {fprintf(stderr, "station x position outside search boundaries %f %f\n",xrmax,xmax); exit(0);}
  yrmin=1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].y<yrmin) yrmin=s[i].p_picks[j].y;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].y<yrmin) yrmin=s[i].s_picks[j].y;}
  if (yrmin<ymin) {fprintf(stderr, "station y position outside search boundaries %f %f\n",yrmin,ymin); exit(0);}
  yrmax=-1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].y>yrmax) yrmax=s[i].p_picks[j].y;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].y>yrmax) yrmax=s[i].s_picks[j].y;}
  if (yrmax>ymax) {fprintf(stderr, "station y position outside search boundaries %f %f\n",yrmax,ymax); exit(0);}
  zrmin=1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].z<zrmin) zrmin=s[i].p_picks[j].z;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].z<zrmin) zrmin=s[i].s_picks[j].z;}
  if (zrmin<zmin) {fprintf(stderr, "station z position outside search boundaries %f %f\n",zrmin,zmin); exit(0);}
  zrmax=-1e30;
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_p; j++) {if (s[i].p_picks[j].z>zrmax) zrmax=s[i].p_picks[j].z;}
  for (i=0; i<old_model.noq; i++) for (j=0; j<s[i].nobs_s; j++) {if (s[i].s_picks[j].z>zrmax) zrmax=s[i].s_picks[j].z;}
  if (zrmax>zmax) {fprintf(stderr, "station z position outside search boundaries %f %f\n",zrmax,zmax); exit(0);}

  fprintf(stderr, "found station x coordinate between %f and %f\n",xrmin,xrmax);
  fprintf(stderr, "found station y coordinate between %f and %f\n",yrmin,yrmax);
  fprintf(stderr, "found station elevations   between %f and %f\n",zrmin,zrmax);

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
 fpout = fopen(argv[2], "w");


// ini counts
 nmod=0;
 nsa=0; nsr=0;
 npa=0; npr=0;
 nqa=0; nqr=0;
 nra=0; nrr=0;
 nba=0; nbr=0;
 nda=0; ndr=0;
 nna=0; nnr=0;
 nma=0; nmr=0;
 nqa=0; nqr=0;

//-------------------------------------------------------------------//

// starting model
  old_model.number=0;
  old_model.p_noise0=start_noise;  old_model.s_noise0=start_noise;
  old_model.p_noise1=start_noise;  old_model.s_noise1=start_noise;
  old_model.p_noise2=start_noise;  old_model.s_noise2=start_noise;
  old_model.p_noise3=start_noise;  old_model.s_noise3=start_noise;

  
// try model until minimum layer thickness 
  fprintf(stderr,"Test starting model "); 

  do
 	{
  	 	if (start_cell_number>1)
 	  	{
	  		old_model.dimension=start_cell_number+(int) rand_gauss_bounded((float) start_cell_number, (float) sdev_start_cell_number, 1.0, (float) gh.nz);
	 	} else old_model.dimension=1;

 		if ( TRIA == 0 )
 		{   	
			for (i=0; i<old_model.dimension; i++) old_model.z[i]=rand_eq_limited(zmin,zmax);
   			for (i=0; i<old_model.dimension; i++)
   			{
   				value=start_vp+(old_model.z[i]-gh.z0)*start_vp_grad;
   				old_model.vp[i]=value+rand_gauss_bounded(value, sdev_start_vp, vpmin, vpmax);
    				old_model.vpvs[i]=start_vpvs+rand_gauss_bounded(start_vpvs, sdev_start_vpvs, vpvsmin, vpvsmax);
   			} 
   		}

   		if ( TRIA == 1 )
   		{	
			old_model.dimension=old_model.dimension+2;
			old_model.z[0]=zmin; 
			value=start_vp+(old_model.z[0]-gh.z0)*start_vp_grad;
   			old_model.vp[0]=value+rand_gauss_bounded(value, sdev_start_vp, vpmin, vpmax);
			old_model.vpvs[0]=start_vpvs+rand_gauss_bounded(start_vpvs, sdev_start_vpvs, vpvsmin, vpvsmax);
			old_model.z[1]=zmax; 
			value=start_vp+(old_model.z[1]-gh.z0)*start_vp_grad;
   			old_model.vp[1]=value+rand_gauss_bounded(value, sdev_start_vp, vpmin, vpmax);			
			old_model.vpvs[1]=start_vpvs+rand_gauss_bounded(start_vpvs, sdev_start_vpvs, vpvsmin, vpvsmax);
			
   			for (i=2; i<old_model.dimension; i++) old_model.z[i]=rand_eq_limited(zmin,zmax);
   			for (i=2; i<old_model.dimension; i++)
   			{
				value=start_vp+(old_model.z[i]-gh.z0)*start_vp_grad;
   				old_model.vp[i]=value+rand_gauss_bounded(value, sdev_start_vp, vpmin, vpmax);		
    				old_model.vpvs[i]=start_vpvs+rand_gauss_bounded(start_vpvs, sdev_start_vpvs, vpvsmin, vpvsmax);
   			} 
   		}		
		
// check model if valid
   		j=model_valid(&old_model, gh.h, zmin, zmax, inv_control);

   fprintf(stderr,". "); 
   		if (j==0) {fprintf(stderr,"\nfound starting model dim=%ld\n",old_model.dimension);} else {/* fprintf(stderr,"failed .. %ld\n",old_model.dimension); */}
  	}
  while (j==1);

// set ini EQ positions
  for (iq=0; iq<old_model.noq; iq++) old_model.eq[iq].x=rand_eq_limited(xmin+(xmax-xmin)/2.0*(1.0-r_start_eqh),xmin+(xmax-xmin)/2.0*(1.0+r_start_eqh));
  for (iq=0; iq<old_model.noq; iq++) old_model.eq[iq].y=rand_eq_limited(ymin+(ymax-ymin)/2.0*(1.0-r_start_eqh),ymin+(ymax-ymin)/2.0*(1.0+r_start_eqh));
  for (iq=0; iq<old_model.noq; iq++) old_model.eq[iq].z=rand_eq_limited(zmin,zmax*r_start_eqv);
  for (iq=0; iq<old_model.noq; iq++) if (s[iq].xfix!=-9999) old_model.eq[iq].x=s[iq].xfix;
  for (iq=0; iq<old_model.noq; iq++) if (s[iq].yfix!=-9999) old_model.eq[iq].y=s[iq].yfix;
  for (iq=0; iq<old_model.noq; iq++) if (s[iq].zfix!=-9999) old_model.eq[iq].z=s[iq].zfix;


// for (iq=0; iq<old_model.noq; iq++) fprintf(stderr,"Starting EQ %f %f %f %f %f \n",old_model.eq[iq].x,old_model.eq[iq].y,old_model.eq[iq].z,xmin,xmax);


// ini ALL station corrections
for (i=0; i<MAX_OBS; i++) old_model.pres[i]=-99999;
for (i=0; i<MAX_OBS; i++) old_model.sres[i]=-99999;

// set ini RES & origin times start_delay, &sdev_start_delay
  for (i=0; i<old_model.nos; i++) old_model.pres[i]=start_delay+rand_gauss_bounded(start_delay, sdev_start_delay, residual_min, residual_max);
  for (i=0; i<old_model.nos; i++) old_model.sres[i]=start_delay+rand_gauss_bounded(start_delay, sdev_start_delay, residual_min, residual_max);
  for (i=0; i<old_model.noq; i++) old_model.origin[i]=0.0;
  
// if fixed reference station correction 
  if (scor_flag==1) old_model.pres[reference_station]=ref_statcor_P;
  if (scor_flag==2) old_model.pres[reference_station]=ref_statcor_P;
  if (scor_flag==2) old_model.sres[reference_station]=ref_statcor_S;


// ----------------------- over write station corrections & model --------------------------------------------------------
// if aflag==3 read model/quakes/st_corrections/noise from data file "model.dat" and assign accordingly to inp_model_switch




  if (mod_from_file==1)
  {	
     if (!(fmodel = fopen ("model.dat", "r"))) { fprintf(stderr, "could not open model file\n"); exit (0);}

// Vp Vp/Vs
// example: STAN  -5.000   4.971   1.256   1.651   0.203   5.081   1.233   1.696   0.162   5.935   1.735 0.00000
//            ^      ^                                       ^               ^
     if (check_string(inp_model_switch,'V')==1)
     {
fprintf(stderr,"reading velocity model \n");

        i=0;
 	while (!feof(fmodel))
  	{
        	read_single_line(fmodel, line); sscanf(line, "%s %f %f %f %f %f %f %f %f %f %f %f %f ", &sd[0], &f1, &df, &df, &df, &df, &f2, &df, &f3, &df, &df, &df, &df); 
                if (strcmp(sd, "STAN") == 0) {old_model.z[i]=f1; old_model.vp[i]=f2; old_model.vpvs[i]=f3; i=i+1;}
		if (i>max_dim) {fprintf(stderr, "model larger than reserved space, increase 'max # of cells/layers' in config file\n"); exit (0);}
  	}
	old_model.dimension=i;

	fseek(fmodel, 0, SEEK_SET);
     }

     
// Quakes
// EQ    1   326.616   -63.163    10.448     0.635     0.881     8.944 1482436219.310  -2.440   1.780   0.75318
//  ^    ^      ^         ^         ^                    
     if (check_string(inp_model_switch,'Q')==1)
 
     {  
fprintf(stderr,"reading quakes \n");
      
        i=0;
 	while (!feof(fmodel))
  	{
        	read_single_line(fmodel, line); sscanf(line, "%s %d %f %f %f %f %f %f %f %f %f %f ", &sd[0], &di, &f1, &f2, &f3, &df, &df, &df, &df, &df, &df, &df); 
                if (strcmp(sd, "EQ") == 0) {old_model.eq[di].x=f1; old_model.eq[di].y=f2; old_model.eq[di].z=f3; i=i+1;}
		if (i>MAX_NOQ) {fprintf(stderr, "model larger than reserved space! increase MAX_NOQ in mc.h\n"); exit (0);}
  	}
	if (old_model.noq!=i) {fprintf(stderr, "number of quakes does not fit pick file!\n"); exit (0);} 
	old_model.noq=i;
	fseek(fmodel, 0, SEEK_SET);
     }

// Station corrections
// RES   44   0.119  -0.512   0.222   0.356
//  ^    ^      ^      ^                
     if (check_string(inp_model_switch,'R')==1)
 
     {  
     
fprintf(stderr,"reading station corrections \n");

        i=0;
 	while (!feof(fmodel))
  	{

        	read_single_line(fmodel, line); sscanf(line, "%s %d %f %f %f %f  ", &sd[0], &di, &f1, &f2, &df, &df); 
                if (strcmp(sd, "RES") == 0) {old_model.pres[di]=f1; old_model.sres[di]=f2; i=i+1;}
		if (i>MAX_STAT) {fprintf(stderr, "model larger than reserved space! increase MAX_STAT in mc.h\n"); exit (0);}
  	}
//	if (old_model.nos!=i) {fprintf(stderr, "number of stations does not fit pick file: is %d should be %ld!\n",i,old_model.nos); exit (0);} 
	old_model.nos=i;
	fseek(fmodel, 0, SEEK_SET);
     }

// Noise
// NOISE   0.303   1.794   1.716   2.379   0.448   1.751   2.691   2.058   0.034   1.220   1.390   1.590   0.097   1.358   1.734   1.805
//  ^        ^       ^       ^       ^        ^      ^       ^       ^            
     if (check_string(inp_model_switch,'N')==1)

     {     
fprintf(stderr,"reading noise \n");
        i=0;
 	while (!feof(fmodel))
  	{
        	read_single_line(fmodel, line); sscanf(line, "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ", &sd[0], &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &df, &df, &df, &df, &df, &df, &df, &df); 
                if (strcmp(sd, "NOISE") == 0) 
		{
			old_model.p_noise0=f1;
			old_model.p_noise1=f2;
			old_model.p_noise2=f3;
			old_model.p_noise3=f4;
			old_model.s_noise0=f5;
			old_model.s_noise1=f6;
			old_model.s_noise2=f7;
			old_model.s_noise3=f8;
		}	
  	}			
	fseek(fmodel, 0, SEEK_SET);
     }
     fclose(fmodel);
  }


//-----------------


// ini chain
 clock_t start = clock(), diff;
 old_misfit = cal_fit_newx(&old_model, s, ne, tttpr, tttsr, gh, 3, &mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0);
 diff = clock() - start;  
 int msec = diff * 1000 / CLOCKS_PER_SEC;
 fprintf(stderr, "Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);



 for (iz=0; iz<gh.nz; iz++) for (ix=0; ix<gh.nz; ix++) for (jx=0; jx<nxmod; jx++) {tttpr_old[iz][ix][jx]=tttpr[iz][ix][jx]; tttsr_old[iz][ix][jx]=tttsr[iz][ix][jx];}


 old_misfit=mfp0/old_model.p_noise0/old_model.p_noise0+mfs0/old_model.s_noise0/old_model.s_noise0;
 old_misfit=old_misfit+mfp1/old_model.p_noise1/old_model.p_noise1+mfs1/old_model.s_noise1/old_model.s_noise1;
 old_misfit=old_misfit+mfp2/old_model.p_noise2/old_model.p_noise2+mfs2/old_model.s_noise2/old_model.s_noise2;
 old_misfit=old_misfit+mfp3/old_model.p_noise3/old_model.p_noise3+mfs3/old_model.s_noise3/old_model.s_noise3;

// fprintf(stderr,"XXX %f %f %d\n",mfp0+mfp1+mfp2+mfp3, mfs0+mfs1+mfs2+mfs3,sum_of_picks);
 old_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);
 old_ll=-old_misfit/2.0;


 best_rms=old_rms;
// best_model=old_model; 
 copy_model(&best_model,&old_model);

 print_model_raw(&old_model,s,old_rms, "sta", "ST");

 fprintf(stderr,"Start model found with loglikelihood %f RMS=%f\n",old_ll,old_rms);

 
// balancing pertubation strings 1st phase
 strncat(pstring_start, "", 1);
 for (i=0; i<strlen(dstring_start); i++)
 {
   ctestcase=dstring_start[i];
   switch(ctestcase)
   {
    case 'Q' :
     for (j=0; j<old_model.noq; j=j+10) strncat(pstring_start, "Q", 2);
    break;
    case 'N' :
     strncat(pstring_start, "N", 2);
    break; 
    case 'R' :
     for (j=0; j<old_model.nos; j=j+10) strncat(pstring_start, "R", 2);
    break;
    case 'M' :
     strncat(pstring_start, "M", 2);
    break; 
    case 'V' :
     strncat(pstring_start, "V", 2);
    break;   
    case 'P' :
     strncat(pstring_start, "P", 2);
    break;   
    case 'B' :
     strncat(pstring_start, "B", 2);
    break; 
    case 'D' :
     strncat(pstring_start, "D", 2);
    break;   
   }
 }
 
// balancing pertubation strings 2nd phase
 strncat(pstring_main, "", 1);
 for (i=0; i<strlen(dstring_main); i++)
 {
   ctestcase=dstring_main[i];
   switch(ctestcase)
   {
    case 'Q' :
     for (j=0; j<old_model.noq; j=j+20) strncat(pstring_main, "Q", 2);
    break;
    case 'N' :
     strncat(pstring_main, "N", 2);
    break; 
    case 'R' :
     for (j=0; j<old_model.nos; j=j+20) strncat(pstring_main, "R", 2);
    break;
    case 'M' :
     strncat(pstring_main, "M", 2);
    break; 
    case 'V' :
     strncat(pstring_main, "V", 2);
    break;   
    case 'P' :
     strncat(pstring_main, "P", 2);
    break;   
    case 'B' :
     strncat(pstring_main, "B", 2);
    break; 
    case 'D' :
     strncat(pstring_main, "D", 2);
    break;   
   }
 }
 
 
 fprintf(stderr,"Start %s %ld\n",pstring_start,strlen(pstring_start));
 fprintf(stderr,"Main  %s %ld\n",pstring_main,strlen(pstring_main));

 revert=(int)(j_max_start+j_max_main/2);


 acce=0; reject=0;
// begin of main loop
 while ( acce < (j_max_start+j_max_main))

 {
  j=acce;  
  if ((j==revert) && (lvz_flag==1))
  {
  	inv_control=-inv_control;
	fprintf(stderr,"start to look for LVZs\n");
  }


  for (iz=0; iz<gh.nz; iz++) for (ix=0; ix<gh.nz; ix++) for (jx=0;jx<nxmod; jx++) {tttpr_old[iz][ix][jx]=tttpr[iz][ix][jx]; tttsr_old[iz][ix][jx]=tttsr[iz][ix][jx];}
  copy_model(&new_model,&old_model);
  new_model.number=j;
// 1st phase: search for epicenters accelerated
   if (j<=j_max_start) {testcase=rand_eq_int(strlen(pstring_start));ctestcase=pstring_start[testcase]; fac=epi_search;}
// 2nd phase: search for model as in config file
   if (j>j_max_start) {testcase=rand_eq_int(strlen(pstring_main));ctestcase=pstring_main[testcase]; fac=1.0;}
  model_not_valid=0;
//fprintf(stderr,"test %c\n",ctestcase);

  switch(ctestcase) 
  {
// change EQ location
   case 'Q'  :
    decision="Q.";
    index=rand_eq_int(new_model.noq);
    dx=rand_gauss_bounded(new_model.eq[index].x, sdevxs*fac, xmin, xmax);
    dy=rand_gauss_bounded(new_model.eq[index].y, sdevys*fac, ymin, ymax);
    dz=rand_gauss_bounded(new_model.eq[index].z, sdevzs*fac, zmin, zmax);
    if (s[index].xfix!=-9999) dx=0.0;
    if (s[index].yfix!=-9999) dy=0.0;
    if (s[index].zfix!=-9999) dz=0.0;

    
    new_model.eq[index].x=new_model.eq[index].x+dx;
    new_model.eq[index].y=new_model.eq[index].y+dy;
    new_model.eq[index].z=new_model.eq[index].z+dz;

    new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 0,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag,eikonal,0); nmod++;
    new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
    new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
    new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
    new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
    new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);




    new_ll=-new_misfit/2.0;
    alpha12=min(1.0,nexp(new_ll-old_ll));
   break; 

// change residual statics
   case 'R'  :
    decision="R.";
    index=rand_eq_int(new_model.nos);
    dx=rand_gauss_bounded(new_model.pres[index], sdevresidual, residual_min,residual_max);
    dy=rand_gauss_bounded(new_model.sres[index], sdevresidual, residual_min,residual_max);
    
// no reference station -> zero mean <delay>==0 for P & S
    if (scor_flag==0)
    {
	new_model.pres[index]=new_model.pres[index]+dx;
	new_model.sres[index]=new_model.sres[index]+dy;
	for (jj=0; jj<new_model.nos; jj++) if (jj!=index) new_model.pres[jj]=new_model.pres[jj]-dx/(new_model.nos-1);
    	for (jj=0; jj<new_model.nos; jj++) if (jj!=index) new_model.sres[jj]=new_model.sres[jj]-dy/(new_model.nos-1);
    }
    
// with reference station
    if (scor_flag!=0)
    {
	if (reference_station==index)
	{
		if (scor_flag==1) {dx=0;}		// P only
    		if (scor_flag==2) {dx=0; dy=0;}		// P&S
	}
	new_model.pres[index]=new_model.pres[index]+dx;
    	new_model.sres[index]=new_model.sres[index]+dy;
    }
    
    new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 0,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++;
    new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
    new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
    new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
    new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
    new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);


    new_ll=-new_misfit/2.0;
    alpha12=min(1.0,nexp(new_ll-old_ll));
   break; 

// change Vp value	
   case 'P'  :
    decision="P.";
    do
    {
     copy_model(&new_model,&old_model);
     index=rand_eq_int(new_model.dimension);
     dvp=rand_gauss_bounded(new_model.vp[index], sdevvp, vpmin,vpmax);
     new_model.vp[index]=new_model.vp[index]+dvp;
    } while (model_valid(&new_model, gh.h, zmin, zmax, inv_control)!=0);

    new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 3,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++;
    new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
    new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
    new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
    new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
    new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);

    new_ll=-new_misfit/2.0;
    alpha12=min(1.0,nexp(new_ll-old_ll));
   break; 

// change Vp/Vs value	
   case 'V'  :
    decision="V.";
    do
    {
     copy_model(&new_model,&old_model);
     index=rand_eq_int(new_model.dimension);    
     dvs=rand_gauss_bounded(new_model.vpvs[index], sdevvpvs, vpvsmin,vpvsmax);
     new_model.vpvs[index]=new_model.vpvs[index]+dvs;
    } while (model_valid(&new_model, gh.h, zmin, zmax, inv_control)!=0);

    new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 2,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++;
    new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
    new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
    new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
    new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
    new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);

    new_ll=-new_misfit/2.0;
    alpha12=min(1.0,nexp(new_ll-old_ll));
   break;

// move
   case 'M'  :
    decision="M.";
    ff=1;
    if (( new_model.dimension > 1 ) && ( TRIA == 0 )) ff=0;
    if (( new_model.dimension > 3 ) && ( TRIA == 1 )) ff=0;
    if (ff==0) // if movable get new model
    {
     do
     {
      copy_model(&new_model,&old_model);
      if (TRIA==0) index=rand_eq_int(new_model.dimension);
      if (TRIA==1) index=2+rand_eq_int(new_model.dimension-2);
      if (index<0) {index=0; ff=1; exit(0);} // emergency exit
       dz=rand_gauss_bounded(new_model.z[index], sdevz, zmin, zmax);
      new_model.z[index]=new_model.z[index]+dz;
     } while (model_valid(&new_model, gh.h, zmin, zmax, inv_control)!=0);
    
     new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 3,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++;
     new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
     new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
     new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
     new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
     new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);

     new_ll=-new_misfit/2.0;
     alpha12=min(1.0,nexp(new_ll-old_ll));
//    not movable
    } else {alpha12=0; model_not_valid=1;}
   break; 

// birth	
   case 'B'  :
// check if ndim>ndimmax !!!
    decision="B.";
    if ((new_model.dimension+1) < (max_dim/(1+sqrt(inv_control*inv_control)))) // if birth possible
    {
     do
     {
      copy_model(&new_model,&old_model);
      newz=rand_eq_limited(zmin,zmax);
// find location of new cell -> index
      index=find_in_cell(&new_model,newz);
// Vp
      dvp=rand_gauss_bounded(new_model.vp[index], sdevvp, vpmin,vpmax);
      dvs=rand_gauss_bounded(new_model.vpvs[index], sdevvpvs, vpvsmin,vpvsmax);
// update model dimension
      new_model.dimension=new_model.dimension+1;
      new_model.vp[new_model.dimension-1]=new_model.vp[index]+dvp;
      new_model.vpvs[new_model.dimension-1]=new_model.vpvs[index]+dvs;
      new_model.z[new_model.dimension-1]=newz;
     } while (model_valid(&new_model, gh.h, zmin, zmax, inv_control)!=0);
     log_fac=log(sdevvp*sqrt(2.0*PI)/(vpmax-vpmin))+(new_model.vp[new_model.dimension-1]-new_model.vp[index])*(new_model.vp[new_model.dimension-1]-new_model.vp[index])/2.0/sdevvp/sdevvp;
     if (sdevvpvs!=0) log_fac=log_fac+log(sdevvpvs*sqrt(2.0*PI)/(vpvsmax-vpvsmin))+(new_model.vpvs[new_model.dimension-1]-new_model.vpvs[index])*(new_model.vpvs[new_model.dimension-1]-new_model.vpvs[index])/2.0/sdevvpvs/sdevvpvs;

     new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 3,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++;
     new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
     new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
     new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
     new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
     new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);
     new_ll=-new_misfit/2.0;
     alpha12=min(1.0,nexp(log_fac+new_ll-old_ll));
// birth was impossible because model too large
    } else {alpha12=0.0; model_not_valid=1;}
   break; 

// death	
   case 'D'  :
// check if mod dimension > 1
    decision="D.";
// in case of interpolational models exclude first 2 cells (top & bottom)
    ff=1;
    if (( new_model.dimension > 1 ) && ( TRIA == 0 )) ff=0;
    if (( new_model.dimension > 3 ) && ( TRIA == 1 )) ff=0;
    if (ff==0) // if deletable
    {
     do
     {
      copy_model(&new_model,&old_model);
      if (TRIA==0) index=rand_eq_int(new_model.dimension);
      if (TRIA==1) index=2+rand_eq_int(new_model.dimension-2);
      if (index<0) {index=0; ff=1; exit(0);} // emergency exit
      ideath=index;
      index=find_neighbor_cell(&new_model,ideath);   log_fac=log((vpmax-vpmin)/sdevvp/sqrt(2.0*PI))-(new_model.vp[ideath]-new_model.vp[index])*(new_model.vp[ideath]-new_model.vp[index])/2.0/sdevvp/sdevvp;
      if (sdevvpvs!=0)log_fac=log_fac+log((vpvsmax-vpvsmin)/sdevvpvs/sqrt(2.0*PI))-(new_model.vpvs[ideath]-new_model.vpvs[index])*(new_model.vpvs[ideath]-new_model.vpvs[index])/2.0/sdevvpvs/sdevvpvs;
// update model dimension  remove cell ideath, replace cell value[index]
      new_model.dimension=new_model.dimension-1;
      for (i=ideath; i<new_model.dimension; i++)
      {
       new_model.z[i]=new_model.z[i+1];
       new_model.vp[i]=new_model.vp[i+1];
       new_model.vpvs[i]=new_model.vpvs[i+1];
      }
     } while (model_valid(&new_model, gh.h, zmin, zmax, inv_control)!=0);

     new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 3,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++;
     new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
     new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
     new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
     new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
     new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);
     new_ll=-new_misfit/2.0;
     alpha12=min(1.0,nexp(log_fac+new_ll-old_ll));
// death was impossible because model too small
    } else {alpha12=0.0; model_not_valid=1;}
   break;

// noise pertubation
   case 'N'  :
    decision="N.";
    dnp=rand_gauss_bounded(old_model.p_noise0, sdevn, noise_min, noise_max);
    dns=rand_gauss_bounded(old_model.s_noise0, sdevn, noise_min, noise_max);
    new_model.p_noise0=old_model.p_noise0+dnp;
    new_model.s_noise0=old_model.s_noise0+dns;
    dnp=rand_gauss_bounded(old_model.p_noise1, sdevn, noise_min, noise_max);
    dns=rand_gauss_bounded(old_model.s_noise1, sdevn, noise_min, noise_max);
    new_model.p_noise1=old_model.p_noise1+dnp;
    new_model.s_noise1=old_model.s_noise1+dns;
    dnp=rand_gauss_bounded(old_model.p_noise2, sdevn, noise_min, noise_max);
    dns=rand_gauss_bounded(old_model.s_noise2, sdevn, noise_min, noise_max);
    new_model.p_noise2=old_model.p_noise2+dnp;
    new_model.s_noise2=old_model.s_noise2+dns;
    dnp=rand_gauss_bounded(old_model.p_noise3, sdevn, noise_min, noise_max);
    dns=rand_gauss_bounded(old_model.s_noise3, sdevn, noise_min, noise_max);
    new_model.p_noise3=old_model.p_noise3+dnp;
    new_model.s_noise3=old_model.s_noise3+dns;

    log_fac=n_ppicks0*log(old_model.p_noise0/new_model.p_noise0)+n_spicks0*log(old_model.s_noise0/new_model.s_noise0);
    log_fac=log_fac+n_ppicks1*log(old_model.p_noise1/new_model.p_noise1)+n_spicks1*log(old_model.s_noise1/new_model.s_noise1);
    log_fac=log_fac+n_ppicks2*log(old_model.p_noise2/new_model.p_noise2)+n_spicks2*log(old_model.s_noise2/new_model.s_noise2);
    log_fac=log_fac+n_ppicks3*log(old_model.p_noise3/new_model.p_noise3)+n_spicks3*log(old_model.s_noise3/new_model.s_noise3);
   
    new_misfit=cal_fit_newx(&new_model, s, ne, tttpr, tttsr, gh, 0,&mfp0,&mfs0,&mfp1,&mfs1, &mfp2, &mfs2, &mfp3, &mfs3, aflag, eikonal,0); nmod++; 

    new_misfit=mfp0/new_model.p_noise0/new_model.p_noise0+mfs0/new_model.s_noise0/new_model.s_noise0;
    new_misfit=new_misfit+mfp1/new_model.p_noise1/new_model.p_noise1+mfs1/new_model.s_noise1/new_model.s_noise1;
    new_misfit=new_misfit+mfp2/new_model.p_noise2/new_model.p_noise2+mfs2/new_model.s_noise2/new_model.s_noise2;
    new_misfit=new_misfit+mfp3/new_model.p_noise3/new_model.p_noise3+mfs3/new_model.s_noise3/new_model.s_noise3;
    new_rms=sqrt((mfp0+mfp1+mfp2+mfp3+mfs0+mfs1+mfs2+mfs3)/sum_of_picks);

    new_ll=-new_misfit/2.0;
    alpha12=min(1.0,nexp(log_fac+new_ll-old_ll));
   break;
  }
// end of primary pertubation


// accept if prior output only
  if ( aflag == 1 ) alpha12=1;
// reject if model not valid
  if (( model_not_valid == 1 ) && (aflag == 0)) alpha12=0;
//fprintf(stderr,"XXXX  %8d %2s %4d RMSnew=%16.10f [s] RMSold=%16.10f alpha=%f nll=%f oll=%f\n",j,decision,new_model.dimension,new_rms,old_rms,alpha12,new_ll,old_ll);

// check if accepted
  if (rand_eq()<alpha12) 
  {
   new_model.number=acce;
   acce=acce+1;
   if (ctestcase=='P') npa++;
   if (ctestcase=='V') nsa++;
   if (ctestcase=='Q') nqa++;
   if (ctestcase=='R') nra++;
   if (ctestcase=='M') nma++;
   if (ctestcase=='B') nba++;
   if (ctestcase=='D') nda++;
   if (ctestcase=='N') nna++;

   fprintf(stderr,"Test  %8d %2s %4ld RMS=%16.10f [s] %f %5.1f accepted\n",j,decision,new_model.dimension,new_rms,alpha12,(100.0*acce/(acce+reject)));
//   old_model=new_model;
   copy_model(&old_model,&new_model);
   old_ll=new_ll;
   old_rms=new_rms;
   old_misfit=new_misfit;
// update tt tables
   for (iz=0; iz<gh.nz; iz++) for (ix=0; ix<gh.nz; ix++) for (jx=0; jx<nxmod; jx++) {tttpr_old[iz][ix][jx]=tttpr[iz][ix][jx]; tttsr_old[iz][ix][jx]=tttsr[iz][ix][jx];}
// output
   if (((int)(acce/deci))*deci==acce) print_model_raw(&old_model,s,old_rms, "mod",decision);
  }
  else
// reject model
  {  
  
//   fprintf(stderr,"Test R %8d %2s %4ld RMS=%16.10f [s] %f %5.1f accepted\n",j,decision,new_model.dimension,new_rms,alpha12,(100.0*acce/(acce+reject)));
// restore old tt tables  
   for (iz=0; iz<gh.nz; iz++) for (ix=0; ix<gh.nz; ix++) for (jx=0; jx<nxmod; jx++) {tttpr[iz][ix][jx]=tttpr_old[iz][ix][jx]; tttsr[iz][ix][jx]=tttsr_old[iz][ix][jx];}
// count rejections
   reject=reject+1;
   if (ctestcase=='P') npr++; 
   if (ctestcase=='V') nsr++; 
   if (ctestcase=='Q') nqr++;
   if (ctestcase=='R') nrr++;
   if (ctestcase=='M') nmr++; 
   if (ctestcase=='B') nbr++;
   if (ctestcase=='D') ndr++; 
   if (ctestcase=='N') nnr++;
  }


// test for BAT, best_likelihood
  if ( old_rms < best_rms ) 
  {
//   best_model=old_model;
   copy_model(&best_model,&old_model);
   best_rms=old_rms;
  }
 }

// end of main loop
// output BAT
 print_model_raw(&best_model,s,best_rms, "bat","BF");

// output diagnostics
 fprintf(fpout,"cnt RMS tested   %8d\n",nmod);
 fprintf(fpout,"cnt noise    a/r %8d %8d\n",nna,  nnr);
 fprintf(fpout,"cnt P-vel    a/r %8d %8d\n",npa,  npr);
 fprintf(fpout,"cnt Vp/Vs    a/r %8d %8d\n",nsa,  nsr);
 fprintf(fpout,"cnt quake    a/r %8d %8d\n",nqa,  nqr);
 fprintf(fpout,"cnt resid    a/r %8d %8d\n",nra,  nrr);
 fprintf(fpout,"cnt move     a/r %8d %8d\n",nma,  nmr);
 fprintf(fpout,"cnt birth    a/r %8d %8d\n",nba,  nbr);
 fprintf(fpout,"cnt death    a/r %8d %8d\n",nda,  ndr);

 fclose(fpout);
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
			if (i==MAX_NOQ) {fprintf(stderr, "number of quakes to large! MAX_NOQ in mc.h!\n"); exit(0);}
			memset (ph, '\0', 2);
			sscanf(buffer, "%s %d %s %f %f %f %lf %d\n", 
				dum, &st_id, ph, &x, &y, &z, &t, &cl);
			if (strchr(ph,'P')!=NULL)
			{
				if (np==MAX_OBS) {fprintf(stderr, "number of picks to large! MAX_OBS in mc.h!\n"); exit(0);}
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
				if (cl>3) {fprintf(stderr, "pick class to large! modify source code!\n"); exit(0);}
				np++;
			}
			else
			{
				if (ns==MAX_OBS) {fprintf(stderr, "number of picks to large! MAX_OBS in mc.h!\n"); exit(0);}
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
				if (cl>3) {fprintf(stderr, "pick class to large! modify source code!\n"); exit(0);}
				ns++;
			}
			k++;
		}
		else
		{	
			if (j!=0) i++;
			d[i].nobs_p0=0; d[i].nobs_p1=0; d[i].nobs_p2=0; d[i].nobs_p3=0;
			d[i].nobs_s0=0; d[i].nobs_s1=0; d[i].nobs_s2=0; d[i].nobs_s3=0;
			d[i].xfix=-9999.0;
			d[i].yfix=-9999.0;
			d[i].zfix=-9999.0;

//			fprintf(stderr, "EVENT %d\n", i);
			sscanf(buffer, "%s %d %d %d %lf %lf %lf %lf", dum, &d[i].eq_id,
					 &d[i].nobs_p, &d[i].nobs_s, &d[i].reftime, &d[i].xfix, &d[i].yfix, &d[i].zfix);
//			fprintf(stderr, "EVENT %d eqx = %f eqy = %f eqz = %lf\n", i, d[i].xfix, d[i].yfix, d[i].zfix);

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



