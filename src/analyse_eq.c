//   Program to analyse output from mcmc_eq 

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

#define NRANSI
int TRIA, DR, aflag;


#define MB 20000			// max # of bins
#define Max_data 210000L	// max number of data points/number of models

FILE *finp;


int read_mcmcdata (FILE *f, struct DATA *d);
void write_mcmcdata (struct DATA *d, int ne);

float cal_fit_new (struct Model *m, struct DATA *d, int ne, float **tttp, float **ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, int flag, int eikonal); 
float traveltimet (float **ttt, int nx, int ny, int nz, float h, float dist, float z, float z0);
float dst (float x1, float x2, float y1, float y2);
void setup_table (struct Model *m, float **ttt, struct GRDHEAD gh, int ps);
int output_ttt (float **ttt, int nx, int nz, float h);
void read_single_line(FILE *fp, char x[]);
void copy_model(struct Model *dest, struct Model *src);
int find_neighbor_cell(struct Model *modx, int n);
int find_in_cell(struct Model *modx, float x);
//int gettimeofday(struct timeval *tv, struct timezone *tz);



void sortn(unsigned long n, float data[])
{
	long i,flag;
	float tmp;
	
	do
	{
		flag=0;
		for (i=1; i<n; i++)
		{
			if (data[i-1]>data[i])
			{
				tmp=data[i-1]; data[i-1]=data[i]; data[i]=tmp; flag=1;
			}
		}
	}
	while (flag==1);

}


long calc_cdf(unsigned long n, float data[], float a[], float b[])
{
	float inc;
	long i,j;

	sortn(n,data);
	inc=1.0/n;
	a[0]=data[0]; b[0]=inc;
	j=0;
	for (i=1; i<n; i++)
	{
		if (data[i]==data[i-1]) b[j]=b[j]+inc;
		else
		{
			j=j+1;
			a[j]=data[i];
			b[j]=b[j-1]+inc;
		}
	}
	j=j+1;
	return(j);
}
#define  pi 3.14159265

float Phi(float x )
{
	float y;
	y=0.5*(1.0+erf(x/sqrt(2.0)));
	return(y);
}

float gauss_cdf(float x, float a, float b, float m, float s)
{
	float phia, phib, phix, y;

	phia=Phi((a-m)/s);
	phib=Phi((b-m)/s);
	phix=Phi((x-m)/s);
	if (x<=a) y=0;
	if ((x>a) && (x<b)) y=(phix-phia)/(phib-phia);
	if (x>=b) y=1;
// test if NaN
	if ( y != y ) y=1.0E20;
	return(y);
}

float misfit(long n, float datax[], float datay[], float lb, float ub, float m, float s)
{
	float rms,z;
	long i;

// out of range
	if (s<=0) return(1.0E20);
//	if (m<=0) return(1.0E20);

	rms=0;
	for (i=0; i<n; i++)
	{
		z=gauss_cdf(datax[i], lb, ub, m, s)-datay[i];
		rms=rms+z*z;
	}
	
	return(rms);
}

void gsearch(float data[Max_data], long nsamp, float *mean, float *sdev, float dx, float *misfit1, float *misfit2)
// call as gsearch(data,nrp, &m, &s, 0.01);
{
	float	datax[Max_data], datay[Max_data];
	long	ncdf, i, flag;
	float	m,s;
	float	avr0, sdev0, udb, ldb;
	float	xx,mftc,mms0,mps0,m0sp,m0sm,mmsm,mmsp,mpsm,mpsp;

// get avr/sdev classical
	avr0=0;	for (i=0; i<nsamp; i++) avr0=avr0+data[i]; avr0=avr0/nsamp;
	sdev0=0; for (i=0; i<nsamp; i++) sdev0=sdev0+(data[i]-avr0)*(data[i]-avr0); sdev0=sqrt(sdev0/nsamp);
	
// convert to CDF
	ncdf=calc_cdf(nsamp,data, datax, datay);

// get upper/lower bounds from data
	
	ldb=datax[0];
	udb=datax[ncdf-1];

// starting at classical position on grid	
	m=dx*(int)(avr0/dx);
	s=dx*(int)(sdev0/dx);	
	mftc=misfit(ncdf, datax, datay, ldb,udb,avr0,sdev0);	
	*misfit1=sqrt(mftc);
	
// search for (local!) minimum on grid

	mftc=misfit(ncdf, datax, datay, ldb,udb,m,s);
	i=0;
	do
	{
		mps0=misfit(ncdf, datax, datay, ldb,udb,m+dx,s);
		mms0=misfit(ncdf, datax, datay, ldb,udb,m-dx,s);
		m0sp=misfit(ncdf, datax, datay, ldb,udb,m,s+dx);
		m0sm=misfit(ncdf, datax, datay, ldb,udb,m,s-dx);

		mmsm=misfit(ncdf, datax, datay, ldb,udb,m-dx,s-dx);
		mpsm=misfit(ncdf, datax, datay, ldb,udb,m+dx,s-dx);
		mpsp=misfit(ncdf, datax, datay, ldb,udb,m+dx,s+dx);
		mmsp=misfit(ncdf, datax, datay, ldb,udb,m-dx,s+dx);				
			
		flag=1;
			   
		xx=mps0; if ((xx<mftc)              && (xx<mms0) && (xx<m0sp) && (xx<m0sm) && (xx<mmsm) && (xx<mpsm) && (xx<mpsp) && (xx<mmsp)) {flag=0; m=m+dx; mftc=xx;}
		xx=mms0; if ((xx<mftc) && (xx<mps0)              && (xx<m0sp) && (xx<m0sm) && (xx<mmsm) && (xx<mpsm) && (xx<mpsp) && (xx<mmsp)) {flag=0; m=m-dx; mftc=xx;}	
		xx=m0sp; if ((xx<mftc) && (xx<mps0) && (xx<mms0)              && (xx<m0sm) && (xx<mmsm) && (xx<mpsm) && (xx<mpsp) && (xx<mmsp)) {flag=0; s=s+dx; mftc=xx;}
		xx=m0sm; if ((xx<mftc) && (xx<mps0) && (xx<mms0) && (xx<m0sp)              && (xx<mmsm) && (xx<mpsm) && (xx<mpsp) && (xx<mmsp)) {flag=0; s=s-dx; mftc=xx;}
		xx=mmsm; if ((xx<mftc) && (xx<mps0) && (xx<mms0) && (xx<m0sp) && (xx<m0sm)              && (xx<mpsm) && (xx<mpsp) && (xx<mmsp)) {flag=0; m=m-dx; s=s-dx; mftc=xx;}
		xx=mpsm; if ((xx<mftc) && (xx<mps0) && (xx<mms0) && (xx<m0sp) && (xx<m0sm) && (xx<mmsm)              && (xx<mpsp) && (xx<mmsp)) {flag=0; m=m+dx; s=s-dx; mftc=xx;}	
		xx=mpsp; if ((xx<mftc) && (xx<mps0) && (xx<mms0) && (xx<m0sp) && (xx<m0sm) && (xx<mmsm) && (xx<mpsm)              && (xx<mmsp)) {flag=0; m=m+dx; s=s+dx; mftc=xx;}
		xx=mmsp; if ((xx<mftc) && (xx<mps0) && (xx<mms0) && (xx<m0sp) && (xx<m0sm) && (xx<mmsm) && (xx<mpsm) && (xx<mpsp)             ) {flag=0; m=m-dx; s=s+dx; mftc=xx;}
		
		i++;
	}
	while (flag==0);
	*mean=m;
	*sdev=s;
	mftc=misfit(ncdf, datax, datay, ldb,udb,m,s);
	*misfit2=sqrt(mftc);

//	fprintf(stderr,"m %f s %f\n",m,s);
}

void map_search(float data[Max_data], long nsamp, float *map)
{
	long	i,j,nob,maxb;
	float	min, max, bin_width, m;
	int bdata[MB];
	
	
// get min/max
	min=data[0];
	max=data[0];
	for (i=0; i<nsamp; i++) if (data[i]>max) max=data[i];
	for (i=0; i<nsamp; i++) if (data[i]<min) min=data[i];
	bin_width=(max-min)/sqrt(nsamp);
	
	nob=(int) sqrt(nsamp) +1;
 	if (nob>MB) {fprintf(stderr,"Number of required bins too large, modify MB and compile\n"); exit(0);}
	for (i=0; i<nob; i++) bdata[i]=0;
 	for (i=0; i<nsamp; i++) bdata[(int)((data[i]-min)/bin_width)]++;

	maxb=bdata[i];
	for (i=0; i<nob; i++) if (bdata[i]>maxb) {maxb=bdata[i];j=i;}
	m=j*bin_width+min;

	*map=m;
//	fprintf(stderr,"m %f s %f\n",m,s);
}


void stats(float data[Max_data],float data2[Max_data],int ndata,int *ndata2,float vmin,float vmax,float dv,float *mean,float *sdev,float *mean2,float *sdev2) 
{ 
 int bdata[MB],bdata2[MB]; 
 int i,nob;
 float sum,x,y;

 nob=(int) ((vmax-vmin)/dv) +1;
 if (nob>MB) {fprintf(stderr,"Number of required bins too large, modify MB and compile\n"); exit(0);}
 for (i=0; i<nob; i++) bdata[i]=0;
 for (i=0; i<ndata; i++) bdata[(int)((data[i]-vmin)/dv)]++;
// first values
 sum=0.0;
 for (i=0; i<ndata; i++) sum=sum+data[i];
 *mean=sum/ndata;
 sum=0;
 for (i=0; i<ndata; i++) sum=sum+(data[i]-*mean)*(data[i]-*mean);
 *sdev=sqrt(sum/ndata);

// second values
// prior distribution for voronoi cells
 if (TRIA==0) for (i=0; i<nob; i++) bdata2[i]=(int)(-1.0*ndata/nob);
// empirical prior distribution for triangles

 if (TRIA==1) 
 {
  for (i=0; i<nob; i++) 
  {
   x=i*dv/(vmax-vmin);
   y=dv/0.001*ndata/1000000.0/(vmax-vmin)*(-15.483392084+2063.40295127*x+22054.6992247*x*x-65489.6918778*x*x*x+72439.7863684*x*x*x*x-40425.0269607*x*x*x*x*x+4758.48553348*x*x*x*x*x*x+4595.49354021*x*x*x*x*x*x*x);
   bdata2[i]=(int)(-y);
  }
 }

 
//fprintf(stderr, "XXX %d\n",bdata2[0]);
 *ndata2=0;
 for (i=0; i<ndata; i++) 
 {
   bdata2[(int)((data[i]-vmin)/dv)]++;
   if (bdata2[(int)((data[i]-vmin)/dv)]>0) 
   {
    data2[*ndata2]=data[i];
    *ndata2=*ndata2+1;
   }
 }
//fprintf(stderr, "XXX %d of %d\n",*ndata2, ndata);


 sum=0.0;
 for (i=0; i<*ndata2; i++) sum=sum+data2[i];
 *mean2=sum/(*ndata2);
 sum=0;
 for (i=0; i<*ndata2; i++) sum=sum+(data2[i]-*mean2)*(data2[i]-*mean2);
 *sdev2=sqrt(sum/(*ndata2));
}


// --------------main-------------------------------------------------------------


int main(int argc, char *argv[])
{
 int i,j,k;

 FILE *finp;
 struct Model model,temp_mod;

 struct GRDHEAD  gh;			/* FD grid header structure 	*/
 float rms;
 char *buf;
 char *c,*c2;
 int l, mcount, ii;
 float vmin, vmax, dv,  vv, vvx;
 int ndata2;
 float fdummy,xx,yy,zz;
 int idummy;
 char   sdummy[1000];
 FILE 		*config_file;		/* config file	*/

 char		    line[1024];

 int    pflag;
 int **hcountp,**hcounts,ndv, ndvpvs;
 double tt,dt,dtp,dts;
 int iip,iis,kkp,kks, noq, eqi, nos, sti;
 float pmean, pmean2, psdev, psdev2;
 float smean, smean2, ssdev, ssdev2;

 float data[Max_data], data2[Max_data];

 double eqx[MAX_NOQ], eqy[MAX_NOQ], eqz[MAX_NOQ], seqx[MAX_NOQ], seqy[MAX_NOQ], seqz[MAX_NOQ];
 double eqz2[MAX_NOQ],seqz2[MAX_NOQ],misfit1[MAX_NOQ], misfit2[MAX_NOQ],eqx3[MAX_NOQ], eqy3[MAX_NOQ],eqz3[MAX_NOQ];
 double eqt[MAX_NOQ], eqdt[MAX_NOQ], seqdt[MAX_NOQ];
 double resdtp[MAX_OBS], sresdtp[MAX_OBS];
 double resdts[MAX_OBS], sresdts[MAX_OBS];
 double np0, np1, np2, np3, ns0, ns1, ns2, ns3;
 double snp0, snp1, snp2, snp3, sns0, sns1, sns2, sns3;

 float **vs, **vp;
 float **eq_depth, **eq_x, **eq_y;


 float vpvsmin, vpvsmax, dvpvs,z0;

 int si, sj, sk, inv_flag;
 float sa, sb, xxx;
 float mm,ss,mfit1, mfit2;
 float boundary[1000];

//*  --- usage */
        if (argc<5)
        {
                fprintf(stderr, "usage: xxx config.dat input_file flag [x z x z x z ...]\n");
                fprintf(stderr, "Ch. Haberland GFZ Potsdam July 2016, Lena Delta\n");
                fprintf(stderr, "V. %04.1f\n", VERSION);
                exit (0);
        }



        pflag=atoi(argv[3]);
/* ---  read config file */
	if (!(config_file = fopen (argv[1], "r"))) { fprintf(stderr, "could not open config file %s\n",argv[1]); exit (0);}
// read config file
        read_single_line(config_file, line); sscanf(line, "%f ", &gh.h); 
        read_single_line(config_file, line); sscanf(line, "%d ", &gh.nx); 
        read_single_line(config_file, line); sscanf(line, "%d ", &gh.ny); 
        read_single_line(config_file, line); sscanf(line, "%d ", &gh.nz); 
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &z0); 
        read_single_line(config_file, line); sscanf(line, "%d ", &idummy); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vmin); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vmax); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vpvsmin); 
        read_single_line(config_file, line); sscanf(line, "%f ", &vpvsmax); 
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f %d", &fdummy, &idummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy); //res max
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
	read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%d ", &idummy);  //DR
        read_single_line(config_file, line); sscanf(line, "%f ", &fdummy);
        read_single_line(config_file, line); sscanf(line, "%d ", &TRIA);
        read_single_line(config_file, line); sscanf(line, "%d ", &idummy);  
        read_single_line(config_file, line); sscanf(line, "%d ", &idummy);  
        read_single_line(config_file, line); sscanf(line, "%d ", &idummy);  
        read_single_line(config_file, line); sscanf(line, "%s ", sdummy); // string
        read_single_line(config_file, line); sscanf(line, "%d ", &inv_flag);  



	fclose (config_file);


 dv=atof(argv[3]);
 ndv=(int) ((vmax-vmin)/dv)+1;
 dvpvs=atof(argv[4]);
 ndvpvs=(int) ((vpvsmax-vpvsmin)/dvpvs)+1;

fprintf(stderr, "ndv %d %f %f\n",ndv,vmin,vmax); 
fprintf(stderr, "ndvpvs %d %f %f\n",ndvpvs,vpvsmin,vpvsmax); 
fprintf(stderr, "nx ny nz %d %d %d\n",gh.nx, gh.ny, gh.nz); 

fprintf(stderr, "%f %f TRIA %d\n",vpvsmin, vpvsmax, TRIA); 



// allocate arrays
 if (!(hcountp = malloc(sizeof(int *) * ndv))) {fprintf(stderr, "hcountp 1 malloc failed\n"); exit(0);}
 for (i = 0; i < ndv; i++) {if (!(hcountp[i] = malloc(sizeof(int) * gh.nz))) {fprintf(stderr, "hcountp 2 malloc failed\n"); exit(0);}} 
 if (!(hcounts = malloc(sizeof(int *) * ndvpvs))) {fprintf(stderr, "hcounts malloc failed\n"); exit(0);}
 for (i = 0; i < ndvpvs; i++) {if (!(hcounts[i] = malloc(sizeof(int) * gh.nz))) {fprintf(stderr, "hcounts malloc failed\n"); exit(0);}} 

// ini array
 for (i=0; i<ndv; i++) for (j=0; j<gh.nz; j++) hcountp[i][j]=0;
 for (i=0; i<ndvpvs; i++) for (j=0; j<gh.nz; j++) hcounts[i][j]=0;
 for (i=0; i<gh.nz; i++) boundary[i]=0;

// ini mask



//	for ( i=0; i<pflag; i++ ) fprintf(stderr, "XXX test %d at x=%f z=%f\n",i,xp[i],zp[i]);

 if (!(buf=malloc(5000*3*20*sizeof(char))))
  {
  	fprintf(stderr, "malloc failed for buf\n");
	exit(0);
  }

// allocate arrays


 if (!(vp = malloc(sizeof(float *) * Max_data))) {fprintf(stderr, "malloc failed\n"); exit(0);}
 if (!(vs = malloc(sizeof(float *) * Max_data))) {fprintf(stderr, "malloc failed\n"); exit(0);}
 if (!(eq_depth = malloc(sizeof(float *) * Max_data))) {fprintf(stderr, "malloc failed\n"); exit(0);}
 if (!(eq_x = malloc(sizeof(float *) * Max_data))) {fprintf(stderr, "malloc failed\n"); exit(0);}
 if (!(eq_y = malloc(sizeof(float *) * Max_data))) {fprintf(stderr, "malloc failed\n"); exit(0);}

// ini model

   for (i=0; i<MAX_NOQ; i++) eqx[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) eqy[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) eqz[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) seqx[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) seqy[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) seqz[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) eqt[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) eqdt[i]=0.0;
   for (i=0; i<MAX_NOQ; i++) seqdt[i]=0.0;
   for (i=0; i<MAX_OBS; i++) resdtp[i]=0.0;
   for (i=0; i<MAX_OBS; i++) sresdtp[i]=0.0;
   for (i=0; i<MAX_OBS; i++) resdts[i]=0.0;
   for (i=0; i<MAX_OBS; i++) sresdts[i]=0.0;
   np0=0.0; np1=0.0; np2=0.0; np3=0.0; ns0=0.0; ns1=0.0; ns2=0.0; ns3=0.0;
   snp0=0.0; snp1=0.0; snp2=0.0; snp0=0.0; sns0=0.0; sns1=0.0; sns2=0.0; sns3=0.0;

   noq=0;
   nos=0;

// first read
   fprintf(stderr, "start reading phase one\n");
   mcount=0;
   finp = fopen(argv[2], "r");
   while (!feof(finp))
   {
    if (fgets(buf, 5000*3*20*sizeof(char), finp)==NULL)
    {
     fprintf(stderr, "\nEOF reached\n");
    } 
    else
    {
     if (strncmp (buf, "mod", 3)==0)
     {
      if (mcount>Max_data) {fprintf(stderr, "increase Max_data and recompile\n"); exit(0);}
      if (!(vp[mcount] = malloc(sizeof(float *) * gh.nz))) {fprintf(stderr, "malloc failed\n"); exit(0);}
      if (!(vs[mcount] = malloc(sizeof(float *) * gh.nz))) {fprintf(stderr, "malloc failed\n"); exit(0);}
      if (!(eq_depth[mcount] = malloc(sizeof(float *) * MAX_NOQ))) {fprintf(stderr, "malloc failed\n"); exit(0);}
      if (!(eq_x[mcount] = malloc(sizeof(float *) * MAX_NOQ))) {fprintf(stderr, "malloc failed\n"); exit(0);}
      if (!(eq_y[mcount] = malloc(sizeof(float *) * MAX_NOQ))) {fprintf(stderr, "malloc failed\n"); exit(0);}

      
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      model.dimension = atoi(strtok (NULL, " "));/* read dimension */	
      rms = atof(strtok (NULL, " "));/* read rms */	
      model.p_noise0 = atof(strtok (NULL, " "));/* read model noise */ np0=np0+model.p_noise0;
      model.p_noise1 = atof(strtok (NULL, " "));/* read model noise */ np1=np1+model.p_noise1;
      model.p_noise2 = atof(strtok (NULL, " "));/* read model noise */ np2=np2+model.p_noise2;
      model.p_noise3 = atof(strtok (NULL, " "));/* read model noise */ np3=np3+model.p_noise3;
      model.s_noise0 = atof(strtok (NULL, " "));/* read model noise */ ns0=ns0+model.s_noise0;
      model.s_noise1 = atof(strtok (NULL, " "));/* read model noise */ ns1=ns1+model.s_noise1;
      model.s_noise2 = atof(strtok (NULL, " "));/* read model noise */ ns2=ns2+model.s_noise2;
      model.s_noise3 = atof(strtok (NULL, " "));/* read model noise */ ns3=ns3+model.s_noise3;

      
      
      l=0; 
      while (l<model.dimension)
      {
       model.z[l] = atof(strtok (NULL, " "));
       model.vp[l] = atof(strtok (NULL, " ")); 
       model.vpvs[l] = atof(strtok (NULL, " ")); 
//       fprintf(stderr, "XXX evaluate model %d at z %f with vp %f vs %f\n", model.number,model.z[l],model.vp[l],model.vpvs[l]);
       l++; 
//fprintf(stderr, "model %ld %d z %f v %f ps %f\n", model.number, l, model.z[l], model.vp[l], model.vpvs[l] );
      }

      fprintf(stderr, "\revaluate model %10d", mcount+1);

// sort mod
	copy_model(&temp_mod,&model);

	sj=1;	
	while ( sj != 0 )
	{
		sj=0;
		for (si=0; si<temp_mod.dimension-1; si++)
		{
			if (temp_mod.z[si]>temp_mod.z[si+1])
			{
				sj=1;
				xxx=temp_mod.z[si+1]; temp_mod.z[si+1]=temp_mod.z[si]; temp_mod.z[si]=xxx;
				xxx=temp_mod.vp[si+1]; temp_mod.vp[si+1]=temp_mod.vp[si]; temp_mod.vp[si]=xxx;
				xxx=temp_mod.vpvs[si+1]; temp_mod.vpvs[si+1]=temp_mod.vpvs[si]; temp_mod.vpvs[si]=xxx;
			}
		}
   	}


      k=0;
      for (i=0; i<gh.nz; i++) 
      {
	zz=(float)i*gh.h+z0;
	if (TRIA==0) 
	{
		vv=model.vp[find_in_cell(&model,zz)];
		vvx=model.vp[find_in_cell(&model,zz-gh.h)];

	}
	if (TRIA==1)
        {
		for (si=0; si<model.dimension-1; si++) if ((zz>=temp_mod.z[si]) && (zz<temp_mod.z[si+1])) sk=si;
		sa=(temp_mod.vp[sk+1]-temp_mod.vp[sk])/(temp_mod.z[sk+1]-temp_mod.z[sk]);
		sb=temp_mod.vp[sk]-sa*temp_mod.z[sk];
		vv= sa*zz+sb;
		vvx= vv;	

	}
// check if boundary present
	if (vv!=vvx) boundary[i]=boundary[i]+1;

// check if value inrange
	if (vv>vmax) vv=vmax; 
	if (vv<vmin) vv=vmin;
        vp[mcount][i]=vv;
	j=(int) ((vv-vmin)/dv);
	hcountp[j][i]=hcountp[j][i]+1;

	if (TRIA==0) vv=model.vpvs[find_in_cell(&model,zz)];
	if (TRIA==1)
        {
		for (si=0; si<model.dimension-1; si++) if ((zz>=temp_mod.z[si]) && (zz<temp_mod.z[si+1])) sk=si;
		sa=(temp_mod.vpvs[sk+1]-temp_mod.vpvs[sk])/(temp_mod.z[sk+1]-temp_mod.z[sk]);
		sb=temp_mod.vpvs[sk]-sa*temp_mod.z[sk];
		vv= sa*zz+sb;
	}

// check if value in range
	if (vv>vpvsmax) vv=vpvsmax; 
	if (vv<vpvsmin) vv=vpvsmin;
        vs[mcount][i]=vv;
	j=(int) ((vv-vpvsmin)/dvpvs);	
	hcounts[j][i]=hcounts[j][i]+1;
      }

      mcount=mcount+1;
     }
// quake
     if (strncmp (buf, "EQ", 2)==0)
     {
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      eqi = atoi(strtok (NULL, " "));/* EQ number */ if (eqi>noq) noq=eqi;
      rms = atof(strtok (NULL, " "));/* read rms */	
      xx = atof(strtok (NULL, " ")); eqx[eqi]=eqx[eqi]+xx;
      yy = atof(strtok (NULL, " ")); eqy[eqi]=eqy[eqi]+yy;
      zz = atof(strtok (NULL, " ")); eqz[eqi]=eqz[eqi]+zz; 
      tt = atof(strtok (NULL, " ")); eqt[eqi]=tt;
      dt = atof(strtok (NULL, " ")); eqdt[eqi]=eqdt[eqi]+dt;       
      eq_x[mcount-1][eqi]=xx;  
      eq_y[mcount-1][eqi]=yy;  
      eq_depth[mcount-1][eqi]=zz;  

     }
// RES
     if (strncmp (buf, "RES", 3)==0)
     {
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      sti = atoi(strtok (NULL, " "));/* station number */ if (sti>nos) nos=sti;
      rms = atof(strtok (NULL, " "));/* read rms */	
      dtp = atof(strtok (NULL, " "));/* read P-res */	resdtp[sti]=resdtp[sti]+dtp;
      dts = atof(strtok (NULL, " "));/* read P-res */	resdts[sti]=resdts[sti]+dts;

     }

    }
   }
   for (i=0; i<gh.nz; i++) boundary[i]=boundary[i]/mcount;
   
   fprintf(stderr, "done phase one with %d models and %d quakes\n",mcount,noq);
   noq=noq+1;
   nos=nos+1;
   fclose(finp);

   for (i=0; i<noq; i++) eqx[i]=eqx[i]/mcount;
   for (i=0; i<noq; i++) eqy[i]=eqy[i]/mcount;
   for (i=0; i<noq; i++) eqz[i]=eqz[i]/mcount;
   for (i=0; i<noq; i++) eqdt[i]=eqdt[i]/mcount;
   for (i=0; i<nos; i++) resdtp[i]=resdtp[i]/mcount;
   for (i=0; i<nos; i++) resdts[i]=resdts[i]/mcount;
   np0=np0/mcount; np1=np1/mcount; np2=np2/mcount; np3=np3/mcount;
   ns0=ns0/mcount; ns1=ns1/mcount; ns2=ns2/mcount; ns3=ns3/mcount;


// second read
   fprintf(stderr, "start reading phase two\n");	
   mcount=0;
   finp = fopen(argv[2], "r");
   while (!feof(finp))
   {
    if (fgets(buf, 5000*3*20*sizeof(char), finp)==NULL)
    {
     fprintf(stderr, "\nEOF reached\n");
    } 
    else
    {
// model
     if (strncmp (buf, "mod", 3)==0)
     {
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      model.dimension = atoi(strtok (NULL, " "));/* read dimension */	
      rms = atof(strtok (NULL, " "));/* read rms */	
      model.p_noise0 = atof(strtok (NULL, " "));/* read model noise */ snp0=snp0+(model.p_noise0-np0)*(model.p_noise0-np0);
      model.p_noise1 = atof(strtok (NULL, " "));/* read model noise */ snp1=snp1+(model.p_noise1-np1)*(model.p_noise1-np1);
      model.p_noise2 = atof(strtok (NULL, " "));/* read model noise */ snp2=snp2+(model.p_noise2-np2)*(model.p_noise2-np2);
      model.p_noise3 = atof(strtok (NULL, " "));/* read model noise */ snp3=snp3+(model.p_noise3-np3)*(model.p_noise3-np3);
      model.s_noise0 = atof(strtok (NULL, " "));/* read model noise */ sns0=sns0+(model.s_noise0-ns0)*(model.s_noise0-ns0);
      model.s_noise1 = atof(strtok (NULL, " "));/* read model noise */ sns1=sns1+(model.s_noise1-ns1)*(model.s_noise1-ns1);
      model.s_noise2 = atof(strtok (NULL, " "));/* read model noise */ sns2=sns2+(model.s_noise2-ns2)*(model.s_noise2-ns2);
      model.s_noise3 = atof(strtok (NULL, " "));/* read model noise */ sns3=sns3+(model.s_noise3-ns3)*(model.s_noise3-ns3);   
      l=0; 
      while (l<model.dimension)
      {
       model.z[l] = atof(strtok (NULL, " "));
       model.vp[l] = atof(strtok (NULL, " "));
       model.vpvs[l] = atof(strtok (NULL, " "));
       l++;
      }
      fprintf(stderr, "\revaluate model %10d", mcount+1);
      mcount=mcount+1;
     }
// quake
     if (strncmp (buf, "EQ", 2)==0)
     {
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      eqi = atoi(strtok (NULL, " "));/* EQ number */
      rms = atof(strtok (NULL, " "));/* read rms */	
      xx = atof(strtok (NULL, " ")); seqx[eqi]=seqx[eqi]+(xx-eqx[eqi])*(xx-eqx[eqi]);
      yy = atof(strtok (NULL, " ")); seqy[eqi]=seqy[eqi]+(yy-eqy[eqi])*(yy-eqy[eqi]);
      zz = atof(strtok (NULL, " ")); seqz[eqi]=seqz[eqi]+(zz-eqz[eqi])*(zz-eqz[eqi]); 
      tt = atof(strtok (NULL, " ")); eqt[eqi]=tt;
      dt = atof(strtok (NULL, " ")); seqdt[eqi]=seqdt[eqi]+(dt-eqdt[eqi])*(dt-eqdt[eqi]);
//      fprintf(stderr, "evaluate EQ %d %d\n", model.number, eqi);
     }
// RES
     if (strncmp (buf, "RES", 3)==0)
     {
      c = strtok (buf, " ");	/* read "mod" */
      c2 = strtok (NULL, " ");	/* read type */
      model.number= atoi(strtok (NULL, " "));	/* read index of model */ 
      sti = atoi(strtok (NULL, " "));/* station number */
      rms = atof(strtok (NULL, " "));/* read rms */	
      dtp = atof(strtok (NULL, " "));/* read P-res */	sresdtp[sti]=sresdtp[sti]+(dtp-resdtp[sti])*(dtp-resdtp[sti]);
      dts = atof(strtok (NULL, " "));/* read P-res */	sresdts[sti]=sresdts[sti]+(dts-resdts[sti])*(dts-resdts[sti]);
//      fprintf(stderr, "evaluate RES\n");

     }
    }
   }
   
   fprintf(stderr, "done phase two %d models\n",mcount);	
   fclose(finp);

   for (i=0; i<noq; i++) seqx[i]=sqrt(seqx[i]/mcount);
   for (i=0; i<noq; i++) seqy[i]=sqrt(seqy[i]/mcount);
   for (i=0; i<noq; i++) seqz[i]=sqrt(seqz[i]/mcount);
   for (i=0; i<noq; i++) seqdt[i]=sqrt(seqdt[i]/mcount);
   for (i=0; i<nos; i++) sresdtp[i]=sqrt(sresdtp[i]/mcount);
   for (i=0; i<nos; i++) sresdts[i]=sqrt(sresdts[i]/mcount);
   snp0=sqrt(snp0/mcount); snp1=sqrt(snp1/mcount); snp2=sqrt(snp2/mcount); snp3=sqrt(snp3/mcount);
   sns0=sqrt(sns0/mcount); sns1=sqrt(sns1/mcount); sns2=sqrt(sns2/mcount); sns3=sqrt(sns3/mcount);

// get mean/sdev for EQ z by cdf
   fprintf(stderr,"analyse quake depth\n");
   for (i=0; i<noq; i++)
   {
      fprintf(stderr,"\rEQz depth analysis of EQ %5d of %5d",i,noq-1);

      for (k=0; k<mcount; k++) data[k]=eq_depth[k][i];
      mm=eqz[i]; ss=seqz[i];
      if ( inv_flag != 1 ) gsearch(data,mcount, &mm, &ss, 0.01, &mfit1, &mfit2);
      eqz2[i]=mm;
      seqz2[i]=ss;
      misfit1[i]=mfit1;
      misfit2[i]=mfit2;
   }
   fprintf(stderr,"\n");

   
// get MAP for EQ z auto binning Square-root choice
   for (i=0; i<noq; i++)
   {
      fprintf(stderr,"\rMAP depth analysis of EQ %5d of %5d",i,noq-1);

      for (k=0; k<mcount; k++) data[k]=eq_x[k][i];
      map_search(data,mcount, &mm);
      eqx3[i]=mm;
      for (k=0; k<mcount; k++) data[k]=eq_y[k][i];
      map_search(data,mcount, &mm);
      eqy3[i]=mm;     
      for (k=0; k<mcount; k++) data[k]=eq_depth[k][i];
      map_search(data,mcount, &mm);
      eqz3[i]=mm;     
   }

   
// output
    fprintf(stderr, "\nOutput results for %d models\n",mcount);

    for (i=0; i<gh.nz; i++) 
    {

      for (k=0; k<mcount; k++) data[k]=vp[k][i];
      stats(data, data2, mcount, &ndata2, vmin,vmax,dv,&pmean,&psdev,&pmean2,&psdev2);


      for (k=0; k<mcount; k++) data[k]=vs[k][i];
      stats(data, data2, mcount, &ndata2, vpvsmin,vpvsmax,dvpvs,&smean,&ssdev,&smean2,&ssdev2);
// get max prob
      kkp=hcountp[0][i]; iip=0;
      for (ii=0; ii<ndv; ii++) {if (hcountp[ii][i]>kkp) {kkp=hcountp[ii][i]; iip=ii;}}
      kks=hcounts[0][i]; iis=0;
      for (ii=0; ii<ndvpvs; ii++) {if (hcounts[ii][i]>kks) {kks=hcounts[ii][i]; iis=ii;}}

      fprintf(stdout, "STAN %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.5f\n", ((float)i*gh.h)+z0,pmean,psdev,smean,ssdev,  pmean2,psdev2,smean2,ssdev2, vmin+(iip+0.5)*dv,vpvsmin+(iis+0.5)*dvpvs, boundary[i]);
    }


    for (i=0; i<noq; i++) fprintf(stdout, "EQ %4d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %14.3lf %7.3lf %7.3lf %9.5f\n", i, eqx[i], eqy[i], eqz[i], seqx[i], seqy[i], seqz[i],eqt[i],eqdt[i],seqdt[i],misfit1[i]);

    for (i=0; i<noq; i++) fprintf(stdout, "EZ %4d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %14.3lf %7.3lf %7.3lf %9.5f\n", i, eqx[i], eqy[i], eqz2[i], seqx[i], seqy[i], seqz2[i],eqt[i],eqdt[i],seqdt[i],misfit2[i]);
    
    for (i=0; i<noq; i++) fprintf(stdout, "EM %4d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %14.3lf %7.3lf %7.3lf %9.5f\n", i, eqx3[i], eqy3[i], eqz3[i], seqx[i], seqy[i], 0.0,eqt[i],eqdt[i],seqdt[i],0.0);
    
    for (i=0; i<nos; i++) fprintf(stdout, "RES %4d %7.3f %7.3f %7.3f %7.3f\n", i,resdtp[i],resdts[i],sresdtp[i],sresdts[i]);
    fprintf(stdout, "NOISE %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", np0, np1, np2, np3, ns0, ns1, ns2, ns3, snp0, snp1, snp2, snp3, sns0, sns1, sns2, sns3);
    for (i=0; i<ndv; i++) for (j=0; j<gh.nz; j++) fprintf(stdout, "BINP %7.3f %9.3f %5d\n", vmin+i*dv, j*gh.h+z0,hcountp[i][j]);
    for (i=0; i<ndvpvs; i++) for (j=0; j<gh.nz; j++) fprintf(stdout, "BINV %7.3f %9.3f %5d\n", vpvsmin+i*dvpvs, j*gh.h+z0,hcounts[i][j]);


 free (buf);
 
 exit(0);
}
