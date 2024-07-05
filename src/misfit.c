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
//-----------------------------------------------------
// CCCCCCCCCCCCCCCCCCCCCCCCC
// HHHHHHHHHHHHHHHHHHHHHHHHH

/* ============================================= */
float cal_fit_newx (struct Model *m, struct DATA *d, int ne, float ***tttp, float ***ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3, int flag, int eikonal, int out) 
{
	int i, j;
	float tp, ts, dist;

        float tp_diff[MAX_OBS],ts_diff[MAX_OBS],sum,true_ttp[MAX_OBS],true_tts[MAX_OBS];
	float rmsp0, rmss0, rmsp1, rmss1, rmsp2, rmss2, rmsp3, rmss3; /* rms value (data fit) for this model */
	float station_correction;

	rmsp0=0.0; rmss0=0.0;
	rmsp1=0.0; rmss1=0.0;
	rmsp2=0.0; rmss2=0.0;
	rmsp3=0.0; rmss3=0.0;

// return if no tt calculations

	if (flag==1) {*mfp0=rmsp0; *mfs0=rmss0; *mfp1=rmsp1; *mfs1=rmss1; *mfp2=rmsp2; *mfs2=rmss2; *mfp3=rmsp3; *mfs3=rmss3; return(1.0);}

// update tt fields if needed
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


// loop over all quakes
	for (i=0; i<m->noq; i++)
	{

// P travel times
		for (j=0; j<d[i].nobs_p; j++)
		{
			dist = dst(d[i].p_picks[j].x, m->eq[i].x, d[i].p_picks[j].y, m->eq[i].y);
			if (eikonal==0) tp = sqrt(dist*dist+m->eq[i].z*m->eq[i].z)/m->vp[find_in_cell(m,0.0)];
 			if (eikonal==1) tp = traveltimet(tttp[d[i].p_picks[j].layer], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].p_picks[j].w1+traveltimet(tttp[d[i].p_picks[j].layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].p_picks[j].w2;
			station_correction=m->pres[d[i].p_picks[j].st_id];
			if (station_correction<-1000) {fprintf(stderr, "ERROR points to invalid P station correction\n"); exit(0);}
			tp+=station_correction;	// add station correction 
			tp_diff[j]=tp-d[i].p_picks[j].t; // P diff
			true_ttp[j]=tp; // true P tt

		}

// sum of dtp
		sum=0;
		for (j=0; j<d[i].nobs_p; j++) sum=sum+tp_diff[j];

// S travel times
		for (j=0; j<d[i].nobs_s; j++)
		{
			dist = dst(d[i].s_picks[j].x, m->eq[i].x, d[i].s_picks[j].y, m->eq[i].y);
			if (eikonal==0) ts = sqrt(dist*dist+m->eq[i].z*m->eq[i].z)/(m->vp[find_in_cell(m,0.0)]/m->vpvs[find_in_cell(m,0.0)]);
 			if (eikonal==1) ts = traveltimet (ttts[d[i].s_picks[j].layer], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z,gh.z0)*d[i].s_picks[j].w1+traveltimet(ttts[d[i].s_picks[j].layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].s_picks[j].w2;
			station_correction=m->sres[d[i].s_picks[j].st_id];
			if (station_correction<-1000) {fprintf(stderr, "ERROR points to invalid S station correction\n"); exit(0);}
			ts+=station_correction;	// add station correction
			ts_diff[j]=ts-d[i].s_picks[j].t; // S diff
			true_tts[j]=ts; // true S tt

		}

// sum of dts
		for (j=0; j<d[i].nobs_s; j++) sum=sum+ts_diff[j];
		
// average diff -> origin time - coupling between P & S model/tt's
		sum=sum/(d[i].nobs_s+d[i].nobs_p);
		m->origin[i]=-sum;
		
// de-mean
		for (j=0; j<d[i].nobs_p; j++) tp_diff[j]=tp_diff[j]-sum;
		for (j=0; j<d[i].nobs_s; j++) ts_diff[j]=ts_diff[j]-sum;		

// output			
		if (out==1) 
		{
			fprintf(stdout, "EVENT %d  %lf %f %f %f %f\n", i, d[i].reftime, m->eq[i].x, m->eq[i].y, m->eq[i].z, m->origin[i]);
			for (j=0; j<d[i].nobs_p; j++)
			{
				dist = dst(d[i].p_picks[j].x, m->eq[i].x, d[i].p_picks[j].y, m->eq[i].y);
				fprintf(stdout, "%f %f %f %f %f %f P\n", tp_diff[j],dist,m->eq[i].z,m->origin[i],d[i].p_picks[j].t,true_ttp[j]);
			}
			for (j=0; j<d[i].nobs_s; j++)
			{	
				dist = dst(d[i].s_picks[j].x, m->eq[i].x, d[i].s_picks[j].y, m->eq[i].y);
				fprintf(stdout, "%f %f %f %f %f %f S\n", ts_diff[j],dist,m->eq[i].z,m->origin[i],d[i].s_picks[j].t,true_tts[j]);
			}
		}

// rms calculations
		for (j=0; j<d[i].nobs_p; j++) if (d[i].p_picks[j].cl==0) rmsp0=rmsp0+tp_diff[j]*tp_diff[j];
		for (j=0; j<d[i].nobs_s; j++) if (d[i].s_picks[j].cl==0) rmss0=rmss0+ts_diff[j]*ts_diff[j];
		for (j=0; j<d[i].nobs_p; j++) if (d[i].p_picks[j].cl==1) rmsp1=rmsp1+tp_diff[j]*tp_diff[j];
		for (j=0; j<d[i].nobs_s; j++) if (d[i].s_picks[j].cl==1) rmss1=rmss1+ts_diff[j]*ts_diff[j];
		for (j=0; j<d[i].nobs_p; j++) if (d[i].p_picks[j].cl==2) rmsp2=rmsp2+tp_diff[j]*tp_diff[j];
		for (j=0; j<d[i].nobs_s; j++) if (d[i].s_picks[j].cl==2) rmss2=rmss2+ts_diff[j]*ts_diff[j];
		for (j=0; j<d[i].nobs_p; j++) if (d[i].p_picks[j].cl==3) rmsp3=rmsp3+tp_diff[j]*tp_diff[j];
		for (j=0; j<d[i].nobs_s; j++) if (d[i].s_picks[j].cl==3) rmss3=rmss3+ts_diff[j]*ts_diff[j];
			
	}
	*mfp0=rmsp0; *mfs0=rmss0;
	*mfp1=rmsp1; *mfs1=rmss1;
	*mfp2=rmsp2; *mfs2=rmss2;
	*mfp3=rmsp3; *mfs3=rmss3; 
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


// loop over all source depths
	for (iz=0; iz<gh.nz; iz++)
	{

		for (i=0; i<nxmod*gh.nz; i++) tbuf[i] = 0.;
		xs=0.0;
		zs=(float)iz;	
		
// eikonal solver call
		time_2d(hsbuf, tbuf, nxmod, gh.nz, xs, zs, 0.001, 0);

// loop over all receiver elevations
		for (j=0; j<gh.nz; j++)
		{
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
