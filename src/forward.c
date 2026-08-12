
float cal_fit_newx (struct Model *m, struct DATA *d, float ***tttp, float ***ttts, struct GRDHEAD gh, int calct, float *mfp0, float *mfs0, float *mfp1, float *mfs1, float *mfp2, float *mfs2, float *mfp3, float *mfs3, int flag, int eikonal, int out) 
{
	int i, j, k;
	float tp, ts, dist;

        float sum;
	float rmsp0, rmss0, rmsp1, rmss1, rmsp2, rmss2, rmsp3, rmss3; /* rms value (data fit) for this model */
	float station_correction;


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
//		fprintf(stderr, "EVENT %d: x: %f\ty: %f\tz: %f\n", i, m->eq[i].x, m->eq[i].y, m->eq[i].z);
		for (j=0; j<d[i].nobs_p; j++)
		{
			dist = dst(d[i].p_picks[j].x, m->eq[i].x, d[i].p_picks[j].y, m->eq[i].y);
			tp=0;
			if (eikonal==0) tp = sqrt(dist*dist+m->eq[i].z*m->eq[i].z)/m->vp[find_in_cell(m,0.0)];
 			if (eikonal==1) tp = traveltimet(tttp[d[i].p_picks[j].layer], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].p_picks[j].w1+traveltimet(tttp[d[i].p_picks[j].layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].p_picks[j].w2;
			station_correction=m->pres[d[i].p_picks[j].st_id];
//fprintf(stderr,"XXX P %d %f\n",d[i].p_picks[j].st_id,station_correction);
			if (station_correction<-1000) {fprintf(stderr, "ERROR points to invalid P station correction\n"); exit(0);}
			tp+=station_correction;	/* add residual static correction */
//			fprintf(stderr, "EQ=%d station %03d P dist: %f tpcalc: %f tpobs: %f vstat: %f\n", i, d[i].p_picks[j].st_id, dist, tp, d[i].p_picks[j].t,vpstat);
			d[i].p_picks[j].t_pred=tp;

		}
// diff
// sum
		sum=0;
		for (k=0; k<d[i].nobs_p; k++) sum=sum+d[i].p_picks[k].t_pred-d[i].p_picks[k].t;

// S travel times
		for (j=0; j<d[i].nobs_s; j++)
		{
			dist = dst(d[i].s_picks[j].x, m->eq[i].x, d[i].s_picks[j].y, m->eq[i].y);
			ts=0;
			if (eikonal==0) ts = sqrt(dist*dist+m->eq[i].z*m->eq[i].z)/(m->vp[find_in_cell(m,0.0)]/m->vpvs[find_in_cell(m,0.0)]);
 			if (eikonal==1) ts = traveltimet (ttts[d[i].s_picks[j].layer], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].s_picks[j].w1+traveltimet(ttts[d[i].s_picks[j].layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, m->eq[i].z, gh.z0)*d[i].s_picks[j].w2;

			station_correction=m->sres[d[i].s_picks[j].st_id];
			if (station_correction<-1000) {fprintf(stderr, "ERROR points to invalid S station correction\n"); exit(0);}
			ts+=station_correction;	/* add residual static correction */
//			fprintf(stderr, "EQ=%d station %03d S dist: %f tscalc: %f tsobs: %f vstat: %f\n", i, d[i].s_picks[j].st_id, dist, ts, d[i].s_picks[j].t,vsstat);
			d[i].s_picks[j].t_pred=ts;
		}
// diff
		for (k=0; k<d[i].nobs_s; k++) sum=sum+d[i].s_picks[k].t_pred-d[i].s_picks[k].t;

		sum=sum/(d[i].nobs_s+d[i].nobs_p);

// origin time
		m->origin[i]=-sum;
		if (out==1) fprintf(stdout, "EVENT %d  %lf %f %f %f %f\n", i, d[i].reftime, m->eq[i].x, m->eq[i].y, m->eq[i].z, m->origin[i]);


//fprintf(stderr, "XXX picktime=%f elstat=%f vp=%f\n", d[i].p_picks[0].t,d[i].p_picks[0].z/vpstat,vpstat);



//		fprintf(stderr, "QQQ quake %d nopicks %d\n", i, d[i].nobs_p);
		if (out==1)
		{
			for (k=0; k<d[i].nobs_p; k++)
			{
				dist = dst(d[i].p_picks[k].x, m->eq[i].x, d[i].p_picks[k].y, m->eq[i].y);
				fprintf(stdout, "%f %f %f %f %f %f P\n", d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum,dist,m->eq[i].z,m->origin[i],d[i].p_picks[k].t,d[i].p_picks[k].t_pred);
			}

			for (k=0; k<d[i].nobs_s; k++)
			{	
				dist = dst(d[i].s_picks[k].x, m->eq[i].x, d[i].s_picks[k].y, m->eq[i].y);
				fprintf(stdout, "%f %f %f %f %f %f S\n", d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum,dist,m->eq[i].z,m->origin[i],d[i].s_picks[k].t,d[i].s_picks[k].t_pred);
			}
		}



// rms
		
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==0) rmsp0=rmsp0+(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum)*(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum);
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==0) rmss0=rmss0+(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum)*(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum);
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==1) rmsp1=rmsp1+(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum)*(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum);
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==1) rmss1=rmss1+(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum)*(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum);
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==2) rmsp2=rmsp2+(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum)*(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum);
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==2) rmss2=rmss2+(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum)*(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum);
		for (k=0; k<d[i].nobs_p; k++) if (d[i].p_picks[k].cl==3) rmsp3=rmsp3+(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum)*(d[i].p_picks[k].t_pred-d[i].p_picks[k].t-sum);
		for (k=0; k<d[i].nobs_s; k++) if (d[i].s_picks[k].cl==3) rmss3=rmss3+(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum)*(d[i].s_picks[k].t_pred-d[i].s_picks[k].t-sum);
				
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

			k=0;
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
			v= vp[iz];
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



