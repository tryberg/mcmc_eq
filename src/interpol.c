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
