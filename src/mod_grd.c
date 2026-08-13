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

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>

#include <string.h>
#include "mc.h"

/* procedure to read one text file line */
void read_single_line(FILE *fp, char x[])
{
 long i;
 char c;
 i=0;
 do 
  {
   if ((c=getc(fp))!='\n' && !feof(fp)) x[i++]=c;
  }
 while (c!='\n' && !feof(fp));
 while (i<80) x[i++]=' ';
 x[i]='\0';
}

void copy_model(struct Model *dest, struct Model *src)
{
 memcpy (dest, src, sizeof (struct Model));

}

// find cell number closest to cell n
int find_neighbor_cell(struct Model *model, int n)
{
 int i,j;
 float y;
 j=0; 
 y=FLT_MAX;
 for (i=0; i<model->dimension; i++)
 {
  if (i!=n)
  {
  if ((model->z[i]-model->z[n])*(model->z[i]-model->z[n])<=y)
  {
   y=(model->z[i]-model->z[n])*(model->z[i]-model->z[n]);
   j=i;
  }
  }
 }
 return (j);
}

// find cell number closest to z
int find_in_cell(struct Model *model, float z)
{
 int i,j;
 float y;
// fprintf(stdout,"IN     %f %f %d\n",x,z,model->dimension);
 j=0; 
 y=FLT_MAX;
 for (i=0; i<model->dimension; i++)
 {
  if ((model->z[i]-z)*(model->z[i]-z)<=y)
  {
   y=(model->z[i]-z)*(model->z[i]-z);
   j=i;
  }
 }
//fprintf(stdout,"XXXXXX %d %f\n",j,z);
 return (j);
}

