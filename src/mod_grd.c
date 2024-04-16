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

