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
