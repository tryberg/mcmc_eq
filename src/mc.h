#ifdef SINGLE			/* somehow needed for Triangle */
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */
#include "fdtimes.h"

#define MD 1000L		// max dimension of model (== layers)
#define PI   3.141592653	// pi
#define VERSION 1.0
#define MAX_OBS 1000		/* maximum number of picks per station   */
#define MAX_STAT 1000		/* maximum number of stations   */
#define MAX_NOQ 3500		/* maximum number of quakes   */
#define    min(x,y)    ( ((x)<(y))? (x):(y) )
#define DBUG 0


extern int TRIA, DR, aflag;

struct QUAKE			/* quake 	*/
{
	float x;		/* x-coord (km) 		*/
	float y;		/* y-coord (km) 		*/
	float z;		/* z-coord (km) 		*/
};

struct Model 
{
   long number;
   long dimension;
   long noq;
   long nos;
   float pres[MAX_STAT];
   float sres[MAX_STAT];
   float origin[MAX_NOQ];		// p origin time
   float p_noise0;
   float s_noise0;
   float p_noise1;
   float s_noise1;
   float p_noise2;
   float s_noise2;
   float p_noise3;
   float s_noise3;
   float z[MD];
   float vp[MD];
   float vpvs[MD];
   struct QUAKE eq[MAX_NOQ];	/* quakes 	*/
};

struct GRDHEAD			/* specs of the FD grid		*/
{
	int nx;			/* number of nodes in x-direction */
	int ny;			/* number of nodes in x-direction */
	int nz;			/* number of nodes in y-direction */
	float h;		/* mesh spacing 		  */	
	float x0;		/* mesh x origin 		  */
	float y0;		/* mesh x origin 		  */
	float z0;		/* mesh y origin 		  */
};

struct OBS			/* receiver/observation struct 	*/
{    
	int   st_id;
	float x;		/* x-coord (km) 		*/
	float y;		/* y-coord (km) 		*/
	float z;		/* z-coord (km) 		*/
	float t;		/* travel time (s) 		*/
	int cl;			/* pick class 			*/
	int layer;		// layer number
	float w1;		// pre-calculated weight1
	float w2;		// pre-calculated weigth2
};

struct DATA			/* shot structure		*/
{
	int	eq_id;	  	/* EQ ID */
	double reftime;
	double xfix;            // if !=-99999  event x fixed
	double yfix;            // if !=-99999  event y fixed
	double zfix;            // if !=-99999  event z fixed
	int   nobs_p; 		/* number of observations 	*/
	int   nobs_s; 		/* number of observations 	*/
	int   nobs_p0; 		/* number of observations 	*/
	int   nobs_s0; 		/* number of observations 	*/
	int   nobs_p1; 		/* number of observations 	*/
	int   nobs_s1; 		/* number of observations 	*/
	int   nobs_p2; 		/* number of observations 	*/
	int   nobs_s2; 		/* number of observations 	*/
	int   nobs_p3;
	int   nobs_s3;
	struct OBS p_picks[MAX_OBS];	/* array with obs information 	*/
	struct OBS s_picks[MAX_OBS];	/* array with obs information 	*/
};





