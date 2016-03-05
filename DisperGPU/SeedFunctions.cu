#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cmath>
#include <fstream>
#include <netcdf.h>
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <math.h>
#include "Header.cuh"


extern "C" void readseedfile(char seedfile[],int npart,int nx,int ny, float *xcoord,float *ycoord, float4* &partpos)
{
	int seedtype;
	//Read seed file
	//first line is a string header
	// Second line is a seed type code
	//0: point seed: 1 line x,y,z,t
	//1: square release 4x,y; 1t 1z
	//2: polygon release nx,y; 1t 1z
	//3: explicit n x,y,z,t

	float pxo, pyo, pzo, pto;
	float pxi, pyj;

	FILE * fseed;
	fseed = fopen(seedfile, "r");
	fscanf(fseed, "%*s %d\t%*s", &seedtype);

	// define varible to be used in polygon seed
	float xmin, xmax, ymin, ymax;
	int nvertx;
	float * vertx, *verty;
	float Rndx, Rndy;


	switch (seedtype){
	case 0:
		printf("Point seed\n");
		//Next line should be pxo,pyo,pzo,pto as x release, yrelease, zrelease(Not used in 2d), t release, (t is the particle age (should be negative for delayed release))
		fscanf(fseed, "%f,%f,%f,%f\t%*s", &pxo, &pyo, &pzo, &pto);
		xyz2ijk(nx, ny, xcoord, ycoord, pxi, pyj, pxo, pyo, pzo);
		for (int np = 0; np < npart; np++)
		{
			partpos[np] = make_float4(pxi, pyj, pzo, pto);
		}
		break;
	case 1:
		printf("square release\n");
		// xcenter,ycenter,width,height,zo,to
		//float xcentre, ycentre, width, height;
		//fscanf(fseed, "%f,%f,%f,%f,%f,%f\t%*s", &xcentre, &ycentre,&width, &height, &pzo, &pto);
		//int npx = round(sqrtf(npart*width / height));
		//float incx = (width) / npx;
		//pxo = xcentre - 0.5f*width;
		//pyo = ycentre - 0.5f*height;
		//Need to be finished...

		break;
	case 2:
		printf("polygon release\n");
		 //read polygon length
		fscanf(fseed, "%d",&nvertx);
		// create array to store polygon x and y
		
		vertx = (float *)malloc(nvertx*sizeof(float));
		verty = (float *)malloc(nvertx*sizeof(float));
		// read Polygon x and y
		for (int v = 0; v < nvertx; v++)
		{
			fscanf(fseed,"%f,%f",&vertx[v],&verty[v]);

		}
		//printf("Vx[0]=%f; Vy[0]=%f\n", vertx[0], verty[0]);
		//printf("Vx[2]=%f; Vy[2]=%f\n", vertx[2], verty[2]);
		//printf("Vx[5]=%f; Vy[5]=%f\n", vertx[5], verty[5]);
		//printf("Vx[10]=%f; Vy[10]=%f\n", vertx[10], verty[10]);
		//printf("Vx[50]=%f; Vy[50]=%f\n", vertx[50], verty[50]);
		
		//find min and max of polygon;
		xmin = vertx[0];
		xmax = vertx[0];
		ymin = verty[0];
		ymax = verty[0];
		for (int v = 0; v < nvertx; v++)
		{
			xmin = min(xmin, vertx[v]);
			xmax = max(xmax, vertx[v]);
			ymin = min(ymin, verty[v]);
			ymax = max(ymax, verty[v]);
		}

		printf("xmin=%f; xmax=%f\t", xmin, xmax);
		printf("ymin=%f; ymax=%f\n", ymin, ymax);


		//Seed particles
		
		
		int inpol;
		for (int np = 0; np < npart; np++)
		{
			inpol = 0;
			while (inpol != 1)
			{
				Rndx = ((float)rand() / (float)(RAND_MAX));
				Rndy = ((float)rand() / (float)(RAND_MAX));
				//printf("Rndx=%f; Rndy=%f\n", Rndx, Rndy);
				// place the particle randomly in the boundbox
				pxo = xmin + Rndx*(xmax - xmin);
				pyo = ymin + Rndy*(ymax - ymin);
				pzo = 0.0;
				pto = 0.0;
				// check whether it is in the polygon
				inpol = wn_PnPoly(pxo, pyo, vertx, verty, nvertx-1);
				//printf("Rndx=%f; Rndy=%f\n", Rndx, Rndy);
				//if not try again
			}
			//printf("pxo=%f; pyo=%f\t", pxo, pyo);
			xyz2ijk(nx, ny, xcoord, ycoord, pxi, pyj, pxo, pyo, pzo);
			//printf("pxi=%f; pyi=%f\n", pxi, pyj);
			partpos[np] = make_float4(pxi, pyj, pzo, pto);
			
		}

		

		break;
	case 3:
		printf("explicit release\n");
		
		for (int np = 0; np < npart; np++)
		{
			fscanf(fseed, "%f,%f,%f,%f", &pxo, &pyo, &pzo, &pto);
			xyz2ijk(nx,ny,xcoord,ycoord,pxi,pyj,pxo,pyo, pzo);
			partpos[np] = make_float4(pxi, pyj, pzo, pto);
		}

		break;
	case 4:
		printf("Reuse release\n");
		for (int np = 0; np < npart; np++)
		{
			fscanf(fseed, "%f,%f,%f,%f,%f,%f", &pxo, &pyo, &pzo, &pto, &pxi, &pyj);
			//xyz2ijk(nx, ny, xcoord, ycoord, pxi, pyj, pxo, pyo, pzo);
			partpos[np] = make_float4(pxi, pyj, pzo, pto);
		}
		
	}
	fclose(fseed);
}

