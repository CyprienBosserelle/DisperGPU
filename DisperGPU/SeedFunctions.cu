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
		break;
	case 3:
		printf("explicit release\n");
		
		for (int np = 0; np < npart; np++)
		{
			fscanf(fseed, "%f,%f,%f,%f\t%*s", &pxo, &pyo, &pzo, &pto);
			xyz2ijk(nx,ny,xcoord,ycoord,pxi,pyj,pxo,pyo, pzo);
			partpos[np] = make_float4(pxi, pyj, pzo, pto);
		}

		break;
		
	}
	fclose(fseed);
}

