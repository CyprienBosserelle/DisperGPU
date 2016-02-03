//////////////////////////////////////////////////////////////////////////////////
//DisperGPU   v0.0                                                                 //
//Copyright (C) 2015 Bosserelle                                                 //
//                                                                              //
//This program is free software: you can redistribute it and/or modify          //
//it under the terms of the GNU General Public License as published by          //
//the Free Software Foundation.                                                 //
//                                                                              //
//This program is distributed in the hope that it will be useful,               //
//but WITHOUT ANY WARRANTY; without even the implied warranty of                //    
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//GNU General Public License for more details.                                  //
//                                                                              //
//You should have received a copy of the GNU General Public License             //
//along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//////////////////////////////////////////////////////////////////////////////////


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

#include "Disper_kernel.cu"

char ncfile[256];
char Uvarname[256];
char Vvarname[256];
char hhvarname[256];
char ncoutfile[256];
char seedfile[256];

FILE * logfile;

int np;//Numer of particles
int partmode; // Particle model type: 0:test set up the model but do not run any loop
              // 1: 2D, passive particle no buyancy taken into acount
		      // 2: Quasi 3D, buoyant/sinking particle advected with 2d depth averaged hydrodynamics, Log profile is assumed
              // 3D: 3D model buoyant/sinking particle advected with 3d hydrodynamcis

float4 * partpos,*partpos_g; //Particule position x,y,z,t

int nx, ny, nz, nt; //HD input may have a more complex structure with staggered grid 

float *Uo, *Un; //U velocity, Step 0 and step n
float *Vo, *Vn; //V velocity, Step 0 and step n
float *hho, *hhn; // Water depth, Step 0 and step n
float *Uo_g, *Un_g, *Ux_g; //Same on GPU plus at t particle step
float *Vo_g, *Vn_g, *Vx_g; // Same on GPU plus at t particle step
float *hho_g, *hhn_g, *hhx_g;// Same on GPU plus at t particle step

int hdstep, hdstart, hdend; // HD model step, HD step start and HD step end
float hddt; // HD model tme step
int lev; //Level for 3D HD but 2D/Q3D particle model
int geocoord; //Geographic coordinate system switch 0 is metric 1 is degrees

float *Nincel, *cNincel, *cTincel; // Number of particle in cell, Cumulative Nincel, Cumulative time in cell CPU
float *Nincel_g, *cNincel_g, *cTincel_g; // Number of particle in cell, Cumulative Nincel, Cumulative time in cell on GPU

float *distX, *distY; // Distance calculated between cells
float *xcoord, *ycoord; // REal world coordinates

int stp, outstep, nextoutstep, outtype; // Model step, output step, next output step, output file type

int backswitch; // 0 run HD model forward 1 run the model backward
float dt; // particle model time step
float Eh, Ev; // Eddy viscosity horizontale, vertical
float minrwdepth; // Minimum depth for using Eddy viscosity

int GPUDEV = 0; // GPU device in use  (default is 0, aka first available device from device query) negative value means force use of cpu and other positive value a dpecific GPU in a multi GPU system

int SEED = 777; //Seed for random number generator
float * d_Rand; //GPU random number array
curandGenerator_t gen; // Random number generator using curand

cudaError CUDerr; // Catching CUDA errors 

cudaArray* Ux_gp; // Cuda array to pre-store HD vel data before converting to textures
cudaArray* Vx_gp; // Cuda array to pre-store HD vel data before converting to textures
cudaArray* hhx_gp; // Cuda array to pre-store HD depth data before converting to textures
cudaArray* distX_gp; // Cuda array to pre-store HD distance before converting to textures
cudaArray* distY_gp; // Cuda array to pre-store HD distance before converting to textures

cudaArray* xcoord_gp; // Cuda array to pre-store HD coordinates before converting to textures
cudaArray* ycoord_gp; // Cuda array to pre-store HD coordinates before converting to textures

// Below create channels between cuda arrays (see above) and textures
cudaChannelFormatDesc channelDescU = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc channelDescV = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc channelDeschh = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc channelDescdX = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc channelDescdY = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc channelDescxcoord = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc channelDescycoord = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);


void CUDA_CHECK(cudaError CUDerr)
{


	if (cudaSuccess != CUDerr) {

		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \

			__FILE__, __LINE__, cudaGetErrorString(CUDerr));

		exit(EXIT_FAILURE);

	}
}





int main()
{
	char logfilename[] = "DisperGPU.log";
	logfile = fopen(logfilename, "w");
	fprintf(logfile, "DisperGPU v0.0\n");


	//////////////////////////////////////////////////////
	/////             Read Operational file           /////
	//////////////////////////////////////////////////////
	fprintf(logfile, "Reading DisperGPU.dat...\t");
	printf( "Reading DisperGPU.dat\n");
	char opfile[] = "DisperGPU.dat";
	
	


	FILE * fop;
	fop = fopen(opfile, "r");

	if (fop == NULL)
	{
		fprintf(logfile, "Error opening DisperGPU.dat: %s\n", strerror(errno));
		perror("Error opening DisperGPU.dat: ");

		exit(-1);
	}
	

	
	fscanf(fop, "%*s %s\t%*s", &ncfile); //HD file name should have U V velocity and depth
	fscanf(fop, "%s\t%*s", &Uvarname);
	fscanf(fop, "%s\t%*s", &Vvarname);
	fscanf(fop, "%s\t%*s", &hhvarname);
	fscanf(fop, "%f\t%*s", &hddt);
	fscanf(fop, "%d,%d\t%*s", &hdstart, &hdend);
	fscanf(fop, "%d\t%*s", &lev);
	fscanf(fop, "%d\t%*s", &geocoord);
	fscanf(fop, "%d\t%*s", &backswitch);
	fscanf(fop, "%d\t%*s", &partmode);
	fscanf(fop, "%u\t%*s", &np);
	fscanf(fop, "%f\t%*s", &dt);
	fscanf(fop, "%f\t%*s", &Eh);
	fscanf(fop, "%f\t%*s", &Ev);
	fscanf(fop, "%f\t%*s", &minrwdepth);
	fscanf(fop, "%s\t%*s", &seedfile);
	//fscanf(fop, "%d\t%*s", &GPUDEV);

	//fscanf(fop, "%d\t%*s", &outtype);
	//fscanf(fop, "%d\t%*s", &outstep);
	//fscanf(fop, "%s\t%*s", &ncoutfile);

	fclose(fop);

	fprintf(logfile, "Complete\n");
	fprintf(logfile, "Reading netCDF file : %s...\n", ncfile);
	printf("Reading netCDF file:%s...\n", ncfile);
	readgridsize(ncfile, Uvarname, Vvarname, hhvarname,nt, nx, ny,xcoord,ycoord);


	fprintf(logfile, "\t nx=%d\tny=%d\n",nx,ny);
	printf("\t nx=%d\tny=%d\n", nx, ny);
	fprintf(logfile, "...done\n");
	printf("...done\n");


	//set up CPU mem
	printf("Allocate CPU memory... ");
	//Vel ARRAYS
	Uo = (float *)malloc(nx*ny*sizeof(float));
	Un = (float *)malloc(nx*ny*sizeof(float));
	Vo = (float *)malloc(nx*ny*sizeof(float));
	Vn = (float *)malloc(nx*ny*sizeof(float));
	hho = (float *)malloc(nx*ny*sizeof(float));
	hhn = (float *)malloc(nx*ny*sizeof(float));

	distX = (float *)malloc(nx*ny*sizeof(float));
	distY = (float *)malloc(nx*ny*sizeof(float));

	//xcoord = (float *)malloc(nx*ny*sizeof(float));// Already allocated in readgridsize subroutine
	//ycoord = (float *)malloc(nx*ny*sizeof(float));

	//Nincel
	Nincel = (float *)malloc(nx*ny*sizeof(float));
	cNincel = (float *)malloc(nx*ny*sizeof(float));
	cTincel = (float *)malloc(nx*ny*sizeof(float));


	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			Nincel[i + j*nx] = 0.0f;
		}
	}
	printf("...done\n");

	printf("Calculate distance array... ");
	CalcDistXY(nx, ny, geocoord, xcoord, ycoord, distX,distY);
	printf("...done\n");
	//Calculate first HD step
	//outstep=10;
	stp = 0;//hdstart*hddt/dt;
	hdstep = hdstart;
	nextoutstep = outstep + stp;
	//printf("HD step:%d\n ",hdstep);
	if (hdend == 0)
	{
		hdend = nt - 1;
	}

	int steptoread = hdstep;

	if (backswitch>0)
	{
		steptoread = hdend - hdstep;
	}
	//////////////////////////////
	//Read first step in Hd model
	///////////////////////////////

	readHDstep(ncfile, Uvarname, Vvarname, hhvarname, nx, ny, steptoread, lev, Uo, Vo, hho);

	printf("Allocating CPU memory for particle position... ");
	//Initialise particles on CPU
	partpos = (float4 *)malloc(np*sizeof(float4));
	//partpos[50] = make_float4(0.0f, 1.0f, 5.0f, 0.2);
	printf("...done.\n");
	//printf("partpos.x=%f", partpos[50].z);
	


	//printf("partpos.x=%f", partpos[50].x);
	//Find GPU
	int nDevices;

	CUDA_CHECK(cudaGetDeviceCount(&nDevices));

	if (nDevices > 0)
	{
		printf("(%i) Cuda device(s) found!\n",nDevices);
	}
	else
	{
		printf("No GPU found. Using CPU only\n");
		GPUDEV = -1;
	}

	if (GPUDEV > nDevices && GPUDEV>0)
	{
		printf("Specified GPU Device not found, Using Device %i.\n",0);
		GPUDEV = 0;
	}
	if (GPUDEV >= 0)
	{
		printf("Allocating mem on GPU...");
		CUDA_CHECK(cudaSetDevice(GPUDEV)); //Add error handling
		//If GPU available then copy set up GPU mem
		

		CUDA_CHECK(cudaMalloc((void **)&partpos_g, sizeof(float4)));

		CUDA_CHECK(cudaMalloc((void **)&Uo_g, nx*ny* sizeof(float)));
		CUDA_CHECK(cudaMalloc((void **)&Un_g, nx*ny* sizeof(float)));
		CUDA_CHECK(cudaMalloc((void **)&Ux_g, nx*ny* sizeof(float)));

		CUDA_CHECK(cudaMalloc((void **)&Vo_g, nx*ny* sizeof(float)));
		CUDA_CHECK(cudaMalloc((void **)&Vn_g, nx*ny* sizeof(float)));
		CUDA_CHECK(cudaMalloc((void **)&Vx_g, nx*ny* sizeof(float)));

		CUDA_CHECK(cudaMalloc((void **)&Nincel_g, nx*ny* sizeof(float)));
		CUDA_CHECK(cudaMalloc((void **)&cNincel_g, nx*ny* sizeof(float)));
		CUDA_CHECK(cudaMalloc((void **)&cTincel_g, nx*ny* sizeof(float)));


		printf(" ...done\n");

		// Loading random number generator
		curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);

		CUDA_CHECK(cudaMalloc((void **)&d_Rand, np*sizeof(float)));

		

	}


	//read seed file //calculate seed position on the GPU if available
	readseedfile(seedfile, np, partpos);

	//read input HD model

	//Run CPU/GPU loop

	//Close and clean up
    
	fclose(logfile);
    return 0;
}

