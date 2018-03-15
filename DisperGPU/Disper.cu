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


#include "Header.cuh"

#include "Disper_kernel.cu"

FILE * logfile;

float4 * partpos,*partpos_g; //Particule position x,y,z,t


float *Uo, *Un; //U velocity, Step 0 and step n
float *Vo, *Vn; //V velocity, Step 0 and step n
float *Ux, *Vx, *hhx; // U and V velocity at the dispersal step
float *hho, *hhn; // Water depth, Step 0 and step n
float *Uo_g, *Un_g, *Ux_g; //Same on GPU plus at t particle step
float *Vo_g, *Vn_g, *Vx_g; // Same on GPU plus at t particle step
float *hho_g, *hhn_g, *hhx_g;// Same on GPU plus at t particle step

float *Nincel, *cNincel, *cTincel; // Number of particle in cell, Cumulative Nincel, Cumulative time in cell CPU
float *Nincel_g, *cNincel_g, *cTincel_g; // Number of particle in cell, Cumulative Nincel, Cumulative time in cell on GPU

float *distX, *distY; // Distance calculated between cells
float *xcoord, *ycoord; // REal world coordinates


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
void GPUstep()
{
	dim3 blockDimHD(16, 16, 1);
	//dim3 gridDimHD(ceil(max(netau,netav)*max(nxiu,nxiv) / (float)blockDimHD.x), 1, 1);
	dim3 gridDimHD((int)ceil((float)nx / (float)blockDimHD.x), (int)ceil((float)ny / (float)blockDimHD.y), 1);

	if (totaltime >= hddt*(hdstep - hdstart + 1))//+1 because we only read the next step when time exeed the previous next step
	{
		//Read next step
		printf("Reading Next step\n");
		hdstep++;

		int steptoread = hdstep;

		if (backswitch>0)
		{
			steptoread = hdend - hdstep;
		}

		NextHDstep<<<gridDimHD, blockDimHD, 0 >>>(nx, ny, Uo_g, Un_g);
		CUDA_CHECK(cudaDeviceSynchronize());

		NextHDstep<<<gridDimHD, blockDimHD, 0 >>>(nx, ny, Vo_g, Vn_g);
		CUDA_CHECK(cudaDeviceSynchronize());

		//NextHDstep<<<gridDimHD, blockDimHD, 0 >>>(nx, ny, hho_g, hhn_g);
		//CUDA_CHECK(cudaDeviceSynchronize());

		
		//readHDstepHYCOM(ncfile, Uvarname, Vvarname, nx, ny, steptoread, lev, Un, Vn, hhn);
		readHDstep(ncfile, Uvarname, Vvarname, hhvarname, nx, ny, steptoread, lev, Un, Vn, hhn);

		CUDA_CHECK(cudaMemcpy(Un_g, Un, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Vn_g, Vn, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		//CUDA_CHECK(cudaMemcpy(hhn_g, hhn, nx*ny*sizeof(float), cudaMemcpyHostToDevice));

	}
	//printf("Run GPU step\n");

	//printf("Nincel Reset\n");
	ResetNincel<<<gridDimHD, blockDimHD, 0 >>>(nx, ny, Nincel_g);
	CUDA_CHECK(cudaDeviceSynchronize());

	//printf("interp HD\n");
	//int interpstep = hdstep - hdstart + 1;
	//InterpstepCPU(nx, ny, backswitch, hdstep, totaltime, hddt, Ux, Uo, Un);
	HD_interp <<< gridDimHD, blockDimHD, 0 >>>(nx, ny, backswitch, hdstep - hdstart, totaltime, hddt, Uo_g, Un_g, Ux_g);
	CUDA_CHECK(cudaDeviceSynchronize());

	HD_interp << <gridDimHD, blockDimHD, 0 >>>(nx, ny, backswitch, hdstep - hdstart, totaltime, hddt, Vo_g, Vn_g, Vx_g);
	CUDA_CHECK(cudaDeviceSynchronize());

	//HD_interp <<<gridDimHD, blockDimHD, 0 >>>(nx, ny, backswitch, hdstep, totaltime, hddt, hho_g, hhn_g, hhx_g);
	//CUDA_CHECK(cudaDeviceSynchronize());


	//printf("Mem copy to array\n");
	CUDA_CHECK(cudaMemcpyToArray(Ux_gp, 0, 0, Ux_g, nx*ny* sizeof(float), cudaMemcpyDeviceToDevice));
	CUDA_CHECK(cudaMemcpyToArray(Vx_gp, 0, 0, Vx_g, nx*ny* sizeof(float), cudaMemcpyDeviceToDevice));
	//CUDA_CHECK(cudaMemcpyToArray(hhx_gp, 0, 0, hhx_g, nx*ny* sizeof(float), cudaMemcpyDeviceToDevice));
	//Generate some random numbers
	// Set seed 
	//curandSetPseudoRandomGeneratorSeed(gen, SEED);
	// Generate n floats on device 

	//printf("Rnd gen\n");
	curandGenerateUniform(gen, d_Rand, np);

	//printf("Part position\n");
	//run the model
	//int nbblocks=npart/256;
	dim3 blockDim(256, 1, 1);
	dim3 gridDim(np / blockDim.x, 1, 1);

	//Calculate particle step

	updatepartpos <<<gridDim, blockDim, 0 >>>(np, dt, Eh, d_Rand, partpos_g);
	CUDA_CHECK(cudaDeviceSynchronize());

	//ij2lonlat <<<gridDim, blockDim, 0 >> >(np, xl_g, yl_g, xp_g, yp_g);
	//CUDA_CHECK(cudaThreadSynchronize());
	//printf("Calc Nincel\n");
	CalcNincel <<<gridDim, blockDim, 0 >>>(np, nx, ny, partpos_g, Nincel_g, cNincel_g, cTincel_g);
	CUDA_CHECK(cudaDeviceSynchronize());
}




void CPUstep()
{
	if (totaltime >= hddt*(hdstep - hdstart + 1))//+1 because we only read the next step when time exeed the previous next step
	{
		//Read next step

		hdstep++;

		int steptoread = hdstep;

		if (backswitch>0)
		{
			steptoread = hdend - hdstep;
		}

		NextstepCPU(nx,ny, Uo, Vo, hho, Un, Vn, hhn);
		//readHDstepHYCOM(ncfile, Uvarname, Vvarname, nx, ny, steptoread, lev, Un, Vn, hhn);
		readHDstep(ncfile, Uvarname, Vvarname, hhvarname, nx, ny, steptoread, lev, Un, Vn, hhn);


	}

	//Interpolate U vel
	InterpstepCPU( nx, ny, backswitch, hdstep, totaltime, hddt, Ux, Uo, Un);

	//Interpolate V vel
	InterpstepCPU(nx, ny, backswitch, hdstep, totaltime, hddt, Vx, Vo, Vn);

	//Interpolate Water depth
	InterpstepCPU(nx, ny, 1.0f, hdstep, totaltime, hddt, hhx, hho, hhn);

	// Reseed random number
	//not needed here?

	// Update particle position
	updatepartposCPU(nx, ny, np, dt, Eh, Ux, Vx, hhx, distX, distY, partpos);

	//update Nincel
	calcNincelCPU(np, nx, ny, partpos, Nincel, cNincel, cTincel);


}



int main()
{
	//Model starts Here//

	//The main function setups all the init of the model and then calls the mainloop to actually run the model


	//First part reads the inputs to the model 
	//then allocate memory on GPU and CPU
	//Then prepare and initialise memory and arrays on CPU and GPU
	// Prepare output file
	// Run main loop
	// Clean up and close


	// Start timer to keep track of time 
	clock_t startcputime, endcputime;


	startcputime = clock();
	Param Dparam;
	Control Dcontrol;

	// Initialise totaltime
	Dcontrol.totaltime = 0.0;
	Dcontrol.nextouttime = 0.0;

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
	

	
	fscanf(fop, "%*s %s\t%*s", &Dparam.ncfile); //HD file name should have U V velocity and depth
	fscanf(fop, "%s\t%*s", &Dparam.Uvarname);
	fscanf(fop, "%s\t%*s", &Dparam.Vvarname);
	fscanf(fop, "%s\t%*s", &Dparam.hhvarname);
	fscanf(fop, "%f\t%*s", &Dparam.hddt);
	fscanf(fop, "%d,%d\t%*s", &Dparam.hdstart, &Dparam.hdend);
	fscanf(fop, "%d\t%*s", &Dparam.lev);
	fscanf(fop, "%d\t%*s", &Dparam.geocoord);
	fscanf(fop, "%d\t%*s", &Dparam.backswitch);
	fscanf(fop, "%d\t%*s", &Dparam.partmode);
	fscanf(fop, "%u\t%*s", &Dparam.np);
	fscanf(fop, "%f\t%*s", &Dcontrol.dt);
	fscanf(fop, "%f\t%*s", &Dparam.Eh);
	fscanf(fop, "%f\t%*s", &Dparam.Ev);
	fscanf(fop, "%f\t%*s", &Dparam.minrwdepth);
	fscanf(fop, "%s\t%*s", &Dparam.seedfile);
	//fscanf(fop, "%d\t%*s", &GPUDEV);

	//fscanf(fop, "%d\t%*s", &outtype);
	fscanf(fop, "%f\t%*s", &Dcontrol.outtime);
	fscanf(fop, "%s\t%*s", &Dparam.ncoutfile);

	fclose(fop);

	fprintf(logfile, "Complete\n");
	fprintf(logfile, "Reading netCDF file : %s...\n", Dparam.ncfile);
	printf("Reading netCDF file: %s...\n", Dparam.ncfile);
	//readgridsizeHYCOM(ncfile, Uvarname, Vvarname, nt, nx, ny, xcoord, ycoord);
	readgridsize(Dparam.ncfile, Dparam.Uvarname, Dparam.Vvarname, Dparam.hhvarname, Dparam.nt, Dparam.nx, Dparam.ny, xcoord, ycoord);


	fprintf(logfile, "\t nx=%d\tny=%d\tnt=%d\n", Dparam.nx, Dparam.ny, Dparam.nt);
	printf("\t nx=%d\tny=%d\tnt=%d\n", Dparam.nx, Dparam.ny, Dparam.nt);
	fprintf(logfile, "...done\n");
	printf("...done\n");

	int nx = Dparam.nx;
	int ny = Dparam.ny;
	int nt = Dparam.nt;


	//set up CPU mem
	printf("Allocate CPU memory... ");
	//Vel ARRAYS
	Uo = (float *)malloc(nx*ny*sizeof(float));
	Un = (float *)malloc(nx*ny*sizeof(float));
	Vo = (float *)malloc(nx*ny*sizeof(float));
	Vn = (float *)malloc(nx*ny*sizeof(float));
	hho = (float *)malloc(nx*ny*sizeof(float));
	hhn = (float *)malloc(nx*ny*sizeof(float));

	Ux = (float *)malloc(nx*ny*sizeof(float));
	Vx = (float *)malloc(nx*ny*sizeof(float));
	hhx = (float *)malloc(nx*ny*sizeof(float));

	distX = (float *)malloc(nx*ny*sizeof(float));
	distY = (float *)malloc(nx*ny*sizeof(float));

	/* initialize random seed: */
	srand((unsigned int)time(NULL));


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
			cNincel[i + j*nx] = 0.0f;
			cTincel[i + j*nx] = 0.0f;
		}
	}
	printf("...done\n");

	printf("Calculate distance array... ");
	CalcDistXY(nx, ny, Dparam.geocoord, xcoord, ycoord, distX,distY);
	printf("...done\n");
	//Calculate first HD step
	//outstep=10;
	Dcontrol.stp = 0;//hdstart*hddt/dt;
	Dcontrol.hdstep = Dcontrol.hdstart;


	
	//printf("HD step:%d\n ",hdstep);
	if (Dcontrol.hdend == 0)
	{
		Dcontrol.hdend = Dparam.nt - 1;
	}

	if (Dcontrol.outtime == 0.0f)
	{
		Dcontrol.outtime = Dparam.hddt*Dcontrol.hdend;
	}
	Dcontrol.nextouttime = Dcontrol.outtime;

	int steptoread = Dcontrol.hdstep;

	if (Dparam.backswitch>0)
	{
		steptoread = Dcontrol.hdend - Dcontrol.hdstep;
	}

	
	//////////////////////////////
	//Read first step in Hd model
	///////////////////////////////

	//readHDstepHYCOM(ncfile, Uvarname, Vvarname, nx, ny, steptoread, lev, Uo, Vo, hho);
	
	
	//Also read next step?
	//readHDstepHYCOM(ncfile, Uvarname, Vvarname, nx, ny, steptoread+1, lev, Un, Vn, hhn);

	readHDstep(Dparam.ncfile, Dparam.Uvarname, Dparam.Vvarname, Dparam.hhvarname, Dparam.nx, Dparam.ny, steptoread, Dparam.lev, Uo, Vo, hho);

	//Also read next step?
	readHDstep(Dparam.ncfile, Dparam.Uvarname, Dparam.Vvarname, Dparam.hhvarname, Dparam.nx, Dparam.ny, steptoread + 1, Dparam.lev, Un, Vn, hhn);

	//Calculate best dt
	if (!(Dcontrol.dt > 0.0f))// if dt==0.0
	{
		Calcmaxstep(nx, ny, Dcontrol.dt, Dparam.hddt, Uo, Vo, Un, Vn, distX, distY);
	}
	Dcontrol.olddt = Dcontrol.dt;
	printf("Allocating CPU memory for particle position... ");
	//Initialise particles on CPU
	partpos = (float4 *)malloc(Dparam.np*sizeof(float4));
	d_Rand = (float *)malloc(Dparam.np*sizeof(float));

	//partpos[50] = make_float4(0.0f, 1.0f, 5.0f, 0.2);
	printf("...done.\n");
	//printf("partpos.x=%f", partpos[50].z);
	


	//printf("partpos.x=%f", partpos[50].x);
	//Find GPU
	int nDevices;

	cudaGetDeviceCount(&nDevices);// Crash when using CUDA check?
	//GPUDEV = -1;
	if (nDevices > 0)
	{
		printf("(%i) Cuda device(s) found!\n",nDevices);
	}
	else
	{
		printf("No GPU found. Using CPU only\n");
		Dparam.GPUDEV = -1;
	}

	if (Dparam.GPUDEV > nDevices && Dparam.GPUDEV>0)
	{
		printf("Specified GPU Device not found, Using Device %i.\n",0);
		Dparam.GPUDEV = 0;
	}
	if (Dparam.GPUDEV >= 0)
	{
		printf("Allocating mem on GPU...");
		CUDA_CHECK(cudaSetDevice(Dparam.GPUDEV)); //Add error handling
		//If GPU available then copy set up GPU mem
		

		CUDA_CHECK(cudaMalloc((void **)&partpos_g, np*sizeof(float4)));

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

		printf("Transfert vectors to GPU memory... ");
		CUDA_CHECK(cudaMemcpy(Uo_g, Uo, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Un_g, Uo, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Ux_g, Uo, nx*ny*sizeof(float), cudaMemcpyHostToDevice));

		CUDA_CHECK(cudaMemcpy(Vo_g, Vo, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Vn_g, Vo, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(Vx_g, Vo, nx*ny*sizeof(float), cudaMemcpyHostToDevice));

		CUDA_CHECK(cudaMemcpy(Nincel_g, Nincel, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(cNincel_g, Nincel, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(cTincel_g, Nincel, nx*ny*sizeof(float), cudaMemcpyHostToDevice));

		CUDA_CHECK(cudaMemcpy(partpos_g, partpos, np*sizeof(float4), cudaMemcpyHostToDevice));

		//CUDA_CHECK(cudaMemcpy(partpos_g, partpos, np*sizeof(float4), cudaMemcpyHostToDevice));
		//done later

		// Loading random number generator
		curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);

		CUDA_CHECK(cudaMalloc((void **)&d_Rand, np*sizeof(float)));

		printf(" ...done\n");

		printf("Create textures on GPU memory... ");
		// Copy velocity arrays
		CUDA_CHECK(cudaMallocArray(&Ux_gp, &channelDescU, nx, ny));
		CUDA_CHECK(cudaMallocArray(&Vx_gp, &channelDescV, nx, ny));

		CUDA_CHECK(cudaMemcpyToArray(Ux_gp, 0, 0, Uo, nx*ny* sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpyToArray(Vx_gp, 0, 0, Vo, nx*ny* sizeof(float), cudaMemcpyHostToDevice));

		texU.addressMode[0] = cudaAddressModeWrap;
		texU.addressMode[1] = cudaAddressModeWrap;
		texU.filterMode = cudaFilterModeLinear;
		texU.normalized = false;


		CUDA_CHECK(cudaBindTextureToArray(texU, Ux_gp, channelDescU));

		texV.addressMode[0] = cudaAddressModeWrap;
		texV.addressMode[1] = cudaAddressModeWrap;
		texV.filterMode = cudaFilterModeLinear;
		texV.normalized = false;

		CUDA_CHECK(cudaBindTextureToArray(texV, Vx_gp, channelDescV));

		CUDA_CHECK(cudaMallocArray(&distX_gp, &channelDescdX, nx, ny));
		//CUDA_CHECK( cudaMallocArray( &distXV_gp, &channelDescdXV, netav, nxiv ));
		//CUDA_CHECK( cudaMallocArray( &distYU_gp, &channelDescdYU, netau, nxiu ));
		CUDA_CHECK(cudaMallocArray(&distY_gp, &channelDescdY, nx, ny));

		CUDA_CHECK(cudaMallocArray(&xcoord_gp, &channelDescxcoord, nx, ny));
		CUDA_CHECK(cudaMallocArray(&ycoord_gp, &channelDescycoord, nx, ny));
		//CUDA_CHECK( cudaMallocArray( &lon_vgp, &channelDesclonv, netav, nxiv ));
		//CUDA_CHECK( cudaMallocArray( &lat_vgp, &channelDesclatv, netav, nxiv ));

		CUDA_CHECK(cudaMemcpyToArray(distX_gp, 0, 0, distX, nx*ny* sizeof(float), cudaMemcpyHostToDevice));
		//CUDA_CHECK( cudaMemcpyToArray(distYU_gp, 0, 0, distYU, netau*nxiu* sizeof(float), cudaMemcpyHostToDevice));
		//CUDA_CHECK( cudaMemcpyToArray(distXV_gp, 0, 0, distXV, netav*nxiv* sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpyToArray(distY_gp, 0, 0, distY, nx*ny* sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpyToArray(xcoord_gp, 0, 0, xcoord, nx*ny* sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpyToArray(ycoord_gp, 0, 0, ycoord, nx*ny* sizeof(float), cudaMemcpyHostToDevice));
		//CUDA_CHECK( cudaMemcpyToArray(lon_vgp, 0, 0, lon_v, netav*nxiv* sizeof(float), cudaMemcpyHostToDevice));
		//CUDA_CHECK( cudaMemcpyToArray(lat_vgp, 0, 0, lat_v, netav*nxiv* sizeof(float), cudaMemcpyHostToDevice));

		texlonu.addressMode[0] = cudaAddressModeWrap;
		texlonu.addressMode[1] = cudaAddressModeWrap;
		texlonu.filterMode = cudaFilterModeLinear;
		texlonu.normalized = false;

		CUDA_CHECK(cudaBindTextureToArray(texlonu, xcoord_gp, channelDescxcoord));

		texlatu.addressMode[0] = cudaAddressModeWrap;
		texlatu.addressMode[1] = cudaAddressModeWrap;
		texlatu.filterMode = cudaFilterModeLinear;
		texlatu.normalized = false;

		CUDA_CHECK(cudaBindTextureToArray(texlatu, ycoord_gp, channelDescycoord));

		texdXU.addressMode[0] = cudaAddressModeWrap;
		texdXU.addressMode[1] = cudaAddressModeWrap;
		texdXU.filterMode = cudaFilterModeLinear;
		texdXU.normalized = false;

		CUDA_CHECK(cudaBindTextureToArray(texdXU, distX_gp, channelDescdX));

		texdYV.addressMode[0] = cudaAddressModeWrap;
		texdYV.addressMode[1] = cudaAddressModeWrap;
		texdYV.filterMode = cudaFilterModeLinear;
		texdYV.normalized = false;

		CUDA_CHECK(cudaBindTextureToArray(texdYV, distY_gp, channelDescdY));

		printf(" ...done\n");



	}


	//read seed file //calculate seed position on the GPU if available
	
	readseedfile(seedfile, np, nx, ny, xcoord, ycoord, partpos);
	//Output seed information for sanity checks
	writexyz(np, nx, ny, xcoord, ycoord, partpos, "OutSeed_000T.xyz");
	//writexyz(xp, yp, zp, tp, xl, yl, npart, fileoutn);
	//create netcdf file
	
	creatncfile(ncoutfile, nx, ny, np, xcoord, ycoord, 0.0f, Nincel, cNincel, cTincel, partpos);


	//Run CPU/GPU loop
	totaltime = 0.0f;
	
	if (partmode > 0)
	{
		if (GPUDEV < 0) //CPU mainloop
		{
			printf("Model starting using CPU. dt=%f; \n", dt);
			printf("step %f of %f\n", totaltime, hddt*(hdend - hdstart));
			while ((hddt*hdend - totaltime) > 0.0f)
			{
				dt = min(dt, nextouttime - totaltime);
				CPUstep();
				totaltime = totaltime + dt;
				stp++;

				if ((nextouttime - totaltime) < 0.001f) // Round off error checking
				{
					//WriteoutCPU();
					char fileoutn[15];
					sprintf(fileoutn, "Part_%d.xyz", stp);
					writexyz(np, nx, ny, xcoord, ycoord, partpos, fileoutn);
					//writestep2nc(ncoutfile, nx, ny, totaltime, Nincel, cNincel, cTincel);
					writestep2nc(ncoutfile, nx, ny, np, totaltime, xcoord, ycoord, Nincel, cNincel, cTincel, partpos);
					nextouttime = nextouttime + outtime;
					dt = olddt;
					//reset Nincel 
					resetNincelCPU(nx, ny, Nincel);
				}


			}
			printf("Model Completed\n Total Number of step:%d\t total nuber of outputs steps:%d\n", stp, 0);

		}
		else //GPU main loop
		{
			//Initial particle position transfert to GPU
			printf("Copy Particle position to GPU.\n");
			CUDA_CHECK(cudaMemcpy(partpos_g, partpos, np*sizeof(float4), cudaMemcpyHostToDevice));

			printf("Model starting using GPU.\n");
			while ((hddt*hdend - totaltime) > 0.0f)
			{
				dt = min(dt, nextouttime - totaltime);
				//printf("dt=%f.\n",dt);
				GPUstep();
				totaltime = totaltime + dt;
				stp++;

				if ((nextouttime - totaltime) < 0.001f) // Round off error checking
				{
					//WriteoutCPU();
					char fileoutn[15];
					sprintf(fileoutn, "Part_%d.xyz", stp);
					//writexyz(np, nx, ny, xcoord, ycoord, partpos, fileoutn);
					//writestep2nc(ncoutfile, nx, ny, totaltime, Nincel, cNincel, cTincel);
					CUDA_CHECK(cudaMemcpy(partpos, partpos_g, np*sizeof(float4), cudaMemcpyDeviceToHost));
					CUDA_CHECK(cudaMemcpy(Nincel, Nincel_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
					CUDA_CHECK(cudaMemcpy(cNincel, cNincel_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
					CUDA_CHECK(cudaMemcpy(cTincel, cTincel_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));


					writestep2nc(ncoutfile, nx, ny, np, totaltime, xcoord, ycoord, Nincel, cNincel, cTincel, partpos);
					nextouttime = nextouttime + outtime;
					dt = olddt;
					//reset Nincel 
					resetNincelCPU(nx, ny, Nincel);
				}


			}
		}
	}
	//Close and clean up
    
	fclose(logfile);
    return 0;
}

