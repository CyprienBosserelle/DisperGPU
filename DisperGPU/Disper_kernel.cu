//////////////////////////////////////////////////////////////////////////////////
//DisperGPU                                                                    //
//Copyright (C) 2016 Bosserelle                                                 //
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



// declare texture reference for 2D float texture

texture<float, 2, cudaReadModeElementType> texU;
texture<float, 2, cudaReadModeElementType> texV;
texture<float, 2, cudaReadModeElementType> texlonu;
texture<float, 2, cudaReadModeElementType> texlatu;
texture<float, 2, cudaReadModeElementType> texdXU;
texture<float, 2, cudaReadModeElementType> texdYV;



__global__ void HD_interp(int nx, int ny, int stp, int backswitch, int nhdstp, float dt, float hddt/*,float *Umask*/, float * Uold, float * Unew, float * UU)
{
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int tx = threadIdx.x;
	unsigned int ty = threadIdx.y;

	__shared__ float Uxo[16][16];
	__shared__ float Uxn[16][16];
	//	__shared__ float Ums[16];


	float fac = 1.0;
	/*Ums[tx]=Umask[ix];*/


	if (backswitch>0)
	{
		fac = -1.0f;
	}


	if (ix<nx && iy<ny)
	{
		Uxo[tx][ty] = fac*Uold[ix + nx*iy]/**Ums[tx]*/;
		Uxn[tx][ty] = fac*Unew[ix + nx*iy]/**Ums[tx]*/;

		UU[ix + nx*iy] = Uxo[tx][ty] + (stp*dt - hddt*nhdstp)*(Uxn[tx][ty] - Uxo[tx][ty]) / hddt;
	}
}


__global__ void NextHDstep(int nx, int ny, float * Uold, float * Unew)
{
	//int ix = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;


	if (ix<nx && iy<ny)
	{
		Uold[ix + iy*nx] = Unew[ix + iy*nx];
	}
}

__global__ void ResetNincel(int nx, int ny, float *Nincel)
{
	//int ix = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;


	if (ix<nx && iy<ny)
	{

		Nincel[ix + iy*nx] = 0.0f;
	}
}