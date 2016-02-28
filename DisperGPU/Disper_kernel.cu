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

__global__ void updatepartpos(int npart, float dt, float Eh, float * dd_rand, float *xx, float *yy, float *zz, float *tt)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;

	float Ux = 0.0f;
	float Vx = 0.0f;
	float Xd = 0.0f; //Diffusion term
	float Yd = 0.0f;

	float distu, distv;


	float xxx, yyy, ttt;
	xxx = xx[i];
	yyy = yy[i];
	ttt = tt[i];

	if (ttt >= 0.0f)
	{
		//Interpolate wter depth, Uvel Vvel at the particle position

		Ux = tex2D(texU, xxx, yyy);
		Vx = tex2D(texV, xxx + 0.5, yyy - 0.5);// U and V don't have the same coordinates but in the number of nodes it is just off by half a grid node in both dimension
		distu = tex2D(texdXU, xxx, yyy);
		distv = tex2D(texdYV, xxx + 0.5, yyy - 0.5);
		if (distu>0.0001 && distv>0.0001)//Avoid the div by zero which makes i and j #INF
		{
			// old formulation
			//Xd=(dd_rand[i]*2-1)*sqrtf(6*Eh*dt);
			//Yd=(dd_rand[npart-i]*2-1)*sqrtf(6*Eh*dt);

			//formulation used in Viikmae et al.
			Xd = sqrtf(-4 * Eh*dt*logf(1 - dd_rand[i]))*cosf(2 * pi*dd_rand[npart - i]);
			Yd = sqrtf(-4 * Eh*dt*logf(1 - dd_rand[i]))*sinf(2 * pi*dd_rand[npart - i]);

			xx[i] = xxx + (Ux*dt + Xd) / distu; // Need to add the runge kutta scheme here or not
			yy[i] = yyy + (Vx*dt + Yd) / distv;
		}
	}

	tt[i] = ttt + dt;




}

__global__ void ij2lonlat(int npart, float * xx, float *yy, float *xp, float *yp)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;
	float lon;
	float lat;
	float xxx, yyy;
	xxx = xx[i];
	yyy = yy[i];

	lon = tex2D(texlonu, xxx, yyy);
	lat = tex2D(texlatu, xxx, yyy);

	xp[i] = lon;
	yp[i] = lat;

	//

}

__global__ void CalcNincel(int npart, int nx, int ny, float *xl, float * yl, float *tt, float *Nincel, float *cNincel, float *cTincel)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;
	float xxx, yyy, ttt;
	int ix, iy;

	xxx = xl[i];
	yyy = yl[i];
	ttt = tt[i];

	if (ttt >= 0)
	{
		if (xxx>0 && xxx<(nx - 1) && yyy>0 && yyy<ny - 1)
		{
			ix = floor(xxx);
			iy = floor(yyy);
			Nincel[ix + iy*nx] = Nincel[ix + iy*nx] + 1;
			cNincel[ix + iy*nx] = cNincel[ix + iy*nx] + 1;
			cTincel[ix + iy*nx] = cTincel[ix + iy*nx] + ttt;
		}
	}
}
