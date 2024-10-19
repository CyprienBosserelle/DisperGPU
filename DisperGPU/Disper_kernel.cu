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

#define pi 3.14159265f

// declare texture reference for 2D float texture

texture<float, 2, cudaReadModeElementType> texU;
texture<float, 2, cudaReadModeElementType> texV;
texture<float, 2, cudaReadModeElementType> texH;
texture<float, 2, cudaReadModeElementType> texlonu;
texture<float, 2, cudaReadModeElementType> texlatu;
texture<float, 2, cudaReadModeElementType> texdXU;
texture<float, 2, cudaReadModeElementType> texdYV;



__global__ void HD_interp(int nx, int ny, int backswitch, int nhdstp, float totaltime, float hddt, float * Uold, float * Unew, float * UU)
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

		UU[ix + nx*iy] = Uxo[tx][ty] + (totaltime - hddt*nhdstp)*(Uxn[tx][ty] - Uxo[tx][ty]) / hddt;
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

__global__ void updatepartpos(int npart, float dt, float Eh,float mindepth, float * dd_rand, float4 * partpos)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;

	float Ux = 0.0f;
	float Vx = 0.0f;
	float Hx = 1.0f;
	float Xd = 0.0f; //Diffusion term
	float Yd = 0.0f;

	float distu, distv;


	float xxx, yyy, zzz, ttt;

	xxx = partpos[i].x; //should be in i,j
	yyy = partpos[i].y;
	zzz = partpos[i].z;
	ttt = partpos[i].w;
	

	if (ttt >= 0.0f)
	{
		//Interpolate wter depth, Uvel Vvel at the particle position

		Ux = tex2D(texU, xxx, yyy);
		Vx = tex2D(texV, xxx + 0.5, yyy - 0.5);// U and V don't have the same coordinates but in the number of nodes it is just off by half a grid node in both dimension
		Hx = tex2D(texH, xxx, yyy);
		distu = tex2D(texdXU, xxx, yyy);
		distv = tex2D(texdYV, xxx + 0.5, yyy - 0.5);
		
		if (distu>0.0001 && distv>0.0001)//Avoid the div by zero which makes i and j #INF
		{
			if (Hx >= mindepth)
			{
				// old formulation
				//Xd=(dd_rand[i]*2-1)*sqrtf(6*Eh*dt);
				//Yd=(dd_rand[npart-i]*2-1)*sqrtf(6*Eh*dt);

				//formulation used in Viikmae et al.
				Xd = sqrtf(-4.0f * Eh*dt*logf(1 - dd_rand[i]))*cosf(2.0f * pi*dd_rand[npart - i]);
				Yd = sqrtf(-4.0f * Eh*dt*logf(1 - dd_rand[i]))*sinf(2.0f * pi*dd_rand[npart - i]);

				xxx = xxx + (Ux*dt + Xd) / distu; // Need to add the runge kutta scheme here or not if time step is small enough
				yyy = yyy + (Vx*dt + Yd) / distv;
				zzz = zzz;
			}
		}
	}

	ttt = ttt + dt;
	partpos[i] = make_float4(xxx, yyy, zzz, ttt);




}
__global__ void updatewoodpos(int npart, float dt, float Eh, float mindepth, float * dd_rand, float4 * partpos)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;

	float Ux = 0.0f; //Velocity at wood centre of G
	float Vx = 0.0f;
	float Us = 0.0f;//Velocity at wood start
	float Ue = 0.0f;//Velocity at wood end
	float Vs = 0.0f;
	float Ve = 0.0f;
	float Hx = 1.0f;
	float Xd = 0.0f; //Diffusion term
	float Yd = 0.0f;

	float distu, distv;
	float Lw, Aw, Dw, Asub, rhow, Cd, alphaw, thetaw;
	float U,Ulim,Ulog,g,rho,mubed;

	float xxx, yyy, zzz, ttt;

	xxx = partpos[i].x; //should be in i,j
	yyy = partpos[i].y;
	zzz = partpos[i].z;
	ttt = partpos[i].w;
	Lw = 10.0f;
	alphaw = 0.0f;
	thetaw = 0.0f;
	Dw = 0.3f;
	g = 9.81f;
	rhow = 900.0f;
	rho = 1000.0f;
	mubed = 0.1f;
	Ulim = 0.0f;

	if (ttt >= 0.0f)
	{
		//Interpolate wter depth, Uvel Vvel at the particle position

		
		distu = tex2D(texdXU, xxx, yyy);
		distv = tex2D(texdYV, xxx + 0.5, yyy - 0.5);

		Ux = tex2D(texU, xxx, yyy);
		Vx = tex2D(texV, xxx + 0.5, yyy - 0.5);// U and V don't have the same coordinates but in the number of nodes it is just off by half a grid node in both dimension
		Hx = tex2D(texH, xxx, yyy);

		if (Hx <= Dw)
		{
			Aw = pi*Dw*Dw / 4.0f;
			Asub = pi*Hx*Hx / 4.0f;

			Ulim = ((g*rhow*Lw*Aw) - (g*rho*Asub*Lw))*(mubed*cosf(alphaw) - sinf(alphaw)) / (0.5f*Cd*rho*(Lw*Hx*sinf(thetaw) + Aw*cosf(thetaw)));
			Ulim = sqrtf(Ulim);
		}

		U = sqrt(Ux*Ux + Vx*Vx);

		Ulog = max(U - Ulim, 0.0f);

		Ux = Ux*Ulog / max(U, 0.000001f); // Wait all this is ointless if U is zeor then nothing moves...
		Vx = Vx*Ulog / max(U, 0.000001f);

		
		if (distu>0.0001 && distv>0.0001)//Avoid the div by zero which makes i and j #INF
		{
			//if (Hx >= mindepth)
			{
				// old formulation
				//Xd=(dd_rand[i]*2-1)*sqrtf(6*Eh*dt);
				//Yd=(dd_rand[npart-i]*2-1)*sqrtf(6*Eh*dt);

				//formulation used in Viikmae et al.
				Xd = sqrtf(-4.0f * Eh*dt*logf(1 - dd_rand[i]))*cosf(2.0f * pi*dd_rand[npart - i]);
				Yd = sqrtf(-4.0f * Eh*dt*logf(1 - dd_rand[i]))*sinf(2.0f * pi*dd_rand[npart - i]);

				xxx = xxx + (Ux*dt + Xd) / distu;
				yyy = yyy + (Vx*dt + Yd) / distv;
				zzz = zzz;
			}
		}
	}

	ttt = ttt + dt;
	partpos[i] = make_float4(xxx, yyy, zzz, ttt);




}



__global__ void updatepartposQ3D(int npart, float dt, float Eh, float Ev, float mindepth, float ws, float * dd_rand, float4 * partpos)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;

	float Ux = 0.0f;
	float Vx = 0.0f;
	float Hx = 1.0f;
	float Xd = 0.0f; //Diffusion term
	float Yd = 0.0f;
	float Vd = 0.0f;

	float distu, distv;


	float xxx, yyy, zzz, ttt;

	xxx = partpos[i].x; //should be in i,j
	yyy = partpos[i].y;
	zzz = partpos[i].z;
	ttt = partpos[i].w;


	if (ttt >= 0.0f)
	{
		//Interpolate wter depth, Uvel Vvel at the particle position

		Ux = tex2D(texU, xxx, yyy);
		Vx = tex2D(texV, xxx + 0.5, yyy - 0.5);// U and V don't have the same coordinates but in the number of nodes it is just off by half a grid node in both dimension
		Hx = tex2D(texH, xxx, yyy);
		distu = tex2D(texdXU, xxx, yyy);
		distv = tex2D(texdYV, xxx + 0.5, yyy - 0.5);

		if (distu>0.0001 && distv>0.0001)//Avoid the div by zero which makes i and j #INF
		{
			if (Hx >= mindepth)
			{
				// old formulation
				//Xd=(dd_rand[i]*2-1)*sqrtf(6*Eh*dt);
				//Yd=(dd_rand[npart-i]*2-1)*sqrtf(6*Eh*dt);

				//formulation used in Viikmae et al.
				Xd = sqrtf(-4.0f * Eh*dt*logf(1 - dd_rand[i]))*cosf(2.0f * pi*dd_rand[npart - i]);
				Yd = sqrtf(-4.0f * Eh*dt*logf(1 - dd_rand[i]))*sinf(2.0f * pi*dd_rand[npart - i]);
				Vd = sqrtf(-4.0f * Ev*dt*logf(1 - dd_rand[i]))*cosf(2.0f * pi*dd_rand[npart - i]);

				xxx = xxx + (Ux*dt + Xd) / distu; // Need to add the runge kutta scheme here or not if time step is small enough
				yyy = yyy + (Vx*dt + Yd) / distv;
				zzz = min(max((zzz*Hx + (Vd - ws)*dt)/Hx,0.0f),1.0f);
			}
		}
	}

	ttt = ttt + dt;
	partpos[i] = make_float4(xxx, yyy, zzz, ttt);




}

__global__ void updatepartposQ3DCB(int npart, float dt, float Eh, float Ev, float mindepth, float ws, float* dd_rand, float4* partpos)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;

	float Ux = 0.0f;
	float Vx = 0.0f;
	float Hx = 1.0f;
	float Xd = 0.0f; //Diffusion term
	float Yd = 0.0f;
	float Vd = 0.0f;

	float zo = 0.3;

	float zow = 0.3;

	float Wu = 1.0;
	float Wv = 9.0;


	float facvel = 1.0f;

	float facwind = 0.0f;

	float hwind = 1.0;


	float distu, distv;


	float xxx, yyy, zzz, ttt;

	xxx = partpos[i].x; //should be in i,j
	yyy = partpos[i].y;
	zzz = partpos[i].z;
	ttt = partpos[i].w;


	if (ttt >= 0.0f)
	{
		//Interpolate wter depth, Uvel Vvel at the particle position

		Ux = tex2D(texU, xxx, yyy);
		Vx = tex2D(texV, xxx + 0.5, yyy - 0.5);// U and V don't have the same coordinates but in the number of nodes it is just off by half a grid node in both dimension
		Hx = tex2D(texH, xxx, yyy);
		distu = tex2D(texdXU, xxx, yyy);
		distv = tex2D(texdYV, xxx + 0.5, yyy - 0.5);

		if (distu > 0.0001 && distv > 0.0001)//Avoid the div by zero which makes i and j #INF
		{
			if (Hx >= mindepth)
			{
				// old formulation
				//Xd=(dd_rand[i]*2-1)*sqrtf(6*Eh*dt);
				//Yd=(dd_rand[npart-i]*2-1)*sqrtf(6*Eh*dt);

				float zr = zzz * Hx;

				facvel = log10f((zr)/zo) / log10(0.37*Hx/zo);

				if (zr > Hx - hwind)
				{
					facwind = 0.03f * log10f(max(zr-(Hx-hwind),zo) / zo) / log10(0.37 * hwind / zo);
				}

				//formulation used in Viikmae et al.
				Xd = sqrtf(-4.0f * Eh * dt * logf(1 - dd_rand[i])) * cosf(2.0f * pi * dd_rand[npart - i]);
				Yd = sqrtf(-4.0f * Eh * dt * logf(1 - dd_rand[i])) * sinf(2.0f * pi * dd_rand[npart - i]);
				Vd = sqrtf(-4.0f * Ev * dt * logf(1 - dd_rand[i])) * cosf(2.0f * pi * dd_rand[npart - i]);

				xxx = xxx + ((Ux * facvel + Wu*facwind) * dt + Xd) / distu; // Need to add the runge kutta scheme here or not if time step is small enough
				yyy = yyy + ((Vx * facvel + Wv*facwind) * dt + Yd) / distv;
				zzz = min(max((zr + (Vd - ws) * dt) / Hx, 0.0f), 1.0f);
			}
		}
	}

	ttt = ttt + dt;
	partpos[i] = make_float4(xxx, yyy, zzz, ttt);




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

__global__ void CalcNincel(int npart, int nx, int ny, float4* partpos, float *Nincel, float *cNincel, float *cTincel)
{
	int i = blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y + threadIdx.x;
	float xxx, yyy, ttt;
	int ix, iy;

	xxx = partpos[i].x;
	yyy = partpos[i].y;
	ttt = partpos[i].w;

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
