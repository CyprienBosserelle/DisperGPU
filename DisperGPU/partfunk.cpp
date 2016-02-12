#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <fstream>
#include <netcdf.h>
#include "Header.cuh"



void CalcDistXY(int nx, int ny, int geocoord, float *xcoord, float * ycoord, float * &distX,float *&distY)
{
	//Calculate the distance in meter between each cells for u grid 
	float R = 6372797.560856f;

	for (int i = 0; i < nx - 1; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			if (geocoord == 1)
			{
				//calc distance between each i using haversine formula
				float dlat = (ycoord[(i + 1) + j*nx] - ycoord[i + j*nx])*pi / 180.0f;
				float lat1 = ycoord[i + j*nx] * pi / 180.0f;
				float lat2 = ycoord[(i + 1) + j*nx] * pi / 180.0f;
				float dlon = (xcoord[(i + 1) + j*nx] - xcoord[i + j*nx])*pi / 180.0f;

				float a = sin(dlat / 2)*sin(dlat / 2) + cos(lat1)*cos(lat2)*sin(dlon / 2)*sin(dlon / 2);
				float c = 2 * atan2f(sqrtf(a), sqrtf(1 - a));
				distX[i + j*nx] = c*R;
			}
			else{
				distX[i + j*nx] = xcoord[(i + 1) + j*nx] - xcoord[i + j*nx];
			}


		}
	}
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny - 1; j++)
		{
			if (geocoord == 1)
			{
				//calc distance between each j using haversine formula
				float dlat = (ycoord[i + (j + 1)*nx] - ycoord[i + j*nx])*pi / 180.0f;
				float lat1 = ycoord[i + j*nx] * pi / 180.0f;
				float lat2 = ycoord[i + (j + 1)*nx] * pi / 180.0f;
				float dlon = (xcoord[i + (j + 1)*nx] - xcoord[i + j*nx])*pi / 180.0f;

				float a = sin(dlat / 2)*sin(dlat / 2) + cos(lat1)*cos(lat2)*sin(dlon / 2)*sin(dlon / 2);
				float c = 2 * atan2f(sqrtf(a), sqrtf(1 - a));
				distY[i + j*nx] = c*R;
			}
			else{
				distY[i + j*nx] = ycoord[i + (j + 1)*nx] - ycoord[i + j*nx];
			}
		}
	}

	//fill in boundaries
	for (int j = 0; j < ny; j++)
	{
		//
		distX[nx - 1 + j*nx] = distX[nx - 2 + j*nx];
	}
	for (int i = 0; i < nx; i++)
	{
		//
		distY[i + (ny - 1)*nx] = distY[i + (ny - 2)*nx];
	}
}

void xyz2ijk(int nx, int ny,float * xcoord, float * ycoord, float &xi,float &yj, float xreal, float yreal, float zreal)
{
	float dist;
	//for (int p = 0; p < np; p++) // run the function rather than the loop inside it...
	//{
		float px = xreal;
		float py = yreal;
		float mindist = hypot(xcoord[0] - px, ycoord[0] - py);
		int minii = 0;
		int minjj = 0;
		//First find the closest node to the particle
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				dist = hypot(xcoord[i+nx*j] - px, ycoord[i+nx*j] - py);
				if (dist <= mindist)
				{
					mindist = dist;
					minii = i;
					minjj = j;
				}
			}
		}

		// find out which quad contain the particle
		//EV[i] = V[i+1] - V[i], where V[] - vertices in order
		//PV[i] = P - V[i]
		//Cross[i] = CrossProduct(EV[i], PV[i]) = EV[i].X * PV[i].Y - EV[i].Y * PV[i].X
		
		//first quad is the one to the top right
		bool quadtest = false;
		float v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y;
		float pxi, pxj;
		if (minii<nx-1 && minjj<ny-1)
		{
			v1x = xcoord[minii + minjj*nx];
			v1y = ycoord[minii + minjj*nx];
			v2x = xcoord[minii+1 + minjj*nx];
			v2y = ycoord[minii+1 + minjj*nx];
			v3x = xcoord[minii+1 + (minjj+1)*nx];
			v3y = ycoord[minii+1 + (minjj+1)*nx];
			v4x = xcoord[minii + (minjj+1)*nx];
			v4y = ycoord[minii + (minjj+1)*nx];
			pxi = minii*1.0f;//int to float
			pxj = minjj*1.0f;

			quadtest = isinquad(v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, px, py);
		}
		if (quadtest == false && minii < nx - 1 && minjj>0)//bottom right quad
		{
			v4x = xcoord[minii + minjj*nx];
			v4y = ycoord[minii + minjj*nx];
			v3x = xcoord[minii + 1 + minjj*nx];
			v3y = ycoord[minii + 1 + minjj*nx];
			v2x = xcoord[minii + 1 + (minjj - 1)*nx];
			v2y = ycoord[minii + 1 + (minjj - 1)*nx];
			v1x = xcoord[minii + (minjj - 1)*nx];
			v1y= ycoord[minii + (minjj - 1)*nx];

			pxi = minii*1.0f;
			pxj = (minjj-1)*1.0f;
			quadtest = isinquad(v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, px, py);

		}
		if (quadtest == false && minii > 0 && minjj>0)//bottom left quad
		{
			v3x = xcoord[minii + minjj*nx];
			v3y = ycoord[minii + minjj*nx];
			v4x = xcoord[minii - 1 + minjj*nx];
			v4y = ycoord[minii - 1 + minjj*nx];
			v1x = xcoord[minii - 1 + (minjj - 1)*nx];
			v1y = ycoord[minii - 1 + (minjj - 1)*nx];
			v2x = xcoord[minii + (minjj - 1)*nx];
			v2y = ycoord[minii + (minjj - 1)*nx];
			pxi = (minii-1)*1.0f;
			pxj = (minjj-1)*1.0f;
			quadtest = isinquad(v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, px, py);

		}
		if (quadtest == false && minii > 0 && minjj<ny-1)//top left quad
		{
			v2x = xcoord[minii + minjj*nx];
			v2y = ycoord[minii + minjj*nx];
			v1x = xcoord[minii - 1 + minjj*nx];
			v1y = ycoord[minii - 1 + minjj*nx];
			v4x = xcoord[minii - 1 + (minjj + 1)*nx];
			v4y = ycoord[minii - 1 + (minjj + 1)*nx];
			v3x = xcoord[minii + (minjj + 1)*nx];
			v3y = ycoord[minii + (minjj + 1)*nx];
			pxi = (minii-1)*1.0f;
			pxj = (minjj)*1.0f;
			quadtest = isinquad(v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, px, py);

		}
		
		//v1 v2 v3 and v4 are the coordinates of the quad that contain the particle
		//starting with v1 the botom left corner (in i,j coordinates)
		
		// inverse bilinear to obtain the di dj

		float qx = 0.5;
		float qy = 0.5;
		float rx = 1;
		float ry = 1;
		float s, t, Jsx, Jsy, Jtx, Jty, jtjxa, jtjxb, jtjya, jtjyb;
		float jtrx, jtry, Ijtjxa, Ijtjxb, Ijtjya, Ijtjyb,Corrx,Corry;
		while (rx >= 0.000001 && ry >= 0.000001)
		{
			s = qx;
			t = qy;
			//r = p1*(1 - s)*(1 - t) + p2*s*(1 - t) + p3*s*t + p4*(1 - s)*t - p; %residual
			rx = v1x*(1 - s)*(1 - t) + v2x*s*(1 - t) + v3x*s*t + v4x*(1 - s)*t - px;
			ry = v1y*(1 - s)*(1 - t) + v2y*s*(1 - t) + v3y*s*t + v4y*(1 - s)*t - py;

			//Js = -p1*(1 - t) + p2*(1 - t) + p3*t - p4*t; %dr / ds
			Jsx = -v1x*(1 - t) + v2x*(1 - t) + v3x*t - v4x*t;
			Jsy = -v1y*(1 - t) + v2y*(1 - t) + v3y*t - v4y*t;


			//Jt = -p1*(1 - s) - p2*s + p3*s + p4*(1 - s); %dr / dt
			Jtx = -v1x*(1 - s) - v2x*s + v3x*s + v4x*(1 - s);
			Jty = -v1y*(1 - s) - v2y*s + v3y*s + v4y*(1 - s);

			//J = [Js, Jt];
			jtjxa = Jsx*Jsx + Jsy*Jsy;
			jtjxb = Jsx*Jtx + Jsy*Jty;
			jtjya = Jsy*Jty + Jsx*Jtx;
			jtjyb = Jtx*Jtx + Jty*Jty;


			jtrx = Jsx*rx + Jsy*ry;
			jtry = Jtx*rx + Jty*ry;


			Ijtjxa = 1 / (jtjxa*jtjyb - jtjxb*jtjya)*jtjyb;
			Ijtjxb = 1 / (jtjxa*jtjyb - jtjxb*jtjya)*-1 * jtjxb;
			Ijtjya = 1 / (jtjxa*jtjyb - jtjxb*jtjya)*-1 * jtjya;
			Ijtjyb = 1 / (jtjxa*jtjyb - jtjxb*jtjya)*jtjxa;

			Corrx = Ijtjxa*jtrx + Ijtjxb*jtry;
			Corry = Ijtjya*jtrx + Ijtjyb*jtry;

			qx = qx - Corrx;
			qy = qy - Corry;


		}
		xi = pxi + qx;
		yj = pxj + qy;
		// Now deal with z 
		//needs to be sorted out later !!!! for now it is assumed to be all independant of z anyway


	
	
}


bool isinquad(float v1x, float v1y, float v2x, float v2y, float v3x, float v3y, float v4x, float v4y, float px, float py)
{
	float EV1x = v2x - v1x;
	float EV1y = v2y - v1y;

	float PV1x = px - v1x;
	float PV1y = py - v1y;

	float EV2x = v3x - v2x;
	float EV2y = v3y - v2y;

	float PV2x = px - v2x;
	float PV2y = py - v2y;

	float EV3x = v4x - v3x;
	float EV3y = v4y - v3y;

	float PV3x = px - v3x;
	float PV3y = py - v3y;

	float EV4x = v1x - v4x;
	float EV4y = v1y - v4y;

	float PV4x = px - v4x;
	float PV4y = py - v4y;

	float Xs1 = EV1x*PV1y - EV1y*PV1x;
	float Xs2 = EV2y*PV2y - EV2y*PV2x;
	float Xs3 = EV3y*PV3y - EV3y*PV3x;
	float Xs4 = EV4y*PV4y - EV4y*PV4x;


	if (Xs1 < 0.0 && Xs2 < 0.0 && Xs3 < 0.0 && Xs4 < 0.0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//void WriteoutCPU(float4 * partpos )
//{
	//
//}
