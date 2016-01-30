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
