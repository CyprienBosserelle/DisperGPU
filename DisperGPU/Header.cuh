
#ifndef DISPERGPU_H
#define DISPERGPU_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
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


#include <time.h>
#include <string>

#include <sstream>
#include <iterator>
#include <netcdf.h>
#include <algorithm>
#include <vector>
#include <ctime>
#include <random>
#include <chrono>

#define pi 3.14159265f

// Class definitions
class Param{
public:
	//Files
	std::string ncfile; //Should use strings here
	std::string Uvarname;
	std::string Vvarname;
	std::string hhvarname;
	std::string ncoutfile;
	std::string seedfile;


	int np;//Numer of particles
	int partmode; // Particle model type: 0:test set up the model but do not run any loop
	// 1: 2D, passive particle no buyancy taken into acount
	// 2: Quasi 3D, buoyant/sinking particle advected with 2d depth averaged hydrodynamics, Log profile is assumed
	// 3D: 3D model buoyant/sinking particle advected with 3d hydrodynamcis
	
	int nx, ny, nz, nt; //HD input may have a more complex structure with staggered grid 

	
	double hddt; // HD model tme step
	int lev; //Level for 3D HD but 2D/Q3D particle model
	int geocoord; //Geographic coordinate system switch 0 is metric 1 is degrees
	int backswitch; // 0 run HD model forward 1 run the model backward
	double Eh, Ev; // Eddy viscosity horizontale, vertical
	double minrwdepth; // Minimum depth for using Eddy viscosity

	int GPUDEV = 0; // GPU device in use  (default is 0, aka first available device from device query) negative value means force use of cpu and other positive value a dpecific GPU in a multi GPU system



};

class Control{
public:
	float totaltime; // needed to track total time as dt can vary
	float nextouttime;
	float outtime;
	int stp, outstep, outtype; // Model step, output step, next output step, output file type

	float dt, olddt; // particle model time step
	int hdstep, hdstart, hdend; // HD model step, HD step start and HD step end
	int SEED = 777; //Seed for random number generator



};


// Shared functions
template <class T> const T& min(const T& a, const T& b);
template <class T> const T& max(const T& a, const T& b);
template <class T> const T& round(const T& a);

void readgridsize(char ncfile[], char Uvar[], char Vvar[], char hhvar[], int &nt, int &nx, int &ny, float *&xcoord, float *&ycoord);

void readgridsizeHYCOM(char ncfile[], char Uvar[], char Vvar[], int &nt, int &nx, int &ny, float *&xcoord, float *&ycoord);
     

void readHDstep(char ncfile[], char Uvar[], char Vvar[], char hhvar[], int nx, int ny, int hdstep, int lev, float *&Uo, float *&Vo, float *&hho);

void readHDstepHYCOM(char ncfile[], char Uvar[], char Vvar[], int nx, int ny, int hdstep, int lev, float *&Uo, float *&Vo, float *&hho);

void CalcDistXY(int nx, int ny, int geocoord, float *xcoord, float * ycoord, float * &distX, float *&distY);

void Calcmaxstep(int nx, int ny, float &dt, float hddt, float *Uo, float *Vo, float *Un, float *Vn, float * distX, float *distY);
void NextstepCPU(int nx, int ny, float *&Uo, float *&Vo, float *&hho, float *Un, float *Vn, float *hhn);
void InterpstepCPU(int nx, int ny, int backswitch, int hdstep, float totaltime, float hddt, float *&Ux, float *Uo, float *Un);
float interp2posCPU(int nx, int ny, float x, float y, float *Ux);

void updatepartposCPU(int nx, int ny, int np, float dt, float Eh, float *Ux, float *Vx, float *hhx, float *distX, float *distY, float4 *&partpos);
void xyz2ijk( int nx, int ny, float * xcoord, float * ycoord, float &xi, float &yj, float xreal, float yreal, float zreal);

bool isinquad(float v1x, float v1y, float v2x, float v2y, float v3x, float v3y, float v4x, float v4y, float px, float py);
extern "C" void readseedfile(char seedfile[], int npart, int nx, int ny, float *xcoord, float *ycoord, float4* &partpos);
void writexyz(int npart, int nx, int ny, float * xcoord, float * ycoord, float4 * partpos, char outfile[]);
void calcNincelCPU(int np, int nx, int ny, float4 * partpos, float * Nincel, float * cNincel, float *cTincel);
void resetNincelCPU(int nx, int ny, float * Nincel);

void creatncfile(char outfile[], int nx, int ny, int np, float *xval, float *yval, float totaltime, float *Nincel, float *cNincel, float *cTincel, float4 * PartPos);
void writestep2nc(char outfile[], int nx, int ny, int np, float totaltime, float *xval, float *yval, float *Nincel, float *cNincel, float * cTincel, float4 *PartPos);

float isLeft(float P0x, float P0y, float P1x, float P1y, float P2x, float P2y);
int cn_PnPoly(float Px, float Py, float* Vx, float *Vy, int n);
int wn_PnPoly(float Px, float Py, float* Vx, float* Vy, int n);


// End of global definition
#endif
