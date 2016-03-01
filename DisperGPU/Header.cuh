#include <cuda.h>
#include <cuda_runtime.h>

#define pi 3.14159265f
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