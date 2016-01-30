#include <cuda.h>
#include <cuda_runtime.h>

#define pi 3.14159265f
template <class T> const T& min(const T& a, const T& b);
template <class T> const T& max(const T& a, const T& b);
template <class T> const T& round(const T& a);

void readgridsize(char ncfile[], char Uvar[], char Vvar[], char hhvar[], int &nt, int &nx, int &ny, float *&xcoord, float *&ycoord);
void readHDstep(char ncfile[], char Uvar[], char Vvar[], char hhvar[], int nx, int ny, int hdstep, int lev, float *&Uo, float *&Vo, float *&hho);
void CalcDistXY(int nx, int ny, int geocoord, float *xcoord, float * ycoord, float * &distX, float *&distY);
extern "C" void readseedfile(char seedfile[], int npart, float4* &partpos);
