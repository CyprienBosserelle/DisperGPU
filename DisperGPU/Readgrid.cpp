//////////////////////////////////////////////////////////////////////////////////
//DispROMS_GPU                                                                    //
//Copyright (C) 2013 Bosserelle                                                 //
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



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <fstream>
#include <netcdf.h>
#include "Header.cuh"

// Define Global variables



void handle_error(int status) {
     if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
		//fprintf(logfile, "Netcdf: %s\n", nc_strerror(status));
        exit(-1);
        }
     }

template <class T> const T& max (const T& a, const T& b) {
  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

template <class T> const T& round(const T& a)
  {
  return floor( a + 0.5 );
  }

void readgridsize(char ncfile[], char Uvar[], char Vvar[], char hhvar[],int &nt, int &nx, int &ny, float *&xcoord, float *&ycoord)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, ndimsU, ndimsV, ndimshh,ndims;
	
	int varid;
	
	
	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char coordname[NC_MAX_NAME + 1];
	size_t *ddimU, *ddimV,*ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	size_t nlev;


	//Open NC file
	status =nc_open(ncfile,0,&ncid);
	if (status != NC_NOERR) handle_error(status);

	//inquire variable by name
	//printf("nvars=%s\n",Uvar);
	status = nc_inq_varid(ncid, Uvar, &varid);
	if (status != NC_NOERR) 
		handle_error(status);
	
	status = nc_inq_varndims(ncid, varid, &ndimsU);
	if (status != NC_NOERR) handle_error(status);
	

	status = nc_inq_vardimid(ncid, varid, dimids);
	if (status != NC_NOERR) handle_error(status);

	ddimU = (size_t *)malloc(ndimsU*sizeof(size_t));

	//Read dimensions nx_u ny_u 
	for (int iddim = 0; iddim < ndimsU; iddim++)
	{
		status = nc_inq_dimlen(ncid, dimids[iddim], &ddimU[iddim]);
		if (status != NC_NOERR) handle_error(status);

		//printf("dim:%d=%d\n", iddim, ddimU[iddim]);
	}

	status = nc_inq_varid(ncid, Uvar, &varid);
	if (status != NC_NOERR)
		handle_error(status);
	


	status = nc_inq_varndims(ncid, varid, &ndimsV);
	if (status != NC_NOERR) handle_error(status);
	//printf("VVar:%d dims\n", ndimsV);

	status = nc_inq_vardimid(ncid, varid, dimids);
	if (status != NC_NOERR) handle_error(status);

	ddimV = (size_t *)malloc(ndimsV*sizeof(size_t));

	//Read dimensions nx_u ny_u 
	for (int iddim = 0; iddim < ndimsV; iddim++)
	{
		status = nc_inq_dimlen(ncid, dimids[iddim], &ddimV[iddim]);
		if (status != NC_NOERR) handle_error(status);

		//printf("dim:%d=%d\n", iddim, ddimV[iddim]);
	}



	status = nc_inq_varid(ncid, hhvar, &varid);
	if (status != NC_NOERR)
		handle_error(status);
	


	status = nc_inq_varndims(ncid, varid, &ndimshh);
	if (status != NC_NOERR) handle_error(status);
	//printf("hhVar:%d dims\n", ndimshh);

	status = nc_inq_vardimid(ncid, varid, dimids);
	if (status != NC_NOERR) handle_error(status);

	ddimhh = (size_t *)malloc(ndimshh*sizeof(size_t));

	//Read dimensions nx_u ny_u 
	for (int iddim = 0; iddim < ndimshh; iddim++)
	{
		status = nc_inq_dimlen(ncid, dimids[iddim], &ddimhh[iddim]);
		if (status != NC_NOERR) handle_error(status);

		//printf("dim:%d=%d\n", iddim, ddimhh[iddim]);
	}

	if (ndimsU != ndimsV || ndimsU<3){
		printf("Variable dimension problems\n");
	}
	

	nt = ddimV[0];
	if (ndimsU > 3){
		printf("U Variable is 4D\n");
		nlev = ddimV[1];
	}

	if (ndimshh > 2)
	{
		ny = ddimhh[1];
		nx = ddimhh[2];
	}
	else
	{
		ny = ddimhh[0];
		nx = ddimhh[1];
	}
	
	//allocate
	xcoord = (float *)malloc(nx*ny*sizeof(float));
	ycoord = (float *)malloc(nx*ny*sizeof(float));
	
	//inquire variable name for x dimension
	//aka x dim of hh
	int ycovar,xcovar;
	
	if (ndimshh > 2)
	{
		ycovar = dimids[1];
		xcovar = dimids[2];
	}
	else
	{
		ycovar = dimids[0];
		xcovar = dimids[1];
	}

	//ycoord
	status = nc_inq_dimname(ncid, ycovar, coordname);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, coordname, &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varndims(ncid, varid, &ndims);
	if (status != NC_NOERR) handle_error(status);

	if (ndims < 2)
	{
		float * ytempvar;
		ytempvar = (float *)malloc(ny*sizeof(float));
		size_t start[] = { 0};
		size_t count[] = { ny};
		status = nc_get_vara_float(ncid, varid, start, count, ytempvar);
		if (status != NC_NOERR) handle_error(status);

		for (int i = 0; i<nx; i++)
		{
			for (int j = 0; j<ny; j++)
			{
				
				ycoord[i + j*nx] = ytempvar[j];
				
			}
		}
	}
	else
	{
		size_t start[] = { 0, 0 };
		size_t count[] = { ny, nx };
		status = nc_get_vara_float(ncid, varid, start, count, ycoord);
		if (status != NC_NOERR) handle_error(status);
			
	}
	//xcoord
	status = nc_inq_dimname(ncid, xcovar, coordname);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, coordname, &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varndims(ncid, varid, &ndims);
	if (status != NC_NOERR) handle_error(status);

	if (ndims < 2)
	{
		float * xtempvar;
		xtempvar = (float *)malloc(nx*sizeof(float));
		size_t start[] = { 0 };
		size_t count[] = { nx };
		status = nc_get_vara_float(ncid, varid, start, count, xtempvar);
		if (status != NC_NOERR) handle_error(status);

		for (int i = 0; i<nx; i++)
		{
			for (int j = 0; j<ny; j++)
			{

				xcoord[i + j*nx] = xtempvar[i];

			}
		}
	}
	else
	{
		size_t start[] = { 0, 0 };
		size_t count[] = { ny, nx };
		status = nc_get_vara_float(ncid, varid, start, count, xcoord);
		if (status != NC_NOERR) handle_error(status);

	}





	status = nc_close(ncid);


}


void readHDstep(char ncfile[], char Uvar[], char Vvar[], char hhvar[], int nx, int ny, int hdstep, int lev, float *&Uo, float *&Vo,float *&hho)
{
	//
	int status;
	int ncid;
	
	int uu_id, vv_id,hh_id;

	printf("Reading HD step: %d ...", hdstep);
	//size_t startl[]={hdstep-1,lev,0,0};
	//size_t countlu[]={1,1,netau,nxiu};
	//size_t countlv[]={1,1,netav,nxiv};
	size_t startl[] = { hdstep, 0, 0 };
	size_t countlu[] = { 1, ny, nx };
	size_t countlv[] = { 1, ny, nx };

	//static ptrdiff_t stridel[]={1,1,1,1};
	static ptrdiff_t stridel[] = { 1, 1, 1 };

	//Open NC file
	status = nc_open(ncfile, 0, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//status = nc_inq_varid (ncid, "u", &uu_id);
	status = nc_inq_varid(ncid, Uvar, &uu_id);
	if (status != NC_NOERR) handle_error(status);
	//status = nc_inq_varid (ncid, "v", &vv_id);
	status = nc_inq_varid(ncid, Vvar, &vv_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, hhvar, &hh_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_vara_float(ncid, uu_id, startl, countlu, Uo);
	if (status != NC_NOERR) handle_error(status);
	status = nc_get_vara_float(ncid, vv_id, startl, countlv, Vo);
	if (status != NC_NOERR) handle_error(status);
	status = nc_get_vara_float(ncid, hh_id, startl, countlv, hho);
	if (status != NC_NOERR) handle_error(status);

	//Set land flag to 0.0m/s to allow particle to stick to the coast

	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			if (Uo[i + j*nx]>10.0f)
			{
				Uo[i + j*nx] = 0.0f;
			}
		}
	}

	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			if (Vo[i + j*nx]>10.0f)
			{
				Vo[i + j*nx] = 0.0f;
			}
		}
	}



	

	status = nc_close(ncid);
	printf("...done\n");
}
