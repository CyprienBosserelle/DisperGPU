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




#include "Header.cuh"

// Define Global variables



void handle_error(int status) {
     if (status != NC_NOERR) {
        fprintf(stderr, "Netcdf %s\n", nc_strerror(status));
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


HDParam readgridsize(HDParam HD, float *&xcoord, float *&ycoord)
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
	printf("Open file\n");

	status =nc_open(HD.ncfileU.c_str(),0,&ncid);
	if (status != NC_NOERR) handle_error(status);

	//inquire variable by name
	printf("Reading information about %s...",HD.Uvarname.c_str());
	status = nc_inq_varid(ncid, HD.Uvarname.c_str(), &varid);

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


	status = nc_close(ncid);


	printf(" %s...", HD.Vvarname.c_str());

	status = nc_open(HD.ncfileV.c_str(), 0, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, HD.Vvarname.c_str(), &varid);

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

	status = nc_close(ncid);

	status = nc_open(HD.ncfileH.c_str(), 0, &ncid);
	if (status != NC_NOERR) handle_error(status);


	printf(" %s...\n", HD.Hvarname.c_str());
	status = nc_inq_varid(ncid, HD.Hvarname.c_str(), &varid);

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
	

	HD.nt = ddimV[0];
	if (ndimsU > 3){
		printf("U Variable is 4D\n");
		nlev = ddimV[1];
	}

	if (ndimshh > 2)
	{
		HD.ny = ddimhh[1];
		HD.nx = ddimhh[2];
		HD.ndim = 4;
	}
	else
	{
		HD.ny = ddimhh[0];
		HD.nx = ddimhh[1];
		HD.ndim = 3;
	}
	
	//allocate
	xcoord = (float *)malloc(HD.nx*HD.ny*sizeof(float));
	ycoord = (float *)malloc(HD.nx*HD.ny*sizeof(float));
	
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
		ytempvar = (float *)malloc(HD.ny*sizeof(float));
		size_t start[] = { 0};
		size_t count[] = { HD.ny};
		status = nc_get_vara_float(ncid, varid, start, count, ytempvar);
		if (status != NC_NOERR) handle_error(status);

		for (int i = 0; i<HD.nx; i++)
		{
			for (int j = 0; j<HD.ny; j++)
			{
				
				ycoord[i + j*HD.nx] = ytempvar[j];
				
			}
		}
	}
	else
	{
		size_t start[] = { 0, 0 };
		size_t count[] = { HD.ny, HD.nx };
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
		xtempvar = (float *)malloc(HD.nx*sizeof(float));
		size_t start[] = { 0 };
		size_t count[] = { HD.nx };
		status = nc_get_vara_float(ncid, varid, start, count, xtempvar);
		if (status != NC_NOERR) handle_error(status);

		for (int i = 0; i<HD.nx; i++)
		{
			for (int j = 0; j<HD.ny; j++)
			{

				xcoord[i + j*HD.nx] = xtempvar[i];

			}
		}
	}
	else
	{
		size_t start[] = { 0, 0 };
		size_t count[] = { HD.ny, HD.nx };
		status = nc_get_vara_float(ncid, varid, start, count, xcoord);
		if (status != NC_NOERR) handle_error(status);

	}





	status = nc_close(ncid);

	return HD;

}
void readgridsizeHYCOM(char ncfile[], char Uvar[], char Vvar[], int &nt, int &nx, int &ny, float *&xcoord, float *&ycoord)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, ndimsU, ndimsV, ndimshh, ndims;

	int varid;


	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char coordname[NC_MAX_NAME + 1];
	size_t *ddimU, *ddimV, *ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	size_t nlev;
	char ncfileU[256];
	sprintf(ncfileU, "./uvel/%s001_00_3zu.nc",ncfile);

	//Open NC file
	printf("Open file\n");
	status = nc_open(ncfileU, 0, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//inquire variable by name
	printf("Reading information aboout %s...", Uvar);
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
	status = nc_close(ncid);


	//char ncfileU[256];
	sprintf(ncfileU, "./vvel/%s001_00_3zv.nc",ncfile);

	//Open NC file
	printf("Open file\n");
	status = nc_open(ncfileU, 0, &ncid);
	if (status != NC_NOERR) handle_error(status);



	printf(" %s...", Vvar);
	status = nc_inq_varid(ncid, Vvar, &varid);
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


	

	if (ndimsU != ndimsV || ndimsU<3){
		printf("Variable dimension problems\n");
	}


	nt = ddimV[0];
	printf("ntimes in files: %d",nt);
	if (ndimsU > 3){
		printf("U Variable is 4D\n");
		nlev = ddimV[1];
		printf("nlevel in files: %d", nlev);
	}

	ny = ddimU[2];
	nx = ddimU[3];
	

	//allocate
	xcoord = (float *)malloc(nx*ny*sizeof(float));
	ycoord = (float *)malloc(nx*ny*sizeof(float));

	//inquire variable name for x dimension
	//aka x dim of hh
	int ycovar, xcovar;

	
	ycovar = dimids[2];
	xcovar = dimids[3];
	

	//ycoord
	status = nc_inq_dimname(ncid, ycovar, coordname);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, "Latitude", &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varndims(ncid, varid, &ndims);
	if (status != NC_NOERR) handle_error(status);

	if (ndims < 2)
	{
		float * ytempvar;
		ytempvar = (float *)malloc(ny*sizeof(float));
		size_t start[] = { 0 };
		size_t count[] = { ny };
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

	status = nc_inq_varid(ncid, "Longitude", &varid);
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


void readHDstepHYCOM(char ncfile[], char Uvar[], char Vvar[], int nx, int ny, int hdstep, int lev, float *&Uo, float *&Vo,float *&hho)
{
	//
	int status;
	int ncid;
	
	int uu_id, vv_id,hh_id;

	printf("Reading HD step: %d ...", hdstep);
	//size_t startl[]={hdstep-1,lev,0,0};
	//size_t countlu[]={1,1,netau,nxiu};
	//size_t countlv[]={1,1,netav,nxiv};
	size_t startl[] = { 0, lev, 0, 0 };
	size_t countlu[] = { 1, 1, ny, nx };
	size_t countlv[] = { 1, 1, ny, nx };

	//static ptrdiff_t stridel[]={1,1,1,1};
	static ptrdiff_t stridel[] = { 1, 1, 1, 1 };

	char ncfileU[256];
	sprintf(ncfileU, "./uvel/%s%3.3d_00_3zu.nc",ncfile,hdstep);

	//Open NC file
	status = nc_open(ncfileU, 0, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//status = nc_inq_varid (ncid, "u", &uu_id);
	status = nc_inq_varid(ncid, Uvar, &uu_id);
	if (status != NC_NOERR) handle_error(status);
	//status = nc_inq_varid (ncid, "v", &vv_id);
	

	status = nc_get_vara_float(ncid, uu_id, startl, countlu, Uo);
	if (status != NC_NOERR) handle_error(status);

	status = nc_close(ncid);

	
	sprintf(ncfileU, "./vvel/%s%3.3d_00_3zv.nc", ncfile, hdstep);

	//Open NC file
	status = nc_open(ncfileU, 0, &ncid);
	if (status != NC_NOERR) handle_error(status);
	status = nc_inq_varid(ncid, Vvar, &vv_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_vara_float(ncid, vv_id, startl, countlv, Vo);
	if (status != NC_NOERR) handle_error(status);
	status = nc_close(ncid);

	//Set land flag to 0.0m/s to allow particle to stick to the coast

	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			if (abs(Uo[i + j*nx])>10.0f)
			{
				Uo[i + j*nx] = 0.0f;
			}
		}
	}

	for (int i = 0; i<nx; i++)
	{
		for (int j = 0; j<ny; j++)
		{
			if (abs(Vo[i + j*nx])>10.0f)
			{
				Vo[i + j*nx] = 0.0f;
			}
		}
	}



	

	
	printf("...done\n");
}


void readHDstep(HDParam HD, int steptoread, float *&Uo, float *&Vo, float *&hho)
{
	//
	int status;
	int ncid;
	float NanValU, NanValV, NanValH;
	int uu_id, vv_id, hh_id;

	size_t *startl, *countlu, *countlv;
	ptrdiff_t *stridel;
	// step to read should be adjusted in each variables so that it keeps using the last output and teh model keeps on going
	// right now the model will catch anexception 
	printf("Reading HD step: %d ...", steptoread);
	//size_t startl[]={hdstep-1,lev,0,0};
	//size_t countlu[]={1,1,netau,nxiu};
	//size_t countlv[]={1,1,netav,nxiv};

	size_t startlhh[] = { steptoread, 0, 0 };
	size_t countlhh[] = { 1, HD.ny, HD.nx };
	
	static ptrdiff_t stridelhh[] = { 1, 1, 1 };

	if (HD.ndim == 3)
	{
		//(float *)malloc(nx*ny*sizeof(float));
		startl = (size_t *)malloc(3 * sizeof(size_t));
		countlu = (size_t *)malloc(3 * sizeof(size_t));
		countlv = (size_t *)malloc(3 * sizeof(size_t));
		stridel = (ptrdiff_t*)malloc(3 * sizeof(ptrdiff_t));

		startl[0] = steptoread;
		startl[1] = 0;
		startl[2] = 0;

		countlu[0] = 1;
		countlu[1] = HD.ny;
		countlu[2] = HD.nx;

		countlv[0] = 1;
		countlv[1] = HD.ny;
		countlv[2] = HD.nx;

		stridel[0] = 1;
		stridel[1] = 1;
		stridel[2] = 1;
		//size_t startl[] = { steptoread, 0, 0 };
		//size_t countlu[] = { 1, HD.ny, HD.nx };
		//size_t countlv[] = { 1, HD.ny, HD.nx };
		//static ptrdiff_t stridel[] = { 1, 1, 1 };
	}
	else
	{

		startl = (size_t *)malloc(4 * sizeof(size_t));
		countlu = (size_t *)malloc(4 * sizeof(size_t));
		countlv = (size_t *)malloc(4 * sizeof(size_t));
		stridel = (ptrdiff_t*)malloc(4 * sizeof(ptrdiff_t));

		startl[0] = steptoread;
		startl[1] = HD.lev;
		startl[2] = 0;
		startl[3] = 0;

		countlu[0] = 1;
		countlu[1] = 1;
		countlu[2] = HD.ny;
		countlu[3] = HD.nx;

		countlv[0] = 1;
		countlv[1] = 1;
		countlv[2] = HD.ny;
		countlv[3] = HD.nx;

		stridel[0] = 1;
		stridel[1] = 1;
		stridel[2] = 1;
		stridel[3] = 1;
		//size_t startl[] = { steptoread, HD.lev, 0, 0 };
		//size_t countlu[] = { 1, 1, HD.ny, HD.nx };
		//size_t countlv[] = { 1, 1, HD.ny, HD.nx };
		//static ptrdiff_t stridel[] = { 1, 1, 1 ,1};
	}
	

	//static ptrdiff_t stridel[]={1,1,1,1};
	

	//Open NC file

	status = nc_open(HD.ncfileU.c_str(), 0, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//status = nc_inq_varid (ncid, "u", &uu_id);
	status = nc_inq_varid(ncid, HD.Uvarname.c_str(), &uu_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_vara_float(ncid, uu_id, startl, countlu, Uo);
	if (status != NC_NOERR) handle_error(status);
	//
	//status = nc_inq_att(ncid, uu_id, "_FillValue")
	status = nc_get_att_float(ncid, uu_id, "_FillValue", &NanValU);
	if (status != NC_NOERR)
	{
		NanValU = 0.0f;
	}

	status = nc_close(ncid);


	status = nc_open(HD.ncfileV.c_str(), 0, &ncid);
	if (status != NC_NOERR) handle_error(status);
	//status = nc_inq_varid (ncid, "v", &vv_id);
	status = nc_inq_varid(ncid, HD.Vvarname.c_str(), &vv_id);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_vara_float(ncid, vv_id, startl, countlv, Vo);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_att_float(ncid, vv_id, "_FillValue", &NanValV);
	if (status != NC_NOERR)
	{
		NanValV = 0.0f;
	}

	status = nc_close(ncid);


	status = nc_open(HD.ncfileH.c_str(), 0, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, HD.Hvarname.c_str(), &hh_id);
	if (status != NC_NOERR) handle_error(status);

	
	
	status = nc_get_vara_float(ncid, hh_id, startlhh, countlhh, hho);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_att_float(ncid, hh_id, "_FillValue", &NanValH);
	if (status != NC_NOERR)
	{
		NanValH = 0.0f;
	}

	//printf("hho=%f\n", hho[10 + 330 * HD.nx]);
	status = nc_close(ncid);
	
	for (int i = 0; i<HD.nx; i++)
	{
		for (int j = 0; j<HD.ny; j++)
		{
			if (Uo[i + j*HD.nx] == NanValU)
			{
				Uo[i + j*HD.nx] = 0.0f;
			}
		}
	}

	for (int i = 0; i<HD.nx; i++)
	{
		for (int j = 0; j<HD.ny; j++)
		{
			if (Vo[i + j*HD.nx] == NanValV)
			{
				Vo[i + j*HD.nx] = 0.0f;
			}
		}
	}

	

	// Apply scale factor and offset
	//if Vscale!=1?
	for (int i = 0; i < HD.nx; i++)
	{
		for (int j = 0; j < HD.ny; j++)
		{

			Uo[i + j*HD.nx] = Uo[i + j*HD.nx] * HD.Vscale + HD.Voffset;
			Vo[i + j*HD.nx] = Vo[i + j*HD.nx] * HD.Vscale + HD.Voffset;
			//hho[i + j*HD.nx] = hho[i + j*HD.nx] * HD.Hscale + HD.Hoffset;
		}
	}

	//printf("Uo=%f\n", Uo[10 + 330 * HD.nx]);

	//Set land flag to 0.0m/s to allow particle to stick to the coast



	

	if (HD.zs2hh > 0)
	{
		// do zb
		//WARNING
		// Here I assume zb is only 2d but it could be 3d this needs to be much more flexible
		//Also reading zb every time is a bit much
		float * zb = (float *)malloc(HD.nx*HD.ny * sizeof(float));
		//size_t countlv[]={1,1,netav,nxiv};
		size_t startl[] = { 0, 0 };
		size_t countlu[] = { HD.ny, HD.nx };
		size_t countlv[] = { HD.ny, HD.nx };
		status = nc_open(HD.ncfileZB.c_str(), 0, &ncid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, HD.ZBvarname.c_str(), &hh_id);
		if (status != NC_NOERR) handle_error(status);



		status = nc_get_vara_float(ncid, hh_id, startl, countlv, zb);
		if (status != NC_NOERR) handle_error(status);

		status = nc_close(ncid);
		//hho is in fact zso so we correct this using zb
		// hh=zs-zb;
		//printf("hho=%f\tzb=%f\n", hho[10 + 330 * HD.nx], zb[10 + 330 * HD.nx]);
		for (int i = 0; i < HD.nx; i++)
		{
			for (int j = 0; j < HD.ny; j++)
			{
				if (hho[i + j*HD.nx] != NanValH)
				{
					hho[i + j*HD.nx] = (hho[i + j*HD.nx] * HD.Hscale + HD.Hoffset) - (zb[i + j*HD.nx] * HD.ZBscale + HD.ZBoffset);
					//hho[i + j*HD.nx] = max(hho[i + j*HD.nx], 0.0f);
					//printf("hho=%f\n", hho[i + j*HD.nx]);
				}
				else
				{
					hho[i + j*HD.nx] = 0.0f;
				}

			}
		}


		//printf("hho=%f\n", hho[10 + 330 * HD.nx]);
		//free zb
		free(zb);
	}
	else
	{
		for (int i = 0; i < HD.nx; i++)
		{
			for (int j = 0; j < HD.ny; j++)
			{
				if (hho[i + j*HD.nx] == NanValH)
				{
					hho[i + j*HD.nx] = 0.0f;
				}
			}
		}
	}
	
	printf("...done\n");
}

void Calcmaxstep(int nx, int ny, double &dt, double hddt, float *Uo, float *Vo, float *Un, float *Vn, float * distX, float *distY)
{
	float tmin = 99999999999.0f;
	float Velmin = 0.0000001f; // Needed to avoid dividing by zero
	float Umax;
	float Vmax;
	float Uttc,Vttc;
	float CFL = 0.9f;
	
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Umax = max(max(Uo[i + j*nx], Un[i + j*nx]), Velmin);
			Vmax = max(max(Uo[i + j*nx], Un[i + j*nx]), Velmin);

			Uttc = distX[i + j*nx] / Umax;
			Vttc = distY[i + j*nx] / Vmax;
			tmin = min(tmin, Uttc);
			tmin = min(tmin, Vttc);
		}
	}
	dt = min(max(CFL*tmin,Velmin), (float) hddt/4.0f); // make sure dt is above round off error and make sure there are going to be at least 4 step in between each HDstep.
	printf("New timestep: dt = %f\n", dt);
}


void NextstepCPU(int nx, int ny, float *&Uo, float *&Vo, float *&hho, float *Un, float *Vn,  float *hhn)
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Uo[i + j*nx] = Un[i + j*nx];
			Vo[i + j*nx] = Vn[i + j*nx];
			hho[i + j*nx] = hhn[i + j*nx];
		}
	}
}

void InterpstepCPU(int nx, int ny, int backswitch, int hdstep, float totaltime,float hddt, float *&Ux, float *Uo, float *Un)
{
	float fac = 1.0;
	float Uxo, Uxn;
	
	/*Ums[tx]=Umask[ix];*/


	if (backswitch>0)
	{
		fac = -1.0f;
	}


	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Uxo = fac*Uo[i + nx*j];
			Uxn = fac*Un[i + nx*j];

			Ux[i + nx*j] = Uxo + (totaltime - hddt*hdstep)*(Uxn - Uxo) / hddt;
		}
	}
}

float interp2posCPU(int nx, int ny, float x, float y, float *Ux)
{
	int x1,x2;
	int y1,y2;

	x1 = floor(x);
	x2 = x1+1;

	y1 = floor(y);
	y2 = y1+1;

	float Q11 = Ux[x1 + y1*nx];
	float Q12 = Ux[x1 + y2*nx];
	float Q21 = Ux[x2 + y1*nx];
	float Q22 = Ux[x2 + y2*nx];

	float R1 = ((x2 - x) / (x2 - x1))*Q11 + ((x - x1)/(x2 - x1))*Q21;


	float R2 = ((x2 - x) / (x2 - x1))*Q12 + ((x - x1) / (x2 - x1))*Q22;

	//	After the two R values are calculated, the value of P can finally be calculated by a weighted average of R1 and R2.

	float P = ((y2 - y) / (y2 - y1))*R1 + ((y - y1) / (y2 - y1))*R2;

	return P;
}

void updatepartposCPU(int nx, int ny, int np, float dt, float Eh, float *Ux, float *Vx, float *hhx,float *distX, float *distY, float4 *&partpos)
{
	float xxx,yyy,zzz,ttt;
	float Up, Vp, hhp;
	float ddx, ddy;
	float Xd, Yd;
	float rna, rnb, rnc, rnd; // to store random numbers

	//float xxn, yyn;
	

	for (int p = 0; p < np; p++)
	{
		xxx = partpos[p].x; //should be in i,j
		yyy = partpos[p].y;
		zzz = partpos[p].z;
		ttt = partpos[p].w;
		if (ttt >= 0.0)
		{
			Up = interp2posCPU(nx, ny, xxx, yyy, Ux);
			Vp = interp2posCPU(nx, ny, xxx, yyy, Vx);
			hhp = interp2posCPU(nx, ny, xxx, yyy, hhx);
			ddx = interp2posCPU(nx, ny, xxx, yyy, distX);
			ddy = interp2posCPU(nx, ny, xxx, yyy, distY);

			rna = ((float)rand() / (float)(RAND_MAX));
			rnb = ((float)rand() / (float)(RAND_MAX));
			rnc = ((float)rand() / (float)(RAND_MAX));
			rnd = ((float)rand() / (float)(RAND_MAX));
			

			Xd = sqrtf(-4.0f * Eh*dt*logf(1.0f - rna))*cosf(2.0f * pi*(rnb));
			Yd = sqrtf(-4.0f * Eh*dt*logf(1.0f - rnc))*sinf(2.0f * pi*(rnd));

			xxx = xxx + (Up*dt + Xd) / ddx; // Need to add the runge kutta scheme here or not
			yyy = yyy + (Vp*dt + Yd) / ddy;
			
			
		}
		partpos[p] = make_float4(xxx, yyy, zzz, ttt);


	}
}

void calcNincelCPU(int np, int nx,int ny, float4 * partpos,float * Nincel, float * cNincel, float *cTincel)
{
	//
	for (int p = 0; p < np; p++)
	{
		int xi = floor(partpos[p].x);
		int yj = floor(partpos[p].y);
		Nincel[xi + nx*yj] = Nincel[xi + nx*yj] + 1;
		cNincel[xi + nx*yj] = cNincel[xi + nx*yj] + 1;
		cTincel[xi + nx*yj] = cTincel[xi + nx*yj] + partpos[p].w;
	}
}

void resetNincelCPU(int nx, int ny, float * Nincel)
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			Nincel[i + nx*j] = 0.0f;
		}
	}
}