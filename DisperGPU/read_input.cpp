//////////////////////////////////////////////////////////////////////////////////

//DispROMS_GPU                                                                    //
//Copyright (C) 2018 Bosserelle                                                 //

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




Param readparamstr(std::string line, Param param)
{


	std::string parameterstr, parametervalue;

	///////////////////////////////////////////////////////
	// General parameters

	parameterstr = "ncoutfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncoutfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}
	
	//

	parameterstr = "seedfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{

		param.seedfile = parametervalue;
	}
		
	
	//
	parameterstr = "np";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.np = std::stoi(parametervalue);
	}

	//
	parameterstr = "partmode";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.partmode = std::stoi(parametervalue);
	}

	//
	parameterstr = "backswitch";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.backswitch = std::stoi(parametervalue);
	}

	//
	parameterstr = "GPUDEV";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.GPUDEV = std::stoi(parametervalue);
	}

	//
	parameterstr = "gpudevice";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.GPUDEV = std::stoi(parametervalue);
	}

	parameterstr = "gpu";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.GPUDEV = std::stoi(parametervalue);
	}

	///////////////////////////////////////////////////////
	// Flow parameters
	//
	parameterstr = "Eh";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Eh = std::stod(parametervalue);
	}
	
	parameterstr = "Ev";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Ev = std::stod(parametervalue);
	}

	parameterstr = "minrwdepth";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.minrwdepth = std::stod(parametervalue);
	}

	parameterstr = "outtime";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.outtime = std::stod(parametervalue);
	}

	return param;
}



HDParam readHDparamstr(std::string line, HDParam param)
{


	std::string parameterstr, parametervalue;

	///////////////////////////////////////////////////////
	// General parameters
	parameterstr = "ncfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncfile = parametervalue;
		param.ncfileU = parametervalue;
		param.ncfileV = parametervalue;
		param.ncfileH = parametervalue;
		param.ncfileZB = parametervalue;
		
	}

	parameterstr = "ncfileU";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncfileU = parametervalue;

		
	}

	parameterstr = "ncfileV";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncfileV = parametervalue;

		
	}

	parameterstr = "ncfileH";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncfileH = parametervalue;

		
	}

	parameterstr = "ncfileZB";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncfileZB = parametervalue;


	}

	//
	parameterstr = "Uvarname";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Uvarname = parametervalue;
	}

	parameterstr = "Vvarname";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Vvarname = parametervalue;
	}

	parameterstr = "Hvarname";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Hvarname = parametervalue;
	}

	parameterstr = "ZBvarname";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ZBvarname = parametervalue;
	}

	//
	parameterstr = "hdstart";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.hdstart = std::stoi(parametervalue);
	}

	//
	parameterstr = "hdend";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.hdend = std::stoi(parametervalue);
	}

	//
	parameterstr = "lev";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.lev = std::stoi(parametervalue);
	}

	//
	parameterstr = "geocoord";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.geocoord = std::stoi(parametervalue);
	}

	//
	
	parameterstr = "hddt";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.hddt = std::stod(parametervalue);
	}


	parameterstr = "Vscale";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Vscale = std::stod(parametervalue);
	}

	parameterstr = "Voffset";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Voffset = std::stod(parametervalue);
	}

	parameterstr = "Hscale";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Hscale = std::stod(parametervalue);
	}
	
	parameterstr = "Hoffset";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.Hoffset = std::stod(parametervalue);
	}

	parameterstr = "ZBscale";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ZBscale = std::stod(parametervalue);
	}

	parameterstr = "ZBoffset";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ZBoffset = std::stod(parametervalue);
	}

	parameterstr = "zs2hh";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.geocoord = std::stoi(parametervalue);
	}


	return param;
}






std::string findparameter(std::string parameterstr, std::string line)
{
	std::size_t found, Numberstart, Numberend;
	std::string parameternumber,left,right;
	std::vector<std::string> splittedstr;
	
	// first look fo an equal sign
	// No equal sign mean not a valid line so skip
	splittedstr=split(line, '=' );
	if (splittedstr.size()>1)
	{
		left = trim(splittedstr[0]," ");
		right = splittedstr[1]; // if there are more than one equal sign in the line the second one is ignored
		found = left.compare(parameterstr);// it needs to strictly compare
		if (found == 0) // found the parameter
		{
			//std::cout <<"found LonMin at : "<< found << std::endl;
			//Numberstart = found + parameterstr.length();
			splittedstr = split(right, ';');
			if (splittedstr.size() >= 1)
			{
				parameternumber = splittedstr[0];
			}
			//std::cout << parameternumber << std::endl;

		}
	}
	return trim(parameternumber, " ");
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (!item.empty())//skip empty tokens
		{
			elems.push_back(item);
		}

	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

std::string trim(const std::string& str, const std::string& whitespace)
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return ""; // no content

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}



void write_text_to_log_file(std::string text)
{
	std::ofstream log_file(
		"DGPU_log.txt", std::ios_base::out | std::ios_base::app);
	log_file << text << std::endl;
	log_file.close(); //destructor implicitly does it
}

void SaveParamtolog(Param Dparam, HDParam HD)
{
	write_text_to_log_file("#################################");
	write_text_to_log_file("# Output file");
	write_text_to_log_file("ncoutfile = " + Dparam.ncoutfile + ";");
		write_text_to_log_file("\n");
	write_text_to_log_file("# Model controls");
	
}



