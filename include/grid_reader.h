/********************************************/
/* This file is distributed under the GNU   */
/* LPGL-2.1-or-later Open Source License.   */
/* See LICENSE file for details.            */
/*                                          */
/* Copyright (c) 2023 Andreas F. Tillack    */
/*                    Althea A. Hansel      */
/*                    Matthew Holcomb       */
/*           Forli Lab @ Scripps Research   */
/********************************************/


#ifndef INCLUDED_GRID_READER
#define INCLUDED_GRID_READER

#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <atomic>
#include <algorithm>

using namespace std;

inline float map2float(const char* c)
// This function converts what we typically find in an autogrid map file into a
// floating point number - just a bit quicker than the usual sscanf()
// -> due to using 32-bit integers this function is limited to 9 digits for both
//    the whole number and the fractional part - a safety check is performed with
//    sscanf() used as the fallback
{
	float result;
	bool negative = false;                       // example: -123.456
	if(*c == '-'){                               // *c = '-'
		negative = true;                     // => negative = true
		c++;
	}
	// safety check
	int len = strlen(c);
	if(len>9){ // no potential issues at or below 9 digits in total
		const char* dp = strchr(c,'.');
		if(dp){
			int d = dp-c;
			if((d>9) || (len-d>9)){ // fall back to sscanf() if numbers are going to be too big for integers
				sscanf(c, "%f", &result);
				if(negative) return -result;
				return result;
			}
		}
	}
	int number = 0;                              // 1. *c = '1': number = 0*10  + 1 = 1
	while((*c >= '0') && (*c <= '9')){           // 2. *c = '2': number = 1*10  + 2 = 12
		number = number * 10 + (*c - '0');   // 3. *c = '3': number = 12*10 + 3 = 123
		c++;                                 // 4. *c = ','
	}
	if(*c == '.') c++; // jump over decimal point
	int decimal = 0;
	int denom = 1;
	while((*c >= '0') && (*c <= '9')){           // 1. *c = '4': decimal = 0*10  + 4 = 4,   denom = 10
		decimal = decimal * 10 + (*c - '0'); // 2. *c = '5': decimal = 4*10  + 5 = 45,  denom = 100
		denom *= 10;                         // 3. *c = '6': decimal = 45*10 + 6 = 456, denom = 1000
		c++;
	}
	// use more expensive division only once
	result = (float)number + (float)decimal/((float)denom);
	if(negative) return -result;
	return result;
}

inline fp_num* read_grid_map(string filename, int &sizeX, int &sizeY, int &sizeZ, fp_num* map_storage = NULL)
{
	unsigned long size;
	ifstream file(filename.c_str(),ifstream::in);
	if (file.fail()==true){
		cout << "Could not open grid map file \"" << filename << "\".\n";
		exit(1);
	}
	// Get file size
	file.seekg(0,ifstream::end);
	size=file.tellg();
	file.seekg(0);
	// Sanity checks
	if (size==0){
		cout << "File \"" << filename << "\" is empty.\n";
		exit(1);
	}
	// Read content
	string line;
	int xnr, ynr, znr;
	unsigned int found_all = 0;
	fp_num center_x, center_y, center_z, spacing;
	do
	{
		std::getline(file, line);
		// capturing number of grid points
		if (line.find("NELEMENTS") == 0)
		{
			sscanf(&line.c_str()[10], "%d %d %d", &xnr, &ynr, &znr);
			if(sizeX==0) sizeX=xnr;
			if(xnr!=sizeX){
				cout << "ERROR: Map dimensions between maps need to be the same.\n";
				exit(3);
			}
			if(sizeY==0) sizeY=ynr;
			if(ynr!=sizeY){
				cout << "ERROR: Map dimensions between maps need to be the same.\n";
				exit(3);
			}
			if(sizeZ==0) sizeZ=znr;
			if(znr!=sizeZ){
				cout << "ERROR: Map dimensions between maps need to be the same.\n";
				exit(3);
			}
			found_all++;
		}
		if(line.find("CENTER") == 0){
			sscanf(&line.c_str()[7], "%f %f %f", &center_x, &center_y, &center_z);
			found_all++;
		}
		if(line.find("SPACING") == 0){
			sscanf(&line.c_str()[8], "%f", &spacing);
			found_all++;
		}
	} while (line.find("CENTER") == std::string::npos);
	if(found_all != 3){
		cout << "ERROR: Grid map file does not NELEMENTS, SPACING, and CENTER fields.";
		exit(7);
	}
	fp_num* data = map_storage;
	if(data == NULL) data = (fp_num*)malloc(sizeof(fp_num)*(sizeX+1)*(sizeY+1)*(sizeZ+1) + 9);
	data[0] = sizeX;
	data[1] = sizeY;
	data[2] = sizeZ;
	data[3] = center_x;
	data[4] = center_y;
	data[5] = center_z;
	data[6] = spacing;
	fp_num* mypoi = data + 9;
	//reading values
	fp_num d;
	fp_num val_min = 1e80;
	fp_num val_max = 0;
	for(int z=0; z<=sizeZ; z++)
		for(int y=0; y<=sizeY; y++)
			for(int x=0; x<=sizeX; x++)
			{
				std::getline(file, line);
				d = map2float(line.c_str());
				val_min = std::min(val_min, d);
				val_max = std::max(val_max, d);
				*(mypoi++) = d;
			}
	file.close();
	data[7] = val_min;
	data[8] = val_max;
	
	return data;
}

#endif // INCLUDED_GRID_READER

