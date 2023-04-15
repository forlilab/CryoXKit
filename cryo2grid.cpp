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

#define USE_SINGLE_PRECISION

#ifdef USE_SINGLE_PRECISION
typedef float fp_num;
#define twoTo23 8388608.0f
inline fp_num fastfloor(fp_num f)
{
	fp_num c = (f >= 0.0 ? -twoTo23 : twoTo23);
	fp_num result = (f - c) + c;
	if(f < result) result -= 1.0;
	return result;
}
#else
typedef double fp_num;
#define twoTo52 4503599627370496.0
inline fp_num fastfloor(fp_num f)
{
	fp_num c = (f >= 0.0 ? -twoTo52 : twoTo52);
	fp_num result = (f - c) + c;
	if(f < result) result -= 1.0;
	return result;
}
#endif

#include "cryo2grid.h"

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <iterator>

int main(int argc, const char* argv[])
{
	string map_file="";
	fp_num X_center, Y_center, Z_center;
	int X_dim, Y_dim, Z_dim;
	fp_num grid_spacing = 0.375;
	// Check for command line parameters
	if(argc>7){ // yes, there are some -- parameter required are map filename, grid center x,y,z, grid x,y,z dimensions, and grid spacing
		map_file = argv[1];
		X_center = atof(argv[2]);
		Y_center = atof(argv[3]);
		Z_center = atof(argv[4]);
		X_dim    = atoi(argv[5]);
		Y_dim    = atoi(argv[6]);
		Z_dim    = atoi(argv[7]);
		if((X_dim <= 0) || (Y_dim <= 0) || (Z_dim <= 0)){
			cout << "ERROR: Please ensure grid dimensions are each greater than 1.\n";
			exit(1);
		}
		if(argc>8) grid_spacing = atoi(argv[8]);
	} else{
		cout << "Syntax: " << argv[0] << " <map file> <center x> <center y> <center z> <x dim> <y dim> <z dim> <optional: grid spacing (default: " << grid_spacing << ")>\n"; // argv[0] is program name
		exit(1);
	}
	
	std::vector<fp_num> density = read_dsn6(
	                                        map_file,
	                                        X_center,
	                                        Y_center,
	                                        Z_center,
	                                        X_dim,
	                                        Y_dim,
	                                        Z_dim,
	                                        grid_spacing,
	                                        true
	                                       );
	
	return 0;
}

