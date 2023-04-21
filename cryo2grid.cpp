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

#include "config.h"
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
	int write_type = write_ad4map;
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
		if(argc>8) grid_spacing = atof(argv[8]);
		if(argc>9) write_type = atoi(argv[9]);
	} else{
		cout << "Syntax: " << argv[0] << " mapfile center_x center_y center_z x_dim y_dim z_dim (spacing [" << grid_spacing << "]) (write [" << write_type << " = AD4 map])\n"; // argv[0] is program name
		exit(1);
	}
	
	std::vector<fp_num> density = read_map(
	                                       map_file,
	                                       automatic,
	                                       X_center,
	                                       Y_center,
	                                       Z_center,
	                                       X_dim,
	                                       Y_dim,
	                                       Z_dim,
	                                       grid_spacing,
	                                       write_type
	                                      );
	return 0;
}

