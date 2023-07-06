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


#ifndef INCLUDED_CRYO2GRID
#define INCLUDED_CRYO2GRID

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <iterator>

inline bool grid_filter(std::string name)
{
	std::size_t ext = name.find_last_of(".");
	if(ext == std::string::npos) return false;
	return !((name.substr(ext).compare(".e")==0) || (name.substr(ext).compare(".d")==0) || (name.substr(ext).compare(".H")==0) || (name.substr(ext).compare(".H")==1));
}

inline std::vector<std::string> filter_grid_files(std::vector<std::string> grid_files)
{
	std::vector<std::string> filtered;
	for(unsigned int i=0; i<grid_files.size(); i++){
		std::size_t ext  = grid_files[i].find_last_of(".");
		if((grid_files[i].substr(ext).compare(".map")==0) && grid_filter(grid_files[i].substr(0,ext)))
			filtered.push_back(grid_files[i]);
	}
	if(filtered.size() == 0){
		cout << "ERROR: Could not find grid map files or only e, d, or H* maps were specified.\n";
		exit(1);
	}
	return filtered;
}

std::vector<GridMap> read_grid_maps(
                                    std::vector<std::string> grid_files,
                                    std::string              rec_name = ""
                                   );

std::string get_grid_receptor_filename(
                                       std::vector<GridMap> grid_maps,
                                       std::vector<std::string> grid_files
                                      );

void write_grid_maps(
                     std::vector<fp_num> density,
                     std::vector<GridMap> grid_maps,
                     std::vector<std::string> grid_files,
                     int write_type = write_grid_ad4
                    );

#endif // INCLUDED_CRYO2GRID

