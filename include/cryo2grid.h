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

#include "cryo2grid.h"

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

std::vector<GridMap> read_grid_maps(std::vector<std::string> grid_files);

void write_grid_maps(
                     std::vector<fp_num> density,
                     std::vector<GridMap> grid_maps,
                     std::vector<std::string> grid_files,
                     int write_type = write_grid_ad4
                    );

#endif // INCLUDED_CRYO2GRID

