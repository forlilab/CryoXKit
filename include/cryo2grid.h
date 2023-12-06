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

void print_version_info();

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

std::vector<fp_num> create_mask(
                                std::vector<fp_num> &grid_or_mask,
                                std::string          mask_pdb,
                                fp_num               rT          = 2,
                                bool                 subtractive = true,
                                bool                 create_new  = true
                               );

std::vector<fp_num> apply_mask(
                               std::vector<fp_num> density,
                               std::vector<fp_num> mask
                              );

void write_density(
                   std::vector<fp_num> density,
                   std::string basename,
                   int write_type = write_grid_mrc
                  );

void write_grid_maps(
                     std::vector<fp_num> density,
                     std::vector<GridMap> grid_maps,
                     std::vector<std::string> grid_files,
                     int write_type = write_grid_ad4
                    );

std::vector<fp_num> average_densities_to_grid(
                                              std::vector<std::string> map_files,
                                              std::vector<std::string> map_receptors,
                                                          std::string  align_rec,
                                                       int map_type,
                                              unsigned int map_x_dim,
                                              unsigned int map_y_dim,
                                              unsigned int map_z_dim,
                                              fp_num       map_x_center,
                                              fp_num       map_y_center,
                                              fp_num       map_z_center,
                                              fp_num       grid_spacing,
                                              bool         repeat_unit_cell = true,
                                              bool         output_align_rec = false
                                             );

#endif // INCLUDED_CRYO2GRID

