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


#ifndef INCLUDED_MAP_WRITER
#define INCLUDED_MAP_WRITER

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

inline std::string num2str(fp_num num)
{
	unsigned int l = fabs(num);
	unsigned int decimals = (fabs(num)-l)*1000 + 0.5;
	return (num<0?"-":"") + to_string(l) + "." + (decimals<100?"0":"") + (decimals<10?"0":"") + to_string(decimals);
}

inline void write_grid_map_ad4(
                               fp_num*      grid_map,
                               std::string &filename,
                               unsigned int map_x_dim,
                               unsigned int map_y_dim,
                               unsigned int map_z_dim,
                               fp_num       map_x_center,
                               fp_num       map_y_center,
                               fp_num       map_z_center,
                               fp_num       grid_spacing,
                               bool         set_extension = true
                              )
{
	if(set_extension){
		std::size_t ext = filename.find_last_of(".");
		filename = filename.substr(0, ext) + ".map";
	}
	cout << "Writing AD4 grid map file [" << filename << "]\n";
	std::ofstream grid_file(filename);
	if(grid_file.fail()){
		cout << "Error: Can't open grid map output file " << filename << ".\n";
		exit(1);
	}
	grid_file << "GRID_PARAMETER_FILE none\n";
	grid_file << "GRID_DATA_FILE none\n";
	grid_file << "MACROMOLECULE none\n";
	grid_file.precision(3);
	grid_file.setf(ios::fixed, ios::floatfield);
	grid_file << "SPACING " << grid_spacing << "\n";
	grid_file << "NELEMENTS " << map_x_dim << " " << map_y_dim << " " << map_z_dim << "\n";
	grid_file << "CENTER " << map_x_center << " " << map_y_center << " " << map_z_center << "\n";
	
	unsigned int grid_points = (map_x_dim + 1) * (map_y_dim + 1) * (map_z_dim + 1);
	std::string data_block;
	data_block.reserve(6.5*grid_points); // at least zero + point + 3 decimals + linebreak = 6
	for(unsigned int i=0; i < grid_points; i++)
		data_block += num2str(grid_map[i]) + "\n";
	
	grid_file.write(data_block.c_str(), data_block.size());
	
	grid_file.close();
}

inline void write_grid_map_mrc(
                               fp_num*      grid_map,
                               std::string &filename,
                               unsigned int map_x_dim,
                               unsigned int map_y_dim,
                               unsigned int map_z_dim,
                               fp_num       map_x_center,
                               fp_num       map_y_center,
                               fp_num       map_z_center,
                               fp_num       grid_spacing,
                               fp_num       rho_min       = -1,
                               fp_num       rho_max       =  1,
                               bool         set_extension = true
                              )
{
	if(set_extension){
		std::size_t ext = filename.find_last_of(".");
		filename = filename.substr(0, ext) + ".grid.mrc";
	}
	cout << "Writing MRC grid map file [" << filename << "]\n";
	std::ofstream map_file(filename, std::ifstream::binary);
	if(map_file.fail()){
		cout << "Error: Can't open grid map output file " << filename << ".\n";
		exit(1);
	}
	struct mrc_header{
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int mode       = 2;
		unsigned int x_start    = 0;
		unsigned int y_start    = 0;
		unsigned int z_start    = 0;
		unsigned int mx;
		unsigned int my;
		unsigned int mz;
		float        cell_a;
		float        cell_b;
		float        cell_c;
		float        alpha      = 90.0f;
		float        beta       = 90.0f;
		float        gamma      = 90.0f;
		unsigned int map_x      = 1;
		unsigned int map_y      = 2;
		unsigned int map_z      = 3;
		float        val_min    = -1.0f;
		float        val_max    = 1.0f;
		float        val_avg    = 0;
		unsigned int spacegroup = 1;
		unsigned int ext_header = 0;
		char extra[100];
		float        x_origin;
		float        y_origin;
		float        z_origin;
		char map_str[4]         = {'M', 'A', 'P', ' '};
		unsigned int mach_str;
		float        val_std    = 1.0f;
		unsigned int nlabel     = 1;
		char labels[800];
	} header;
	memset(header.extra, 0,100);
	if(HOST_LITTLE_ENDIAN){
		header.extra[12] = 0xAD;
		header.extra[13] = 0x4E;
	} else{
		header.extra[14] = 0x4E;
		header.extra[15] = 0xAD;
	}
	memset(header.labels,0,800);
	strncpy(header.labels, "Cryo2Grid MRC grid map", 23);
	header.nx       = map_x_dim + 1;
	header.mx       = map_x_dim + 1;
	header.ny       = map_y_dim + 1;
	header.my       = map_y_dim + 1;
	header.nz       = map_z_dim + 1;
	header.mz       = map_z_dim + 1;
	header.cell_a   = (map_x_dim + 1) * grid_spacing;
	header.cell_b   = (map_y_dim + 1) * grid_spacing;
	header.cell_c   = (map_z_dim + 1) * grid_spacing;
	header.x_origin = map_x_center - map_x_dim * grid_spacing * 0.5;
	header.y_origin = map_y_center - map_y_dim * grid_spacing * 0.5;
	header.z_origin = map_z_center - map_z_dim * grid_spacing * 0.5;
	if((rho_min <= rho_max) && (0 >= std::min(rho_min, rho_max))){
		header.val_min = rho_min;
		header.val_max = rho_max;
	}
	header.mach_str = HOST_LITTLE_ENDIAN ? 17476 : 4369;
	
	map_file.write(reinterpret_cast<char*>(&header),sizeof(header));
	map_file.write(reinterpret_cast<char*>(grid_map), (map_x_dim + 1) * (map_y_dim + 1) * (map_z_dim + 1) * sizeof(float));
	
	map_file.close();
}

inline void write_grid(
                       std::vector<fp_num> grid_map,
                       std::string        &filename,
                       int                 write_mode = write_grid_ad4
                      )
{
	timeval runtime;
	start_timer(runtime);
	switch(write_mode){
		case write_grid_ad4: write_grid_map_ad4(
		                                        grid_map.data() + 9,
		                                        filename,
		                                        grid_map[0],
		                                        grid_map[1],
		                                        grid_map[2],
		                                        grid_map[3],
		                                        grid_map[4],
		                                        grid_map[5],
		                                        grid_map[6],
		                                        true
		                                       );
		                     cout << "<- Finished writing, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";;
		                     break;
		case write_grid_mrc: write_grid_map_mrc(
		                                        grid_map.data() + 9,
		                                        filename,
		                                        grid_map[0],
		                                        grid_map[1],
		                                        grid_map[2],
		                                        grid_map[3],
		                                        grid_map[4],
		                                        grid_map[5],
		                                        grid_map[6],
		                                        grid_map[7],
		                                        grid_map[8],
		                                        true
		                                       );
		                     cout << "<- Finished writing, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";;
		default:             break;
	}
}

#endif // INCLUDED_MAP_WRITER

