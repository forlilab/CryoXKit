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

#ifdef PARALLELIZE
#include <omp.h>
#endif
#include "include/Config.h"
#include "include/grid_reader.h"
#include "include/pdb_reader.h"
#include "include/map_reader.h"
#include "include/map_writer.h"
#include "include/map_modifier.h"
#include "include/cryo2grid.h"
#ifndef _WIN32
// libgen.h contains basename() and dirname() from a fullpath name
// Specific: to open correctly grid map field fiels and associated files
// http://ask.systutorials.com/681/get-the-directory-path-and-file-name-from-absolute-path-linux
#include <libgen.h>
#endif

bool has_absolute_path(const char* filename)
{
	#ifndef _WIN32
	return (filename[0]=='/');
	#else
	char drive_tmp[_MAX_DRIVE];
	char path_tmp[_MAX_DIR];
	_splitpath(filename, drive_tmp, path_tmp, NULL, NULL);
	return ((strlen(drive_tmp)>0) || (path_tmp[0]=='\\') || (path_tmp[0]=='/'));
	#endif
}

std::string get_filepath(const char* filename)
{
	#ifndef _WIN32
	char* ts1 = strdup(filename);
	std::string result = dirname(ts1);
	free(ts1);
	return result;
	#else
	char drive_tmp[_MAX_DRIVE];
	char path_tmp[_MAX_DIR];
	_splitpath(filename, drive_tmp, path_tmp, NULL, NULL);
	return drive_tmp + path_tmp;
	#endif
}

std::vector<GridMap> read_grid_maps(
                                    std::vector<std::string> grid_files,
                                    std::string              rec_name
                                   )
{
	std::vector<GridMap> grid_maps;
	if(grid_files.size() > 0){
		timeval runtime;
		start_timer(runtime);
		int X_dim = 0;
		int Y_dim = 0;
		int Z_dim = 0;
		std::string receptor_file = rec_name;
		cout << "Reading grid map files:\n";
		cout << "\t-> " << grid_files[0] << "\n";
		grid_maps.push_back(read_grid_map(grid_files[0], X_dim, Y_dim, Z_dim, receptor_file));
		grid_maps.resize(grid_files.size());
		#pragma omp parallel for
		for(unsigned int i=1; i<grid_files.size(); i++){
			#pragma omp critical
			cout << "\t-> " << grid_files[i] << "\n";
			grid_maps[i] = read_grid_map(grid_files[i], X_dim, Y_dim, Z_dim, receptor_file);
		}
		cout << "<- Done, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";
	}
	return grid_maps;
}

std::string get_grid_receptor_filename(
                                       std::vector<GridMap> grid_maps,
                                       std::vector<std::string> grid_files
                                      )
{
	std::string rec_name = "";
	if(grid_maps.size() == 0) return rec_name;
	if((unsigned int)(grid_maps[0])[0]<=10) return rec_name;
	rec_name.assign(reinterpret_cast<char*>(grid_maps[0].data() + 10));
	std::string grid_path = get_filepath(grid_files[0].c_str());
	if(grid_path==".") grid_path="";
	if(grid_path.size()>0){
		grid_path  += "/";
		if(!has_absolute_path(rec_name.c_str()))
			rec_name = grid_path + rec_name;
	}
	return rec_name;
}

void write_grid_maps(
                     std::vector<fp_num> density,
                     std::vector<GridMap> grid_maps,
                     std::vector<std::string> grid_files,
                     int write_type
                    )
{
	if(grid_files.size() > 0){
		int X_dim = (grid_maps[0])[1];
		int Y_dim = (grid_maps[0])[2];
		int Z_dim = (grid_maps[0])[3];
		#pragma omp parallel for
		for(unsigned int i=0; i<grid_files.size(); i++){
			int grid_points = (X_dim + 1) * (Y_dim + 1) * (Z_dim + 1) + (unsigned int)((grid_maps[i])[0]);
			int offset = (unsigned int)density[0] - (unsigned int)((grid_maps[i])[0]);
			for(int j=(unsigned int)((grid_maps[i])[0]); j<grid_points; j++)
				(grid_maps[i])[j] += density[j + offset];
			write_grid(
			           grid_maps[i].data(),
			           grid_files[i],
			           write_type,
			           0,
			           false
			          );
		}
		cout << "\n";
	}
}

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
                                              bool         repeat_unit_cell
                                             )
{
	if(map_files.size() < 1){
		cout << "ERROR: No density map file(s) specified, nothing to do.\n";
		exit(1);
	}
	std::vector<std::vector<fp_num>> densities;
	densities.resize(map_files.size());
	#pragma omp parallel for
	for(unsigned int i=0; i<map_files.size(); i++){
		fp_num* grid_align = NULL;
		if(i < map_receptors.size())
			if(map_receptors[i].size() > 4) // i.e. longer than ".pdb"
				grid_align = align_pdb_atoms(
				                             map_receptors[i],
				                             align_rec,
				                             map_x_dim,
				                             map_y_dim,
				                             map_z_dim,
				                             map_x_center,
				                             map_y_center,
				                             map_z_center,
				                             grid_spacing
				                            );
		densities[i] = read_map_to_grid(
		                                map_files[i],
		                                map_type,
		                                map_x_dim,
		                                map_y_dim,
		                                map_z_dim,
		                                map_x_center,
		                                map_y_center,
		                                map_z_center,
		                                grid_spacing,
		                                repeat_unit_cell,
		                                grid_align
		                               );
		if(grid_align != NULL) delete[] grid_align;
		
		// normalize if more than one map file
		if(map_files.size() > 1){
			fp_num rho_avg = 0;
			fp_num rho_std = 0;
			fp_num rho;
			for(unsigned int j=(unsigned int)(densities[i])[0]; j<densities[i].size(); j++){
				rho = (densities[i])[j];
				rho_avg += rho;
				rho_std += rho*rho;
			}
			rho_avg /= densities[i].size() - (unsigned int)(densities[i])[0];
			rho_std /= densities[i].size() - (unsigned int)(densities[i])[0];
			rho_std -= rho_avg * rho_avg;
			rho_std  = sqrt(rho_std) * map_files.size(); // multiply by number of maps here so adding below gives us an average
			for(unsigned int j=(unsigned int)(densities[i])[0]; j<densities[i].size(); j++){
				(densities[i])[j] -= rho_avg;
				(densities[i])[j] /= rho_std;
			}
		}
	}
	for(unsigned int i=1; i<map_files.size(); i++){
		#pragma omp parallel for
		for(unsigned int j=(unsigned int)(densities[i])[0]; j<densities[i].size(); j++){
			(densities[0])[j] += (densities[i])[j];
		}
	}
	// recalculate rho_min, rho_max
	if(map_files.size() > 1){
		fp_num rho_min = 1e80;
		fp_num rho_max = 0;
		for(unsigned int j=(unsigned int)(densities[0])[0]; j<densities[0].size(); j++){
			rho_min = std::min((densities[0])[j], rho_min);
			rho_max = std::max((densities[0])[j], rho_max);
		}
		(densities[0])[8] = rho_min;
		(densities[0])[9] = rho_max;
	}
	return densities[0];
}

int main(int argc, const char* argv[])
{
	timeval runtime;
	start_timer(runtime);
	
	std::vector<std::string> map_files{""};
	std::vector<std::string> map_receptors{""};
	int X_dim = 0;
	int Y_dim = 0;
	int Z_dim = 0;
	fp_num X_center, Y_center, Z_center;
	fp_num grid_spacing = 0.375;
	int write_type      = write_grid_ad4;
	int mod_type        = log_modifier;
	bool argument_error = true;
	std::vector<std::string> grid_files;
	std::string align_rec = "";
	// Check for command line parameters
	if(argc>2){ // yes, there are some -- parameter required are: (grid filename xor grid center x,y,z, grid x,y,z dimensions, grid spacing, and write type) as well as optionally modifier type
		map_files[0]    = argv[1]; // map filename
		string grid     = argv[2]; // grid filename XOR
		std::size_t ext = grid.find_last_of(".");
		if(grid.substr(ext).compare(".map")==0)
		{
			if(grid_filter(grid.substr(0,ext))) grid_files.push_back(grid);
			int count = 3;
			while(count < argc){
				grid = argv[count];
				ext  = grid.find_last_of(".");
				if(grid.substr(ext).compare(".map")!=0) break; // not a grid map file
				if(grid_filter(grid.substr(0,ext))) grid_files.push_back(grid);
				count++;
			}
			while(count < argc){
				grid = argv[count];
				ext  = grid.find_last_of(".");
				if((grid.substr(ext).compare(".pdb")==0) ||
				   (grid.substr(ext).compare(".pdbqt")==0)){
					if(map_receptors[0].size()==0){
						map_receptors[0] = grid;
					} else align_rec = grid;
					count++;
				} else mod_type = atoi(argv[count++]);
			}
			argument_error = (grid_files.size() == 0);
			if(argument_error){
				cout << "ERROR: Could not find grid map files or only e, d, or H* maps were specified.\n";
				exit(1);
			}
		} else if((grid.substr(ext).compare(".pdb")==0) ||
		          (grid.substr(ext).compare(".pdbqt")==0))
		{
			map_receptors[0] = grid;
			if(argc>3){
				grid = argv[3];
				ext  = grid.find_last_of(".");
				if((grid.substr(ext).compare(".pdb")==0) ||
				   (grid.substr(ext).compare(".pdbqt")==0))
					align_rec = grid;
			}
			argument_error = (argc <= 3);
		} else{
			if(argc>7){
				X_center = atof(argv[2]); // grid center
				Y_center = atof(argv[3]);
				Z_center = atof(argv[4]);
				X_dim    = atoi(argv[5]); // dimensions
				Y_dim    = atoi(argv[6]);
				Z_dim    = atoi(argv[7]);
				if((X_dim <= 0) || (Y_dim <= 0) || (Z_dim <= 0)){
					cout << "ERROR: Please ensure grid dimensions are each greater than 1.\n";
					exit(1);
				}
				if(argc>8) grid_spacing = atof(argv[8]); // grid spacing
				if(argc>9) write_type = atoi(argv[9]); // write type
				if(argc>10) mod_type = atoi(argv[10]); // modifier fxn type
				argument_error = false;
			} else argument_error = true;
		}
	}
	if(argument_error){
		cout << "Syntax:\n";
		cout << argv[0] << " mapfile center_x center_y center_z x_dim y_dim z_dim (spacing [" << grid_spacing << "]) (write [" << write_type << " = AD4 map]) (modifier fxn [" << mod_type << " = logistics])\n"; // argv[0] is program name
		cout << "*or* for map modification (automatically excludes e, d, and H* maps):\n";
		cout << argv[0] << " mapfile gridfile (map ligand) (modifier fxn [" << mod_type << " = logistics])\n"; // argv[0] is program name
		exit(1);
	}
	
	std::vector<GridMap> grid_maps;
	if(grid_files.size()>0){
		grid_maps    = read_grid_maps(grid_files, align_rec);
		align_rec    = get_grid_receptor_filename(grid_maps, grid_files);
		X_dim        = (grid_maps[0])[1];
		Y_dim        = (grid_maps[0])[2];
		Z_dim        = (grid_maps[0])[3];
		X_center     = (grid_maps[0])[4];
		Y_center     = (grid_maps[0])[5];
		Z_center     = (grid_maps[0])[6];
		grid_spacing = (grid_maps[0])[7];
	}
	
	std::vector<fp_num> density = average_densities_to_grid(
	                                                        map_files,
	                                                        map_receptors,
	                                                        align_rec,
	                                                        automatic,
	                                                        X_dim,
	                                                        Y_dim,
	                                                        Z_dim,
	                                                        X_center,
	                                                        Y_center,
	                                                        Z_center,
	                                                        grid_spacing,
	                                                        true
	                                                       );
	
	std::vector<fp_num> modified = modify_densities(
	                                                density,
	                                                mod_type
	                                               );
	
	if(grid_files.size()>0){
		write_grid_maps(modified, grid_maps, grid_files, write_type);
	} else{
		write_grid(
		           modified.data(),
		           map_files[0],
		           write_type
		          );
	}
	cout << "Done. Overall runtime was " << seconds_since(runtime)*1000.0 << " ms.\n";
	return 0;
}

