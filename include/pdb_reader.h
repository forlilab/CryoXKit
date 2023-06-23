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


#ifndef INCLUDED_PDB_READER
#define INCLUDED_PDB_READER

#ifndef INCLUDED_VECTORMAT
#include "VecMat.h"
#endif

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

// structure to store relevant PDB atom data (no charge stored as a) PDB doesn't have it, b) it's not needed)
// ATOM      1  N   SER A   1      -2.367   4.481 -16.909  1.00  1.00     0.185 N
typedef struct
{
	bool         hetatm;       // false (wether first field is "HETATM" or not)
	char         name[5];      // "N"
	char         res_name[4];  // "SER"
	char         chain_id[2];  // "A"
	unsigned int res_id;       // 1
	float        x,y,z;        // -2.367, 4.481, -16.909
	char         atom_type[4]; // "N"
} PDBatom;

inline bool point_in_box(
                         Vec3<fp_num> point,
                         Vec3<fp_num> box_dim
                        )
{
	for(unsigned int i=0; i<3; i++)
		if((point.vec[i]<0) || (point.vec[i]>box_dim.vec[i]))
			return false; // if one dimension is outside then the point is outside
	return true;
}


Mat33<fp_num> align_atoms(
                          std::vector<PDBatom> &map_atoms,
                          std::vector<PDBatom> &grid_atoms,
                          int                   map_x_dim,
                          int                   map_y_dim,
                          int                   map_z_dim,
                          fp_num                map_x_center,
                          fp_num                map_y_center,
                          fp_num                map_z_center,
                          fp_num                grid_spacing,
                          Vec3<fp_num>         &translate
                         )
{
	timeval runtime;
	start_timer(runtime);
	cout << "Aligning map receptor to receptor ...\n";
	bool use_grid_box = (map_x_dim > 0) && (map_y_dim > 0) && (map_z_dim > 0);
	Vec3<fp_num> grid_dims;
	grid_dims.vec[0] = map_x_dim * grid_spacing; grid_dims.vec[1] = map_y_dim * grid_spacing; grid_dims.vec[2] = map_z_dim * grid_spacing;
	Vec3<fp_num> grid_start;
	grid_start.vec[0] = map_x_center - grid_dims.vec[0] * 0.5; grid_start.vec[1] = map_y_center - grid_dims.vec[1] * 0.5; grid_start.vec[2] = map_z_center - grid_dims.vec[2] * 0.5;
	Vec3<fp_num> location;
	Vec3<fp_num> center(0.0);
	std::vector<unsigned int> grid_ids;
	// Calculate heavy atom geometric center of atoms in grid box or of all if there is no grid box
	for(unsigned int i=0; i<grid_atoms.size(); i++){
		if((grid_atoms[i].hetatm == false) && (grid_atoms[i].atom_type[0] != 'H')){ // only focus on large molecules in pdb(qt)
			location.vec[0] = grid_atoms[i].x; location.vec[1] = grid_atoms[i].y; location.vec[2] = grid_atoms[i].z;
			if(!use_grid_box || point_in_box(location - grid_start, grid_dims)){
				center += location;
				grid_ids.push_back(i);
			}
		}
	}
	if(grid_ids.size() == 0){ // shouldn't happen
		cout << "ERROR: No receptor atoms inside grid box. This likely means the wrong receptor is listed in the grid file.\n";
		exit(1);
	}
	center /= grid_ids.size();
	// now find atom closest to center
	location.vec[0] = grid_atoms[grid_ids[0]].x - center.vec[0];
	location.vec[1] = grid_atoms[grid_ids[0]].y - center.vec[1];
	location.vec[2] = grid_atoms[grid_ids[0]].z - center.vec[2];
	fp_num closest_dist2 = location * location;
	int closest = 0;
	for(unsigned int i=1; i<grid_ids.size(); i++){
		location.vec[0] = grid_atoms[grid_ids[i]].x - center.vec[0];
		location.vec[1] = grid_atoms[grid_ids[i]].y - center.vec[1];
		location.vec[2] = grid_atoms[grid_ids[i]].z - center.vec[2];
		fp_num d2 = location*location;
		if(d2 < closest_dist2){
			closest_dist2 = d2;
			closest = i;
		}
	}
	// match atoms from grid_ids to map atoms
	std::vector<unsigned int> map_match;
	for(int i=0; i<(int)grid_ids.size(); i++){
		bool found = false;
		for(unsigned int j=0; j<map_atoms.size(); j++){
			if((map_atoms[j].res_id == grid_atoms[grid_ids[i]].res_id) &&
			   (map_atoms[j].chain_id[0] == grid_atoms[grid_ids[i]].chain_id[0]) &&
			   (strcmp(map_atoms[j].res_name, grid_atoms[grid_ids[i]].res_name) == 0) &&
			   (strcmp(map_atoms[j].name, grid_atoms[grid_ids[i]].name) == 0))
			{
				map_match.push_back(j);
				found = true;
				break;
			}
		}
		if(!found){
			grid_ids.erase(grid_ids.begin() + i);
			if(closest==i)
				closest = -1;
			else if(closest>i) closest--;
			i--;
		}
	}
	Vec3<fp_num> map_center;
	if(closest < 0){ // fallback is to you use residue center (assumption is that it won't be different enough to add enough error
		center.vec[0] = 0; center.vec[1] = 0; center.vec[2] = 0;
		map_center.vec[0] = 0; map_center.vec[1] = 0; map_center.vec[2] = 0;
		unsigned int count = 0;
		for(unsigned int i=0; i<grid_ids.size(); i++){
			if((grid_atoms[grid_ids[i]].res_id == grid_atoms[grid_ids[closest]].res_id) &&
			   (grid_atoms[grid_ids[i]].chain_id[0] == grid_atoms[grid_ids[closest]].chain_id[0]) &&
			   (strcmp(grid_atoms[grid_ids[i]].res_name, grid_atoms[grid_ids[closest]].res_name) == 0))
			{
				center.vec[0] += grid_atoms[grid_ids[i]].x; center.vec[1] += grid_atoms[grid_ids[i]].y; center.vec[2] += grid_atoms[grid_ids[i]].z;
				map_center.vec[0] += map_atoms[map_match[i]].x; map_center.vec[1] += map_atoms[map_match[i]].y; map_center.vec[2] += map_atoms[map_match[i]].z;
				count++;
			}
		}
		center /= count;
		map_center /= count;
	} else{
		// new center to align to based on closest atom to geometric center
		center.vec[0] = grid_atoms[grid_ids[closest]].x; center.vec[1] = grid_atoms[grid_ids[closest]].y; center.vec[2] = grid_atoms[grid_ids[closest]].z;
		map_center.vec[0] = map_atoms[map_match[closest]].x; map_center.vec[1] = map_atoms[map_match[closest]].y; map_center.vec[2] = map_atoms[map_match[closest]].z;
	}
	translate = center - map_center;
	cout << "\t-> Translation vector: (" << translate.V3Str(',') << ")\n";
	Mat33<double> map_S, grid_S;
	map_S.M3Zeros(); grid_S.M3Zeros();
	// Calculate gyration tensors for both
	for(unsigned int i=0; i<grid_ids.size(); i++){
		// align to respective centers
		grid_atoms[grid_ids[i]].x -= center.vec[0];
		grid_atoms[grid_ids[i]].y -= center.vec[1];
		grid_atoms[grid_ids[i]].z -= center.vec[2];
		map_atoms[map_match[i]].x -= map_center.vec[0];
		map_atoms[map_match[i]].y -= map_center.vec[1];
		map_atoms[map_match[i]].z -= map_center.vec[2];
	}
	double norm = 1.0 / grid_ids.size();
	norm *= norm;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		for(unsigned int j=i+1; j<grid_ids.size(); j++){
			// S_mn = 1/N^2 * sum_(i>j) (r_m(i) - r_m(j)) * (r_n(i) - r_n(j))
			grid_S.mat[0][0] += norm * (grid_atoms[grid_ids[i]].x - grid_atoms[grid_ids[j]].x) * (grid_atoms[grid_ids[i]].x - grid_atoms[grid_ids[j]].x);
			grid_S.mat[0][1] += norm * (grid_atoms[grid_ids[i]].x - grid_atoms[grid_ids[j]].x) * (grid_atoms[grid_ids[i]].y - grid_atoms[grid_ids[j]].y);
			grid_S.mat[0][2] += norm * (grid_atoms[grid_ids[i]].x - grid_atoms[grid_ids[j]].x) * (grid_atoms[grid_ids[i]].z - grid_atoms[grid_ids[j]].z);
			grid_S.mat[1][1] += norm * (grid_atoms[grid_ids[i]].y - grid_atoms[grid_ids[j]].y) * (grid_atoms[grid_ids[i]].y - grid_atoms[grid_ids[j]].y);
			grid_S.mat[1][2] += norm * (grid_atoms[grid_ids[i]].y - grid_atoms[grid_ids[j]].y) * (grid_atoms[grid_ids[i]].z - grid_atoms[grid_ids[j]].z);
			grid_S.mat[2][2] += norm * (grid_atoms[grid_ids[i]].z - grid_atoms[grid_ids[j]].z) * (grid_atoms[grid_ids[i]].z - grid_atoms[grid_ids[j]].z);
		
			map_S.mat[0][0] += norm * (map_atoms[map_match[i]].x - map_atoms[map_match[j]].x) * (map_atoms[map_match[i]].x - map_atoms[map_match[j]].x);
			map_S.mat[0][1] += norm * (map_atoms[map_match[i]].x - map_atoms[map_match[j]].x) * (map_atoms[map_match[i]].y - map_atoms[map_match[j]].y);
			map_S.mat[0][2] += norm * (map_atoms[map_match[i]].x - map_atoms[map_match[j]].x) * (map_atoms[map_match[i]].z - map_atoms[map_match[j]].z);
			map_S.mat[1][1] += norm * (map_atoms[map_match[i]].y - map_atoms[map_match[j]].y) * (map_atoms[map_match[i]].y - map_atoms[map_match[j]].y);
			map_S.mat[1][2] += norm * (map_atoms[map_match[i]].y - map_atoms[map_match[j]].y) * (map_atoms[map_match[i]].z - map_atoms[map_match[j]].z);
			map_S.mat[2][2] += norm * (map_atoms[map_match[i]].z - map_atoms[map_match[j]].z) * (map_atoms[map_match[i]].z - map_atoms[map_match[j]].z);
		}
	}
	grid_S.mat[1][0] = grid_S.mat[0][1];
	grid_S.mat[2][0] = grid_S.mat[0][2];
	grid_S.mat[2][1] = grid_S.mat[1][2];
	CVec3<double> cew = grid_S.Eigenvalues();
	if(cew.Im()*cew.Im()>EPS){ // gyration tensor eigenvalues need to be real
		cout << "ERROR: Grid gyration tensor eigenvalues need to be real.\n";
		exit(2);
	}
	Vec3<double> ew = cew.Re();
	Mat33<double> grid_rot = RotFromEigenvectors(grid_S.Eigenvectors(ew));

	map_S.mat[1][0] = map_S.mat[0][1];
	map_S.mat[2][0] = map_S.mat[0][2];
	map_S.mat[2][1] = map_S.mat[1][2];
	cew = map_S.Eigenvalues();
	if(cew.Im()*cew.Im()>EPS){ // gyration tensor eigenvalues need to be real
		cout << "ERROR: Map gyration tensor eigenvalues need to be real.\n";
		exit(3);
	}
	ew = cew.Re();
	Mat33<double> map_rot = RotFromEigenvectors(map_S.Eigenvectors(ew));
//	cout << "grid:\n" << grid_rot.M3Str() << "\n";
//	cout << "map:\n" << map_rot.M3Str() << "\n";
	// X * map_rot = grid_rot
	// X = grid_rot * map_rot^T (inverse is transpose for rotation matrices)
	Mat33<fp_num> result;
	result = grid_rot * map_rot.M3Transpose();
	cout << "\t-> Rotation matrix:\n" << result.M3Str() << "\n";
	// calculate RMSD
	fp_num rmsd = 0;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		location.vec[0] = map_atoms[map_match[i]].x;
		location.vec[1] = map_atoms[map_match[i]].y;
		location.vec[2] = map_atoms[map_match[i]].z;
		location = result * location;
		location.vec[0] -= grid_atoms[grid_ids[i]].x;
		location.vec[1] -= grid_atoms[grid_ids[i]].y;
		location.vec[2] -= grid_atoms[grid_ids[i]].z;
		rmsd += (location * location);
	}
	cout << "\t-> RMSD after alignment (" << grid_ids.size() << " atoms): " << sqrt(rmsd/grid_ids.size()) << " A\n";
	cout << "<- Finished alignment, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";
	return result;
}

// trim input string range -- careful: no range checking
static inline void range_trim_to_char(std::string s, unsigned int start, unsigned int end, char* c) {
	unsigned int count = 0;
	for(unsigned int i = start; i<end; i++)
		if(!std::isspace(s[i])) c[count++] = s[i];
	c[count] = '\0';
}

std::vector<PDBatom> read_pdb_atoms(
                                    std::string filename
                                   )
{
	timeval runtime;
	start_timer(runtime);
	cout << "Reading pdb(qt) file [" << filename << "] ... ";
	std::vector<PDBatom> atoms;
	PDBatom current;
	std::ifstream file(filename);
	if(file.fail()){
		cout << "\nERROR: Can't open file.\n";
		exit(1);
	}
	std::string line;
	char tempstr[256];
	bool hetatm;
	while(std::getline(file, line))
	{
		sscanf(line.c_str(),"%255s",tempstr);
		hetatm = (strcmp(tempstr, "HETATM") == 0);
		if (hetatm || (strcmp(tempstr, "ATOM") == 0))
		{
			current.hetatm = hetatm;
			line.insert(54,1,' '); // add spaces to make reading coordinates easier
			line.insert(46,1,' ');
			line.insert(38,1,' ');
			sscanf(&line.c_str()[30], "%f %f %f", &(current.x), &(current.y), &(current.z));
			range_trim_to_char(line, 12, 16, current.name);
			range_trim_to_char(line, 17, 20, current.res_name);
			range_trim_to_char(line, 21, 22, current.chain_id);
			range_trim_to_char(line, 80, 83, current.atom_type); // reading atom type
			line[26]='\0'; // make sure res_id only 4 digits
			sscanf(&line.c_str()[22], "%d", &(current.res_id));
			atoms.push_back(current);
		}
	}
	file.close();
	cout << "Done, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";
	return atoms;
}


#endif // INCLUDED_PDB_READER

