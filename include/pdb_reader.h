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
	char         alt_id;       // ' '
	char         res_name[4];  // "SER"
	char         chain_id;     // 'A'
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
                          Vec3<fp_num>         &map_center,
                          Vec3<fp_num>         &grid_center
                         )
{
	timeval runtime;
	start_timer(runtime);
	cout << "Aligning map receptor to grid receptor\n";
	bool use_grid_box = (map_x_dim > 0) && (map_y_dim > 0) && (map_z_dim > 0);
	Vec3<fp_num> grid_dims;
	grid_dims.vec[0] = map_x_dim * grid_spacing; grid_dims.vec[1] = map_y_dim * grid_spacing; grid_dims.vec[2] = map_z_dim * grid_spacing;
	Vec3<fp_num> grid_start;
	grid_start.vec[0] = map_x_center - grid_dims.vec[0] * 0.5; grid_start.vec[1] = map_y_center - grid_dims.vec[1] * 0.5; grid_start.vec[2] = map_z_center - grid_dims.vec[2] * 0.5;
	Vec3<fp_num> location;
	std::vector<int> grid_ids;
	unsigned int curr_resid;
	// Calculate heavy atom geometric center of atoms in grid box or of all if there is no grid box
	for(unsigned int i=0; i<grid_atoms.size(); i++){
		curr_resid = grid_atoms[i].res_id;
		// only focus on heavy atoms of large molecules in pdb(qt) with no alternative coordinates
		if((grid_atoms[i].alt_id==' ') && (grid_atoms[i].hetatm == false) && (grid_atoms[i].atom_type[0] != 'H'))
		{
			if(!use_grid_box || point_in_box(location - grid_start, grid_dims)){
				grid_ids.push_back(i);
			} else{ // make sure to exclude whole residue
				// exclude already included atoms
				while((grid_ids.size()>0) && (grid_atoms[grid_ids.back()].res_id == curr_resid))
					grid_ids.pop_back();
				// fast-forward to end of residue
				while((i+1<grid_atoms.size()) && (grid_atoms[i+1].res_id == curr_resid)) i++;
			}
		}
	}
	if(grid_ids.size() == 0){ // shouldn't happen
		cout << "ERROR: No receptor atoms inside grid box. This likely means the wrong receptor is listed in the grid file.\n";
		exit(1);
	}
	// compile a list of grid residue centers (and corresponding ids and how many atoms)
	std::vector<unsigned int> grid_res_start, grid_res_idx, grid_res_num;
	std::vector<Vec3<fp_num>> grid_res_center;
	curr_resid         = grid_atoms[grid_ids[0]].res_id + 1;
	char curr_chain_id = grid_atoms[grid_ids[0]].chain_id + 1;
	for(unsigned int i = 0; i < grid_ids.size(); i++){
		location.vec[0] = grid_atoms[grid_ids[i]].x;
		location.vec[1] = grid_atoms[grid_ids[i]].y;
		location.vec[2] = grid_atoms[grid_ids[i]].z;
		if((curr_resid != grid_atoms[grid_ids[i]].res_id) ||
		   (curr_chain_id != grid_atoms[grid_ids[i]].chain_id)){
			grid_res_start.push_back(i);
			curr_resid = grid_atoms[grid_ids[i]].res_id;
			curr_chain_id = grid_atoms[grid_ids[i]].chain_id;
			grid_res_center.push_back(location);
			grid_res_num.push_back(1);
			grid_res_idx.push_back(grid_res_num.size()-1);
		} else{
			grid_res_idx.push_back(grid_res_idx.back());
			grid_res_center[grid_res_idx.back()] += location;
			grid_res_num[grid_res_idx.back()]++;
		}
	}
	grid_res_start.push_back(grid_ids.size());
	for(unsigned int i=0; i<grid_res_center.size(); i++)
		grid_res_center[i] /= grid_res_num[i];
	std::vector<fp_num> grid_res_center_d2;
	std::vector<int> assignment;
	// match atoms from grid_ids to map atoms
	std::vector<int> map_match;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		location.vec[0] = grid_atoms[grid_ids[i]].x;
		location.vec[1] = grid_atoms[grid_ids[i]].y;
		location.vec[2] = grid_atoms[grid_ids[i]].z;
		location       -= grid_res_center[grid_res_idx[i]];
		grid_res_center_d2.push_back(location * location); // <- distance^2 of current grid atom from its residue's center
		assignment.push_back(0);
		map_match.push_back(-1);
	}
	std::vector<unsigned int> grid_type;
	std::vector<unsigned int> candidates;
	std::vector<fp_num> cross_dist2;
	char* curr_atom_type;
	for(unsigned int r=0; r<grid_res_num.size(); r++){ // go over all residues in grid
		unsigned int i = grid_res_start[r];
		while(i < grid_res_start[r+1]){
			grid_type.clear();
			grid_type.push_back(i);
			curr_atom_type = grid_atoms[grid_ids[i]].atom_type;
			curr_chain_id  = grid_atoms[grid_ids[i]].chain_id;
			curr_resid     = grid_atoms[grid_ids[i]].res_id;
			int next_id = -1;
			while(++i < grid_res_start[r+1]){
				if(assignment[i] == 0){
					if(strcmp(curr_atom_type, grid_atoms[grid_ids[i]].atom_type) == 0){
						grid_type.push_back(i);
						assignment[i] = 1;
					} else next_id = (next_id < 0) ? i : next_id;
				}
			}
			i = (next_id < 0) ? grid_res_start[r+1] : next_id;
			// calculate matching map residue's center on the fly
			Vec3<fp_num> map_res_center(fp_num(0));
			unsigned int map_res_count = 0;
			// container for atoms with the same name is the one we're looking for (should mostly only be one ... but you never know)
			candidates.clear();
			for(unsigned int j=0; j<map_atoms.size(); j++){
				if((map_atoms[j].alt_id   == ' ') &&
				   (map_atoms[j].hetatm   == false) &&
				   (map_atoms[j].atom_type[0] != 'H') &&
				   (map_atoms[j].res_id   == curr_resid) &&
				   (map_atoms[j].chain_id == curr_chain_id)){
					map_res_center.vec[0] += map_atoms[j].x;
					map_res_center.vec[1] += map_atoms[j].y;
					map_res_center.vec[2] += map_atoms[j].z;
					map_res_count++;
					if(strcmp(map_atoms[j].atom_type, curr_atom_type) == 0){
						candidates.push_back(j);
					}
				}
			}
			if(candidates.size() == grid_type.size()){ // make sure we're mapping only residues with the same number of atoms
				map_res_center /= map_res_count;
				cross_dist2.clear();
#if DEBUG_LEVEL>4
				cout << "         ";
				for(unsigned int l=0; l<grid_type.size(); l++)
					cout << std::setw(9) << grid_ids[grid_type[l]]+1;
				cout << "\n";
#endif
				for(unsigned int k=0; k<candidates.size(); k++){
#if DEBUG_LEVEL>4
					cout << std::setw(5) << candidates[k]+1 << "    ";
#endif
					for(unsigned int l=0; l<grid_type.size(); l++){
						location.vec[0] = map_atoms[candidates[k]].x;
						location.vec[1] = map_atoms[candidates[k]].y;
						location.vec[2] = map_atoms[candidates[k]].z;
						location -= map_res_center;
						cross_dist2.push_back(fabs((location * location) - grid_res_center_d2[grid_type[l]]));
#if DEBUG_LEVEL>4
						cout.precision(4);
						cout << std::setw(9) << cross_dist2.back();
#endif
					}
#if DEBUG_LEVEL>4
					cout << "\n";
#endif
				}
				while((candidates.size()>0) && (grid_type.size()>0)){
					unsigned int closest_k, closest_l;
					fp_num closest_d2 = 1e8;
					for(unsigned int k=0; k<candidates.size(); k++){
						for(unsigned int l=0; l<grid_type.size(); l++){
							if(cross_dist2[k*grid_type.size()+l] < closest_d2){
								closest_d2 = cross_dist2[k*grid_type.size()+l];
								closest_k = k;
								closest_l = l;
							}
						}
					}
					// only accept residues with atomic deviation <= 0.5 A
					if(cross_dist2[closest_k*grid_type.size()+closest_l] > 0.5){
						for(i = grid_res_start[r]; i < grid_res_start[r+1]; i++) grid_ids[i] = -1;
						i = grid_res_start[r+1];
						break;
					}
					unsigned int cs = cross_dist2.size();
					for(unsigned int k=0; k<candidates.size(); k++)
						for(unsigned int l=0; l<grid_type.size(); l++)
							if((closest_k!=k) && (closest_l!=l)) cross_dist2.push_back(cross_dist2[k*grid_type.size()+l]);
					cross_dist2.erase(cross_dist2.begin(), cross_dist2.begin() + cs);
					map_match[grid_type[closest_l]] = candidates[closest_k];
#if DEBUG_LEVEL>3
					cout << grid_ids[grid_type[closest_l]]+1 << " -> " << candidates[closest_k]+1 << "\n";
#endif
					grid_type.erase(grid_type.begin() + closest_l);
					candidates.erase(candidates.begin() + closest_k);
				}
			} else{ // if no candidate atoms are found for the given residue (i.e. only alt_id atoms or no atoms) exclude residue
				for(i = grid_res_start[r]; i < grid_res_start[r+1]; i++) grid_ids[i] = -1;
				i = grid_res_start[r+1];
			}
		}
	}
	// align map to residue center closest to geometric center
	Vec3<fp_num> center(0.0);
	unsigned int count = 0;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		if(grid_ids[i]<0) continue;
		location.vec[0] = grid_atoms[i].x; location.vec[1] = grid_atoms[i].y; location.vec[2] = grid_atoms[i].z;
		center += location;
		count++;
	}
	center /= count;
	unsigned int closest = 0;
	fp_num closest_dist2 = 1e8;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		if(grid_ids[i]<0) continue;
		location.vec[0] = grid_atoms[grid_ids[i]].x;
		location.vec[1] = grid_atoms[grid_ids[i]].y;
		location.vec[2] = grid_atoms[grid_ids[i]].z;
		location -= center;
		fp_num d2 = location*location;
		if(d2 < closest_dist2){
			closest_dist2 = d2;
			closest = i;
		}
	}
	grid_center.vec[0] = 0; grid_center.vec[1] = 0; grid_center.vec[2] = 0;
	map_center.vec[0]  = 0; map_center.vec[1]  = 0; map_center.vec[2]  = 0;
	count = 0;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		if(grid_ids[i]<0) continue;
		if((grid_atoms[grid_ids[i]].res_id   == grid_atoms[grid_ids[closest]].res_id) &&
		   (grid_atoms[grid_ids[i]].chain_id == grid_atoms[grid_ids[closest]].chain_id))
		{
			grid_center.vec[0] += grid_atoms[grid_ids[i]].x; grid_center.vec[1] += grid_atoms[grid_ids[i]].y; grid_center.vec[2] += grid_atoms[grid_ids[i]].z;
			map_center.vec[0]  += map_atoms[map_match[i]].x; map_center.vec[1]  += map_atoms[map_match[i]].y; map_center.vec[2]  += map_atoms[map_match[i]].z;
			count++;
		}
	}
	grid_center /= count;
	map_center  /= count;
	cout << "\t-> Grid center: (" << grid_center.V3Str(',',4) << ")\n";
	cout << "\t-> Map center:  (" << map_center.V3Str(',',4) << ")\n";
	Mat33<double> B, BTB, BBT, U, V, M;
	B.M3Zeros();
	// Calculate gyration tensors for both
	count = 0;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		if(grid_ids[i]<0) continue;
		// align to respective centers
		grid_atoms[grid_ids[i]].x -= grid_center.vec[0];
		grid_atoms[grid_ids[i]].y -= grid_center.vec[1];
		grid_atoms[grid_ids[i]].z -= grid_center.vec[2];
		map_atoms[map_match[i]].x -= map_center.vec[0];
		map_atoms[map_match[i]].y -= map_center.vec[1];
		map_atoms[map_match[i]].z -= map_center.vec[2];
		// B = 1/N sum r_grid * r_map^T (outer product aka matrix multiplication)
		B.mat[0][0] += grid_atoms[grid_ids[i]].x * map_atoms[map_match[i]].x;
		B.mat[0][1] += grid_atoms[grid_ids[i]].x * map_atoms[map_match[i]].y;
		B.mat[0][2] += grid_atoms[grid_ids[i]].x * map_atoms[map_match[i]].z;
		
		B.mat[1][0] += grid_atoms[grid_ids[i]].y * map_atoms[map_match[i]].x;
		B.mat[1][1] += grid_atoms[grid_ids[i]].y * map_atoms[map_match[i]].y;
		B.mat[1][2] += grid_atoms[grid_ids[i]].y * map_atoms[map_match[i]].z;
		
		B.mat[2][0] += grid_atoms[grid_ids[i]].z * map_atoms[map_match[i]].x;
		B.mat[2][1] += grid_atoms[grid_ids[i]].z * map_atoms[map_match[i]].y;
		B.mat[2][2] += grid_atoms[grid_ids[i]].z * map_atoms[map_match[i]].z;
		count++;
	}
	B /= count;
	
	BBT = B * B.M3Transpose();
	CVec3<double> cew = BBT.Eigenvalues();
	if(cew.Im()*cew.Im()>EPS){ // shouldn't happen
		cout << "ERROR: BB^T eigenvalues need to be real.\n";
		exit(2);
	}
	Vec3<double> ew = cew.Re();
	U = BBT.Eigenvectors(ew, true); // make sure to normalize eigenvalues

	BTB = B.M3Transpose() * B;
	cew = BTB.Eigenvalues();
	if(cew.Im()*cew.Im()>EPS){ // shouldn't happen
		cout << "ERROR: B^TB eigenvalues need to be real.\n";
		exit(3);
	}
	ew = cew.Re();
	V = BTB.Eigenvectors(ew, true); // make sure to normalize eigenvalues
	M.mat[2][2] = U.M3Det() * V.M3Det();
	Mat33<fp_num> result;
	result = U * (M * V.M3Transpose());
	cout << "\t-> Rotation matrix:\n";
	cout.precision(4);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "\t\t" << std::setw(9) << result.mat[0][0] << " " << std::setw(9) << result.mat[1][0] << " " << std::setw(9) << result.mat[2][0] << "\n";
	cout << "\t\t" << std::setw(9) << result.mat[0][1] << " " << std::setw(9) << result.mat[1][1] << " " << std::setw(9) << result.mat[2][1] << "\n";
	cout << "\t\t" << std::setw(9) << result.mat[0][2] << " " << std::setw(9) << result.mat[1][2] << " " << std::setw(9) << result.mat[2][2] << "\n";
	// calculate RMSD
	fp_num rmsd = 0;
	count = 0;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		if(grid_ids[i]<0) continue;
		count++;
		location.vec[0] = map_atoms[map_match[i]].x;
		location.vec[1] = map_atoms[map_match[i]].y;
		location.vec[2] = map_atoms[map_match[i]].z;
		location = result * location;
		center.vec[0] = grid_atoms[grid_ids[i]].x;
		center.vec[1] = grid_atoms[grid_ids[i]].y;
		center.vec[2] = grid_atoms[grid_ids[i]].z;
#if DEBUG_LEVEL>3
		cout.precision(3);
		cout.fill(' ');
		cout.setf(ios::fixed, ios::floatfield);
		cout << "ATOM  ";
		cout << std::setw(5) << count << "  ";
		std::string str = map_atoms[map_match[i]].name;
		str.resize(4,' ');
		cout << str << std::setw(3) << map_atoms[map_match[i]].res_name << " ";
		cout << map_atoms[map_match[i]].chain_id;
		cout << std::setw(4) << map_atoms[map_match[i]].res_id << "    ";
		cout << std::setw(8) << location.vec[0]+grid_center.vec[0] << std::setw(8) << location.vec[1]+grid_center.vec[1] << std::setw(8) << location.vec[2]+grid_center.vec[2];
		cout << "  1.00  0.00          " << std::setw(2) << map_atoms[map_match[i]].atom_type << "\n";
#endif
		location -= center;
		rmsd += (location * location);
	}
	cout.precision(3);
	cout << "\t-> RMSD after alignment (" << count << " atoms): " << sqrt(rmsd/count) << " A\n";
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
			current.alt_id = line[16];
			range_trim_to_char(line, 17, 20, current.res_name);
			current.chain_id = line[21];
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

