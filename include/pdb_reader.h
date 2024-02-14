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
	bool         hetatm;           // false (wether first field is "HETATM" or not)
	char         name[5];          // "N"
	char         alt_id;           // ' '
	char         res_name[4];      // "SER"
	char         chain_id;         // 'A'
	int          res_id;           // 1
	char         ins_id;           // ' '
	float        x,y,z;            // -2.367, 4.481, -16.909
	char         atom_type[4];     // "N"
	unsigned int atom_type_number; //  7 (recognized types are > 0; unrecognized 0)
} PDBatom;

typedef struct{
	float        x,y,z,w;
} point;

#define AD_RECOGNIZED_TYPES 30
inline unsigned int find_atom_type_number(char* type)
{
	const char* ad_atom_types[AD_RECOGNIZED_TYPES]               = {"C", "A", "CG", "G", "CX", "N", "NA", "NS", "NX", "O", "OA", "OS", "OX", "H", "HD", "HS", "S", "SA", "F", "P", "CL", "CA", "MG", "MN", "FE", "ZN", "BR", "I", "SI", "B"};
	const unsigned int ad_atom_type_numbers[AD_RECOGNIZED_TYPES] = { 6 ,  6 ,   6 ,  6 ,   6 ,  7 ,   7 ,   7 ,   7 ,  8 ,   8 ,   8 ,   8 ,  1 ,   1 ,   1 , 16 ,  16 ,  9 , 15 ,  17 ,  20 ,  12 ,  25 ,  26 ,  29 ,  35 , 53 ,  14 ,  5 };
	unsigned int tlen = strlen(type);
	if(tlen == 0) return 0;
	for(unsigned int i=0; i<AD_RECOGNIZED_TYPES; i++){
		if((type[0] == (ad_atom_types[i])[0]) &&
		   (type[1] == (ad_atom_types[i])[1])){ // even if the type string is one letter, we can compare trailing '\0' character
			return ad_atom_type_numbers[i];
		}
	}
	return 0;
}

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


// trim input string range -- careful: no range checking
static inline void range_trim_to_char(std::string s, unsigned int start, unsigned int end, char* c) {
	unsigned int count = 0;
	for(unsigned int i = start; i<end; i++)
		if(!std::isspace(s[i])) c[count++] = s[i];
	c[count] = '\0';
}

inline std::vector<PDBatom> read_pdb_atoms(
                                           std::string filename
                                          )
{
	timeval runtime;
	start_timer(runtime);
	std::vector<PDBatom> atoms;
	PDBatom current;
	std::ifstream file(filename);
	if(file.fail()){
		#pragma omp critical
		cout << "\nERROR: Can't open pdb(qt) file ["<< filename << "].\n";
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
			current.atom_type_number = find_atom_type_number(current.atom_type);
			current.ins_id = line[26];
			line[26] = '\0'; // make sure res_id only 4 digits
			sscanf(&line.c_str()[22], "%d", &(current.res_id));
			atoms.push_back(current);
		}
	}
	file.close();
	if(atoms.size()==0){
		#pragma omp critical
		cout << "\nERROR: Pdb(qt) file ["<< filename << "] contains no atoms.\n";
		exit(2);
	}
	// sort atoms by residue id
	std::sort(atoms.begin(), atoms.end(), [](PDBatom a, PDBatom b){return (a.res_id < b.res_id);});
	return atoms;
}

inline std::vector<point> map_pdb_points(
                                         std::vector<PDBatom> &map_atoms,
                                         char                  map_chain_id,
                                         std::vector<PDBatom> &align_atoms,
                                         char                  align_chain_id,
                                         int                   map_x_dim,
                                         int                   map_y_dim,
                                         int                   map_z_dim,
                                         fp_num                map_x_center,
                                         fp_num                map_y_center,
                                         fp_num                map_z_center,
                                         fp_num                grid_spacing
                                        )
{
	std::vector<point> point_mapping;
	std::vector<int> grid_ids;
	bool use_grid_box = (map_x_dim > 0) && (map_y_dim > 0) && (map_z_dim > 0);
	Vec3<fp_num> grid_dims;
	grid_dims.vec[0] = map_x_dim * grid_spacing; grid_dims.vec[1] = map_y_dim * grid_spacing; grid_dims.vec[2] = map_z_dim * grid_spacing;
	Vec3<fp_num> grid_start;
	grid_start.vec[0] = map_x_center - grid_dims.vec[0] * 0.5; grid_start.vec[1] = map_y_center - grid_dims.vec[1] * 0.5; grid_start.vec[2] = map_z_center - grid_dims.vec[2] * 0.5;
	Vec3<fp_num> location;
	int curr_resid;
	// Calculate heavy atom geometric center of atoms in grid box or of all if there is no grid box
	for(unsigned int i=0; i<align_atoms.size(); i++){
		curr_resid = align_atoms[i].res_id;
		// only focus on heavy atoms of large molecules in pdb(qt) with no alternative coordinates
		if((align_atoms[i].chain_id == align_chain_id) &&
		   (align_atoms[i].alt_id == ' ') &&
		   (align_atoms[i].hetatm == false) &&
		   (align_atoms[i].atom_type_number > 1))
		{
			location.vec[0] = align_atoms[i].x; location.vec[1] = align_atoms[i].y; location.vec[2] = align_atoms[i].z;
			if(!use_grid_box || point_in_box(location - grid_start, grid_dims)){
				grid_ids.push_back(i);
			} else{ // make sure to exclude whole residue
				// exclude already included atoms
				while((grid_ids.size()>0) && (align_atoms[grid_ids.back()].res_id == curr_resid))
					grid_ids.pop_back();
				// fast-forward to end of residue
				while((i+1<align_atoms.size()) && (align_atoms[i+1].res_id == curr_resid)) i++;
			}
		}
	}
	if(grid_ids.size() == 0){
		return point_mapping;
	}
	// compile a list of grid residue centers (and corresponding ids and how many atoms)
	std::vector<std::string> res_dict;
	std::vector<unsigned int> grid_res_start, grid_res_idx;
	std::vector<int> grid_res_dict, grid_res_dict_ids, grid_res_weight;
	std::vector<Vec3<fp_num>> grid_res_center;
	curr_resid = align_atoms[grid_ids[0]].res_id - 1;
	int found, w;
	for(unsigned int i = 0; i < grid_ids.size(); i++){
		location.vec[0] = align_atoms[grid_ids[i]].x;
		location.vec[1] = align_atoms[grid_ids[i]].y;
		location.vec[2] = align_atoms[grid_ids[i]].z;
		w = align_atoms[grid_ids[i]].atom_type_number;
		if(curr_resid  != align_atoms[grid_ids[i]].res_id){
			found = -1;
			for(unsigned int j=0; j<res_dict.size(); j++){
				if(res_dict[j].compare(align_atoms[grid_ids[i]].res_name)==0){
					found = j;
					break;
				}
			}
			if(found < 0){
				res_dict.push_back(align_atoms[grid_ids[i]].res_name);
				found = res_dict.size()-1;
			}
			grid_res_dict.push_back(found);
			grid_res_start.push_back(i);
			curr_resid = align_atoms[grid_ids[i]].res_id;
			grid_res_dict_ids.push_back(curr_resid);
			grid_res_center.push_back(location * w);
			grid_res_weight.push_back(w);
			grid_res_idx.push_back(grid_res_weight.size()-1);
		} else{
			grid_res_idx.push_back(grid_res_idx.back());
			grid_res_center[grid_res_idx.back()] += location * w;
			grid_res_weight[grid_res_idx.back()] += w;
		}
	}
	std::vector<unsigned int> map_res_start;
	std::vector<int> map_res_dict, map_res_dict_ids, map_res_weight;
	std::vector<Vec3<fp_num>> map_res_center;
	curr_resid = map_atoms[0].res_id - 1;
	for(unsigned int i = 0; i < map_atoms.size(); i++){
		if((map_atoms[i].chain_id == map_chain_id) &&
		   (map_atoms[i].alt_id   == ' ') &&
		   (map_atoms[i].hetatm   == false) &&
		   (map_atoms[i].atom_type_number > 1))
		{
			location.vec[0] = map_atoms[i].x;
			location.vec[1] = map_atoms[i].y;
			location.vec[2] = map_atoms[i].z;
			w = map_atoms[i].atom_type_number;
			if(curr_resid != map_atoms[i].res_id){
				found = -1;
				for(unsigned int j=0; j<res_dict.size(); j++){
					if(res_dict[j].compare(map_atoms[i].res_name)==0){
						found = j;
						break;
					}
				}
				if(map_res_dict.size() > 0){ // fill holes in map residues
					for(int j = 1; j < map_atoms[i].res_id - curr_resid; j++){
						map_res_dict.push_back(-1);
						map_res_dict_ids.push_back(-1);
						map_res_start.push_back(i);
						map_res_center.push_back(location * w);
						map_res_weight.push_back(w);
					}
				}
				map_res_dict.push_back(found);
				map_res_start.push_back(i);
				map_res_center.push_back(location * w);
				map_res_weight.push_back(w);
				curr_resid = map_atoms[i].res_id;
				map_res_dict_ids.push_back(curr_resid);
			} else{
				map_res_center.back() += location * w;
				map_res_weight.back() += w;
			}
		}
	}
	// find best (longest) match of grid residues in map residues
	unsigned int dict_start = 0;
	std::vector<unsigned int> residue_mapping, temp_map;
	for(unsigned int i = 0; i < map_res_dict.size(); i++){
		if(map_res_dict[i] == grid_res_dict[dict_start]){ // found a candidate starting position
			temp_map.clear();
			for(unsigned int j = dict_start; j < grid_res_dict.size(); j++){
				unsigned int rel_res_id = grid_res_dict_ids[j] - grid_res_dict_ids[dict_start];
				if(i + rel_res_id >= map_res_dict.size()) break;
				if(grid_res_dict[j] != map_res_dict[i + rel_res_id]) break;
				temp_map.push_back(j);
				temp_map.push_back(i+rel_res_id);
			}
			if(temp_map.size() > residue_mapping.size())
				residue_mapping.swap(temp_map);
		}
		if(i+1 == map_res_dict.size()){ // if we get to through the map atoms already ...
			if((residue_mapping.size() < grid_res_dict.size() - dict_start - 1)){ // ... and haven't found a long enough match
				i = 0; // keep trying
				dict_start++; // ... with one less grid residues
			}
		}
	}
	for(int i=map_res_dict.size(); i-->0;){ // go in reverse just in case
		if(map_res_dict[i] == grid_res_dict[dict_start]){ // found a candidate starting position
			temp_map.clear();
			for(unsigned int j=dict_start; j<grid_res_dict.size(); j++){
				int rel_res_id = grid_res_dict_ids[j] - grid_res_dict_ids[dict_start];
				if(i - rel_res_id >= 0) break;
				if(grid_res_dict[j] != map_res_dict[i - rel_res_id]) break;
				temp_map.push_back(j);
				temp_map.push_back(i-rel_res_id);
			}
			if(temp_map.size() > residue_mapping.size())
				residue_mapping.swap(temp_map);
		}
		if(i == 0){ // if we get to through the map atoms already ...
			if((residue_mapping.size() < grid_res_dict.size() - dict_start - 1)){ // ... and haven't found a long enough match
				i = map_res_dict.size(); // keep trying
				dict_start++; // ... with one less grid residues
			}
		}
	}
#if DEBUG_LEVEL>3
	for(unsigned int i=0; i<residue_mapping.size(); i+=2)
		cout << align_atoms[grid_ids[grid_res_start[residue_mapping[i]]]].res_name << " #" << align_atoms[grid_ids[grid_res_start[residue_mapping[i]]]].res_id << " -> "
		     << map_atoms[map_res_start[residue_mapping[i+1]]].res_name << " #" << map_atoms[map_res_start[residue_mapping[i+1]]].res_id << "\n";
	cout << "\n";
#endif
	grid_res_start.push_back(grid_ids.size());
	map_res_start.push_back(map_atoms.size());
	for(unsigned int i=0; i<grid_res_center.size(); i++)
		grid_res_center[i] /= grid_res_weight[i];
	for(unsigned int i=0; i<map_res_center.size(); i++)
		map_res_center[i] /= map_res_weight[i];
	std::vector<fp_num> grid_res_center_d2;
	std::vector<int> assignment;
	// match atoms from grid_ids to map atoms
	std::vector<int> map_match;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		location.vec[0] = align_atoms[grid_ids[i]].x;
		location.vec[1] = align_atoms[grid_ids[i]].y;
		location.vec[2] = align_atoms[grid_ids[i]].z;
		location       -= grid_res_center[grid_res_idx[i]];
		grid_res_center_d2.push_back(location * location); // <- distance^2 of current grid atom from its residue's center
		assignment.push_back(0);
		map_match.push_back(-1);
	}
	std::vector<unsigned int> grid_type;
	std::vector<unsigned int> candidates;
	std::vector<fp_num> cross_dist2;
	std::vector<bool> res_matched(grid_res_center.size(), false);
	unsigned curr_atom_type;
	for(unsigned int grid_res_id=0; grid_res_id<residue_mapping.size(); grid_res_id+=2){
		unsigned int r = residue_mapping[grid_res_id];
		unsigned int z = residue_mapping[grid_res_id+1];
		unsigned int i = grid_res_start[r];
		while(i < grid_res_start[r+1]){
			grid_type.clear();
			grid_type.push_back(i);
			curr_atom_type = align_atoms[grid_ids[i]].atom_type_number;
			curr_resid     = align_atoms[grid_ids[i]].res_id;
			int next_id = -1;
			while(++i < grid_res_start[r+1]){
				if(assignment[i] == 0){
					if(curr_atom_type == align_atoms[grid_ids[i]].atom_type_number){
						grid_type.push_back(i);
						assignment[i] = 1;
					} else next_id = (next_id < 0) ? i : next_id;
				}
			}
			i = (next_id < 0) ? grid_res_start[r+1] : next_id;
			// container for atoms with the same name is the one we're looking for (should mostly only be one ... but you never know)
			candidates.clear();
			for(unsigned int j = map_res_start[z]; j < map_res_start[z+1]; j++){
				if((map_atoms[j].alt_id   == ' ') &&
				   (map_atoms[j].hetatm   == false) &&
				   (map_atoms[j].atom_type_number > 1) &&
				   (map_atoms[j].atom_type_number == curr_atom_type))
				{
					candidates.push_back(j);
				}
			}
			if(candidates.size() == grid_type.size()){ // make sure we're mapping only residues with the same number of (heavy) atoms
				cross_dist2.clear();
				for(unsigned int k=0; k<candidates.size(); k++){
					for(unsigned int l=0; l<grid_type.size(); l++){
						location.vec[0] = map_atoms[candidates[k]].x;
						location.vec[1] = map_atoms[candidates[k]].y;
						location.vec[2] = map_atoms[candidates[k]].z;
						location -= map_res_center[z];
						cross_dist2.push_back(fabs((location * location) - grid_res_center_d2[grid_type[l]]));
					}
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
					// only add individual atoms with square atomic deviation <= 0.5 A^2
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
					grid_type.erase(grid_type.begin() + closest_l);
					candidates.erase(candidates.begin() + closest_k);
					res_matched[r] = true;
				}
			} else{ // if no candidate atoms are found for the given residue (i.e. only alt_id atoms or no atoms) exclude residue
				for(i = grid_res_start[r]; i < grid_res_start[r+1]; i++) grid_ids[i] = -1;
				i = grid_res_start[r+1];
			}
		}
	}
	point A, B;
	for(unsigned int i=0; i<grid_ids.size(); i++){
		if((grid_ids[i]<0) || (map_match[i]<0)) continue;
		A.x = align_atoms[grid_ids[i]].x; A.y = align_atoms[grid_ids[i]].y; A.z = align_atoms[grid_ids[i]].z; A.w = 1;//align_atoms[grid_ids[i]].atom_type_number;
		B.x = map_atoms[map_match[i]].x; B.y = map_atoms[map_match[i]].y; B.z = map_atoms[map_match[i]].z; B.w = A.w;
		point_mapping.push_back(A);
		point_mapping.push_back(B);
	}
	unsigned int matched = 0;
	for(unsigned int r=0; r<res_matched.size(); r++)
		matched += (res_matched[r] == true);
	if(matched < 4) point_mapping.clear(); // match atoms from at least three different residues
	return point_mapping;
}

inline fp_num* align_mapping(
                             std::vector<point>   &point_map,
                             fp_num               &rmsd,
                             std::stringstream    &output
                            )
{
	if(point_map.size() == 0){ // this really shouldn't happen but better say something if it were to ...
		rmsd = -1;
		return NULL;
	}
	Vec3<fp_num> location;
	Vec3<fp_num> center, map_center, grid_center;
	Mat33<fp_num> grid_rot, best_grid_rot;
	grid_center.vec[0] = 0; grid_center.vec[1] = 0; grid_center.vec[2] = 0;
	map_center.vec[0]  = 0; map_center.vec[1]  = 0; map_center.vec[2]  = 0;
	fp_num weight = 0;
	for(unsigned int i=0; i<point_map.size(); i+=2){
		grid_center.vec[0] += point_map[i].w * point_map[i].x;   grid_center.vec[1] += point_map[i].w * point_map[i].y;   grid_center.vec[2] += point_map[i].w * point_map[i].z;
		map_center.vec[0]  += point_map[i].w * point_map[i+1].x; map_center.vec[1]  += point_map[i].w * point_map[i+1].y; map_center.vec[2]  += point_map[i].w * point_map[i+1].z;
		weight += point_map[i].w;
	}
	grid_center /= weight;
	map_center  /= weight;
	output << "\t-> Grid center: (" << grid_center.V3Str(',',4) << ")\n";
	output << "\t-> Map center: (" << map_center.V3Str(',',4) << ")\n";
	fp_num best_rmsd = -1;
	for(unsigned int i=0; i<point_map.size(); i+=2){
		// align to respective centers
		point_map[i].x   -= grid_center.vec[0];
		point_map[i].y   -= grid_center.vec[1];
		point_map[i].z   -= grid_center.vec[2];
		point_map[i+1].x -= map_center.vec[0];
		point_map[i+1].y -= map_center.vec[1];
		point_map[i+1].z -= map_center.vec[2];
	}
	Mat33<double> B, BTB, BBT, U, V, M, pre_rot;
	for(unsigned int pre=0; pre<10; pre++){
		switch(pre){
			default:
			case 0: pre_rot.M3Zeros();
				break;
			case 1: // 90 degree rotation around x
				pre_rot.mat[0][0] = 1; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = 0; pre_rot.mat[1][2] = -1;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = 1; pre_rot.mat[2][2] = 0;
				break;
			case 2: // 90 degree rotation around y
				pre_rot.mat[0][0] = 0; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = 1;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = 1; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = -1; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = 0;
				break;
			case 3: // 90 degree rotation around z
				pre_rot.mat[0][0] = 0; pre_rot.mat[0][1] = -1; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = 1; pre_rot.mat[1][1] = 0; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = 1;
				break;
			case 4: // 180 degree rotation around x
				pre_rot.mat[0][0] = 1; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = -1; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = -1;
				break;
			case 5: // 180 degree rotation around y
				pre_rot.mat[0][0] = -1; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = 1; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = -1;
				break;
			case 6: // 180 degree rotation around z
				pre_rot.mat[0][0] = -1; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = -1; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = 1;
				break;
			case 7: // 270 degree rotation around x
				pre_rot.mat[0][0] = 1; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = 0; pre_rot.mat[1][2] = 1;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = -1; pre_rot.mat[2][2] = 0;
				break;
			case 8: // 270 degree rotation around y
				pre_rot.mat[0][0] = 0; pre_rot.mat[0][1] = 0; pre_rot.mat[0][2] = -1;
				pre_rot.mat[1][0] = 0; pre_rot.mat[1][1] = 1; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = 1; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = 0;
				break;
			case 9: // 270 degree rotation around z
				pre_rot.mat[0][0] = 0; pre_rot.mat[0][1] = 1; pre_rot.mat[0][2] = 0;
				pre_rot.mat[1][0] = -1; pre_rot.mat[1][1] = 0; pre_rot.mat[1][2] = 0;
				pre_rot.mat[2][0] = 0; pre_rot.mat[2][1] = 0; pre_rot.mat[2][2] = 1;
				break;
		}
		B.M3Zeros();
		weight = 0;
		for(unsigned int i=0; i<point_map.size(); i+=2){
			// align to respective centers
			location.vec[0] = point_map[i].x;
			location.vec[1] = point_map[i].y;
			location.vec[2] = point_map[i].z;
			location = pre_rot * location;
			center.vec[0] = point_map[i+1].x;
			center.vec[1] = point_map[i+1].y;
			center.vec[2] = point_map[i+1].z;
			// B = 1/N sum r_grid * r_map^T (outer product aka matrix multiplication)
			B.mat[0][0] += point_map[i].w * location.vec[0] * center.vec[0];
			B.mat[0][1] += point_map[i].w * location.vec[0] * center.vec[1];
			B.mat[0][2] += point_map[i].w * location.vec[0] * center.vec[2];
			
			B.mat[1][0] += point_map[i].w * location.vec[1] * center.vec[0];
			B.mat[1][1] += point_map[i].w * location.vec[1] * center.vec[1];
			B.mat[1][2] += point_map[i].w * location.vec[1] * center.vec[2];
			
			B.mat[2][0] += point_map[i].w * location.vec[2] * center.vec[0];
			B.mat[2][1] += point_map[i].w * location.vec[2] * center.vec[1];
			B.mat[2][2] += point_map[i].w * location.vec[2] * center.vec[2];
			weight += point_map[i].w;
		}
		B /= weight;
		BBT = B * B.M3Transpose();
		CVec3<double> cew = BBT.Eigenvalues();
		if(cew.Im()*cew.Im()>EPS){ // shouldn't happen
		#pragma omp critical
			cout << output.str() << "ERROR: BB^T eigenvalues need to be real.\n";
			exit(4);
		}
		Vec3<double> ew = cew.Re();
		if(fabs(ew.vec[0]) < fabs(ew.vec[1])){
			double tmp = ew.vec[0];
			ew.vec[0] = ew.vec[1];
			ew.vec[1] = tmp;
		}
		if(fabs(ew.vec[1]) < fabs(ew.vec[2])){
			double tmp = ew.vec[1];
			ew.vec[1] = ew.vec[2];
			ew.vec[2] = tmp;
		}
		if(fabs(ew.vec[0]) < fabs(ew.vec[1])){
			double tmp = ew.vec[0];
			ew.vec[0] = ew.vec[1];
			ew.vec[1] = tmp;
		}
		U = MGS(BBT.Eigenvectors(ew)); // make sure to ortho-normalize eigenvalues
		
		BTB = B.M3Transpose() * B;
		cew = BTB.Eigenvalues();
		if(cew.Im()*cew.Im()>EPS){ // shouldn't happen
			#pragma omp critical
			cout << output.str() << "ERROR: B^TB eigenvalues need to be real.\n";
			exit(5);
		}
		ew = cew.Re();
		if(fabs(ew.vec[0]) < fabs(ew.vec[1])){
			double tmp = ew.vec[0];
			ew.vec[0] = ew.vec[1];
			ew.vec[1] = tmp;
		}
		if(fabs(ew.vec[1]) < fabs(ew.vec[2])){
			double tmp = ew.vec[1];
			ew.vec[1] = ew.vec[2];
			ew.vec[2] = tmp;
		}
		if(fabs(ew.vec[0]) < fabs(ew.vec[1])){
			double tmp = ew.vec[0];
			ew.vec[0] = ew.vec[1];
			ew.vec[1] = tmp;
		}
		V = MGS(BTB.Eigenvectors(ew)); // make sure to ortho-normalize eigenvalues
		
		M.M3Eye();
		M.mat[2][2] = U.M3Det() * V.M3Det();
		grid_rot = (V * (M * U.M3Transpose())) * pre_rot;
		// calculate RMSD
		rmsd = 0;
		for(unsigned int i=0; i<point_map.size(); i+=2){
			location.vec[0] = point_map[i].x;
			location.vec[1] = point_map[i].y;
			location.vec[2] = point_map[i].z;
			center          = grid_rot * location;
			center.vec[0]  -= point_map[i+1].x;
			center.vec[1]  -= point_map[i+1].y;
			center.vec[2]  -= point_map[i+1].z;
			rmsd           += (center * center);
		}
		rmsd = sqrt(2.0*rmsd/point_map.size());
		if((rmsd < best_rmsd) || (best_rmsd < 0)){
			best_rmsd     = rmsd;
			best_grid_rot = grid_rot;
		}
	}
	grid_rot = best_grid_rot;
	rmsd     = best_rmsd;
	output << "\t-> Rotation matrix:\n";
	output.precision(4);
	output.setf(ios::fixed, ios::floatfield);
	output << "\t\t" << std::setw(9) << grid_rot.mat[0][0] << " " << std::setw(9) << grid_rot.mat[0][1] << " " << std::setw(9) << grid_rot.mat[0][2] << "\n";
	output << "\t\t" << std::setw(9) << grid_rot.mat[1][0] << " " << std::setw(9) << grid_rot.mat[1][1] << " " << std::setw(9) << grid_rot.mat[1][2] << "\n";
	output << "\t\t" << std::setw(9) << grid_rot.mat[2][0] << " " << std::setw(9) << grid_rot.mat[2][1] << " " << std::setw(9) << grid_rot.mat[2][2] << "\n";
	output.precision(3);
	output << "\t-> RMSD after alignment (" << (point_map.size()>>1) << " atoms): " << rmsd << " A\n";
	fp_num* grid_align = new fp_num[9 + 3 + 3];
	memcpy(grid_align, grid_rot.mat, 9 * sizeof(fp_num));
	memcpy(grid_align + 9, map_center.vec, 3 * sizeof(fp_num));
	memcpy(grid_align + 12, grid_center.vec, 3 * sizeof(fp_num));
	
	return grid_align;
}

inline fp_num* align_pdb_atoms(
                               std::string map_ligand,
                               std::string align_lig,
                               int         map_x_dim,
                               int         map_y_dim,
                               int         map_z_dim,
                               fp_num      map_x_center,
                               fp_num      map_y_center,
                               fp_num      map_z_center,
                               fp_num      grid_spacing,
                               bool        output_align_rec = false
                              )
{
	timeval runtime;
	start_timer(runtime);
	stringstream output;
	output << "Aligning density map receptor to grid receptor\n";
	std::vector<PDBatom> map_atoms, align_atoms;
	output << "\t-> Reading density map receptor [" << map_ligand << "]\n";
	map_atoms = read_pdb_atoms(map_ligand);
	if(align_lig.length() != 0){
		output << "\t-> Reading grid map receptor [" << align_lig << "]\n";
		align_atoms = read_pdb_atoms(align_lig);
	} else{
		#pragma omp critical
		cout << output.str() << "ERROR: No receptor specified in grid map files.\n";
		exit(1);
	}
	std::vector<char> map_chain_ids, grid_chain_ids;
	char curr_chain_id = map_atoms[0].chain_id;
	map_chain_ids.push_back(curr_chain_id);
	for(unsigned int i = 0; i < map_atoms.size(); i++){
		if((map_atoms[i].alt_id == ' ') &&
		   (map_atoms[i].hetatm == false) &&
		   (curr_chain_id != map_atoms[i].chain_id)){
			curr_chain_id = map_atoms[i].chain_id;
			bool found    = false;
			for(unsigned int j = 0; j < map_chain_ids.size(); j++){
				if(curr_chain_id == map_chain_ids[j]){
					found = true;
					break;
				}
			}
			if(!found) map_chain_ids.push_back(curr_chain_id);
		}
	}
	bool use_grid_box = (map_x_dim > 0) && (map_y_dim > 0) && (map_z_dim > 0);
	Vec3<fp_num> grid_dims;
	grid_dims.vec[0] = map_x_dim * grid_spacing; grid_dims.vec[1] = map_y_dim * grid_spacing; grid_dims.vec[2] = map_z_dim * grid_spacing;
	Vec3<fp_num> grid_start;
	grid_start.vec[0] = map_x_center - grid_dims.vec[0] * 0.5; grid_start.vec[1] = map_y_center - grid_dims.vec[1] * 0.5; grid_start.vec[2] = map_z_center - grid_dims.vec[2] * 0.5;
	Vec3<fp_num> location;
	curr_chain_id = '\0';
	for(unsigned int i=0; i<align_atoms.size(); i++){
		// only focus on heavy atoms of large molecules in pdb(qt) with no alternative coordinates
		if((align_atoms[i].alt_id==' ') &&
		   (align_atoms[i].hetatm == false) &&
		   (curr_chain_id != align_atoms[i].chain_id))
		{
			location.vec[0] = align_atoms[i].x; location.vec[1] = align_atoms[i].y; location.vec[2] = align_atoms[i].z;
			if(!use_grid_box || point_in_box(location - grid_start, grid_dims)){
				curr_chain_id = align_atoms[i].chain_id;
				bool found    = false;
				for(unsigned int j = 0; j < grid_chain_ids.size(); j++){
					if(curr_chain_id == grid_chain_ids[j]){
						found = true;
						break;
					}
				}
				if(!found) grid_chain_ids.push_back(curr_chain_id);
			}
		}
	}
	fp_num best_rmsd = -1;
	fp_num* best_align = new fp_num[9 + 3 + 3];
	std::stringstream best_out;
	for(unsigned int i=0; i<map_chain_ids.size(); i++){
		for(unsigned int j=0; j<grid_chain_ids.size(); j++){
			std::vector<point> point_map = map_pdb_points(
			                                              map_atoms,
			                                              map_chain_ids[i],
			                                              align_atoms,
			                                              grid_chain_ids[j],
			                                              map_x_dim,
			                                              map_y_dim,
			                                              map_z_dim,
			                                              map_x_center,
			                                              map_y_center,
			                                              map_z_center,
			                                              grid_spacing
			                                             );
			fp_num rmsd = -1;
			std::stringstream align_out;
			fp_num* grid_align = align_mapping(
			                                   point_map,
			                                   rmsd,
			                                   align_out
			                                  );
			if((rmsd >= 0) && ((rmsd < best_rmsd) || (best_rmsd < 0))){
				best_rmsd = rmsd;
				memcpy(best_align, grid_align, (9+3+3)*sizeof(fp_num));
				best_out.swap(align_out);
			}
			delete[] grid_align;
		}
	}
	if(best_rmsd < 0){
		#pragma omp critical
		cout << output.str() << "ERROR: Could not find enough common atoms between density and grid map receptor.\n";
		exit(1);
	}
	output << best_out.str();
	output.precision(3);
	output << "<- Finished alignment, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";
	if(output_align_rec){
		std::size_t ext  = align_lig.find_last_of(".");
		string filename  = align_lig.substr(0, ext) + "_ALIGNED_TO_";
		std::size_t mext = map_ligand.find_last_of(".");
		if(mext==std::string::npos) mext = map_ligand.size();
		for(unsigned int i=0; i<mext; i++) filename += (map_ligand[i]!='.') && (map_ligand[i]!='/') ? map_ligand[i] : '_';
		filename += ".pdb";
		std::ofstream align_file(filename);
		if(align_file.fail()){
			cout << "Error: Can't open aligned grid map receptor output file " << filename << ".\n";
			exit(1);
		}
		Vec3<fp_num> map_center, grid_center;
		Mat33<fp_num> grid_rot;
		memcpy(grid_rot.mat, best_align, 9 * sizeof(fp_num));
		memcpy(map_center.vec, best_align + 9, 3 * sizeof(fp_num));
		memcpy(grid_center.vec, best_align + 12, 3 * sizeof(fp_num));
		for(unsigned int i=0; i<align_atoms.size(); i++){
			location.vec[0]  = align_atoms[i].x - grid_center.vec[0];
			location.vec[1]  = align_atoms[i].y - grid_center.vec[1];
			location.vec[2]  = align_atoms[i].z - grid_center.vec[2];
			location = grid_rot * location;
			align_file.precision(3);
			align_file.fill(' ');
			align_file.setf(ios::fixed, ios::floatfield);
			align_file << "ATOM  ";
			align_file << std::setw(5) << i+1 << "  ";
			std::string str = align_atoms[i].name;
			str.resize(4,' ');
			align_file << str << std::setw(3) << align_atoms[i].res_name << " ";
			align_file << align_atoms[i].chain_id;
			align_file << std::setw(4) << align_atoms[i].res_id << "    ";
			align_file << std::setw(8) << location.vec[0]+map_center.vec[0] << std::setw(8) << location.vec[1]+map_center.vec[1] << std::setw(8) << location.vec[2]+map_center.vec[2];
			align_file << "  1.00  0.00          " << std::setw(2) << align_atoms[i].atom_type << "\n";
		}
		align_file.close();
	}
	#pragma omp critical
	cout << output.str();
	return best_align;
}

#endif // INCLUDED_PDB_READER

