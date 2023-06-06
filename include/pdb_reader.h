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
	char         name[5];      // "N"
	char         res_name[4];  // "SER"
	char         chain_id[2];  // "A"
	unsigned int res_id;       // 1
	float        x,y,z;        // -2.367, 4.481, -16.909
	char         atom_type[4]; // "N"
} PDBatom;


std::vector<ReceptorAtom> read_pdb_atoms(
                                         std::string filename
                                        )
{
	std::vector<PDBatom> atoms;
	PDBatom current;
	std::ifstream file(filename);
	if(file.fail()){
		cout << "Error: Can't open pdb(qt) file " << filename << ".\n";
		exit(1);
	}
	std::string line;
	char tempstr[256];
	while(std::getline(file, line))
	{
		sscanf(line.c_str(),"%255s",tempstr);
		if ((strcmp(tempstr, "HETATM") == 0) || (strcmp(tempstr, "ATOM") == 0))
		{
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
	return atoms;
}


#endif // INCLUDED_PDB_READER

