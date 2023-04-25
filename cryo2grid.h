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


#ifndef INCLUDED_MAP_READER
#define INCLUDED_MAP_READER

#define DSN6_BLOCKSIZE 512

enum map_types_supported{
	automatic = 0,
	dsn6      = 1,
	dsn6_swap = 2,
	brix      = 3,
	mrc       = 4
};

enum grid_map_write_modes{
	write_nothing  = 0,
	write_grid_ad4 = 1,
	write_grid_mrc = 2
};

#define MAPEPS 1e-4

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

//                              start     extent     grid      cell axes    angles       origin     min, max, avg, std
const int mrc_offsets[22] = {16, 20, 24, 0, 4, 8, 28, 32, 36, 40, 44, 48, 52, 56, 60, 196, 200, 204, 76, 80, 84, 216};

static unsigned short static_one = 1;
#define HOST_LITTLE_ENDIAN (*(unsigned char*)&static_one == 1)

#define swap_32bit(byte_data) (HOST_LITTLE_ENDIAN ? (int)((*((unsigned char*)byte_data)<<24) | (*((unsigned char*)byte_data+1)<<16) | (*((unsigned char*)byte_data+2)<<8) | (*((unsigned char*)byte_data+3))) : (int)(*((unsigned char*)byte_data) | (*((unsigned char*)byte_data+1)<<8)  |(*((unsigned char*)byte_data+2)<<16) | (*((unsigned char*)byte_data+3)<<24)))
#define read_32bit(byte_data) (endian_swap ? &(tempval = swap_32bit(byte_data)) : &(tempval = *(reinterpret_cast<int*>(byte_data))))

#define short_swap(byte_data) (HOST_LITTLE_ENDIAN ? (short int)((*((unsigned char*)byte_data)<<8) | *((unsigned char*)byte_data+1)) : (short int)(*((unsigned char*)byte_data) | (*((unsigned char*)byte_data+1)<<8)))
#define read_short(byte_data) (endian_swap ? short_swap(byte_data) : *(reinterpret_cast<short int*>(byte_data)))

#define xyz_idx(x,y,z) ((x) + (y)*x_dim + (z)*xy_stride)

// McKie & McKie: Essentials of Crystallography, Blackwell Scientific, Oxford (1986)
// (r_x)   [ a  b*cos(gamma)              c*cos(beta) ]   (f_x)
// (r_y) = [ 0  b*sin(gamma)                      c*n ] * (f_y)
// (r_z)   [ 0       0        c*sqrt(sin^2(beta)-n^2) ]   (f_z)
//
// n = (cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)
//
// => r_x = f_x*a + f_y*b*cos(gamma) + f_z*c*cos(beta)
//    r_y =         f_y*b*sin(gamma) + f_z*c*n
//    r_z =                            f_z*c*sqrt(sin^2(beta)-n^2)
//
#define x_f2c(af_x, bf_y, cf_z) (af_x + bf_y*cos_gamma + cf_z*cos_beta)
#define y_f2c(bf_y, cf_z)       (bf_y*sin_gamma + cf_z*n)
#define z_f2c(cf_z)             (cf_z*sqrt_factor)
#define f2c(f_x, f_y, f_z) { f_x *= a_unit; f_y *= b_unit; f_z *= c_unit; f_x = x_f2c(f_x, f_y, f_z); f_y = y_f2c(f_y, f_z); f_z = z_f2c(f_z); }
// inverse to go from cartesian to fractional:
// (f_x)   [ 1/a -1/a*cos(gamma)/sin(gamma)  1/a*(n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2) ]   (r_x)
// (f_y) = [ 0             1/b*1/sin(gamma)                      -1/b*1/sin(gamma)*n*1/sqrt(sin^2(beta)-n^2) ] * (r_y)
// (f_z)   [ 0               0                                                   1/c*1/sqrt(sin^2(beta)-n^2) ]   (r_z)
//
// => f_x = [r_x - r_y*cos(gamma)/sin(gamma) + r_z * (n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2)]/a
//    f_y = [      r_y*1/sin(gamma)          - r_z*1/sin(gamma)*n/sqrt(sin^2(beta)-n^2)]/b
//    f_z = [                                  r_z*1/sqrt(sin^2(beta)-n^2)]/c
#define a_c2f(r_x, r_y, r_z) ((r_x - r_y*cos_inv_sin_gamma + r_z*long_inv_term)*inv_a_unit)
#define b_c2f(r_y, r_z)      (r_y*inv_sin_gamma_b - r_z*n_inv_sin_sqrt_b)
#define c_c2f(r_z)           (r_z*inv_sqrt_factor_c)
#define c2f(r_x, r_y, r_z) { r_x = a_c2f(r_x, r_y, r_z); r_y = b_c2f(r_y, r_z); r_z = c_c2f(r_z); }


inline char* find_block(char* &brix_string)
{
	while(*brix_string == ' ') brix_string++; // remove white-space
	char* start = brix_string; // number begins here
	while(*brix_string != ' ') brix_string++; // find next space
	*brix_string++ = '\0'; // end number string and advance to next spot
	return start;
}

inline fp_num brix_number_float(char* &brix_string)
{
	return stof(find_block(brix_string));
}

inline int brix_number_int(char* &brix_string)
{
	return stoi(find_block(brix_string));
}

inline bool brix_entry(char* &brix_string, const char* name, bool report_error = true)
{
	bool success = true;
	while(*brix_string == ' ') brix_string++; // remove white-space
	unsigned int count = 0;
	unsigned int nlen  = strlen(name);
	while(*brix_string != ' '){ // find next space ...
		if(count >= nlen){
			success = false;
			break;
		}
		if(toupper(name[count])!=toupper(*brix_string)){
			success = false;
			break;
		}
		brix_string++;
		count++;
	}
	if(report_error){
		if(!success){
			cout << "\nERROR: Brix \"" << name << "\"entry not found.\n";
			exit(37);
		}
		*brix_string++ = '\0'; // end number string and advance to next spot
	}
	return success;
}

#define next_entry(name) { if(map_type == brix) brix_entry(header, name, true); }

#define read_entry(header,type) ((map_type < 3) ? read_short(header+2*data_count) : ((map_type < 4) ? brix_number_##type(header) : *(reinterpret_cast<type*>(read_32bit(header+mrc_offsets[data_count])))));data_count++

#define read_mrc(result, data) {                                                  \
	switch(mrc_mode){                                                         \
		case 0: /* signed 8-bit  char */                                  \
		        result = *(data);                                         \
		        data_count++;                                             \
		        break;                                                    \
		case 1: /* signed 16-bit short int */                             \
		        result = read_short(data);                                \
		        data_count+=2;                                            \
		        break;                                                    \
		default:                                                          \
		case 2: /* 32-bit float */                                        \
		        result = *(reinterpret_cast<float*>(read_32bit(data)));   \
		        data_count+=4;                                            \
		        break;                                                    \
	}                                                                         \
}

inline int determine_map_type(char* &header)
{
	int map_type = automatic;
	bool BRIX = false;
	if((header[208]=='M') && (header[209]=='A') && (header[210]=='P')){
		map_type = mrc;
		return map_type;
	}
	
	short int norm       = *(reinterpret_cast<short int*>(header+36));
	bool endian_swap     = (norm != 100);
	norm                 = read_short(header+36);
	if(norm == 100){
		map_type     = dsn6;
		if(endian_swap) map_type = dsn6_swap;
	}
	if(brix_entry(header, ":-)",false)){
		BRIX = true;
	}
	if(BRIX && (norm == 100)){ // found the one map file whose origin encodes a smiley but that isn't a BRIX file
		BRIX=false;
	}
	if(BRIX){
		map_type    = brix;
	}
	return map_type;
}

inline std::string num2str(fp_num num)
{
	unsigned int l = fabs(num);
	unsigned int decimals = (fabs(num)-l)*1000 + 0.5;
	return (num<0?"-":"") + to_string(l) + "." + (decimals<100?"0":"") + (decimals<10?"0":"") + to_string(decimals);
}

inline void write_grid_map_ad4(
                               fp_num*       grid_map,
                               std::string  &filename,
                               unsigned int  map_x_dim,
                               unsigned int  map_y_dim,
                               unsigned int  map_z_dim,
                               fp_num        map_x_center,
                               fp_num        map_y_center,
                               fp_num        map_z_center,
                               fp_num        grid_spacing,
                               bool          set_extension = true
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
                               fp_num*       grid_map,
                               std::string  &filename,
                               unsigned int  map_x_dim,
                               unsigned int  map_y_dim,
                               unsigned int  map_z_dim,
                               fp_num        map_x_center,
                               fp_num        map_y_center,
                               fp_num        map_z_center,
                               fp_num        grid_spacing,
                               fp_num        rho_min       = -1,
                               fp_num        rho_max       =  1,
                               bool          set_extension = true
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
	if((rho_min <= rho_max) && (std::min(rho_min, rho_max) >= 0)){
		header.val_min = rho_min;
		header.val_max = rho_max;
	}
	header.mach_str = HOST_LITTLE_ENDIAN ? 17476 : 4369;
	
	map_file.write(reinterpret_cast<char*>(&header),sizeof(header));
	map_file.write(reinterpret_cast<char*>(grid_map), (map_x_dim + 1) * (map_y_dim + 1) * (map_z_dim + 1) * sizeof(float));
	
	map_file.close();
}

inline void write_grid(
                       fp_num*       grid_map,
                       std::string  &filename,
                       int           write_mode = write_grid_ad4
                      )
{
	timeval runtime;
	start_timer(runtime);
	switch(write_mode){
		case write_grid_ad4: write_grid_map_ad4(
		                                        grid_map + 9,
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
		                                        grid_map + 9,
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

inline std::vector<fp_num> read_map_to_grid(
                                            std::string  &filename,
                                                     int  map_type,
                                            unsigned int  map_x_dim,
                                            unsigned int  map_y_dim,
                                            unsigned int  map_z_dim,
                                            fp_num        map_x_center,
                                            fp_num        map_y_center,
                                            fp_num        map_z_center,
                                            fp_num        grid_spacing
                                           )
{
	timeval runtime;
	start_timer(runtime);
	
	std::vector<fp_num> densities;
	
	cout << "Reading map file [" << filename << "]\n";
	std::ifstream map_file(filename, std::ifstream::binary);
	if(map_file.fail()){
		cout << "\nERROR: Can't open map file " << filename << ".\n";
		exit(2);
	}
	std::streamoff filesize  = map_file.tellg();
	                           map_file.seekg(0, std::ios::end);
	filesize                 = map_file.tellg() - filesize;
	                           map_file.seekg(0, std::ios::beg);
	char map_header[DSN6_BLOCKSIZE];
	if(!map_file.read(map_header, DSN6_BLOCKSIZE)){
		cout << "\nERROR: Can't reader header.\n";
		exit(3);
	}
	char* header = map_header;
	if(map_type == automatic)
		map_type         = determine_map_type(header);
	unsigned int norm        = 1;
	unsigned int header_end  = DSN6_BLOCKSIZE;
	bool endian_swap         = false;
	int mrc_mode             = -1;
	int tempval;
	switch(map_type){
		case dsn6_swap: endian_swap = true;
		case dsn6:      cout << "\t-> DSN6";
		                norm = 100;
		                break;
		case brix:      cout << "\t-> BRIX";
		                break;
		case mrc:       cout << "\t-> MRC";
		                header_end  = 1024;
		                if(*((char*)&norm) != (*header || *(header+1))) endian_swap = true;
		                header_end += *(reinterpret_cast<int*>(read_32bit(header+92))); // add bytes of the extended header
		                mrc_mode    = *(reinterpret_cast<int*>(read_32bit(header+12)));
		                break;
		default:
		case automatic: cout << "\nERROR: Unknown map file type.\n";
		                exit(42);
	}
	cout << (endian_swap ?" endian-swapped":"") << ", file size: " << filesize << "\n";
	
	unsigned int data_count  = 0;
	next_entry("origin");
	fp_num x_start           = read_entry(header, int);
	fp_num y_start           = read_entry(header, int);
	fp_num z_start           = read_entry(header, int);
	next_entry("extent");
	unsigned int x_dim       = read_entry(header, int);
	unsigned int y_dim       = read_entry(header, int);
	unsigned int z_dim       = read_entry(header, int);
	unsigned long xy_stride  = x_dim * y_dim;
	unsigned int f_x_stride  = (((x_dim&7)>0) + (x_dim>>3)) << 6;
	unsigned int f_xy_stride = f_x_stride * (((y_dim&7)>0) + (y_dim>>3)) << 3;
	std::streamoff expected  = f_xy_stride * (((z_dim&7)>0) + (z_dim>>3)) + header_end;
	if(mrc_mode>=0){ // reading mrc file
		f_x_stride       = x_dim;
		f_xy_stride      = xy_stride;
		switch(mrc_mode){
			case 0:  break;
			case 1:  f_x_stride  <<= 1;
			         f_xy_stride <<= 1;
			         break;
			case 2:  f_x_stride  <<= 2;
			         f_xy_stride <<= 2;
			         break;
			default: cout << "ERROR: Only mode 0, 1, 2 are supported for CCP4/MRC map files.\n";
			         exit(33);
		}
		expected         = f_xy_stride * z_dim + header_end;
	}
	if(filesize < expected){
		cout << "\nERROR: Map file size is too small based on header information (expected: " << expected << " Bytes).\n";
		exit(5);
	}
	next_entry("grid");
	fp_num inv_x_step        = read_entry(header, int);
	fp_num inv_y_step        = read_entry(header, int);
	fp_num inv_z_step        = read_entry(header, int);
	fp_num x_step            = 1.0 / inv_x_step;
	fp_num y_step            = 1.0 / inv_y_step;
	fp_num z_step            = 1.0 / inv_z_step;
	next_entry("cell");
	fp_num a_unit            = read_entry(header, float);
	fp_num b_unit            = read_entry(header, float);
	fp_num c_unit            = read_entry(header, float);
	fp_num alpha             = read_entry(header, float);
	fp_num beta              = read_entry(header, float);
	fp_num gamma             = read_entry(header, float);
	fp_num rho_scale         = 1;
	fp_num offset            = 0;
	if(map_type<4){
		next_entry("prod");
		rho_scale        = norm / read_entry(header, float);
		next_entry("plus");
		offset           = read_entry(header, float);
	}
	fp_num rho_min, rho_max;
	double rho_avg, rho_std;
	bool calc_rho_stat       = true;
	fp_num x_origin          = 0;
	fp_num y_origin          = 0;
	fp_num z_origin          = 0;
	if(map_type == mrc){ // http://situs.biomachina.org/fmap.pdf
		x_origin         = read_entry(header, float);
		y_origin         = read_entry(header, float);
		z_origin         = read_entry(header, float);
		if(!((x_origin==0) && (y_origin==0) && (z_origin==0))){ // we tried reading origin fields first and if they are empty, fall back on n*start
			x_start  = 0;
			y_start  = 0;
			z_start  = 0;
		}
		rho_min          = read_entry(header, float);
		rho_max          = read_entry(header, float);
		rho_avg          = read_entry(header, float);
		rho_std          = read_entry(header, float);
		// follows note 5 in documentation:
		// https://www.ccpem.ac.uk/mrc_format/mrc2014.php#note5
		calc_rho_stat    = !((rho_min <= rho_max) && (rho_avg >= std::min(rho_min, rho_max)) && (rho_std > 0));
	}
	
	if(map_type<3){
		fp_num scale     = 1.0 / read_entry(header, int);
		a_unit          *= scale;
		b_unit          *= scale;
		c_unit          *= scale;
		alpha           *= scale;
		beta            *= scale;
		gamma           *= scale;
	}
	cout << "\t-> x_dim = " << x_dim << ", y_dim = " << y_dim << ", z_dim = " << z_dim << "\n";
	cout << "\t-> a_unit = " << a_unit << ", b_unit = " << b_unit << ", c_unit = " << c_unit << "\n";
	cout << "\t-> alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << "\n";
	
	alpha                   *= PI / 180.0;
	beta                    *= PI / 180.0;
	gamma                   *= PI / 180.0;
	
	fp_num inv_a_unit        = 1/a_unit;
	fp_num inv_b_unit        = 1/b_unit;
	fp_num inv_c_unit        = 1/c_unit;
	// r_x = f_x*a + f_y*b*cos(gamma) + f_z*c*cos(beta)
	// r_y =         f_y*b*sin(gamma) + f_z*c*n
	// r_z =                            f_z*c*sqrt(sin^2(beta)-n^2)
	fp_num cos_gamma         = cos(gamma);
	fp_num cos_beta          = cos(beta);
	fp_num sin_gamma         = sin(gamma);
	fp_num inv_sin_gamma     = 1/sin_gamma;
	fp_num cos_inv_sin_gamma = cos_gamma*inv_sin_gamma;
	fp_num n                 = (cos(alpha) - cos_gamma*cos_beta)/sin_gamma;
	fp_num sqrt_factor       = sqrt(1.0 - cos_beta*cos_beta - n*n);
	fp_num inv_sqrt_factor   = 1/sqrt_factor;
	// f_x = [r_x - r_y*cos(gamma)/sin(gamma) + r_z * (n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2)]/a
	// f_y = [      r_y*1/sin(gamma)          - r_z*1/sin(gamma)*n/sqrt(sin^2(beta)-n^2)]/b
	// f_z = [                                  r_z*1/sqrt(sin^2(beta)-n^2)]/c
	fp_num long_inv_term     = (n*cos_inv_sin_gamma - cos_beta)*inv_sqrt_factor;
	fp_num inv_sin_gamma_b   = inv_sin_gamma*inv_b_unit;
	fp_num n_inv_sin_sqrt_b  = n*inv_sin_gamma*inv_sqrt_factor*inv_b_unit;
	fp_num inv_sqrt_factor_c = inv_sqrt_factor*inv_c_unit;
	// (r_x)   [ a  b*cos(gamma)              c*cos(beta) ]   (f_x)
	// (r_y) = [ 0  b*sin(gamma)                     c*n  ] * (f_y)
	// (r_z)   [ 0       0        c*sqrt(sin^2(beta)-n^2) ]   (f_z)
	cout << "\t-> Fractional to cartesian conversion matrix:\n";
	cout.precision(4);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "\t\t" << std::setw(9) << a_unit << " " << std::setw(9) << b_unit*cos_gamma << " " << std::setw(9) << c_unit*cos_beta << "\n";
	cout << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << b_unit*sin_gamma << " " << std::setw(9) << c_unit*n << "\n";
	cout << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << 0 << " " << std::setw(9) << c_unit*sqrt_factor << "\n";
	// inverse to go from cartesian to fractional:
	// (f_x)   [ 1/a -1/a*cos(gamma)/sin(gamma)  1/a*(n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2) ]   (r_x)
	// (f_y) = [ 0             1/b*1/sin(gamma)                      -1/b*1/sin(gamma)*n*1/sqrt(sin^2(beta)-n^2) ] * (r_y)
	// (f_z)   [ 0               0                                                   1/c*1/sqrt(sin^2(beta)-n^2) ]   (r_z)
	cout << "\t-> Cartesian to fractional conversion matrix:\n";
	cout.precision(4);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "\t\t" << std::setw(9) << inv_a_unit << " " << std::setw(9) << -cos_inv_sin_gamma*inv_a_unit << " " << std::setw(9) << long_inv_term*inv_a_unit << "\n";
	cout << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << inv_sin_gamma_b << " " << std::setw(9) << n_inv_sin_sqrt_b << "\n";
	cout << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << 0 << " " << std::setw(9) << inv_sqrt_factor_c << "\n";
	// Make sure we have enough data
	fp_num density_x_start = x_start * x_step;
	fp_num density_y_start = y_start * y_step;
	fp_num density_z_start = z_start * z_step;
	fp_num density_x_end   = (x_start + x_dim - 1) * x_step;
	fp_num density_y_end   = (y_start + y_dim - 1) * y_step;
	fp_num density_z_end   = (z_start + z_dim - 1) * z_step;
	f2c(density_x_start, density_y_start, density_z_start);
	f2c(density_x_end, density_y_end, density_z_end);
	density_x_start       += x_origin;
	density_x_end         += x_origin;
	density_y_start       += y_origin;
	density_y_end         += y_origin;
	density_z_start       += z_origin;
	density_z_end         += z_origin;
	cout.precision(4);
	cout << "\t-> density coordinate range: (" << density_x_start << ", " << density_y_start << ", " << density_z_start << ") A to (" << density_x_end << ", " << density_y_end << ", " << density_z_end << ") A\n";
	fp_num map_x_start  = map_x_center - map_x_dim * grid_spacing * 0.5;
	fp_num map_y_start  = map_y_center - map_y_dim * grid_spacing * 0.5;
	fp_num map_z_start  = map_z_center - map_z_dim * grid_spacing * 0.5;
	fp_num map_x_end    = map_x_start  + map_x_dim * grid_spacing;
	fp_num map_y_end    = map_y_start  + map_y_dim * grid_spacing;
	fp_num map_z_end    = map_z_start  + map_z_dim * grid_spacing;
	// test if there is data for our grid box
	if(((map_x_start - density_x_start) < -MAPEPS) || ((map_y_start - density_y_start) < -MAPEPS) || ((map_z_start - density_z_start) < -MAPEPS) ||
	   ((map_x_end   - density_x_end)   >  MAPEPS) || ((map_y_end   - density_y_end)   >  MAPEPS) || ((map_z_end   - density_z_end)   >  MAPEPS))
	{
		cout << "\nERROR: The specified grid box is (partially) outside of the file's density data.\n";
		cout << "       Grid coordinate range: (" << map_x_start << ", " << map_y_start << ", " << map_z_start << ") A to (" << map_x_end << ", " << map_y_end << ", " << map_z_end << ") A\n";
		exit(6);
	}
	densities.resize(xy_stride*z_dim);
	
	fp_num rho;
	if(calc_rho_stat){
		rho_min = 1e80;
		rho_max = 0;
		rho_avg = 0;
		rho_std = 0;
	}
	
	map_file.seekg(header_end);
	char* data_block     = new char[f_xy_stride];
	if(map_type < 4){ // DSN6 and BRIX
		unsigned z_block, y_block, x_block, data_offset;
		for(unsigned int z = 0; z < z_dim; z += 8){ // z slow
			z_block = std::min(z_dim - z, 8u);
			map_file.read(data_block, f_xy_stride);
			data_offset = 0;
			for(unsigned int y = 0; y < y_dim; y += 8){ // y medium
				y_block = std::min(y_dim - y, 8u);
				for(unsigned int x = 0; x < x_dim; x += 8){ // x fast
					x_block = std::min(x_dim - x, 8u);
					for(unsigned int xb = 0; xb < x_block; xb++){
						for(unsigned int yb = 0; yb < y_block; yb++){
							for(unsigned int zb = 0; zb < z_block; zb++){
								rho = (((unsigned char)data_block[data_offset+((xb+(yb<<3)+(zb<<6))^endian_swap)])-offset)*rho_scale;
								densities[xyz_idx(x+xb, y+yb, z+zb)] = rho;
								rho_min  = std::min(rho, rho_min);
								rho_max  = std::max(rho, rho_max);
								rho_avg += rho;
								rho_std += rho*rho;
							}
						}
					}
					data_offset += 512;
				}
			}
		}
	} else{
		for(unsigned int z = 0; z < z_dim; z++){ // z slow
			map_file.read(data_block, f_xy_stride);
			data_count = 0;
			for(unsigned int y = 0; y < y_dim; y++){ // y medium
				for(unsigned int x = 0; x < x_dim; x++){ // x fast
					read_mrc(rho, data_block+data_count);
					densities[xyz_idx(x, y, z)] = rho;
					rho_min = std::min(rho, rho_min);
					rho_max = std::max(rho, rho_max);
					if(calc_rho_stat){
						rho_avg += rho;
						rho_std += rho*rho;
					}
				}
			}
		}
	}
	delete[] data_block;
	if(calc_rho_stat){
		rho_avg /= x_dim * y_dim * z_dim;
		rho_std /= x_dim * y_dim * z_dim;
		rho_std -= rho_avg * rho_avg;
		rho_std  = sqrt(rho_std);
	}
	if(map_type == mrc){
		rho_min -= rho_avg;
		rho_min /= rho_std;
		rho_max -= rho_avg;
		rho_max /= rho_std;
	}
	cout.precision(3);
	cout << "\t-> density range: " << rho_min << " to " << rho_max << std::setprecision(6) << " (average: " << rho_avg << " +/- " << rho_std << ")\n";
	map_file.close();
	double file_reading_ms = seconds_since(runtime)*1000.0;
	cout.precision(3);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "<- Finished reading densities, took " << file_reading_ms << " ms.\n\n";
	
	cout << "Interpolating density data for " << map_x_dim << "x" << map_y_dim << "x" << map_z_dim << " grid (spacing: " << grid_spacing << " A)\n";
	cout << "\t-> grid start:  (" << map_x_start << ", " << map_y_start << ", " << map_z_start << ") A\n";
	cout << "\t-> grid size:   (" << map_x_dim * grid_spacing << ", " << map_y_dim * grid_spacing << ", " << map_z_dim * grid_spacing << ") A\n";
	
	fp_num grid_a, grid_b, grid_c, density;
	std::vector<fp_num> grid_map;
	unsigned int g1 = map_x_dim + 1;
	unsigned int g2 = g1 * (map_y_dim + 1);
	grid_map.resize(g2 * (map_z_dim + 1) + 9);
	grid_map[0]     = map_x_dim;
	grid_map[1]     = map_y_dim;
	grid_map[2]     = map_z_dim;
	grid_map[3]     = map_x_center;
	grid_map[4]     = map_y_center;
	grid_map[5]     = map_z_center;
	grid_map[6]     = grid_spacing;
	
	fp_num* density_data = densities.data();
	fp_num* data_point;
#ifdef MGLTOOLS_MATH_COMPARISON
	fp_num ga, gb, gc;
#endif
	rho_min = 0;
	rho_max = 0;
	for(unsigned int z = 0; z <= map_z_dim; z++){
		for(unsigned int y = 0; y <= map_y_dim; y++){
			for(unsigned int x = 0; x <= map_x_dim; x++){
				grid_a = map_x_start + x * grid_spacing - x_origin;
				grid_b = map_y_start + y * grid_spacing - y_origin;
				grid_c = map_z_start + z * grid_spacing - z_origin;
				c2f(grid_a, grid_b, grid_c);
				grid_a *= inv_x_step;
				grid_b *= inv_y_step;
				grid_c *= inv_z_step;
				grid_a -= x_start;
				grid_b -= y_start;
				grid_c -= z_start;
#ifdef MGLTOOLS_MATH_COMPARISON
				ga     = (grid_a - density_x_start) / (x_step * a_unit);
				gb     = (grid_b - density_y_start) / (y_step * b_unit);
				gc     = (grid_c - density_z_start) / (z_step * c_unit);
				cout << grid_a << ":" << ga << " ; " << grid_b << ":" << gb << " ; " << grid_c << ":" << gc << "\n";
#endif
				// Getting coordinates
				unsigned int x_low  = grid_a; // conversion to integer always truncates float
				unsigned int y_low  = grid_b; // <- this looks dangerous but we made sure that grid_a,b,c can only be >= 0
				unsigned int z_low  = grid_c;
				
				data_point = density_data + (unsigned long)xyz_idx(x_low,y_low,z_low);
				
				fp_num dx   = grid_a - x_low;
				fp_num omdx = 1.0 - dx;
				fp_num dy   = grid_b - y_low;
				fp_num omdy = 1.0 - dy;
				fp_num dz   = grid_c - z_low;
				fp_num omdz = 1.0 - dz;
				
				/*          Z                          */
				/*          '                          */
				/*          2 - - - - 3                */
				/*         /.        /|                */
				/*        6 - - - - 7 |                */
				/*        | '       | |                */
				/*        | 0 - - - + 1 -- Y           */
				/*        '/        |/                 */
				/*        4 - - - - 5                  */
				/*       /                             */
				/*      X                              */
				density = omdz * (omdy * (data_point[0]*omdx + data_point[1]*dx) +
				                    dy * (data_point[x_dim]*omdx + data_point[x_dim + 1]*dx)) +
				          dz   * (omdy * (data_point[xy_stride]*omdx + data_point[xy_stride + 1]*dx) +
				                    dy * (data_point[x_dim + xy_stride]*omdx + data_point[x_dim + xy_stride + 1]*dx));
				if(map_type == mrc){
					density -= rho_avg;
					density /= rho_std;
				}
				grid_map[x + y*g1 + z*g2 + 9] = density;
				rho_min = std::min(density, rho_min);
				rho_max = std::max(density, rho_max);
			}
		}
	}
	grid_map[7] = rho_min;
	grid_map[8] = rho_max;
	cout << "<- Finished interpolating grid map, took " << seconds_since(runtime)*1000.0 - file_reading_ms << " ms.\n\n";
	
	return grid_map;
}

#endif // INCLUDED_MAP_READER

