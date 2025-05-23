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

#define MAPEPS 1e-4
#define MEDIAN_BINS 10000

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

#include "ScalarMat.h"

using namespace std;

//                              start     extent     grid      cell axes    angles       origin     min, max, avg, std
const int mrc_offsets[22] = {16, 20, 24, 0, 4, 8, 28, 32, 36, 40, 44, 48, 52, 56, 60, 196, 200, 204, 76, 80, 84, 216};

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

inline void apply_periodicity(
                              fp_num len,
                              fp_num &r
                             )
{
	if(r >= len){
		r -= len;
		while(r > len) r -= len;
		return;
	}
	while(r < 0) r += len;
}

#define GAUSS_ZIGGURAT_NORM (1.0/3.442619855899)

inline std::vector<fp_num> add_normal_noise(
                                            std::vector<fp_num> density,
                                            fp_num  sigma  = 1
                                           )
{
	unsigned int off     = density[0];
	std::vector<fp_num> result(density.size());
	memcpy(result.data(), density.data(), off * sizeof(fp_num));
	fp_num rho_min       =  1e80;
	fp_num rho_max       = -1e80;
	fp_num rho_avg       = 0;
	fp_num rho_std       = 0;
	__uint32_t seed = time(NULL);
	init_CMWC4096(seed);
	zigset();
	double scale = GAUSS_ZIGGURAT_NORM * sigma;
	for(unsigned int i = off;
	                 i < density.size();
	                 i++)
	{
		fp_num val = density[i] + ran_n()*scale;
		result[i]  = val;
		rho_min    = std::min(val, rho_min);
		rho_max    = std::max(val, rho_max);
		rho_avg   += val;
		rho_std   += val*val;
	}
	rho_avg   /= density.size() - off;
	rho_std   /= density.size() - off;
	rho_std   -= rho_avg * rho_avg;
	rho_std    = sqrt(rho_std);
	result[8]  = rho_min;
	result[9]  = rho_max;
	result[11] = rho_std;
	return result;
}

inline std::vector<fp_num> gaussian_convolution(
                                                std::vector<fp_num> density,
                                                fp_num  sigma  = 2,
                                                fp_num  cutoff = 8
                                               )
{
	unsigned int off     = density[0];
	unsigned int map_x_p = density[1] + 1;
	unsigned int map_y_p = density[2] + 1;
	unsigned int map_z_p = density[3] + 1;
	unsigned int g1      = map_x_p;
	unsigned int g2      = g1 * map_y_p;
	std::vector<fp_num> result(density.size());
	memcpy(result.data(), density.data(), off * sizeof(fp_num));
	fp_num cut           = cutoff / density[7];
	fp_num g_factor      = density[7] * density[7] / (sigma * sigma);
	fp_num cut2          = cut*cut;
	#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int i = 0;
	                 i < g2*map_z_p;
	                 i++)
	{
		// idx = x + g1*y + g2*z (g2 = g1 * (map_y + 1))
		unsigned int gz = i / g2; // idx / g2 = z
		// idx / g1 = y + g2/g1*z = y + (map_y+1)*z
		unsigned int gy = i / g1 - map_y_p * gz;
		unsigned int gx = i - g1*gy - g2*gz; // idx = x + g1*y + g2*z
		int x_start      = std::max(0,(int)floor(gx-cut));
		int x_end        = std::min((int)map_x_p,(int)ceil(gx+cut));
		int y_start      = std::max(0,(int)floor(gy-cut));
		int y_end        = std::min((int)map_y_p,(int)ceil(gy+cut));
		int z_end        = std::min((int)map_z_p,(int)ceil(gz+cut));
		fp_num integral  = 0;
		fp_num wsum      = 0;
		for(int z = std::max(0,(int)floor(gz-cut)); z < z_end; z++){
			fp_num dist_z2 = (gz-z)*(gz-z);
			if(dist_z2 < cut2){ // no need to go through the motions if we're already too far ...
				for(int y = y_start; y < y_end; y++){
					fp_num dist_zy2 = dist_z2 + (gy-y)*(gy-y);
					if(dist_zy2 < cut2){
						for(int x = x_start; x < x_end; x++){
							// calculate (square) distance from current grid point to current location
							fp_num dist2  = (gx-x)*(gx-x) + dist_zy2;
							if(dist2 < cut2){
								fp_num weight = gaussfit(g_factor * dist2);
								integral     += weight * density[off + x  + y*g1  + z*g2];
								wsum         += weight;
							}
						}
					}
				}
			}
		}
		result[i + off] = (wsum > 0) ? integral / wsum : 0;
	}
	fp_num rho_min       =  1e80;
	fp_num rho_max       = -1e80;
	fp_num rho_avg       = 0;
	fp_num rho_std       = 0;
	for(unsigned int i = off;
	                 i < result.size();
	                 i++)
	{
		fp_num val = result[i];
		rho_min    = std::min(val, rho_min);
		rho_max    = std::max(val, rho_max);
		rho_avg   += val;
		rho_std   += val*val;
	}
	rho_avg   /= g2*map_z_p;
	rho_std   /= g2*map_z_p;
	rho_std   -= rho_avg * rho_avg;
	rho_std    = sqrt(rho_std);
	result[8]  = rho_min;
	result[9]  = rho_max;
	result[11] = rho_std;
	return result;
}

inline std::vector<fp_num> read_map_to_grid(
                                            std::string  filename,
                                                     int map_type,
                                            unsigned int map_x_dim,
                                            unsigned int map_y_dim,
                                            unsigned int map_z_dim,
                                            fp_num       map_x_center,
                                            fp_num       map_y_center,
                                            fp_num       map_z_center,
                                            fp_num       grid_spacing,
                                            fp_num       gaussian_filter_sigma = 0,
                                            fp_num       noise_std_range       = 0,
                                            bool         repeat_unit_cell      = true,
                                            fp_num*      grid_align            = NULL
                                           )
{
	timeval runtime;
	start_timer(runtime);
	stringstream output;
	
	std::vector<fp_num> densities;
	
	output << "Reading map file [" << filename << "]\n";
	std::ifstream map_file(filename, std::ifstream::binary);
	if(map_file.fail()){
		#pragma omp critical
		cout << output.str() << "\nERROR: Can't open map file " << filename << ".\n";
		exit(1);
	}
	std::streamoff filesize  = map_file.tellg();
	                           map_file.seekg(0, std::ios::end);
	filesize                 = map_file.tellg() - filesize;
	                           map_file.seekg(0, std::ios::beg);
	char map_header[256];
	if(!map_file.read(map_header, 256)){ // the first 256 Bytes are all that's needed really
		#pragma omp critical
		cout << output.str() << "\nERROR: Can't read header.\n";
		exit(2);
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
		case dsn6:      output << "\t-> DSN6";
		                norm = 100;
		                break;
		case brix:      output << "\t-> BRIX";
		                break;
		case mrc:       output << "\t-> MRC";
		                header_end  = 1024;
		                if(*((char*)&norm) != (*header || *(header+1))) endian_swap = true;
		                header_end += *(reinterpret_cast<int*>(read_32bit(header+92))); // add bytes of the extended header
		                mrc_mode    = *(reinterpret_cast<int*>(read_32bit(header+12)));
		                break;
		default:
		case automatic:
		                #pragma omp critical
		                cout << output.str() << "\nERROR: Unknown map file type.\n";
		                exit(42);
	}
	output << (endian_swap ?" endian-swapped":"") << ", file size: " << filesize << "\n";
	
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
			default: 
			         #pragma omp critical
			         cout << output.str() << "ERROR: Only mode 0, 1, 2 are supported for CCP4/MRC map files.\n";
			         exit(33);
		}
		expected         = f_xy_stride * z_dim + header_end;
	}
	if(filesize < expected){
		#pragma omp critical
		cout << output.str() << "\nERROR: Map file size is too small based on header information (expected: " << expected << " Bytes).\n";
		exit(3);
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
	output << "\t-> x_dim = " << x_dim << ", y_dim = " << y_dim << ", z_dim = " << z_dim << "\n";
	output << "\t-> a_unit = " << a_unit << ", b_unit = " << b_unit << ", c_unit = " << c_unit << "\n";
	output << "\t-> alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << "\n";
	
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
	output << "\t-> Fractional to cartesian conversion matrix:\n";
	output.precision(4);
	output.setf(ios::fixed, ios::floatfield);
	output << "\t\t" << std::setw(9) << a_unit << " " << std::setw(9) << b_unit*cos_gamma << " " << std::setw(9) << c_unit*cos_beta << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << b_unit*sin_gamma << " " << std::setw(9) << c_unit*n << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << 0 << " " << std::setw(9) << c_unit*sqrt_factor << "\n";
	// inverse to go from cartesian to fractional:
	// (f_x)   [ 1/a -1/a*cos(gamma)/sin(gamma)  1/a*(n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2) ]   (r_x)
	// (f_y) = [ 0             1/b*1/sin(gamma)                      -1/b*1/sin(gamma)*n*1/sqrt(sin^2(beta)-n^2) ] * (r_y)
	// (f_z)   [ 0               0                                                   1/c*1/sqrt(sin^2(beta)-n^2) ]   (r_z)
	output << "\t-> Cartesian to fractional conversion matrix:\n";
	output.precision(4);
	output.setf(ios::fixed, ios::floatfield);
	output << "\t\t" << std::setw(9) << inv_a_unit << " " << std::setw(9) << -cos_inv_sin_gamma*inv_a_unit << " " << std::setw(9) << long_inv_term*inv_a_unit << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << inv_sin_gamma_b << " " << std::setw(9) << n_inv_sin_sqrt_b << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << 0 << " " << std::setw(9) << inv_sqrt_factor_c << "\n";
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
	output.precision(4);
	output << "\t-> density unit cell range: (" << density_x_start << ", " << density_y_start << ", " << density_z_start << ") A to (" << density_x_end << ", " << density_y_end << ", " << density_z_end << ") A\n";
	fp_num map_x_start  = map_x_center - map_x_dim * grid_spacing * 0.5;
	fp_num map_y_start  = map_y_center - map_y_dim * grid_spacing * 0.5;
	fp_num map_z_start  = map_z_center - map_z_dim * grid_spacing * 0.5;
	fp_num map_x_end    = map_x_start  + map_x_dim * grid_spacing;
	fp_num map_y_end    = map_y_start  + map_y_dim * grid_spacing;
	fp_num map_z_end    = map_z_start  + map_z_dim * grid_spacing;
	fp_num bounding_x_start = map_x_start;
	fp_num bounding_y_start = map_y_start;
	fp_num bounding_z_start = map_z_start;
	fp_num bounding_x_end = map_x_end;
	fp_num bounding_y_end = map_y_end;
	fp_num bounding_z_end = map_z_end;
	fp_num gx, gy, gz;
	fp_num ga, gb, gc;
	if(grid_align != NULL){
		for(unsigned int i=0; i<8; i++){ // go over all corners of the grid box
			// subtract grid center
			gx = map_x_center + map_x_dim * grid_spacing * ((i&1) ? 0.5 : -0.5) - grid_align[12];
			gy = map_y_center + map_y_dim * grid_spacing * ((i&2) ? 0.5 : -0.5) - grid_align[13];
			gz = map_z_center + map_x_dim * grid_spacing * ((i&4) ? 0.5 : -0.5) - grid_align[14];
			// rotate
			ga = gx * grid_align[0] + gy * grid_align[1] + gz * grid_align[2];
			gb = gx * grid_align[3] + gy * grid_align[4] + gz * grid_align[5];
			gc = gx * grid_align[6] + gy * grid_align[7] + gz * grid_align[8];
			// move to map center
			ga += grid_align[9];
			gb += grid_align[10];
			gc += grid_align[11];
			// move bounding box if needed
			bounding_x_start = (bounding_x_start > ga) ? ga : bounding_x_start;
			bounding_y_start = (bounding_y_start > gb) ? gb : bounding_y_start;
			bounding_z_start = (bounding_z_start > gc) ? gc : bounding_z_start;
			bounding_x_end = (bounding_x_end < ga) ? ga : bounding_x_end;
			bounding_y_end = (bounding_y_end < gb) ? gb : bounding_y_end;
			bounding_z_end = (bounding_z_end < gc) ? gc : bounding_z_end;
		}
	}
	// test if there is data for our grid box
	if(!repeat_unit_cell &&
	   (((bounding_x_start - density_x_start) < -MAPEPS) || ((bounding_y_start - density_y_start) < -MAPEPS) || ((bounding_z_start - density_z_start) < -MAPEPS) ||
	   ((bounding_x_end   - density_x_end)   >  MAPEPS) || ((bounding_y_end   - density_y_end)   >  MAPEPS) || ((bounding_z_end   - density_z_end)   >  MAPEPS)))
	{
		#pragma omp critical
		{
			cout << output.str() << "\nERROR: The specified grid box is (partially) outside of the file's density data.\n";
			if(grid_align == NULL)
				cout << "       Grid coordinate range: (" << map_x_start << ", " << map_y_start << ", " << map_z_start << ") A to (" << map_x_end << ", " << map_y_end << ", " << map_z_end << ") A\n";
			else
				cout << "       Grid bounding box: (" << bounding_x_start << ", " << bounding_y_start << ", " << bounding_z_start << ") A to (" << bounding_x_end << ", " << bounding_y_end << ", " << bounding_z_end << ") A\n";
		}
		exit(4);
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
	} else{ // MRC
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
	// calculate median
	std::vector<unsigned int> density_hist(MEDIAN_BINS+1, 0);
	fp_num inv_binwidth = MEDIAN_BINS / (rho_max - rho_min);
	data_count = densities.size();
	for(unsigned int i=0; i<data_count; i++)
		density_hist[(unsigned int)floor((densities[i]-rho_min) * inv_binwidth)]++;
	unsigned int half_count = data_count >> 1; // find median == find bin number with just more than half the points
	unsigned int median_idx = 0;
	data_count = 0;
	while(data_count < half_count)
		data_count += density_hist[median_idx++];
	if(map_type == mrc){
		rho_min -= rho_avg;
		rho_min /= rho_std;
		rho_max -= rho_avg;
		rho_max /= rho_std;
	}
	output.precision(3);
	output << "\t-> density range: " << rho_min << " to " << rho_max << std::setprecision(6) << " (average: " << rho_avg << " +/- " << rho_std << "; median: " << median_idx / inv_binwidth + rho_min << ")\n";
	map_file.close();
	double file_reading_ms = seconds_since(runtime)*1000.0;
	output.precision(3);
	output.setf(ios::fixed, ios::floatfield);
	output << "<- Finished reading densities, took " << file_reading_ms << " ms.\n\n";
	
	output << "Interpolating density data for " << map_x_dim << "x" << map_y_dim << "x" << map_z_dim << " grid (spacing: " << grid_spacing << " A)\n";
	output << "\t-> grid start:  (" << map_x_start << ", " << map_y_start << ", " << map_z_start << ") A\n";
	output << "\t-> grid size:   (" << map_x_dim * grid_spacing << ", " << map_y_dim * grid_spacing << ", " << map_z_dim * grid_spacing << ") A\n";
	if(grid_align != NULL)
		output << "\t-> grid bounding box: (" << bounding_x_start << ", " << bounding_y_start << ", " << bounding_z_start << ") A to (" << bounding_x_end << ", " << bounding_y_end << ", " << bounding_z_end << ") A\n";
	
	fp_num grid_a, grid_b, grid_c, density;
	unsigned int g1 = map_x_dim + 1;
	unsigned int g2 = g1 * (map_y_dim + 1);
	std::vector<fp_num> grid_map(g2 * (map_z_dim + 1) + 12);
	grid_map[0]     = 12;
	grid_map[1]     = map_x_dim;
	grid_map[2]     = map_y_dim;
	grid_map[3]     = map_z_dim;
	grid_map[4]     = map_x_center;
	grid_map[5]     = map_y_center;
	grid_map[6]     = map_z_center;
	grid_map[7]     = grid_spacing;
	
	fp_num* density_data = densities.data();
	fp_num* data_point;
	rho_min = 1e80;
	rho_max = 0;
	for(unsigned int z = 0; z <= map_z_dim; z++){
		for(unsigned int y = 0; y <= map_y_dim; y++){
			for(unsigned int x = 0; x <= map_x_dim; x++){
				grid_a = map_x_start + x * grid_spacing - x_origin;
				grid_b = map_y_start + y * grid_spacing - y_origin;
				grid_c = map_z_start + z * grid_spacing - z_origin;
				if(grid_align != NULL){ // align grid coordinates to density map
					// subtract grid center
					gx = grid_a - grid_align[12];
					gy = grid_b - grid_align[13];
					gz = grid_c - grid_align[14];
					// rotate
					grid_a = gx * grid_align[0] + gy * grid_align[1] + gz * grid_align[2];
					grid_b = gx * grid_align[3] + gy * grid_align[4] + gz * grid_align[5];
					grid_c = gx * grid_align[6] + gy * grid_align[7] + gz * grid_align[8];
					// move to map center
					grid_a += grid_align[9];
					grid_b += grid_align[10];
					grid_c += grid_align[11];
				}
				c2f(grid_a, grid_b, grid_c);
				grid_a *= inv_x_step;
				grid_b *= inv_y_step;
				grid_c *= inv_z_step;
				grid_a -= x_start;
				grid_b -= y_start;
				grid_c -= z_start;
				if(repeat_unit_cell){
					apply_periodicity(x_dim, grid_a);
					apply_periodicity(y_dim, grid_b);
					apply_periodicity(z_dim, grid_c);
				}
#ifdef MGLTOOLS_MATH_COMPARISON
				ga     = (grid_a - density_x_start) / (x_step * a_unit);
				gb     = (grid_b - density_y_start) / (y_step * b_unit);
				gc     = (grid_c - density_z_start) / (z_step * c_unit);
				output << grid_a << ":" << ga << " ; " << grid_b << ":" << gb << " ; " << grid_c << ":" << gc << "\n";
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
				int x_high = x_dim - 1;
				    x_high = ((int)x_low < x_high) ? 1 : -x_high;
				int y_high = y_dim - 1;
				    y_high = ((int)y_low < y_high) ? x_dim : -y_high*x_dim;
				int z_high = z_dim - 1;
				    z_high = ((int)z_low < z_high) ? xy_stride : -z_high*xy_stride;
				
				density = omdz * (omdy * (data_point[0]               * omdx + data_point[x_high]                   * dx) +
				                    dy * (data_point[y_high]          * omdx + data_point[x_high + y_high]          * dx)) +
				          dz   * (omdy * (data_point[z_high]          * omdx + data_point[x_high          + z_high] * dx) +
				                    dy * (data_point[y_high + z_high] * omdx + data_point[x_high + y_high + z_high] * dx));
				if(map_type == mrc){
					density -= rho_avg;
					density /= rho_std;
				}
				grid_map[x + y*g1 + z*g2 + 12] = density;
				rho_min = std::min(density, rho_min);
				rho_max = std::max(density, rho_max);
			}
		}
	}
	double current_ms = seconds_since(runtime)*1000.0;
	output << "<- Finished interpolating densities to grid map points, took " << current_ms - file_reading_ms << " ms.\n\n";
	grid_map[8]  = rho_min;
	grid_map[9]  = rho_max;
	grid_map[11] = rho_std;
	if(noise_std_range < -EPS){
		output << "Applying Normal-distributed noise with width of " << fabs(noise_std_range*rho_std) << " A\n";
		grid_map = add_normal_noise(
		                            grid_map,
		                            fabs(noise_std_range*rho_std)
		                           );
		output << "<- Done, took " << seconds_since(runtime)*1000.0 - current_ms << " ms.\n\n";
		current_ms = seconds_since(runtime)*1000.0;
	}
	if(gaussian_filter_sigma > EPS){
		output << "Applying Gaussian filter with width of " << gaussian_filter_sigma << " A\n";
		grid_map = gaussian_convolution(
		                                grid_map,
		                                gaussian_filter_sigma
		                               );
		output << "<- Done, took " << seconds_since(runtime)*1000.0 - current_ms << " ms.\n\n";
		current_ms = seconds_since(runtime)*1000.0;
	}
	if(noise_std_range > EPS){
		output << "Applying Normal-distributed noise with width of " << noise_std_range*rho_std << " A\n";
		grid_map = add_normal_noise(
		                            grid_map,
		                            noise_std_range*rho_std
		                           );
		output << "<- Done, took " << seconds_since(runtime)*1000.0 - current_ms << " ms.\n\n";
		current_ms = seconds_since(runtime)*1000.0;
	}
	// calculate median
	output << "Calculating median\n";
	rho_min      = grid_map[8]; // in case the gaussian filter changed them
	rho_max      = grid_map[9];
	memset(density_hist.data(), 0, MEDIAN_BINS*sizeof(unsigned int));
	inv_binwidth = MEDIAN_BINS / (rho_max - rho_min);
	for(unsigned int i=grid_map[0]; i<grid_map.size(); i++)
		density_hist[(unsigned int)floor((grid_map[i]-rho_min) * inv_binwidth)]++;
	half_count   = (grid_map.size() - (int)grid_map[0]) >> 1; // find median == find bin number with just more than half the points
	data_count   = 0;
	for(median_idx = 0; (data_count < half_count) && (median_idx < MEDIAN_BINS); median_idx++)
		data_count += density_hist[median_idx];
	grid_map[10] = (fp_num)median_idx / MEDIAN_BINS;
	output << "\t-> range: " << grid_map[8] << " to " << grid_map[9] << " (median: " << grid_map[10] * (grid_map[9] - grid_map[8]) + grid_map[8] << ")\n";
	output << "<- Done, took " << seconds_since(runtime)*1000.0 - current_ms << " ms.\n\n";
	#pragma omp critical
	cout << output.str();
	
	return grid_map;
}

inline void convert_map_to_mrc(std::string  filename)
{
	timeval runtime;
	start_timer(runtime);
	stringstream output;
	
	std::vector<float> densities;
	
	output << "Reading map file [" << filename << "]\n";
	std::ifstream map_file(filename, std::ifstream::binary);
	if(map_file.fail()){
		#pragma omp critical
		cout << output.str() << "\nERROR: Can't open map file " << filename << ".\n";
		exit(1);
	}
	std::streamoff filesize  = map_file.tellg();
	                           map_file.seekg(0, std::ios::end);
	filesize                 = map_file.tellg() - filesize;
	                           map_file.seekg(0, std::ios::beg);
	char map_header[256];
	if(!map_file.read(map_header, 256)){ // the first 256 Bytes are all that's needed really
		#pragma omp critical
		cout << output.str() << "\nERROR: Can't read header.\n";
		exit(2);
	}
	char* header             = map_header;
	int map_type             = determine_map_type(header);
	unsigned int norm        = 1;
	unsigned int header_end  = DSN6_BLOCKSIZE;
	bool endian_swap         = false;
	int mrc_mode             = -1;
	int tempval;
	switch(map_type){
		case dsn6_swap: endian_swap = true;
		case dsn6:      output << "\t-> DSN6";
		                norm = 100;
		                break;
		case brix:      output << "\t-> BRIX";
		                break;
		case mrc:       output << "\t-> MRC";
		                header_end  = 1024;
		                if(*((char*)&norm) != (*header || *(header+1))) endian_swap = true;
		                header_end += *(reinterpret_cast<int*>(read_32bit(header+92))); // add bytes of the extended header
		                mrc_mode    = *(reinterpret_cast<int*>(read_32bit(header+12)));
		                break;
		default:
		case automatic:
		                #pragma omp critical
		                cout << output.str() << "\nERROR: Unknown map file type.\n";
		                exit(42);
	}
	output << (endian_swap ?" endian-swapped":"") << ", file size: " << filesize << "\n";
	
	unsigned int data_count  = 0;
	next_entry("origin");
	float x_start           = read_entry(header, int);
	float y_start           = read_entry(header, int);
	float z_start           = read_entry(header, int);
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
			default: 
			         #pragma omp critical
			         cout << output.str() << "ERROR: Only mode 0, 1, 2 are supported for CCP4/MRC map files.\n";
			         exit(33);
		}
		expected         = f_xy_stride * z_dim + header_end;
	}
	if(filesize < expected){
		#pragma omp critical
		cout << output.str() << "\nERROR: Map file size is too small based on header information (expected: " << expected << " Bytes).\n";
		exit(3);
	}
	next_entry("grid");
	float inv_x_step        = read_entry(header, int);
	float inv_y_step        = read_entry(header, int);
	float inv_z_step        = read_entry(header, int);
	float x_step            = 1.0 / inv_x_step;
	float y_step            = 1.0 / inv_y_step;
	float z_step            = 1.0 / inv_z_step;
	next_entry("cell");
	float a_unit            = read_entry(header, float);
	float b_unit            = read_entry(header, float);
	float c_unit            = read_entry(header, float);
	float alpha             = read_entry(header, float);
	float beta              = read_entry(header, float);
	float gamma             = read_entry(header, float);
	float rho_scale         = 1;
	float offset            = 0;
	if(map_type<4){
		next_entry("prod");
		rho_scale        = norm / read_entry(header, float);
		next_entry("plus");
		offset           = read_entry(header, float);
	}
	float rho_min, rho_max;
	double rho_avg, rho_std;
	bool calc_rho_stat       = true;
	float x_origin          = 0;
	float y_origin          = 0;
	float z_origin          = 0;
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
		float scale     = 1.0 / read_entry(header, int);
		a_unit          *= scale;
		b_unit          *= scale;
		c_unit          *= scale;
		alpha           *= scale;
		beta            *= scale;
		gamma           *= scale;
	}
	output << "\t-> x_dim = " << x_dim << ", y_dim = " << y_dim << ", z_dim = " << z_dim << "\n";
	output << "\t-> a_unit = " << a_unit << ", b_unit = " << b_unit << ", c_unit = " << c_unit << "\n";
	output << "\t-> alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << "\n";
	
	alpha                   *= PI / 180.0;
	beta                    *= PI / 180.0;
	gamma                   *= PI / 180.0;
	
	float inv_a_unit        = 1/a_unit;
	float inv_b_unit        = 1/b_unit;
	float inv_c_unit        = 1/c_unit;
	// r_x = f_x*a + f_y*b*cos(gamma) + f_z*c*cos(beta)
	// r_y =         f_y*b*sin(gamma) + f_z*c*n
	// r_z =                            f_z*c*sqrt(sin^2(beta)-n^2)
	float cos_gamma         = cos(gamma);
	float cos_beta          = cos(beta);
	float sin_gamma         = sin(gamma);
	float inv_sin_gamma     = 1/sin_gamma;
	float cos_inv_sin_gamma = cos_gamma*inv_sin_gamma;
	float n                 = (cos(alpha) - cos_gamma*cos_beta)/sin_gamma;
	float sqrt_factor       = sqrt(1.0 - cos_beta*cos_beta - n*n);
	float inv_sqrt_factor   = 1/sqrt_factor;
	// f_x = [r_x - r_y*cos(gamma)/sin(gamma) + r_z * (n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2)]/a
	// f_y = [      r_y*1/sin(gamma)          - r_z*1/sin(gamma)*n/sqrt(sin^2(beta)-n^2)]/b
	// f_z = [                                  r_z*1/sqrt(sin^2(beta)-n^2)]/c
	float long_inv_term     = (n*cos_inv_sin_gamma - cos_beta)*inv_sqrt_factor;
	float inv_sin_gamma_b   = inv_sin_gamma*inv_b_unit;
	float n_inv_sin_sqrt_b  = n*inv_sin_gamma*inv_sqrt_factor*inv_b_unit;
	float inv_sqrt_factor_c = inv_sqrt_factor*inv_c_unit;
	// (r_x)   [ a  b*cos(gamma)              c*cos(beta) ]   (f_x)
	// (r_y) = [ 0  b*sin(gamma)                     c*n  ] * (f_y)
	// (r_z)   [ 0       0        c*sqrt(sin^2(beta)-n^2) ]   (f_z)
	output << "\t-> Fractional to cartesian conversion matrix:\n";
	output.precision(4);
	output.setf(ios::fixed, ios::floatfield);
	output << "\t\t" << std::setw(9) << a_unit << " " << std::setw(9) << b_unit*cos_gamma << " " << std::setw(9) << c_unit*cos_beta << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << b_unit*sin_gamma << " " << std::setw(9) << c_unit*n << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << 0 << " " << std::setw(9) << c_unit*sqrt_factor << "\n";
	// inverse to go from cartesian to fractional:
	// (f_x)   [ 1/a -1/a*cos(gamma)/sin(gamma)  1/a*(n*cos(gamma)/sin(gamma)-cos(beta))*1/sqrt(sin^2(beta)-n^2) ]   (r_x)
	// (f_y) = [ 0             1/b*1/sin(gamma)                      -1/b*1/sin(gamma)*n*1/sqrt(sin^2(beta)-n^2) ] * (r_y)
	// (f_z)   [ 0               0                                                   1/c*1/sqrt(sin^2(beta)-n^2) ]   (r_z)
	output << "\t-> Cartesian to fractional conversion matrix:\n";
	output.precision(4);
	output.setf(ios::fixed, ios::floatfield);
	output << "\t\t" << std::setw(9) << inv_a_unit << " " << std::setw(9) << -cos_inv_sin_gamma*inv_a_unit << " " << std::setw(9) << long_inv_term*inv_a_unit << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << inv_sin_gamma_b << " " << std::setw(9) << n_inv_sin_sqrt_b << "\n";
	output << "\t\t" << std::setw(9) << 0 << " " << std::setw(9) << 0 << " " << std::setw(9) << inv_sqrt_factor_c << "\n";
	// Make sure we have enough data
	float density_x_start = x_start * x_step;
	float density_y_start = y_start * y_step;
	float density_z_start = z_start * z_step;
	float density_x_end   = (x_start + x_dim - 1) * x_step;
	float density_y_end   = (y_start + y_dim - 1) * y_step;
	float density_z_end   = (z_start + z_dim - 1) * z_step;
	f2c(density_x_start, density_y_start, density_z_start);
	f2c(density_x_end, density_y_end, density_z_end);
	density_x_start       += x_origin;
	density_x_end         += x_origin;
	density_y_start       += y_origin;
	density_y_end         += y_origin;
	density_z_start       += z_origin;
	density_z_end         += z_origin;
	output.precision(4);
	output << "\t-> density unit cell range: (" << density_x_start << ", " << density_y_start << ", " << density_z_start << ") A to (" << density_x_end << ", " << density_y_end << ", " << density_z_end << ") A\n";
	densities.resize(xy_stride*z_dim);
	
	float rho;
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
	} else{ // MRC
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
	// calculate median
	std::vector<unsigned int> density_hist(MEDIAN_BINS+1, 0);
	float inv_binwidth = MEDIAN_BINS / (rho_max - rho_min);
	data_count = densities.size();
	for(unsigned int i=0; i<data_count; i++)
		density_hist[(unsigned int)floor((densities[i]-rho_min) * inv_binwidth)]++;
	unsigned int half_count = data_count >> 1; // find median == find bin number with just more than half the points
	unsigned int median_idx = 0;
	data_count = 0;
	while(data_count < half_count)
		data_count += density_hist[median_idx++];
	if(map_type == mrc){
		rho_min -= rho_avg;
		rho_min /= rho_std;
		rho_max -= rho_avg;
		rho_max /= rho_std;
	}
	output.precision(3);
	output << "\t-> density range: " << rho_min << " to " << rho_max << std::setprecision(6) << " (average: " << rho_avg << " +/- " << rho_std << "; median: " << median_idx / inv_binwidth + rho_min << ")\n";
	map_file.close();
	double file_reading_ms = seconds_since(runtime)*1000.0;
	output.precision(3);
	output.setf(ios::fixed, ios::floatfield);
	output << "<- Finished reading densities, took " << file_reading_ms << " ms.\n\n";
	#pragma omp critical
	cout << output.str();
	
	std::size_t ext = filename.find_last_of(".");
	filename = filename.substr(0, ext) + ".convert.mrc";
#ifdef PARALLELIZE
	#pragma omp critical
#endif
	cout << "Writing MRC grid map file [" << filename << "]\n";
	std::ofstream out_file(filename, std::ifstream::binary);
	if(out_file.fail()){
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
	} out_header;
	memset(out_header.extra, 0,100);
	if(HOST_LITTLE_ENDIAN){
		out_header.extra[12] = 0xAD;
		out_header.extra[13] = 0x4E;
	} else{
		out_header.extra[14] = 0x4E;
		out_header.extra[15] = 0xAD;
	}
	memset(out_header.labels,0,800);
	strncpy(out_header.labels, "CryoXKit MRC grid map", 22);
	out_header.nx       = x_dim;
	out_header.mx       = inv_x_step;
	out_header.ny       = y_dim;
	out_header.my       = inv_y_step;
	out_header.nz       = z_dim;
	out_header.mz       = inv_z_step;
	out_header.x_start  = x_start;
	out_header.y_start  = y_start;
	out_header.z_start  = z_start;
	out_header.cell_a   = a_unit;
	out_header.cell_b   = b_unit;
	out_header.cell_c   = c_unit;
	out_header.alpha    = alpha * 180.0/PI;
	out_header.beta     = beta * 180.0/PI;
	out_header.gamma    = gamma * 180.0/PI;
	out_header.x_origin = x_origin;
	out_header.y_origin = y_origin;
	out_header.z_origin = z_origin;
	if((rho_min <= rho_max) && (0 >= std::min(rho_min, rho_max))){
		out_header.val_min = rho_min;
		out_header.val_max = rho_max;
	}
	out_header.mach_str = HOST_LITTLE_ENDIAN ? 17476 : 4369;
	
	out_file.write(reinterpret_cast<char*>(&out_header),sizeof(out_header));
	out_file.write(reinterpret_cast<char*>(densities.data()), densities.size() * sizeof(float));
	
	out_file.close();
	cout << "<- Finished writing, took " << seconds_since(runtime)*1000.0-file_reading_ms << " ms.\n\n";
	cout << "Done. Overall runtime was " << seconds_since(runtime)*1000.0 << " ms.\n";
}

#endif // INCLUDED_MAP_READER

