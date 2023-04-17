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

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <atomic>
#include <algorithm>

#ifdef PARALLELIZE
#include <omp.h>
#endif

using namespace std;

#include <time.h>
#ifndef _WIN32
// Time measurement
#include <sys/time.h>
#endif

template<typename T>
inline double seconds_since(T& time_start)
{
#ifndef _WIN32
	timeval time_end;
	gettimeofday(&time_end,NULL);
        double num_sec     = time_end.tv_sec  - time_start.tv_sec;
        double num_usec    = time_end.tv_usec - time_start.tv_usec;
        return (num_sec + (num_usec/1000000));
#else
	return 0.0;
#endif
}

template<typename T>
inline void start_timer(T& time_start)
{
#ifndef _WIN32
	gettimeofday(&time_start,NULL);
#endif
}

#define short_swap(byte_data) ((short int)((*((unsigned char*)byte_data)<<8) | *((unsigned char*)byte_data+1)))
#define read_short(byte_data) (endian_swap ? short_swap(byte_data) : *(reinterpret_cast<short int*>(byte_data)))

#define PI 3.14159265358979323846
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

inline std::string num2str(fp_num num)
{
	long l = (long)fastfloor(num);
	return to_string(l) + "." + to_string((long)fastfloor((num-l)*1000+0.5));
}

inline char* find_block(char* &brix_string)
{
	while(*brix_string == ' ') brix_string++; // remove white-space
	char* start = brix_string; // number begins here
	while(*brix_string != ' ') brix_string++; // find next space
	*brix_string++ = '\0'; // end number string and advance to next spot
	return start;
}

inline int brix_number(char* &brix_string)
{
	return stoi(find_block(brix_string));
}

inline fp_num brix_float(char* &brix_string)
{
	return stof(find_block(brix_string));
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

#define next_entry(name) { if(BRIX) brix_entry(brix_header, name, true); }

inline std::vector<fp_num> read_dsn6(
                                     std::string   filename,
                                     fp_num        map_x_center,
                                     fp_num        map_y_center,
                                     fp_num        map_z_center,
                                     unsigned int  map_x_dim,
                                     unsigned int  map_y_dim,
                                     unsigned int  map_z_dim,
                                     fp_num        grid_spacing,
                                     bool          write_map = false
                                    )
{
	timeval runtime;
	start_timer(runtime);
	
	std::vector<fp_num> densities;
	
	cout << "Reading map file " << filename << "\n";
	std::ifstream map_file(filename, std::ifstream::binary);
	if(map_file.fail()){
		cout << "\nERROR: Can't open map file " << filename << ".\n";
		exit(2);
	}
	std::streamoff filesize = map_file.tellg();
	map_file.seekg(0, std::ios::end);
	filesize = map_file.tellg() - filesize;
	map_file.seekg(0, std::ios::beg);
	cout << "\t-> file size: " << filesize << "\n";
	char header[DSN6_BLOCKSIZE];
	char* brix_header = header;
	bool BRIX = false;
	if(!map_file.read(header, DSN6_BLOCKSIZE)){
		cout << "\nERROR: Can't reader header.\n";
		exit(3);
	}
	if(brix_entry(brix_header, ":-)",false)){
		BRIX = true;
	}
	short int norm       = *(reinterpret_cast<short int*>(header+36));
	bool endian_swap     = (norm != 100);
	norm                 = read_short(header+36);
	if((norm != 100) && !BRIX){
		cout << "\nERROR: File is not a valid DSN6 file.\n";
		exit(4);
	}
	if(BRIX && (norm == 100)){ // found the one DSN6 file whose origin encodes a smiley
		BRIX=false;
	}
	if(BRIX){
		endian_swap = false;
		norm = 1;
		cout << "\t-> BRIX format\n";
	} else cout << "\t-> endian swap: " << endian_swap << ", Norm: " << norm << "\n";
	
	next_entry("origin");
	fp_num x_start       = BRIX ? brix_number(brix_header) : read_short(header);
	fp_num y_start       = BRIX ? brix_number(brix_header) : read_short(header+2);
	fp_num z_start       = BRIX ? brix_number(brix_header) : read_short(header+4);
	next_entry("extent");
	unsigned int x_dim   = BRIX ? brix_number(brix_header) : read_short(header+6);
	unsigned int y_dim   = BRIX ? brix_number(brix_header) : read_short(header+8);
	unsigned int z_dim   = BRIX ? brix_number(brix_header) : read_short(header+10);
	unsigned long xy_stride      = x_dim * y_dim;
	unsigned int x_file_stride   = (((x_dim&7)>0) + (x_dim>>3)) << 6;
	unsigned int xy_file_stride  = x_file_stride * (((y_dim&7)>0) + (y_dim>>3)) << 3;
	std::streamoff size_expected = xy_file_stride * (((z_dim&7)>0) + (z_dim>>3)) + 512;
	if(filesize < size_expected){
		cout << "\nERROR: Not enough data blocks in provided DSN6 file.\n";
		exit(5);
	}
	next_entry("grid");
	fp_num inv_x_step    = BRIX ? brix_number(brix_header) : read_short(header+12);
	fp_num inv_y_step    = BRIX ? brix_number(brix_header) : read_short(header+14);
	fp_num inv_z_step    = BRIX ? brix_number(brix_header) : read_short(header+16);
	fp_num x_step        = 1.0 / inv_x_step;
	fp_num y_step        = 1.0 / inv_y_step;
	fp_num z_step        = 1.0 / inv_z_step;
	next_entry("cell");
	fp_num a_unit        = BRIX ? brix_float(brix_header) : read_short(header+18);
	fp_num b_unit        = BRIX ? brix_float(brix_header) : read_short(header+20);
	fp_num c_unit        = BRIX ? brix_float(brix_header) : read_short(header+22);
	fp_num alpha         = BRIX ? brix_float(brix_header) : read_short(header+24);
	fp_num beta          = BRIX ? brix_float(brix_header) : read_short(header+26);
	fp_num gamma         = BRIX ? brix_float(brix_header) : read_short(header+28);
	next_entry("prod");
	fp_num rho_scale     = norm / (BRIX ? brix_float(brix_header) : (fp_num)read_short(header+30));
	next_entry("plus");
	fp_num offset        = BRIX ? brix_float(brix_header) : read_short(header+32);
	
	if(!BRIX){
		fp_num unit_scale    = 1.0 / read_short(header+34);
		a_unit              *= unit_scale;
		b_unit              *= unit_scale;
		c_unit              *= unit_scale;
		alpha               *= unit_scale;
		beta                *= unit_scale;
		gamma               *= unit_scale;
	}
	cout << "\t-> x_dim = " << x_dim << ", y_dim = " << y_dim << ", z_dim = " << z_dim << "\n";
	cout << "\t-> a_unit = " << a_unit << ", b_unit = " << b_unit << ", c_unit = " << c_unit << "\n";
	cout << "\t-> alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << "\n";
	
	alpha               *= PI / 180.0;
	beta                *= PI / 180.0;
	gamma               *= PI / 180.0;
	
	fp_num inv_a_unit    = 1/a_unit;
	fp_num inv_b_unit    = 1/b_unit;
	fp_num inv_c_unit    = 1/c_unit;
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
	fp_num density_y_end   = (y_start + y_dim - 1)  * y_step;
	fp_num density_z_end   = (z_start + z_dim - 1)  * z_step;
	f2c(density_x_start, density_y_start, density_z_start);
	f2c(density_x_end, density_y_end, density_z_end);
	cout << "\t-> density coordinate range: (" << density_x_start << ", " << density_y_start << ", " << density_z_start << ") A to (" << density_x_end << ", " << density_y_end << ", " << density_z_end << ") A\n";
	fp_num map_x_start  = map_x_center - map_x_dim * grid_spacing * 0.5;
	fp_num map_y_start  = map_y_center - map_y_dim * grid_spacing * 0.5;
	fp_num map_z_start  = map_z_center - map_z_dim * grid_spacing * 0.5;
	fp_num map_x_end    = map_x_start + map_x_dim * grid_spacing;
	fp_num map_y_end    = map_y_start + map_y_dim * grid_spacing;
	fp_num map_z_end    = map_z_start + map_z_dim * grid_spacing;
	// test if there is data for our grid box
	if((map_x_start < density_x_start) || (map_y_start < density_y_start) || (map_z_start < density_z_start) ||
	   (map_x_end   > density_x_end)   || (map_y_end   > density_y_end)   || (map_z_end   > density_z_end))
	{
		cout << "\nERROR: The specified grid box is (partially) outside of the file's density data.\n";
		exit(6);
	}
	densities.resize(xy_stride*z_dim);
	
	fp_num rho;
	fp_num rho_min = 255*rho_scale;
	fp_num rho_max = 0;
	unsigned z_block, y_block, x_block, data_offset;
	char* data_block     = new char[xy_file_stride];
	for(unsigned int z = 0; z < z_dim; z += 8){ // z slow
		z_block = std::min(z_dim - z, 8u);
		map_file.read(data_block, xy_file_stride);
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
							rho_min = std::min(rho, rho_min);
							rho_max = std::max(rho, rho_max);
						}
					}
				}
				data_offset += 512;
			}
		}
	}
	cout << "\t-> density range: " << rho_min << " to " << rho_max << "\n";
	delete[] data_block;
	map_file.close();
	double file_reading_ms = seconds_since(runtime)*1000.0;
	cout << "<- Finished reading densities, took " << file_reading_ms << " ms.\n\n";
	
	cout.precision(3);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Interpolating density data for " << map_x_dim << "x" << map_y_dim << "x" << map_z_dim << " grid (spacing: " << grid_spacing << " A)\n";
	cout << "\t-> grid start:  (" << map_x_start << ", " << map_y_start << ", " << map_z_start << ") A\n";
	cout << "\t-> grid size:   (" << map_x_dim * grid_spacing << ", " << map_y_dim * grid_spacing << ", " << map_z_dim * grid_spacing << ") A\n";
	
	fp_num grid_a, grid_b, grid_c, density;
	std::vector<fp_num> grid_map;
	grid_map.resize((map_x_dim + 1) * (map_y_dim + 1) * (map_z_dim + 1));
	fp_num* density_data = densities.data();
	fp_num* data_point;
#ifdef MGLTOOLS_MATH_COMPARISON
	fp_num ga, gb, gc;
#endif
	std::ofstream grid_file;
	if(write_map){
		std::size_t ext = filename.find_last_of(".");
		std::string map_name = filename.substr(0, ext) + ".map";
		grid_file.open(map_name);
		if(grid_file.fail()){
			cout << "Error: Can't open grid map output file " << map_name << ".\n";
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
	}
	for(unsigned int z = 0; z <= map_z_dim; z++){
		for(unsigned int y = 0; y <= map_y_dim; y++){
			for(unsigned int x = 0; x <= map_x_dim; x++){
				grid_a = map_x_start + x * grid_spacing;
				grid_b = map_y_start + y * grid_spacing;
				grid_c = map_z_start + z * grid_spacing;
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
				fp_num x_low  = fastfloor(grid_a);
				fp_num y_low  = fastfloor(grid_b);
				fp_num z_low  = fastfloor(grid_c);
				
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
				if(write_map) grid_file << num2str(density) << "\n";
			}
		}
	}
	cout << "<- Finished interpolating ";
	if(write_map){
		grid_file.close();
		cout << "and writing ";
	}
	cout << "grid map, took " << seconds_since(runtime)*1000.0 - file_reading_ms << " ms.\n\n";
	cout << "Done. Overall time was " << seconds_since(runtime)*1000.0 << " ms.\n";
	return grid_map;
}

#endif // INCLUDED_MAP_READER

