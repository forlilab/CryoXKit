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


#ifndef INCLUDED_MAP_MODIFIER
#define INCLUDED_MAP_MODIFIER

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

inline void modify_densities(
                             std::vector<fp_num> &densities,
                             const int            mod_fxn    = no_modifier,
                             const fp_num*        fxn_params = NULL
                            )
{
	if(mod_fxn == no_modifier) return;
	timeval runtime;
	start_timer(runtime);
	cout << "Adjusting density values";
	const unsigned int nr_points = (densities[0] + 1) * (densities[1] + 1) * (densities[2] + 1) + 9;
#ifdef MODIFY_NORMALIZED_DENSITIES
	fp_num norm  = 10.0 / (std::max(fabs(densities[7]), fabs(densities[8])));
	const fp_num dens_min        = densities[7] * norm;
	const fp_num dens_max        = densities[8] * norm;
#else
	fp_num norm  = 1;
	const fp_num dens_min        = densities[7];
	const fp_num dens_max        = densities[8];
#endif
	fp_num rho_min               = 1e80;
	fp_num rho_max               = 0;
	fp_num density;
	if(mod_fxn == log_modifier){
		cout << " using logistics function modifier\n";
		const fp_num log_max    = fxn_params[0];
		const fp_num rate       = fxn_params[1];
		const fp_num exp_shift  = dens_max - (dens_max - dens_min) * fxn_params[2]; // rho_max - x0*(rho_max - rho_min)
		for(unsigned i=9; i<nr_points; i++){
			density = log_max / (1.0 + exp(rate * (exp_shift - densities[i] * norm)));
			densities[i] = density;
			rho_min = std::min(density, rho_min);
			rho_max = std::max(density, rho_max);
		}
	}
	densities[7] = rho_min;
	densities[8] = rho_max;
	cout << "<- Finished adjusting, took " << seconds_since(runtime)*1000.0 << " ms.\n\n";;
}

#endif // INCLUDED_MAP_MODIFIER

