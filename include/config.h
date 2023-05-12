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


#ifndef INCLUDED_CONFIG
#define INCLUDED_CONFIG

#define USE_SINGLE_PRECISION
#define MODIFY_NORMALIZED_DENSITIES
#define PARALLELIZE

enum map_types_supported{
	automatic = 0,
	dsn6      = 1,
	dsn6_swap = 2,
	brix      = 3,
	mrc       = 4
};

enum density_modifier_fxns{
	no_modifier  = 0,
	log_modifier = 1
};

enum grid_map_write_modes{
	write_nothing  = 0,
	write_grid_ad4 = 1,
	write_grid_mrc = 2
};

#ifdef PARALLELIZE
#include <omp.h>
#endif

static unsigned short static_one = 1;
#define HOST_LITTLE_ENDIAN (*(unsigned char*)&static_one == 1)

#ifdef USE_SINGLE_PRECISION
typedef float fp_num;
#else
typedef double fp_num;
#endif

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

#define PI      3.14159265359f

#endif

