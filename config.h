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

// Special care needs to be taken when compiling under Windows
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define IN_WINDOWS
#include <stdint.h>
typedef uint32_t __uint32_t;
typedef uint64_t __uint64_t;
typedef int32_t __int32_t;
typedef int64_t __int64_t;
#endif

//Define constants here
#define PI      3.14159265359f
#define pi2     6.28318530718f
#define invpi   0.318309886184f
#define invpi2  0.159154943092f
#define pi_half 1.57079632679f
#define EPS 1.2e-7f
#define NA 6.02214179E23f
#define two_to_onesixth 1.12246204831f
#define e_in_esu 4.80320427f // 1 e in 10^(-10) stat coulomb (= 10^-(10) esu)
// 1 atm = 101.325 kPa = 1.01325x10^5 J/m^3 ; 1 pErg/Angström^3 = 1E-19 J / (1E-30 m^3) = 1E11 J/m^3 = 1E11 Pa => 1 Pa = 1E-11 pErg/Angström^3 => 1 atm = 1.01325E-6 pErg/Angström^3
#define atm_in_pErg_per_Ang3 1.01325E-6f
#define perg_to_kJ_per_mol (NA*1E-19/1000)
#define perg_to_kcal_per_mol (perg_to_kJ_per_mol/4.184f)
#define kB 1.3806488e-4f
// 1 Debye = 10^-21 C*m^2/s / c 
// => perg/Debye = 10^-19 J/(10^-21 C*m^2/s)*2.99792458*10^8*m/s = 2.99792458x10^10 N/C
// 1 MV/m = 10^6 V/m = 10^6 N/C (CV = J => CV/m = N => N/C = V/m)
// => 1 MV/m = 10^6/(2.99792458*10^10) perg/Debye = 3.335646x10^-5 perg/D
#define MV_per_m_to_perg_per_Debye 3.335646E-5f
// hbar = 1.05457E-34 Js
#define hbar_perg_ps (1.05457E-34f/(1E-19*1E-12))
// 1 amu = 1.66053886E-27 J/(m/s)^2 -> Js^2/m^2
#define amu_to_perg_ps2_per_Ang2 (1.66053886E-27f/(1E-19*1E-24)*(1E-20))

#endif

