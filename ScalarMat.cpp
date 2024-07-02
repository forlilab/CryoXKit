/**************************************************************/
/* Parts of this code are from the FIGSiM project distributed */
/* under the University of Illinois/NCSA Open Source License. */
/* See LICENSE_UIUC file for details.                         */
/*                                                            */
/* Copyright (c) 2016 FIGSiM developers                       */
/*                                                            */
/* As part of Cryo2Grid this file is distributed under the    */
/* GNU LPGL-2.1-or-later Open Source License.                 */
/* See LICENSE file for details.                              */
/*                                                            */
/* Copyright (c) Andreas F. Tillack                           */
/*               Forli Lab @ Scripps Research                 */
/**************************************************************/


#include "include/ScalarMat.h"

#define PHI 0x9e3779b9

static uint32_t Q[4096], c=362436; // CMWC4096 variables
static __int32_t hz; // Ziggurat variables
static __uint32_t ic, iz, kn[128];
static double wn[128],fn[128];
static __uint64_t rancount=0;

void init_CMWC4096(__uint32_t x)
{
	init_CMWC4096(x, Q, ic, c);
}

void init_CMWC4096(__uint32_t x, uint32_t* Q, uint32_t &ic, uint32_t &c)
{
	ic=4095;
	c=362436;
	hz=8;
	Q[0]=x;
	Q[1]=Q[0]+PHI;
	Q[2]=Q[1]+PHI;
	for(unsigned int i=3; i<4096; i++) Q[i]=Q[i-3] ^ Q[i-2] ^ PHI ^ i;
	for(unsigned int i=0; i<4096; i++) CMWC4096(ic,c,Q);
	rancount=0;
}

//see G. Marsaglia, JMASM May 2003, Vol 2, No 1, 2-13
__uint32_t CMWC4096()
{
	return CMWC4096(ic,c,Q);
}
__uint32_t CMWC4096(__uint32_t &ic, __uint32_t &c, uint32_t* Q)
{
	__uint32_t x, r = 0xfffffffe;
	__uint64_t t, a = 18782LL;
	ic = (ic+1)&4095;
	t = a*Q[ic]+c;
	c = (t>>32);
	x = t+c;
	if(x<c){
		x++;
		c++;
	}
	rancount++;
	return(Q[ic] = r-x);
}

__uint64_t ran_count()
{
	return rancount;
}

double ranQ(__uint32_t &ic, __uint32_t &c, uint32_t* Q)
{
	return CMWC4096(ic,c,Q)/4294967296.0;
}

double ranQ()
{
	return CMWC4096()/4294967296.0;
}

/* nfix() generates variates from the residue when rejection in RNOR occurs. */
inline double nfix_wtail()
{
	const double r = 3.442619855899; /* The start of the right tail */
	const double inv_r = 1/r;
	static double x, y;
	for(;;){
		x=hz*wn[iz];
		/* iz==0, handles the base strip */
		if(iz==0){
			do{
				x=-log(ranQ())*inv_r;
				y=-log(ranQ());
			} while(y+y<x*x);
			if(hz>0) return r+x; else return -r-x;
		}
		/* iz>0, handle the wedges of other strips */
		if(fn[iz]+ranQ()*(fn[iz-1]-fn[iz]) < exp(-.5*x*x)) return x;
		/* initiate, try to exit for(;;) for loop*/
		hz=CMWC4096();
		iz=hz&127;
		if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
	}
}

/* nfix() generates variates from the residue when rejection in RNOR occurs. */
inline double nfix_notail()
{
	static double x;
	for(;;){
		x=hz*wn[iz];
		if((iz>0) && (fn[iz]+ranQ()*(fn[iz-1]-fn[iz]) < exp(-.5*x*x))) return x;
		/* initiate, try to exit for(;;) for loop*/
		hz=CMWC4096();
		iz=hz&127;
		if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
	}
}

double ran_n()
{
	hz=(__int32_t)CMWC4096();
	iz=hz&127;
	if(fabs(hz)<kn[iz]) return hz*wn[iz]; else return nfix_wtail();
}

double ran_n2()
{
	hz=(__int32_t)CMWC4096();
	iz=hz&127;
	if(fabs(hz)<kn[iz]) return hz*wn[iz]; else return nfix_notail();
}

/*--------This procedure sets the seed and creates the tables------*/
void zigset()
{
	const double m1 = 2147483648.0;
	double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
	rancount=0;
	
	/* Set up tables for RNOR */
	q=vn/exp(-.5*dn*dn);
	kn[0]=(__uint32_t)((dn/q)*m1);
	kn[1]=0;
	
	wn[0]=q/m1;
	wn[127]=dn/m1;
	
	fn[0]=1.;
	fn[127]=exp(-.5*dn*dn);
	
	for(unsigned int i=126;i>=1;i--){
		dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
		kn[i+1]=(__uint32_t)((dn/tn)*m1);
		tn=dn;
		fn[i]=exp(-.5*dn*dn);
		wn[i]=dn/m1;
	}
}

/// Averages a 1D array of length numel
double average(double *A, int numel)
{
	double temp = 0;
	
	for(int i = 0; i<numel; i++)
	{
		temp += A[i];
	}
	
	double avg = temp/(double)numel;
	return avg;
}

/// Averages a range of a 1D array
double average_range(double *A, int first, int last)
{
	double temp = 0;
	
	for(int i = first; i<last; i++)
	{
		temp += A[i];
	}
	
	double avg = temp/(double)(last-first);
	return avg;
}

/// Averages a 2D array
void average2x(double **A, double *B, int numelx, int numely)
{
	double temp;
	
	for(int i = 0; i<numely; i++)
	{
		temp = 0;
		B[i] = 0;
		
		for(int j = 0; j<numelx; j++)
		{
			temp +=A[i][j];
		}
		
		B[i] += temp/(double)numelx;
	}
}

/// Averages a range of a 2D array
void average_range2x(double **A, double *B, int firstx, int lastx, int firsty, int lasty)
{
	double temp;
	
	for(int i = firsty; i<lasty; i++)
	{
		temp = 0;
		B[i] = 0;
		
		for(int j = firstx; j<lastx; j++)
		{
			temp +=A[i][j];
		}
		
		B[i] += temp/(double)(lastx-firstx);
	}
}

/// Sums an array of ints
int totalcounti(int *A, int numel)
{
	int tc = 0;
	
	for(int i = 0; i<numel; i++)
	{
		tc += A[i];
	}
	
	return tc;
}

/// sums an array of doubles
double totalcountd(double *A, int numel)
{
	double tc = 0;
	
	for(int i = 0; i<numel; i++)
	{
		tc += A[i];
	}
	
	return tc;
}

/// Advance line function -- moves string pointer to after next \n
void advanceline(fstream &file, int nline)
{
	int i;
	char spt;
	
	for(i=0; i<nline;i++)
	{
		do file.get(spt);
		while (spt != '\n');
	}
}

/// Get tab-delimited int function -- reads in chars until encounters delimiter, then converts to int
int gettdint(fstream &file)
{
	char ch;
	string buffer;
	while (file.get(ch) && ch != '\t' && ch != '\n')
	{
		buffer += ch;
	}
	
	int i = atoi(buffer.c_str());
	return i;
}

/// Get tab-delimited double function -- reads in chars until encounters delimiter, then converts to double
double gettddoub(fstream &file)
{
	char ch;
	string buffer;
	while (file.get(ch) && ch != '\t' && ch != '\n')
	{
		buffer += ch;
	}
	
	double d = atof(buffer.c_str());
	return d;
}

