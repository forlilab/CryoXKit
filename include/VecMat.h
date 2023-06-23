/**************************************************************/
/* Parts of this code are from the FIGSiM project distributed */
/* under the University of Illinois/NCSA Open Source License. */
/* See LICENSE_UIUC file for details.                         */
/*                                                            */
/* Copyright (c) 2016 FIGSiM developers                       */
/*                                                            */
/* As part of AD-NG this file is distributed under the        */
/* GNU LPGL-2.1-or-later Open Source License.                 */
/* See LICENSE file in top directory for details.             */
/*                                                            */
/* Copyright (c) 2020 ADsandbox developers:                   */
/*           Dr. Andreas F. Tillack                           */
/*           Dr. Diogo M. Santos-Martins                      */
/*           Dr. Stefano Forli                                */
/*           Forli Lab @ Scripps Research                     */
/**************************************************************/


#ifndef INCLUDED_VECTORMAT
#define INCLUDED_VECTORMAT

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif

#include <complex>
#include <cstdio>
#include <ctime>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

const double sqrt3 = sqrt(3.0);

template<class T = double>
class Mat33; // Forward declaration of class Mat33 - AT Feb 5, 2011

/// Vec3 class - dedicated class for 3-dimensional vectors
template<class T = double>
class Vec3
{
	public:
		T vec[3];
		friend class Mat33<T>;
		
	// Constructors and destructors
		Vec3();
		Vec3(const T, const T, const T);
		Vec3(const T);
		Vec3(const T*); // construct from array
		template<typename U>
		Vec3(const Vec3<U> &);
		~Vec3() {};
	
	// Utility functions
		template<typename U>
		Vec3& operator= (const Vec3<U> &);
		void V3Zeros();
		void V3Swap(Vec3 &);
		
	//Comparison operators
		bool operator== (const Vec3 &);
		bool operator!= (const Vec3 &);

	//Level 1 BLAS
		Vec3 operator+ (const T);
		Vec3& operator+= (const T);
		Vec3 operator- (const T);
		Vec3& operator-= (const T);
		Vec3 operator* (const T);
		Vec3& operator*= (const T);
		Vec3 operator/ (const T);
		Vec3& operator/= (const T);
		
		Vec3 operator+ (const Vec3 &);
		template<typename U>
		Vec3& operator+= (const Vec3<U> &);
		Vec3 operator- (const Vec3 &);
		Vec3& operator-= (const Vec3 &);
		template<typename U>
		T operator* (const Vec3<U> &) const;

	//Level 2 BLAS
		//These do not compile properly. Ask someone who knows C++ better than I do. - LEJ 01/07/10
		// -- solved by AT Feb 5, 2011 (needed forward declaration of Mat33 class)
		Vec3 operator* (const Mat33<T> &);
		Vec3& operator*= (const Mat33<T> &);

	//Special functions
		Vec3 V3Cross(const Vec3 &);
		Mat33<T> V3TensProd(const Vec3 &); //Does not compile properly. Ask someone who knows C++ better than I do. - LEJ 01/07/10 - does now (AT, Feb 5, 2011)
		T V3Norm() const;
		T V3Sum() const;
		T V3Prod() const;
		std::string V3Str() const
		{
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed <<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2];
			
			std::string vecstring = converter.str();
			return vecstring;
		};
		std::string V3Str(const char &dlm, const unsigned int prec = 6) const
		{
			//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
			char delim;
			if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))) {
				delim = '\t';
			}
			else delim = dlm;
			
			std::stringstream converter;
			converter << std::fixed << std::setprecision(prec) <<vec[0]<< delim <<vec[1]<< delim <<vec[2];
			
			std::string vecstring = converter.str();
			return vecstring;
		};
};

/// complex class - 3d complex vector class, only minimal function set defined needed for eigenvalue problem
template<class T = double>
class CVec3{
	public:
		complex<T> cvec[3];
	// Constructors and destructor
		CVec3();
		CVec3(const complex<T>, const complex<T>, const complex<T>);
		CVec3(const Vec3<T> &);
		CVec3(const CVec3 &);
		~CVec3(){};
	// Operators
		bool operator== (const CVec3 &);
		bool operator!= (const CVec3 &);
		bool operator== (const Vec3<T> &);
		bool operator!= (const Vec3<T> &);
		CVec3& operator = (const CVec3 &);
		CVec3 operator- (const T); 
		CVec3& operator-= (const T);
	// Functions
		Vec3<T> Re();
		Vec3<T> Im();
		Vec3<T> Abs();
		std::string CV3Str() const
		{
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed <<cvec[0]<<"\t"<<cvec[1]<<"\t"<<cvec[2];
			
			return converter.str();
		};
		std::string CV3Str(const char &dlm) const
		{
			//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
			char delim;
			
			if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))){
				delim = '\t';
			} else delim = dlm;
			
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed << cvec[0] << delim << cvec[1] << delim <<cvec[2];
			
			return converter.str();
		};
};

///Four-element vector (needed for X3D rotations)
template<class T = double>
class Vec4 {
	public:
		T vec[4];
		friend class Vec3<T>;
		friend class Mat33<T>;
	//Constructors and destructors
		Vec4();
		template<typename U>
		Vec4(const U, const U, const U, const U);
		Vec4(const Vec3<T> &, const T);
		Vec4(const Vec4<T> &);
		~Vec4() {};
	//Operators
		bool operator== (const Vec4 &);
		bool operator!= (const Vec4 &);
		Vec4 operator+ (const Vec4 &);
		Vec4 operator* (const T);
		Vec4 operator/ (const T);
	//Functions
		std::string V4Str() const
		{
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed <<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<"\t"<<vec[3];
			
			std::string vecstring = converter.str();
			return vecstring;
		};
		std::string V4Str(const char &dlm) const
		{
			//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
			char delim;
			if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))) {
				delim = '\t';
			}
			else delim = dlm;
			
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed <<vec[0]<< delim <<vec[1]<< delim <<vec[2] << delim << vec[3];
			
			std::string vecstring = converter.str();
			return vecstring;
		};
};

/// Mat33 class - for when you need to do transformations on 3-vectors…
template<class T>
class Mat33
{
	public:
		T mat[3][3];
		friend class Vec3<T>;
		
		//Constructors and destructors
		inline Mat33();
		Mat33(const Mat33 &);
		Mat33(const Vec3<T> &, const Vec3<T> &);
		Mat33(const T);
		Mat33(const T a, const T b, const T c); // diagonal constructor
		Mat33(const T*);
		~Mat33(){};
		
		//Utility functions
		template<typename U>
		Mat33& operator= (const Mat33<U> &);
		void M3Eye();
		inline void M3Zeros();
		
		//Comparison operators
		bool operator== (const Mat33 &);
		bool operator!= (const Mat33 &);
		
		//Matrix-scalar operations
		Mat33 operator+ (const T);
		Mat33& operator+= (const T);
		Mat33 operator- (const T);
		Mat33& operator-= (const T);
		Mat33 operator* (const T);
		Mat33& operator*= (const T);
		Mat33 operator/ (const T);
		Mat33& operator/= (const T);
		
		//Level 2 BLAS
		template<typename U>
		Vec3<T> operator* (const Vec3<U> &);
		//No *= operator due to datatype mismatch.
		Mat33 operator- (const Vec3<T> &); // subtract vector from diagonal
		Mat33& operator-= (const Vec3<T> &); // subtract vector from diagonal
		
		//Level 3 BLAS
		Mat33 operator+ (const Mat33 &);
		Mat33& operator+= (const Mat33 &);
		Mat33 operator- (const Mat33 &);
		Mat33& operator -= (const Mat33 &);
		template<typename U>
		Mat33 operator* (const Mat33<U> &);
		template<typename U>
		Mat33& operator*= (const Mat33<U> &);
		//These do not compile properly. Ask someone who knows C++ better than I do. - LEJ 01/07/10
		Mat33 operator/ (const Mat33 &);
		Mat33& operator/= (const Mat33 &);
		
		//Special functions
		CVec3<T> Eigenvalues();
		Mat33 Eigenvectors(Vec3<T> &ew){ bool* multiples=new bool[3]; Mat33 A=Eigenvectors(ew,multiples,false); delete[] multiples; return A; };
		Mat33 Eigenvectors(Vec3<T> &ew, bool normalize){ bool* multiples=new bool[3]; Mat33 A=Eigenvectors(ew,multiples,normalize); delete[] multiples; return A; };
		Mat33 Eigenvectors(Vec3<T> &ew, bool* multiples, bool normalize);
		inline Mat33 M3MulDiag(const T a, const T b, const T c);
		Mat33 M3Transpose();
		Mat33 M3RowSwap(const int, const int);
		inline Vec3<T> ColumnVec3(const int);
		inline Mat33 M3Inv(bool &sing);
		inline Mat33 MulSymM3(const Mat33 &C);
		inline Mat33 TransMulM3(const Mat33 &C);
		inline Vec3<T> TransMulVec(const Vec3<T> &u);
		inline Mat33 SymMulM3(const Mat33 &C);
		inline Mat33 SymMulSymM3(const Mat33 &C);
		inline Mat33 SymM3Inv(bool &sing);
		inline Vec3<T> SymM3InvMult(bool &sing, Vec3<T> &z);
		T M3Trace();
		inline T M3Det();
		Vec3<T> M3Diag();
		void M3GEPP(Vec3<T> &x);
		Vec3<T> M3LinSolve(Vec3<T> &u,bool &multiples);
		void LUDecomposition();
		bool M3BackSub(Vec3<T> &x); // returns true if multiples of vector are also solutions
		std::string M3Str()
		{
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed <<mat[0][0]<<"\t"<<mat[0][1]<<"\t"<<mat[0][2]<<"\n";
			converter << std::fixed <<mat[1][0]<<"\t"<<mat[1][1]<<"\t"<<mat[1][2]<<"\n";
			converter << std::fixed <<mat[2][0]<<"\t"<<mat[2][1]<<"\t"<<mat[2][2];
			
			std::string matstring = converter.str();
			return matstring;
		};
		std::string M3RowStr(const int i)
		{
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed <<mat[i][0]<<"\t"<<mat[i][1]<<"\t"<<mat[i][2];
			
			std::string matstring = converter.str();
			return matstring;
		};
		std::string M3RowStr(const int i, const char &dlm)
		{
			//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
			char delim;
			
			if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))) {
				delim = '\t';
			} else delim = dlm;
			
			std::stringstream converter;
			std::cout.precision(6);
			
			converter << std::fixed << mat[i][0] << delim << mat[i][1] << delim <<mat[i][2];
			
			std::string matstring = converter.str();
			return matstring;
		};
};

/*	Function implementations are below. Math functions and constructors are inlined for speed.
 *	Non-math functions are kept in .cpp file. The vast majority of the class is contained in
 *	This file, organized in the same order that the prototypes are listed
 */

/// get solutions to third order polynomial in reduced form: z^3 + p*z + q = 0 -- AT
template<typename T>
CVec3<T> SolvePolynomial3(const T p, const T q)
{
	T D=q*q/4.0+p*p*p/27.0;
#if DEBUG_LEVEL>2
	cout << "Start SolvePolynomial3\n";
	cout << "p = " << p << ", q = " << q << " => D = " << D << "\n";
#endif
	T u, v;
	T minusqhalf=-0.5*q;
	CVec3<T> cv;
	
	if(D>EPS*EPS){ // one real, two complex solutions
		T sqrtD=sqrt(D);
		u=cbrt(sqrtD+minusqhalf);
		v=cbrt(-sqrtD+minusqhalf);
#if DEBUG_LEVEL>2
	cout << "u = " << u << ", v = " << v << "\n";
#endif
		cv.cvec[0]=u+v;
		cv.cvec[1]=complex<T>(-0.5*(u+v),0.5*sqrt3*(u-v));
		cv.cvec[2]=conj(cv.cvec[1]);
	} else{
		if(D<-EPS*EPS){ // three different real solutions, now things get complex to calculate ;-)
			// for any complex number z = Re+Im*i = A*e^i*theta ; A = sqrt(Re^2+Im^2) ; theta = atan(Im/Re)
			// sqrt(D) = i sqrt(-D)
			T sqrtD=sqrt(-D); // is imaginary now ...
			// u = cbrt(minusqhalf+i*sqrt(-D)) = cbrt(sqrt(minusqhalf^2-D)*e^(i*atan(sqrt(-D)/minusqhalf)))
			// => u = cbrt(sqrt(minusqhalf^2-D)*e^(i*1/3*atan(sqrt(-D)/minusqhalf))
			// similar for v: v = cbrt(sqrt(minusqhalf^2-D))*e^(-i*1/3*atan(sqrt(-D)/minusqhalf));
			T uv_length=cbrt(sqrt(minusqhalf*minusqhalf-D)); // cubicroot(sqrt(Re^2 + Im^2) (Im^2 = sqrt(-D)^2)
			T uv_phase=atan2(sqrtD,minusqhalf)/3.0;
#if DEBUG_LEVEL>2
			cout << "magnitude = " << uv_length << ", phase = " << uv_phase << "\n";
#endif
			u=uv_length*cos(uv_phase); // here, u is real part
			v=uv_length*sin(uv_phase); // here, v is imaginary part
#if DEBUG_LEVEL>2
			cout << "u = v* = (" << u << ", " << v << ")\n";
#endif
			cv.cvec[0]=u+u;
			cv.cvec[1]=sqrt3*v-u;
			cv.cvec[2]=-sqrt3*v-u;
		} else{ // D=0: either threefold real solution, or real and two-fold real solution (two distinct solutions)
#if DEBUG_LEVEL>2
			cout << "u = v = " << cbrt(minusqhalf) << "\n";
#endif
			cv.cvec[0]=0.0;
			if(fabs(p)>EPS) cv.cvec[0]=3*q/p;
			cv.cvec[1]=cv.cvec[0];
			cv.cvec[1]*=-0.5;
			cv.cvec[2]=cv.cvec[1];
		}
	}
#if DEBUG_LEVEL>2
	cout << "Finished SolvePolynomial3.\n";
#endif
	return cv;
}

/*	***Vec3***	*/

//Constructors and Destructors

/// Default constructor (zeros)
template<typename T>
inline Vec3<T>::Vec3()
{
	vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;
}

/// Construct from doubles
template<typename T>
inline Vec3<T>::Vec3(const T a, const T b, const T c)
{
	vec[0] = a;
	vec[1] = b;
	vec[2] = c;
}

template<typename T>
inline Vec3<T>::Vec3(const T d)
{
	vec[0] = d; vec[1] = d; vec[2] = d;
}
/// Construct from an array of 3 doubles
template<typename T>
inline Vec3<T>::Vec3(const T* A)
{
	vec[0] = A[0];
	vec[1] = A[1];
	vec[2] = A[2];
}

/// Copy constructor
template<typename T>
template<typename U>
inline Vec3<T>::Vec3(const Vec3<U> &u)
{
	vec[0] = (T)u.vec[0];
	vec[1] = (T)u.vec[1];
	vec[2] = (T)u.vec[2];
}


//Non-BLAS utility functions

/// Zero vector v = [0 0 0]
template<typename T>
inline void Vec3<T>::V3Zeros()
{
	vec[0] = 0.0;
	vec[1] = 0.0;
	vec[2] = 0.0;
}

/// Assignment operator
template<typename T>
template<typename U>
inline Vec3<T>& Vec3<T>::operator=(const Vec3<U> &u)
{
	vec[0] = (T)u.vec[0];
	vec[1] = (T)u.vec[1];
	vec[2] = (T)u.vec[2];
	return *this;
}

/// Swap values in two vectors u <-> v
template<typename T>
inline void Vec3<T>::V3Swap(Vec3<T> &u)
{
	T temp;
	for (int i = 0; i < 3; i ++)
	{
		temp = vec[i];
		vec[i] = u.vec[i];
		u.vec[i] = temp;
	}
}


// Comparison operators
template<typename T>
inline bool Vec3<T>::operator==(const Vec3<T> &u)
{
	if (fabs(u.vec[0]-vec[0])>EPS) return false;
	if (fabs(u.vec[1]-vec[1])>EPS) return false;
	if (fabs(u.vec[2]-vec[2])>EPS) return false;
	return true;
}

template<typename T>
inline bool Vec3<T>::operator!=(const Vec3<T> &u)
{
	if (fabs(u.vec[0]-vec[0])>EPS) return true;
	if (fabs(u.vec[1]-vec[1])>EPS) return true;
	if (fabs(u.vec[2]-vec[2])>EPS) return true;
	return false;
}


//Level 1 BLAS

/// Vector scalar addition v = u + a
template<typename T>
inline Vec3<T> Vec3<T>::operator+ (const T a)
{
	Vec3<T> v(a+vec[0], a+vec[1], a+vec[2]);
	return v;
}

/// Vector scalar addition u = u + a
template<typename T>
inline Vec3<T>& Vec3<T>::operator += (const T a)
{
	vec[0] +=a;
	vec[1] +=a;
	vec[2] +=a;
	return *this;
}

/// Vector scalar subtraction v = u - a
template<typename T>
inline Vec3<T> Vec3<T>::operator- (const T a)
{
	Vec3<T> v(vec[0]-a, vec[1]-a, vec[2]-a);
	return v;
}

/// Vector scalar subtraction u = u - a
template<typename T>
inline Vec3<T>& Vec3<T>::operator -= (const T a)
{
	vec[0] -=a;
	vec[1] -=a;
	vec[2] -=a;
	return *this;
}

/// Vector scalar multiplication v = u*a
template<typename T>
inline Vec3<T> Vec3<T>::operator* (const T a)
{
	Vec3<T> v((T)a*vec[0], (T)a*vec[1], (T)a*vec[2]);
	return v;
}


/// Vector scalar multiplication u = u*a
template<typename T>
inline Vec3<T>& Vec3<T>::operator*= (const T a)
{
	vec[0] *=a;
	vec[1] *=a;
	vec[2] *=a;
	return *this;
}

/// Vector scalar division v = u/a
template<typename T>
inline Vec3<T> Vec3<T>::operator/ (const T a)
{
	T b = 1.0/a;
	Vec3<T> v(vec[0]*b, vec[1]*b, vec[2]*b);
	return v;
}

/// Vector scalar division u = u/a
template<typename T>
inline Vec3<T>& Vec3<T>::operator/= (const T a)
{
	T b = 1.0/a;
	vec[0] *=b;
	vec[1] *=b;
	vec[2] *=b;
	return *this;
}

/// Vector addition v = u + w
template<typename T>
inline Vec3<T> Vec3<T>::operator+ (const Vec3<T> &w)
{
	Vec3<T> v(vec[0]+w.vec[0], vec[1]+w.vec[1], vec[2] +w.vec[2]);
	return v;
}

/// Vector addition u = u + w
template<typename T>
template<typename U>
inline Vec3<T>& Vec3<T>::operator+= (const Vec3<U> &w)
{
	vec[0] += (T)w.vec[0];
	vec[1] += (T)w.vec[1];
	vec[2] += (T)w.vec[2];
	return *this;
}

/// Vector subtraction v = u - w
template<typename T>
inline Vec3<T> Vec3<T>::operator- (const Vec3<T> &w)
{
	Vec3<T> v(vec[0]-w.vec[0], vec[1]-w.vec[1], vec[2] -w.vec[2]);
	return v;
}

/// Vector subtraction u = u - w
template<typename T>
inline Vec3<T>& Vec3<T>::operator-= (const Vec3<T> &w)
{
	vec[0] -= w.vec[0];
	vec[1] -= w.vec[1];
	vec[2] -= w.vec[2];
	return *this;
}

/// Dot product a = v'*u
template<typename T>
template<typename U>
inline T Vec3<T>::operator* (const Vec3<U> &B) const
{
	T a = vec[0]*(T)B.vec[0] + vec[1]*(T)B.vec[1] + vec[2]*(T)B.vec[2];
	return a;
}

//Level 2 BLAS

/// Vector-Matrix Multiply v = u'*A
template<typename T>
inline Vec3<T> Vec3<T>::operator* (const Mat33<T> &A)
{
	Vec3<T> v(vec[0]*A.mat[0][0]+vec[1]*A.mat[0][1]+vec[2]*A.mat[0][2],
	          vec[0]*A.mat[1][0]+vec[1]*A.mat[1][1]+vec[2]*A.mat[1][2],
	          vec[0]*A.mat[2][0]+vec[1]*A.mat[2][1]+vec[2]*A.mat[2][2]);
	return v;
}

/// Vector-Matrix multiply u' = u'*A
template<typename T>
inline Vec3<T>& Vec3<T>::operator*= (const Mat33<T> &A)
{
	Vec3<T> u(*this);
	vec[0] = u.vec[0]*A.mat[0][0]+u.vec[1]*A.mat[0][1]+u.vec[2]*A.mat[0][2];
	vec[1] = u.vec[0]*A.mat[1][0]+u.vec[1]*A.mat[1][1]+u.vec[2]*A.mat[1][2];
	vec[2] = u.vec[0]*A.mat[2][0]+u.vec[1]*A.mat[2][1]+u.vec[2]*A.mat[2][2];
	return *this;
}


//Special functions

/// Vector norm |v|
template<typename T>
inline T Vec3<T>::V3Norm() const
{
	return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

/// Cross product v = u x w
template<typename T>
inline Vec3<T> Vec3<T>::V3Cross(const Vec3<T> &w)
{
	Vec3<T> v(vec[1]*w.vec[2] - vec[2]*w.vec[1],
	          vec[2]*w.vec[0] - vec[0]*w.vec[2],
	          vec[0]*w.vec[1] - vec[1]*w.vec[0]);
	vec[0] = v.vec[0];
	vec[1] = v.vec[1];
	vec[2] = v.vec[2];
	return v;
}

/// Tensor product A = u*w'
template<typename T>
inline Mat33<T> Vec3<T>::V3TensProd(const Vec3<T> &w)
{
	Mat33<T> A;
	A.M3Zeros();
	
	for(int i = 0; i < 3; i++){
		A.mat[i][0] += vec[i]*w.vec[0];
		A.mat[i][1] += vec[i]*w.vec[1];
		A.mat[i][2] += vec[i]*w.vec[2];
	}
	
	return A;
}

/// Sum of vector elements
template<typename T>
inline T Vec3<T>::V3Sum() const
{
	T a = vec[0]+vec[1]+vec[2];
	return a;
}

/// Product of vector elements
template<typename T>
inline T Vec3<T>::V3Prod() const
{
	T a = vec[0]*vec[1]*vec[2];
	return a;
}

/*	***CVec3***	*/

// Constructors

/// Default constructor (everything is zero)
template<typename T>
inline CVec3<T>::CVec3(){
	cvec[0] = complex<T>(0,0);
	cvec[1] = complex<T>(0,0);
	cvec[2] = complex<T>(0,0);
}

/// construct from complex<T>
template<typename T>
inline CVec3<T>::CVec3(const complex<T> a, const complex<T> b, const complex<T> c){
	cvec[0] = a;
	cvec[1] = b;
	cvec[2] = c;
}

/// construct from Vec3
template<typename T>
inline CVec3<T>::CVec3(const Vec3<T> &u){
	cvec[0] = u.vec[0];
	cvec[1] = u.vec[1];
	cvec[2] = u.vec[2];
}

/// Copy constructor
template<typename T>
inline CVec3<T>::CVec3(const CVec3<T> &u){
	cvec[0] = u.cvec[0];
	cvec[1] = u.cvec[1];
	cvec[2] = u.cvec[2];
}

// Comparison operators
template<typename T>
inline bool CVec3<T>::operator==(const CVec3<T> &u)
{
	if (u.cvec[0] != cvec[0]) return false;
	if (u.cvec[1] != cvec[1]) return false;
	if (u.cvec[2] != cvec[2]) return false;
	return true;
}

template<typename T>
inline bool CVec3<T>::operator!=(const CVec3<T> &u)
{
	if (u.cvec[0] != cvec[0]) return true;
	if (u.cvec[1] != cvec[1]) return true;
	if (u.cvec[2] != cvec[2]) return true;
	return false;
}

template<typename T>
inline bool CVec3<T>::operator==(const Vec3<T> &u)
{
	if (cvec[0] != complex<T>(u.vec[0])) return false;
	if (cvec[1] != complex<T>(u.vec[1])) return false;
	if (cvec[2] != complex<T>(u.vec[2])) return false;
	return true;
}

template<typename T>
inline bool CVec3<T>::operator!=(const Vec3<T> &u)
{
	if (cvec[0] != complex<T>(u.vec[0])) return true;
	if (cvec[1] != complex<T>(u.vec[1])) return true;
	if (cvec[2] != complex<T>(u.vec[2])) return true;
	return true;
}

/// Assignment operator
template<typename T>
inline CVec3<T>& CVec3<T>::operator=(const CVec3<T> &u){
	cvec[0] = u.cvec[0];
	cvec[1] = u.cvec[1];
	cvec[2] = u.cvec[2];
	return *this;
}

/// Scalar subtraction z = u - d
template<typename T>
inline CVec3<T> CVec3<T>::operator- (const T d)
{
	CVec3<T> cv(cvec[0]-d, cvec[1]-d, cvec[2]-d);
	return cv;
}

/// Scalar subtraction z = u - d
template<typename T>
inline CVec3<T>& CVec3<T>::operator -= (const T d)
{
	cvec[0] -=d;
	cvec[1] -=d;
	cvec[2] -=d;
	return *this;
}

// Functions

template<typename T>
inline Vec3<T> CVec3<T>::Re(){
	Vec3<T> v(real(cvec[0]),real(cvec[1]),real(cvec[2]));
	return v;
}

template<typename T>
inline Vec3<T> CVec3<T>::Im(){
	Vec3<T> v(imag(cvec[0]),imag(cvec[1]),imag(cvec[2]));
	return v;
}

template<typename T>
inline Vec3<T> CVec3<T>::Abs(){
	Vec3<T> v(abs(cvec[0]),abs(cvec[1]),abs(cvec[2]));
	return v;
}

/*	***Mat33***	*/

// Constructors and destructors

/// Default constructor (identity matrix)
template<typename T>
inline Mat33<T>::Mat33()
{
	mat[0][0] = 1.0; mat[0][1] = 0.0; mat[0][2] = 0.0;
	mat[1][0] = 0.0; mat[1][1] = 1.0; mat[1][2] = 0.0;
	mat[2][0] = 0.0; mat[2][1] = 0.0; mat[2][2] = 1.0;
}

/// Copy constructor
template<typename T>
inline Mat33<T>::Mat33(const Mat33<T> &A)
{
	for (int i = 0; i < 3; i++)
	{
		mat[i][0] = A.mat[i][0];
		mat[i][1] = A.mat[i][1];
		mat[i][2] = A.mat[i][2];
	}
}

/// Tensor constructor
template<typename T>
inline Mat33<T>::Mat33(const Vec3<T> &u, const Vec3<T> &w)
{
	for(int i = 0; i < 3; i++){
		mat[i][0] = u.vec[i]*w.vec[0];
		mat[i][1] = u.vec[i]*w.vec[1];
		mat[i][2] = u.vec[i]*w.vec[2];
	}
}

/// Constant constructor
template<typename T>
inline Mat33<T>::Mat33(const T a)
{
	mat[0][0] = a;
	mat[0][1] = a;
	mat[0][2] = a;
	
	mat[1][0] = a;
	mat[1][1] = a;
	mat[1][2] = a;
	
	mat[2][0] = a;
	mat[2][1] = a;
	mat[2][2] = a;
}

/// Diagonal constructor
template<typename T>
inline Mat33<T>::Mat33(const T a, const T b, const T c)
{
	mat[0][0] = a;
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	
	mat[1][0] = 0.0;
	mat[1][1] = b;
	mat[1][2] = 0.0;
	
	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = c;
}

/// Constant constructor from array
template<typename T>
inline Mat33<T>::Mat33(const T* a)
{
	mat[0][0] = a[0];
	mat[0][1] = a[1];
	mat[0][2] = a[2];
	
	mat[1][0] = a[3];
	mat[1][1] = a[4];
	mat[1][2] = a[5];
	
	mat[2][0] = a[6];
	mat[2][1] = a[7];
	mat[2][2] = a[8];
}

//Utility functions

/// Assignment operator operator B = A
template<typename T>
template<typename U>
inline Mat33<T>& Mat33<T>::operator= (const Mat33<U> &A)
{
	for (int i = 0; i < 3; i++){
		mat[i][0] = (T)A.mat[i][0];
		mat[i][1] = (T)A.mat[i][1];
		mat[i][2] = (T)A.mat[i][2];
	}
	return *this;
}

/// Identity matrix
template<typename T>
inline void Mat33<T>::M3Eye()
{
	mat[0][0] = 1.0; mat[0][1] = 0.0; mat[0][2] = 0.0;
	mat[1][0] = 0.0; mat[1][1] = 1.0; mat[1][2] = 0.0;
	mat[2][0] = 0.0; mat[2][1] = 0.0; mat[2][2] = 1.0;
}

/// Zero matrix
template<typename T>
inline void Mat33<T>::M3Zeros()
{
	mat[0][0] = 0.0; mat[0][1] = 0.0; mat[0][2] = 0.0;
	mat[1][0] = 0.0; mat[1][1] = 0.0; mat[1][2] = 0.0;
	mat[2][0] = 0.0; mat[2][1] = 0.0; mat[2][2] = 0.0;
}


/// Comparison operators
template<typename T>
inline bool Mat33<T>::operator== (const Mat33<T> &A)
{
	for(int i = 0; i < 3; i++){
		if(fabs(A.mat[i][0]-mat[i][0])>EPS) return false;
		if(fabs(A.mat[i][1]-mat[i][1])>EPS) return false;
		if(fabs(A.mat[i][2]-mat[i][2])>EPS) return false;
	}
	return true;
}

template<typename T>
inline bool Mat33<T>::operator!= (const Mat33<T> &A)
{
	for(int i = 0; i < 3; i++){
		if(fabs(A.mat[i][0]-mat[i][0])>EPS) return true;
		if(fabs(A.mat[i][1]-mat[i][1])>EPS) return true;
		if(fabs(A.mat[i][2]-mat[i][2])>EPS) return true;
	}
	return false;
}


//Matrix-scalar operations

/// Matrix-scalar addition B = A + a
template<typename T>
inline Mat33<T> Mat33<T>::operator+ (const T a)
{
	Mat33<T> B;
	B.mat[0][0] = mat[0][0]+a;
	B.mat[0][1] = mat[0][1]+a;
	B.mat[0][2] = mat[0][2]+a;
	B.mat[1][0] = mat[1][0]+a;
	B.mat[1][1] = mat[1][1]+a;
	B.mat[1][2] = mat[1][2]+a;
	B.mat[2][0] = mat[2][0]+a;
	B.mat[2][1] = mat[2][1]+a;
	B.mat[2][2] = mat[2][2]+a;
	return B;
}

/// Matrix-scalar addition A = A + a
template<typename T>
inline Mat33<T>& Mat33<T>::operator+= (const T a)
{
	mat[0][0] += a;
	mat[0][1] += a;
	mat[0][2] += a;
	mat[1][0] += a;
	mat[1][1] += a;
	mat[1][2] += a;
	mat[2][0] += a;
	mat[2][1] += a;
	mat[2][2] += a;
	return *this;
}

/// Matrix-scalar subtraction B = A - a
template<typename T>
inline Mat33<T> Mat33<T>::operator- (const T a)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]-a;
	B.mat[0][1] = mat[0][1]-a;
	B.mat[0][2] = mat[0][2]-a;
	B.mat[1][0] = mat[1][0]-a;
	B.mat[1][1] = mat[1][1]-a;
	B.mat[1][2] = mat[1][2]-a;
	B.mat[2][0] = mat[2][0]-a;
	B.mat[2][1] = mat[2][1]-a;
	B.mat[2][2] = mat[2][2]-a;
	
	return B;
}

/// Matrix-scalar subtraction A = A - a
template<typename T>
inline Mat33<T>& Mat33<T>::operator-= (const T a)
{
	mat[0][0] -= a;
	mat[0][1] -= a;
	mat[0][2] -= a;
	mat[1][0] -= a;
	mat[1][1] -= a;
	mat[1][2] -= a;
	mat[2][0] -= a;
	mat[2][1] -= a;
	mat[2][2] -= a;
	return *this;
}

/// Matrix-vector subtraction from diagonal
template<typename T>
inline Mat33<T> Mat33<T>::operator- (const Vec3<T> &u)
{
	Mat33<T> B;
	B.mat[0][0] = mat[0][0]-u.vec[0];
	B.mat[1][1] = mat[1][1]-u.vec[1];
	B.mat[2][2] = mat[2][2]-u.vec[2];
	return B;
}

/// Matrix-vector subtraction from diagonal
template<typename T>
inline Mat33<T>& Mat33<T>::operator-= (const Vec3<T> &u)
{
	mat[0][0] -= u.vec[0];
	mat[1][1] -= u.vec[1];
	mat[2][2] -= u.vec[2];
	return *this;
}

/// Matrix-scalar multiplication B = A*a
template<typename T>
inline Mat33<T> Mat33<T>::operator* (const T a)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*a;
	B.mat[0][1] = mat[0][1]*a;
	B.mat[0][2] = mat[0][2]*a;
	
	B.mat[1][0] = mat[1][0]*a;
	B.mat[1][1] = mat[1][1]*a;
	B.mat[1][2] = mat[1][2]*a;
	
	B.mat[2][0] = mat[2][0]*a;
	B.mat[2][1] = mat[2][1]*a;
	B.mat[2][2] = mat[2][2]*a;
	
	return B;
}

/// Matrix-scalar multiplication A = A*a
template<typename T>
inline Mat33<T>& Mat33<T>::operator*= (const T a)
{
	mat[0][0] *= a;
	mat[0][1] *= a;
	mat[0][2] *= a;
	mat[1][0] *= a;
	mat[1][1] *= a;
	mat[1][2] *= a;
	mat[2][0] *= a;
	mat[2][1] *= a;
	mat[2][2] *= a;
	return *this;
}

/// Matrix-scalar division B = A/a
template<typename T>
inline Mat33<T> Mat33<T>::operator/ (const T a)
{
	Mat33<T> B(*this);
	T inva = 1.0/a;
	B.mat[0][0] = mat[0][0]*inva;
	B.mat[0][1] = mat[0][1]*inva;
	B.mat[0][2] = mat[0][2]*inva;
	
	B.mat[1][0] = mat[1][0]*inva;
	B.mat[1][1] = mat[1][1]*inva;
	B.mat[1][2] = mat[1][2]*inva;
	
	B.mat[2][0] = mat[2][0]*inva;
	B.mat[2][1] = mat[2][1]*inva;
	B.mat[2][2] = mat[2][2]*inva;
	return B;
}

/// Matrix-scalar division A = A/a
template<typename T>
inline Mat33<T>& Mat33<T>::operator/= (const T a)
{
	T inva = 1.0/a;
	mat[0][0] *= inva;
	mat[0][1] *= inva;
	mat[0][2] *= inva;
	mat[1][0] *= inva;
	mat[1][1] *= inva;
	mat[1][2] *= inva;
	mat[2][0] *= inva;
	mat[2][1] *= inva;
	mat[2][2] *= inva;
	return *this;
}


//Level 2 BLAS

/// Matrix-vector multiplication v = A*u
template<typename T>
template<typename U>
inline Vec3<T> Mat33<T>::operator* (const Vec3<U> &u)
{
	return Vec3<T>(mat[0][0]*u.vec[0]+mat[0][1]*u.vec[1]+mat[0][2]*u.vec[2],
	               mat[1][0]*u.vec[0]+mat[1][1]*u.vec[1]+mat[1][2]*u.vec[2],
	               mat[2][0]*u.vec[0]+mat[2][1]*u.vec[1]+mat[2][2]*u.vec[2]);
}


//Level 3 BLAS

/// Matrix addition B = A + C
template<typename T>
inline Mat33<T> Mat33<T>::operator+ (const Mat33<T> &C)
{
	Mat33<T> B(*this);
	for(int i = 0; i < 3; i++){
		B.mat[i][0] += C.mat[i][0];
		B.mat[i][1] += C.mat[i][1];
		B.mat[i][2] += C.mat[i][2];
	}
	return B;
}

/// Matrix addition C = A + C
template<typename T>
inline Mat33<T>& Mat33<T>::operator+= (const Mat33<T> &C)
{
	for(int i = 0; i < 3; i++){
		mat[i][0] += C.mat[i][0];
		mat[i][1] += C.mat[i][1];
		mat[i][2] += C.mat[i][2];
	}
	return *this;
}

/// Matrix subtraction B = A - C
template<typename T>
inline Mat33<T> Mat33<T>::operator- (const Mat33<T> &C)
{
	Mat33<T> B(*this);
	for(int i = 0; i < 3; i++){
		B.mat[i][0] -= C.mat[i][0];
		B.mat[i][1] -= C.mat[i][1];
		B.mat[i][2] -= C.mat[i][2];
	}
	return B;
}

/// Matrix subtraction A = A - C
template<typename T>
inline Mat33<T>& Mat33<T>::operator -= (const Mat33<T> &C)
{
	for(int i = 0; i < 3; i++){
		mat[i][0] -= C.mat[i][0];
		mat[i][1] -= C.mat[i][1];
		mat[i][2] -= C.mat[i][2];
	}
	return *this;
}

/// Matrix multiplication B = A*C
template<typename T>
template<typename U>
inline Mat33<T> Mat33<T>::operator* (const Mat33<U> &C)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[1][0] + mat[0][2]*C.mat[2][0];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[2][1];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[1][0]*C.mat[0][0] + mat[1][1]*C.mat[1][0] + mat[1][2]*C.mat[2][0];
	B.mat[1][1] = mat[1][0]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[2][1];
	B.mat[1][2] = mat[1][0]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[2][0]*C.mat[0][0] + mat[2][1]*C.mat[1][0] + mat[2][2]*C.mat[2][0];
	B.mat[2][1] = mat[2][0]*C.mat[0][1] + mat[2][1]*C.mat[1][1] + mat[2][2]*C.mat[2][1];
	B.mat[2][2] = mat[2][0]*C.mat[0][2] + mat[2][1]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

/// Matrix multiplication A = A*C
template<typename T>
template<typename U>
inline Mat33<T>& Mat33<T>::operator*= (const Mat33<U> &C)
{
	Mat33<T> A(*this);
	
	mat[0][0] = A.mat[0][0]*C.mat[0][0] + A.mat[0][1]*C.mat[1][0] + A.mat[0][2]*C.mat[2][0];
	mat[0][1] = A.mat[0][0]*C.mat[0][1] + A.mat[0][1]*C.mat[1][1] + A.mat[0][2]*C.mat[2][1];
	mat[0][2] = A.mat[0][0]*C.mat[0][2] + A.mat[0][1]*C.mat[1][2] + A.mat[0][2]*C.mat[2][2];
	mat[1][0] = A.mat[1][0]*C.mat[0][0] + A.mat[1][1]*C.mat[1][0] + A.mat[1][2]*C.mat[2][0];
	mat[1][1] = A.mat[1][0]*C.mat[0][1] + A.mat[1][1]*C.mat[1][1] + A.mat[1][2]*C.mat[2][1];
	mat[1][2] = A.mat[1][0]*C.mat[0][2] + A.mat[1][1]*C.mat[1][2] + A.mat[1][2]*C.mat[2][2];
	mat[2][0] = A.mat[2][0]*C.mat[0][0] + A.mat[2][1]*C.mat[1][0] + A.mat[2][2]*C.mat[2][0];
	mat[2][1] = A.mat[2][0]*C.mat[0][1] + A.mat[2][1]*C.mat[1][1] + A.mat[2][2]*C.mat[2][1];
	mat[2][2] = A.mat[2][0]*C.mat[0][2] + A.mat[2][1]*C.mat[1][2] + A.mat[2][2]*C.mat[2][2];
	
	return *this;
}

//Special functions

template<typename T>
inline Mat33<T> Mat33<T>::TransMulM3(const Mat33<T> &C)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[1][0]*C.mat[1][0] + mat[2][0]*C.mat[2][0];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[1][0]*C.mat[1][1] + mat[2][0]*C.mat[2][1];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[1][0]*C.mat[1][2] + mat[2][0]*C.mat[2][2];
	
	B.mat[1][0] = mat[0][1]*C.mat[0][0] + mat[1][1]*C.mat[1][0] + mat[2][1]*C.mat[2][0];
	B.mat[1][1] = mat[0][1]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[2][1]*C.mat[2][1];
	B.mat[1][2] = mat[0][1]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[2][1]*C.mat[2][2];
	
	B.mat[2][0] = mat[0][2]*C.mat[0][0] + mat[1][2]*C.mat[1][0] + mat[2][2]*C.mat[2][0];
	B.mat[2][1] = mat[0][2]*C.mat[0][1] + mat[1][2]*C.mat[1][1] + mat[2][2]*C.mat[2][1];
	B.mat[2][2] = mat[0][2]*C.mat[0][2] + mat[1][2]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

/// Matrix-vector multiplication v = A*u
template<typename T>
inline Vec3<T> Mat33<T>::TransMulVec(const Vec3<T> &u)
{
	return Vec3<T>(mat[0][0]*u.vec[0]+mat[1][0]*u.vec[1]+mat[2][0]*u.vec[2],
	               mat[0][1]*u.vec[0]+mat[1][1]*u.vec[1]+mat[2][1]*u.vec[2],
	               mat[0][2]*u.vec[0]+mat[1][2]*u.vec[1]+mat[2][2]*u.vec[2]);
}

template<typename T>
inline Vec3<T> Mat33<T>::ColumnVec3(const int col)
{
	Vec3<T> v(mat[0][col],mat[1][col],mat[2][col]);
	return v;
}

/// Matrix determinant d = Det(A)
template<typename T>
inline T Mat33<T>::M3Det()
{
	T det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
			+ mat[1][0]*(mat[2][1]*mat[0][2]-mat[0][1]*mat[2][2])
			+ mat[2][0]*(mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1]);
	return det;
}

// Multiplication of non-symmetric matrix with symmetric matrix C
template<typename T>
inline Mat33<T> Mat33<T>::MulSymM3(const Mat33<T> &C)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[0][1] + mat[0][2]*C.mat[0][2];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[1][2];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[1][0]*C.mat[0][0] + mat[1][1]*C.mat[0][1] + mat[1][2]*C.mat[0][2];
	B.mat[1][1] = mat[1][0]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[1][2];
	B.mat[1][2] = mat[1][0]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[2][0]*C.mat[0][0] + mat[2][1]*C.mat[0][1] + mat[2][2]*C.mat[0][2];
	B.mat[2][1] = mat[2][0]*C.mat[0][1] + mat[2][1]*C.mat[1][1] + mat[2][2]*C.mat[1][2];
	B.mat[2][2] = mat[2][0]*C.mat[0][2] + mat[2][1]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

// Multiplication of symmetric matrix with non-symmetric matrix C
template<typename T>
inline Mat33<T> Mat33<T>::SymMulM3(const Mat33<T> &C)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[1][0] + mat[0][2]*C.mat[2][0];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[2][1];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[0][1]*C.mat[0][0] + mat[1][1]*C.mat[1][0] + mat[1][2]*C.mat[2][0];
	B.mat[1][1] = mat[0][1]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[2][1];
	B.mat[1][2] = mat[0][1]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[0][2]*C.mat[0][0] + mat[1][2]*C.mat[1][0] + mat[2][2]*C.mat[2][0];
	B.mat[2][1] = mat[0][2]*C.mat[0][1] + mat[1][2]*C.mat[1][1] + mat[2][2]*C.mat[2][1];
	B.mat[2][2] = mat[0][2]*C.mat[0][2] + mat[1][2]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

// Multiplication of two symmetric matrices
template<typename T>
inline Mat33<T> Mat33<T>::SymMulSymM3(const Mat33<T> &C)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[0][1] + mat[0][2]*C.mat[0][2];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[1][2];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[0][1]*C.mat[0][0] + mat[1][1]*C.mat[0][1] + mat[1][2]*C.mat[0][2];
	B.mat[1][1] = mat[0][1]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[1][2];
	B.mat[1][2] = mat[0][1]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[0][2]*C.mat[0][0] + mat[1][2]*C.mat[0][1] + mat[2][2]*C.mat[0][2];
	B.mat[2][1] = mat[0][2]*C.mat[0][1] + mat[1][2]*C.mat[1][1] + mat[2][2]*C.mat[1][2];
	B.mat[2][2] = mat[0][2]*C.mat[0][2] + mat[1][2]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

/*!
 * Cramer's rule inversion for 3x3 matrix
 */
template<typename T>
inline Mat33<T> Mat33<T>::M3Inv(bool &sing)
{
	Mat33<T> A(*this);
	
	A.mat[0][0] = mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1];
	A.mat[0][1] = mat[2][1]*mat[0][2]-mat[0][1]*mat[2][2];
	A.mat[0][2] = mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1];
	
	T det = mat[0][0]*A.mat[0][0] + mat[1][0]*A.mat[0][1] + mat[2][0]*A.mat[0][2];
	
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		sing=true;
	}
	
	T invdet = 1.0/det;
	
	A.mat[0][0] *= invdet;
	A.mat[0][1] *= invdet;
	A.mat[0][2] *= invdet;
	
	A.mat[1][0] = (mat[2][0]*mat[1][2]-mat[1][0]*mat[2][2])*invdet;
	A.mat[1][1] = (mat[0][0]*mat[2][2]-mat[2][0]*mat[0][2])*invdet;
	A.mat[1][2] = (mat[1][0]*mat[0][2]-mat[0][0]*mat[1][2])*invdet;
	
	A.mat[2][0] = (mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0])*invdet;
	A.mat[2][1] = (mat[2][0]*mat[0][1]-mat[0][0]*mat[2][1])*invdet;
	A.mat[2][2] = (mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1])*invdet;
	
	return A;
}

/*!
 * Cramer's rule inversion for *symmetric* 3x3 matrix (expects upper half to be filled, WARNING: no checking if matrix is actually symmetric)
 */
template<typename T>
inline Mat33<T> Mat33<T>::SymM3Inv(bool &sing)
{
	Mat33<T> A;
	
	A.mat[0][0] = mat[1][1]*mat[2][2]-mat[1][2]*mat[1][2];
	A.mat[0][1] = mat[1][2]*mat[0][2]-mat[0][1]*mat[2][2];
	A.mat[0][2] = mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1];
	T det = mat[0][0]*A.mat[0][0] + mat[0][1]*A.mat[0][1] + mat[0][2]*A.mat[0][2];
	
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		sing=true;
	}
	
	T invdet = 1.0/det;
	
	A.mat[0][0] *= invdet;
	A.mat[0][1] *= invdet;
	A.mat[0][2] *= invdet;
	
	A.mat[1][0] = A.mat[0][1];
	A.mat[1][1] = (mat[0][0]*mat[2][2]-mat[0][2]*mat[0][2])*invdet;
	A.mat[1][2] = (mat[0][1]*mat[0][2]-mat[0][0]*mat[1][2])*invdet;
	
	A.mat[2][0] = A.mat[0][2];
	A.mat[2][1] = A.mat[1][2];
	A.mat[2][2] = (mat[0][0]*mat[1][1]-mat[0][1]*mat[0][1])*invdet;
	
	return A;
}

/*!
 * Cramer's rule inversion for *symmetric* 3x3 matrix (expects upper half to be filled, WARNING: no checking if matrix is actually symmetric)
 * followed by multiplication with vector
 */
template<typename T>
inline Vec3<T> Mat33<T>::SymM3InvMult(bool &sing, Vec3<T> &z)
{
	Vec3<T> v, r;
	
	v.vec[0]=mat[1][1]*mat[2][2]-mat[1][2]*mat[1][2];
	v.vec[1]=mat[1][2]*mat[0][2]-mat[0][1]*mat[2][2];
	v.vec[2]=mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1];
	
	T det = mat[0][0]*v.vec[0] + mat[0][1]*v.vec[1] + mat[0][2]*v.vec[2];
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		sing=true;
	}
	
	r.vec[0] = v.vec[0]*z.vec[0]+v.vec[1]*z.vec[1]+v.vec[2]*z.vec[2];
	v.vec[0] = (mat[0][1]*mat[0][2]-mat[0][0]*mat[1][2]); // a_12
	r.vec[1] = v.vec[1]*z.vec[0]+(mat[0][0]*mat[2][2]-mat[0][2]*mat[0][2])*z.vec[1]+v.vec[0]*z.vec[2];
	r.vec[2] = v.vec[2]*z.vec[0]+v.vec[0]*z.vec[1]+(mat[0][0]*mat[1][1]-mat[0][1]*mat[0][1])*z.vec[2];
	r/=det;
	
	return r;
}

/// Multiply matrix with diagonal matrix specified through a,b,c C = A*[a,0,0 ; 0,b,0 ; 0,0,c ]
template<typename T>
inline Mat33<T> Mat33<T>::M3MulDiag(const T a, const T b, const T c)
{
	Mat33<T> B;
	
	B.mat[0][0] = mat[0][0]*a;
	B.mat[0][1] = mat[0][1]*b;
	B.mat[0][2] = mat[0][2]*c;
	
	B.mat[1][0] = mat[1][0]*a;
	B.mat[1][1] = mat[1][1]*b;
	B.mat[1][2] = mat[1][2]*c;
	
	B.mat[2][0] = mat[2][0]*a;
	B.mat[2][1] = mat[2][1]*b;
	B.mat[2][2] = mat[2][2]*c;
	
	return B;
}

/// Matrix transpose C = A'
template<typename T>
inline Mat33<T> Mat33<T>::M3Transpose()
{
	Mat33<T> C;
	C.mat[0][0] = mat[0][0];
	C.mat[0][1] = mat[1][0];
	C.mat[0][2] = mat[2][0];
	C.mat[1][0] = mat[0][1];
	C.mat[1][1] = mat[1][1];
	C.mat[1][2] = mat[2][1];
	C.mat[2][0] = mat[0][2];
	C.mat[2][1] = mat[1][2];
	C.mat[2][2] = mat[2][2];
	return C;
}

/// Matrix row swap
template<typename T>
inline Mat33<T> Mat33<T>::M3RowSwap(const int a, const int b)
{
	Mat33<T> A;
	T temp;
	temp = mat[b][0];
	A.mat[b][0] = mat[a][0];
	A.mat[a][0] = temp;
	temp = mat[b][1];
	A.mat[b][1] = mat[a][1];
	A.mat[a][1] = temp;
	temp = mat[b][2];
	A.mat[b][2] = mat[a][2];
	A.mat[a][2] = temp;
	return A;
}

/// Extract diagonal elements of matrix v = Aii -- goes in Mat33 class
template<typename T>
inline Vec3<T> Mat33<T>::M3Diag()
{
	Vec3<T> v(mat[0][0],mat[1][1],mat[2][2]);
	return v;
}

/// Trace of a matrix a = Tr(A)
template<typename T>
inline T Mat33<T>::M3Trace()
{
	T a = mat[0][0]+mat[1][1]+mat[2][2];
	return a;
}

/// Back-substitution on 3x3 triangular system
template<typename T>
inline bool Mat33<T>::M3BackSub(Vec3<T> &x)
{
	bool multiples=false;
	if (fabs(mat[2][2])>EPS){
		x.vec[2] = x.vec[2]/mat[2][2];
	} else{
		x.vec[2]=1.0;
		multiples=true;
	}
	if (fabs(mat[1][1])>EPS){
		x.vec[1] = (x.vec[1] - mat[1][2]*x.vec[2])/mat[1][1];
	} else{
		x.vec[1]=1.0;
		multiples=true;
	}
	if (fabs(mat[0][0])>EPS){
		x.vec[0] = (x.vec[0] - mat[0][1]*x.vec[1] - mat[0][2]*x.vec[2])/mat[0][0];
	} else{
		x.vec[0]=1.0;
		multiples=true;
	}
#if DEBUG_LEVEL>2
	cout << M3Str() << "\n---\n";
#endif
	return multiples;
}

/// Linear Solver for 3x3 systems using GEPP + back-substitution
template<typename T>
inline Vec3<T> Mat33<T>::M3LinSolve(Vec3<T> &u, bool &multiples)
{
	Vec3<T> v(u);
	Mat33<T> A(*this);
	A.M3GEPP(v);
	multiples=A.M3BackSub(v);
	return v;
}

template<typename T>
inline Mat33<T> Vec2Rot(Vec3<T> &rot_vec, const unsigned int a0) // a0 is reference axis {0,1,2} = {x,y,z}
{
	Mat33<T> a,b;
	rot_vec/=rot_vec.V3Norm(); // normalize (just in case)
	unsigned int a1=(a0+1)%3;
	unsigned int a2=(a0+2)%3;
	T sqr_argument=rot_vec.vec[a1]*rot_vec.vec[a1]+rot_vec.vec[a2]*rot_vec.vec[a2];
	if(sqr_argument>EPS){
		T invhyp=1.0/sqrt(sqr_argument);
		T cp=rot_vec.vec[a2]*invhyp;
		T sp=rot_vec.vec[a1]*invhyp;
		T st=sqrt(1.0-rot_vec.vec[a0]*rot_vec.vec[a0]);
		
		// we're doing axis rotation here, that's why the sign of the sin term sp is opposite
		a.mat[a0][a0]=1.0;	a.mat[a0][a1]=0.0;	a.mat[a0][a2]=0.0;
		a.mat[a1][a0]=0.0;	a.mat[a1][a1]=cp;	a.mat[a1][a2]=sp;
		a.mat[a2][a0]=0.0;	a.mat[a2][a1]=-1.0*sp;	a.mat[a2][a2]=cp;
		
		b.mat[a0][a0]=rot_vec.vec[a0];	b.mat[a0][a1]=0.0;	b.mat[a0][a2]=-1.0*st;
		b.mat[a1][a0]=0.0;		b.mat[a1][a1]=1.0;	b.mat[a1][a2]=0.0;
		b.mat[a2][a0]=st;		b.mat[a2][a1]=0.0;	b.mat[a2][a2]=rot_vec.vec[a0];
		
		return a*b; // return rotation matrix if vector not pointing in z-direction
	}
	// if we get here then the vector is in axis-direction
	if(rot_vec.vec[a0]<0.0){ // rotate 180 degrees around next axis
		a.mat[a0][a0]=-1.0;
		a.mat[a2][a2]=-1.0;
	}
	return a; // return unit matrix if no rotation necessary (also applies to zero-vector)
}

/*!
 * Gaussian elimination with partial pivoting for solving Ax = b
 * Used within M3LinSolve, but not inlined due to size
 */
template<typename T>
inline void Mat33<T>::M3GEPP(Vec3<T> &x)
{
	// Perform elimination
	for(unsigned int i=0; i<3; i++){
		// Find largest pivot (to be on diagonal - AT), with index pidx
		unsigned int pidx = i;
		for(unsigned int j=i+1; j<3; j++){
			if(fabs(mat[j][i]) > fabs(mat[pidx][i])) pidx = j; // was compared to mat[pidx][j], but no need to care for the diagonal down the road ...
		}
		
		// Swap rows -- replace when figure out how to use class function within each other LEJ 01/07/10
		if(pidx != i){
			for(unsigned int j=0; j<3; j++){
				T tempij = mat[pidx][j];
				mat[pidx][j] = mat[i][j];
				mat[i][j] = tempij;
			}
			T tempij =  x.vec[pidx];
			x.vec[pidx] = x.vec[i];
			x.vec[i] = tempij;
		}
		// Eliminate
		for(unsigned k=i+1; k<3; k++){
			x.vec[k] -= mat[k][i]/mat[i][i]*x.vec[i];
			for(int j=2; j>(int)i; j--) mat[k][(unsigned int)j] -= (mat[k][i]/mat[i][i])*mat[i][(unsigned int)j];
			mat[k][i]=0.0;
		}
	}
}

/// LU decomposition routine. WARNING: No pivoting at the moment -- AT
template<typename T>
inline void Mat33<T>::LUDecomposition()
{
	int i,j,k;
	for(i=0; i<3; i++){
		for(j=i; j<3; j++){
			for(k=0; k<i-1; k++){
				mat[i][j] -= mat[i][k]*mat[k][j];
			}
		}
		for(j=i+1; j<3; j++){
			for(k=0; k<i-1; k++){
				mat[j][i] -= mat[j][k]*mat[k][i];
			}
#if DEBUG_LEVEL>2
			cout << "a[" << j << "][" << i << "] = " << mat[j][i] << ", a[" << i << "][" << i << "] = " << mat[i][i]<< "\n";
#endif
			mat[j][i] /= mat[i][i];
		}
	}
}

/*!
 * Determine eigenvalues of 3x3 matrix -- AT
 * characteristic equation: det(A-lambda*E)=0 (lambda ... eigenvalues; E=3x3 unit diagonal matrix)
 *
 * |	a11-lambda	a12	a13	|
 * |	a21	a22-lambda	a23	| = 0 = lambda^3 + alpha*lambda^2 + beta*lambda + gamma
 * |	a31	a32	a33-lambda	|
 *
 * alpha = -a11 - a22 - a33 ;
 * beta = a11*(a22+a33) - alphax - a13*a31 - a12*a21
 * gamma = a11*alphax - a12*(a23*a31-a21*a33) - a13*(a21*a32-a31*a22)
 * alphax = a23*a32 - a22*a33
 *
 * reduced form after substition with lambda = z - alpha/3
 * z^3 + p*z + q = 0
 *
 * p = beta - alpha*alpha/3
 * q = (2*alpha^3 - 9*alpha*beta + 27*gamma)/27 = gamma + alpha/3*(2*(alpha/3)^2 - beta)
 */
template<typename T>
inline CVec3<T> Mat33<T>::Eigenvalues()
{
#if DEBUG_LEVEL>2
	cout << "Start Mat33::Eigenvalue\n";
#endif
	T alpha = -mat[0][0]-mat[1][1]-mat[2][2];
	T alphax = mat[1][2]*mat[2][1]-mat[1][1]*mat[2][2];
	T beta = mat[0][0]*(mat[1][1]+mat[2][2])-alphax-mat[0][2]*mat[2][0]-mat[0][1]*mat[1][0];
	T gamma = mat[0][0]*alphax-mat[0][1]*(mat[1][2]*mat[2][0]-mat[1][0]*mat[2][2])-mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
#if DEBUG_LEVEL>2
	cout << "alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << ", alphax = " << alphax << "\n";
#endif
	CVec3<T> cv;
	alphax=alpha/3.0; // alphax redefined
	cv=SolvePolynomial3(beta-alpha*alphax,gamma+alphax*((T)2.0*alphax*alphax-beta));
	cv-=alphax; // lambda = z - alpha/3
#if DEBUG_LEVEL>2
	cout << "Finished Mat33::Eigenvalues\n";
#endif
	return cv;
}

template<typename T>
inline Mat33<T> Mat33<T>::Eigenvectors(Vec3<T> &ew, bool* multiples, bool normalize)
{
	Mat33<T> result;
	for(unsigned int i=0; i<3; i++){
		Mat33<T> A(*this);
		A.mat[0][0]-=ew.vec[i]; A.mat[1][1]-=ew.vec[i]; A.mat[2][2]-=ew.vec[i];
		Vec3<T> zero(0.0);
		Vec3<T> v=A.M3LinSolve(zero,multiples[i]);
		if(normalize) v/=v.V3Norm();
		result.mat[0][i]=v.vec[0]; result.mat[1][i]=v.vec[1]; result.mat[2][i]=v.vec[2];
	}
	return result;
}

/*	***Vec4***	*/

/// Default constructor (zeros)
template<typename T>
inline Vec4<T>::Vec4()
{
	vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0; vec[3] = 0.0;
}

/// Construct from individual numbers
template<typename T>
template<typename U>
inline Vec4<T>::Vec4(const U a, const U b, const U c, const U d)
{
	vec[0] = (T)a;
	vec[1] = (T)b;
	vec[2] = (T)c;
	vec[3] = (T)d;
}

/// Construct from Vec3 and T
template<typename T>
inline Vec4<T>::Vec4(const Vec3<T> &u, const T d)
{
	vec[0] = u.vec[0];
	vec[1] = u.vec[1];
	vec[2] = u.vec[2];
	vec[3] = d;
}

/// Copy constructor
template<typename T>
inline Vec4<T>::Vec4(const Vec4<T> &u)
{
	vec[0] = u.vec[0];
	vec[1] = u.vec[1];
	vec[2] = u.vec[2];
	vec[3] = u.vec[3];
}

/// Vector addition v = u + w
template<typename T>
inline Vec4<T> Vec4<T>::operator+ (const Vec4<T> &w)
{
	Vec4<T> v(vec[0]+w.vec[0], vec[1]+w.vec[1], vec[2] +w.vec[2], vec[3] +w.vec[3]);
	return v;
}

/// Vector scalar multiplication v = u*a
template<typename T>
inline Vec4<T> Vec4<T>::operator* (const T a)
{
	Vec4<T> v(a*vec[0], a*vec[1], a*vec[2], a*vec[3]);
	return v;
}

/// Vector scalar division v = u/a
template<typename T>
inline Vec4<T> Vec4<T>::operator/ (const T a)
{
	T b = 1.0/a;
	Vec4<T> v(b*vec[0], b*vec[1], b*vec[2], b*vec[3]);
	return v;
}

/// Comparison operators
template<typename T>
inline bool Vec4<T>::operator==(const Vec4<T> &u)
{
	if(fabs(u.vec[0]-vec[0])>EPS) return false;
	if(fabs(u.vec[1]-vec[1])>EPS) return false;
	if(fabs(u.vec[2]-vec[2])>EPS) return false;
	if(fabs(u.vec[3]-vec[3])>EPS) return false;
	return true;
}

template<typename T>
inline bool Vec4<T>::operator!=(const Vec4<T> &u)
{
	if(fabs(u.vec[0]-vec[0])>EPS) return true;
	if(fabs(u.vec[1]-vec[1])>EPS) return true;
	if(fabs(u.vec[2]-vec[2])>EPS) return true;
	if(fabs(u.vec[3]-vec[3])>EPS) return true;
	return false;
}

/// Rotation axis and angle from rotation matrix
template<typename T>
inline Vec4<T> Rot2AxisAngle(Mat33<T> R)
{
	Vec4<T> AaA(1.0,0.0,0.0,0.0); // using axis=(1,0,0), theta=0 as failsafe b/c X3D animation would fail with (0,0,0,0) ...
	T cost=(R.M3Trace()-1.0)/2.0;
	if(cost>1.0) cost=1.0; // safety first
	if(cost<-1.0) cost=-1.0;
	T theta = acos(cost);
	T sintheta = sin(theta);
	T n_norm=1.0;
	if(fabs(sintheta)>EPS*EPS){
		n_norm = 1.0/(2.0*sintheta);
	}
	AaA.vec[0] = n_norm*(R.mat[2][1]-R.mat[1][2]);
	AaA.vec[1] = n_norm*(R.mat[0][2]-R.mat[2][0]);
	AaA.vec[2] = n_norm*(R.mat[1][0]-R.mat[0][1]);
	if(AaA.vec[0]*AaA.vec[0]+AaA.vec[1]*AaA.vec[1]+AaA.vec[2]*AaA.vec[2]<EPS) return Vec4<T>(1.0,0.0,0.0,0.0);
	AaA.vec[3] = theta;
	return AaA;
}


/// Rotation matrix from axis and angle
template<typename T>
inline Mat33<T> AxisAngle2Rot(Vec4<T> AxisAngle)
{
	// normalize axis part, just in case
	T invnorm=AxisAngle.vec[0]*AxisAngle.vec[0]+AxisAngle.vec[1]*AxisAngle.vec[1]+AxisAngle.vec[2]*AxisAngle.vec[2];
	Mat33<T> rot;
	if(invnorm>=EPS){ // safety first
		invnorm=1.0/sqrt(invnorm);
		AxisAngle.vec[0]*=invnorm; AxisAngle.vec[1]*=invnorm; AxisAngle.vec[2]*=invnorm;
		
		T x=AxisAngle.vec[0]; T y=AxisAngle.vec[1]; T z=AxisAngle.vec[2];
		T oneminuscos=1.0-cos(AxisAngle.vec[3]);
		T sine=sin(AxisAngle.vec[3]);
		
		rot.mat[0][0]+=oneminuscos*(x*x-1.0);	rot.mat[0][1]=-z*sine+oneminuscos*x*y;	rot.mat[0][2]=y*sine+oneminuscos*x*z;
		rot.mat[1][0]=z*sine+oneminuscos*x*y;	rot.mat[1][1]+=oneminuscos*(y*y-1.0);	rot.mat[1][2]=-x*sine+oneminuscos*y*z;
		rot.mat[2][0]=-y*sine+oneminuscos*x*z;	rot.mat[2][1]=x*sine+oneminuscos*y*z;	rot.mat[2][2]+=oneminuscos*(z*z-1.0);
	}
	return rot;
}

/// Rotation matrix from axis and angle
template<typename T>
inline Mat33<T> AxisAngle2Rot(Vec3<T> &axis, T &angle)
{
	// normalize axis part, just in case
	T invnorm=axis.vec[0]*axis.vec[0]+axis.vec[1]*axis.vec[1]+axis.vec[2]*axis.vec[2];
	Mat33<T> rot;
	if(invnorm>=EPS){ // safety first
		invnorm=1.0/sqrt(invnorm);
		axis.vec[0]*=invnorm; axis.vec[1]*=invnorm; axis.vec[2]*=invnorm;
		
		T x=axis.vec[0]; T y=axis.vec[1]; T z=axis.vec[2];
		T oneminuscos=1.0-cos(angle);
		T sine=sin(angle);
		
		rot.mat[0][0]+=oneminuscos*(x*x-1.0);	rot.mat[0][1]=-z*sine+oneminuscos*x*y;	rot.mat[0][2]=y*sine+oneminuscos*x*z;
		rot.mat[1][0]=z*sine+oneminuscos*x*y;	rot.mat[1][1]+=oneminuscos*(y*y-1.0);	rot.mat[1][2]=-x*sine+oneminuscos*y*z;
		rot.mat[2][0]=-y*sine+oneminuscos*x*z;	rot.mat[2][1]=x*sine+oneminuscos*y*z;	rot.mat[2][2]+=oneminuscos*(z*z-1.0);
	}
	return rot;
}

/// Rotation matrix from axis and angle
template<typename T>
inline Mat33<T> AC2Rot(Vec3<T> axis, T cos_angle)
{
	// normalize axis part, just in case
	T x=axis.vec[0]; T y=axis.vec[1]; T z=axis.vec[2];
	T invnorm=x*x+y*y+z*z;
	Mat33<T> rot;
	if(invnorm>=EPS){ // safety first
		invnorm=1.0/sqrt(invnorm);
		x*=invnorm; y*=invnorm; z*=invnorm;
		
		T oneminuscos = (cos_angle > 1.0) ? 0.0 : ((cos_angle < -1.0) ? 2.0 : 1.0-cos_angle);
		T sine = (cos_angle > 1.0) ? 0.0 : ((cos_angle < -1.0) ? 0.0 : sqrt(1.0-cos_angle*cos_angle));
		
		T zs = z*sine+oneminuscos*x*y;
		T ys = y*sine+oneminuscos*x*z;
		T xs = x*sine+oneminuscos*y*z;
		sine *= 2.0;
		
		rot.mat[0][0] += oneminuscos*(x*x-1.0);	rot.mat[0][1] = zs - z*sine;		rot.mat[0][2] = ys;
		rot.mat[1][0] = zs;			rot.mat[1][1] += oneminuscos*(y*y-1.0);	rot.mat[1][2] = xs - x*sine;
		rot.mat[2][0] = ys - y*sine;		rot.mat[2][1] = xs;			rot.mat[2][2] += oneminuscos*(z*z-1.0);
		// 12 flops diagonal + 18 flops (vs. 12+15=27) => total 12+18=30 flops (vs. 39)
	}
	return rot;
}


template<typename T>
inline bool Point_in_Ellipsoid(Vec3<T> &point, Vec3<T> &saxes, Vec3<T> &center, Mat33<T> &rot, unsigned int sphere_surface_pm = 0, bool* in_surface = NULL)
{
	Vec3<T> dist=(point-center);
	T d2;
	if((fabs(saxes.vec[0]-saxes.vec[1])<EPS) && (fabs(saxes.vec[2]-saxes.vec[1])<EPS)){ // sphere here
		d2 = dist*dist;
		T s2 = saxes.vec[0]*saxes.vec[0];
		bool pie = (d2<=s2);
		if(sphere_surface_pm){
			s2  = saxes.vec[0] - sphere_surface_pm*0.01; // 1 pm = 0.01 Angstrom
			s2 *= s2;
			*in_surface = pie && (d2>=s2);
		}
		return pie;
	}
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	dist=(rot.M3Transpose()*dist);
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	d2 = dist*dist;
	
	return (d2<=1.0); // point is inside if (x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1)
}

#define atan_a -0.012299380859105
#define atan_b 0.054082655552459
#define atan_c -0.11769677376706
#define atan_d 0.19402227554937
#define atan_e -0.33269718723178
#define atan_f 0.99998657415361

template<typename T>
inline T fastatan(T x)
{
	T arg=fabs(x);
	T arg2=x*x;
	if(arg<=1.0){
		return copysign((((((atan_a*arg2+atan_b)*arg2+atan_c)*arg2+atan_d)*arg2+atan_e)*arg2+atan_f)*arg,x);
	} else{
		arg=1.0/arg;
		arg2=arg*arg;
		return copysign(PI/2-(((((atan_a*arg2+atan_b)*arg2+atan_c)*arg2+atan_d)*arg2+atan_e)*arg2+atan_f)*arg,x);
	}
}

template<typename T>
inline T fastatan2(T y, T x)
{
	if(x>0.0) return fastatan(y/x);
	if(y>=0.0) return fastatan(y/x)+PI;
	return fastatan(y/x)-PI;
}

template<typename T>
inline void UnitVec2ThetaPhi(Vec3<T> &v, T &theta, T &phi)
{
	if(v.vec[2]>1.0) v.vec[2]=1.0; // safety first
	if(v.vec[2]<-1.0) v.vec[2]=-1.0;
	theta=acos(v.vec[2]);
	phi=atan(v.vec[1]/v.vec[0]);
}

template<typename T>
inline void Vec2ThetaPhi(Vec3<T> &v, T &theta, T &phi)
{
	Vec3<T> r=v/v.V3Norm();
	if(r.vec[2]>1.0) r.vec[2]=1.0; // safety first
	if(r.vec[2]<-1.0) r.vec[2]=-1.0;
	theta=acos(r.vec[2]);
	phi=atan(r.vec[1]/r.vec[0]); // atan is 1/0 safe ;-)
}

template<typename T>
inline void Vec2ThetaPhi(Vec3<T> &v, T &theta, T &phi, T &r)
{
	r=v.V3Norm();
	Vec3<T> rvec=v/r;
	if(rvec.vec[2]>1.0) rvec.vec[2]=1.0; // safety first
	if(rvec.vec[2]<-1.0) rvec.vec[2]=-1.0;
	theta=acos(rvec.vec[2]);
	phi=atan(rvec.vec[1]/rvec.vec[0]);
}

template<typename T>
inline T VecDist2ThetaPhi(Vec3<T> &v, T r, T &theta, T &phi)
{
	Vec3<T> rvec=v;
	rvec.vec[2]/=r;
	if(rvec.vec[2]>1.0) rvec.vec[2]=1.0; // safety first
	if(rvec.vec[2]<-1.0) rvec.vec[2]=-1.0;
	theta=acos(rvec.vec[2]);
	phi=atan(rvec.vec[1]/rvec.vec[0]);
	return rvec.vec[2];
}

template<typename T>
inline T VecDist2Phi(Vec3<T> &v, T r, T &phi)
{
	T cost=v.vec[2]/r;
	if(cost>1.0) cost=1.0; // safety first
	if(cost<-1.0) cost=-1.0;
	phi=fastatan2(v.vec[1],v.vec[0]);
	return cost;
}

template<typename T>
inline T EllipsoidRmin(T &theta, T &phi, Vec3<T> &saxes)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	T sint=sin(theta);
	Vec3<T> dist(sint*cos(phi),sint*sin(phi),cos(theta));
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dist*dist);
}

template<typename T>
inline T EllipsoidRmin(T &theta, T &phi, Vec3<T> &saxes, Mat33<T> &rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	T sint=sin(theta);
	Vec3<T> direction(sint*cos(phi),sint*sin(phi),cos(theta));
	Vec3<T> dist=rot.M3Transpose()*direction;
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dist*dist);
}

template<typename T>
inline T EllipsoidRmin(Vec3<T> &direction, Vec3<T> &saxes, Mat33<T> rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	Vec3<T> dist=rot.M3Transpose()*direction/direction.V3Norm();
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dist*dist);
}

// Based on P.P. Klein, "On the Ellipsoid and Plane Intersection Equation", Applied Mathematics 3, 1634 (2012)
// see page 1639
template<typename T>
inline T Ellipsoid_Cross_Section(Vec3<T> &invsaxes2, Vec3<T> direction)
{
	T a=direction*direction;
	if(a>EPS*EPS){
		a=1.0/a;
		T c=direction.vec[0]*direction.vec[0]*(invsaxes2.vec[1]*invsaxes2.vec[2])+direction.vec[1]*direction.vec[1]*(invsaxes2.vec[0]*invsaxes2.vec[2])+direction.vec[2]*direction.vec[2]*(invsaxes2.vec[0]*invsaxes2.vec[1]);
		// a*beta^2 + b*beta + c = 0
		// beta_+/- = -b/2a +/- sqrt(b^2 - 4*a*c)/2a
		// beta+ * beta_ = (-b + sqrt(b^2 - 4*a*c)) * (-b - sqrt(b^2 - 4*a*c))/(2a)^2 = (b^2 - (b^2 - 4ac))/4a^2 = 4ac/4a^2 = c/a
		return PI/sqrt(c*a);
	} else return 0.0;
}

template<typename T>
inline T VectorAngle(Vec3<T> &a, Vec3<T> &b)
{
	T theta;
	T cosine=(a*b)/sqrt((a*a)*(b*b));
	if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
		theta=0.0;
	} else{
		if(cosine<-1.0){
			theta=PI;
		} else theta=acos(cosine);
	}
	return theta;
}

template<typename T>
inline T VectorCos(Vec3<T> &a, Vec3<T> &b)
{
	T cosine=(a*b)/sqrt((a*a)*(b*b));
	if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
		cosine=1.0;
	} else{
		if(cosine<-1.0) cosine=-1.0;
	}
	return cosine;
}

/// Obtain rotation matrix mapping vector a to b direction
template<typename T>
inline Mat33<T> RotAtoB(Vec3<T> &a, Vec3<T> &b)
{
	Vec3<T> axis=a;
	axis.V3Cross(b);
	T cos_angle=VectorCos(a,b);
	if(axis*axis<EPS){ // vectors are colinear, two solutions here: everything stays as is (angle=0) or rotate PI around vector perpendicular to a (and b)
		if(cos_angle-EPS < -1.0){
			// find non zero component of a
			if(fabs(a.vec[0])>EPS){
				axis.vec[1]=1.0;
				axis.vec[2]=1.0;
				// solve for axis*a = 0
				// => a_x*axis_x+a_y*axis_y+a_z*axis_z
				// => axis_x = -(a_y+a_z)/a_x
				axis.vec[0]=-(a.vec[1]+a.vec[2])/a.vec[0];
			} else{
				if(fabs(a.vec[1])>EPS){
					axis.vec[0]=1.0;
					axis.vec[2]=1.0;
					axis.vec[1]=-(a.vec[0]+a.vec[2])/a.vec[1];
				} else{
					if(fabs(a.vec[2])>EPS){
						axis.vec[0]=1.0;
						axis.vec[1]=1.0;
						axis.vec[2]=-(a.vec[0]+a.vec[1])/a.vec[2];
					}
				}
			}
		}
	}
	return AC2Rot(axis,cos_angle);
}

template<typename T>
inline Vec3<T> Rot2AlphaBetaGamma(Mat33<T> &rot)
{
	Vec3<T> result;
	/*    Z(a)            Y(b)           Z(c)
	 * ( ca  sa  0) * ( cb  0  sb) * ( cc  sc  0)   ( ca  sa  0) * ( cb*cc  cb*sc sb)   ( ca*cb*cc+sa*sc  ca*cb*sc+sa*cc  ca*sb)
	 * (-sa  ca  0) * ( 0   1   0) * (-sc  cc  0) = (-sa  ca  0) * (  -sc     cc  0 ) = (-sa*cb*cc-ca*sc -sa*cb*sc+ca*cc -sa*sb)
	 * ( 0   0   1) * (-sb  0  cb) * ( 0   0   1)   ( 0   0   1) * (-sb*cc -sb*sc cb)   (    -sb*cc          -sb*sc         cb )
	 */
	
	T cb=rot.mat[2][2];
	if(cb>1.0) cb=1.0;
	if(cb<-1.0) cb=-1.0;
	T y=-rot.mat[1][2];
	T x=rot.mat[0][2];
	
	result.vec[0] = atan2(y,x);
	result.vec[1] = acos(cb);
	
	y=-rot.mat[2][1];
	x=-rot.mat[2][0];
	
	result.vec[2] = atan2(y,x);
	return result;
}

struct double4
{
	double x,y,z,w;
	double4(){}
	template<typename T>
	double4(T v){ x = y = z = w = v; }
	
	double4 operator*(const double4& other)
	{
		double4 tmp;
		tmp.x = x*other.x;
		tmp.y = y*other.y;
		tmp.z = z*other.z;
		tmp.w = w*other.w;
		return tmp;
	}
	
	template<typename T>
	double4 operator*(const T& other)
	{
		double4 tmp;
		tmp.x = x*other;
		tmp.y = y*other;
		tmp.z = z*other;
		tmp.w = w*other;
		return tmp;
	}
	
	double4& operator+=(const double4& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
		return *this;
	}
	
	double4& operator-=(const double4& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
		return *this;
	}
	
	template<typename T>
	double4& operator*=(T scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
		return (*this);
	}
	
	template<typename T>
	double4& operator/=(T scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
		return (*this);
	}
};

inline double4 fabs(const double4& a)
{
	double4 tmp;
	tmp.x = a.x < 0.f ? 0.f  : a.x;
	tmp.y = a.y < 0.f ? 0.f  : a.y;
	tmp.z = a.z < 0.f ? 0.f  : a.z;
	tmp.w = a.w < 0.f ? 0.f  : a.w;
	return tmp;
}

inline double4 operator+(const double4& a,const double4& b)
{
	double4 tmp;
	tmp.x = a.x + b.x;
	tmp.y = a.y + b.y;
	tmp.z = a.z + b.z;
	tmp.w = a.w + b.w;
	return tmp;
}

inline double4 operator-(const double4& a,const double4& b)
{
	double4 tmp;
	tmp.x = a.x - b.x;
	tmp.y = a.y - b.y;
	tmp.z = a.z - b.z;
	tmp.w = a.w - b.w;
	return tmp;
}

inline double4 operator*(const double4& a,const double& s)
{
	double4 tmp;
	tmp.x = a.x*s;
	tmp.y = a.y*s;
	tmp.z = a.z*s;
	tmp.w = a.w*s;
	return tmp;
}

inline double4 operator/(const double4& a,const double& s)
{
	double4 tmp;
	tmp.x = a.x/s;
	tmp.y = a.y/s;
	tmp.z = a.z/s;
	tmp.w = a.w/s;
	return tmp;
}

inline double4 cross(const double4& p0, const double4& p1)
{
	double4 result;
	result.w=0.0;
	
	result.x=p0.y*p1.z-p0.z*p1.y;
	result.y=p0.z*p1.x-p0.x*p1.z;
	result.z=p0.x*p1.y-p0.y*p1.x;
	
	return result;
}

inline double dot(const double4& p0, const double4& p1)
{
	return p0.x*p1.x+p0.y*p1.y+p0.z*p1.z+p0.w*p1.w;
}

inline double4 normalize(const double4& a)
{
	double norm=dot(a,a);
	norm=1.0/sqrt(norm);
	return a*norm;
}

typedef double double16[16];

template<typename T>
inline double4 create_double4(T a, T b, T c, T d)
{
	double4 result;
	result.x = a;
	result.y = b;
	result.z = c;
	result.w = d;
	return result;
}

template<typename T>
inline double4 create_double4(Vec3<T> v)
{
	double4 result;
	result.x = v.vec[0];
	result.y = v.vec[1];
	result.z = v.vec[2];
	result.w = 0.0;
	return result;
}

template<typename T>
inline double4 create_double4(T d)
{
	double4 result;
	result.x = d;
	result.y = d;
	result.z = d;
	result.w = d;
	return result;
}

template<typename T>
inline double4 AxisAngle2Quaternion(Vec4<T>& aa)
{
	double4 result;
	result.w=cos(aa.vec[3]/2.0);
	T norm=sin(aa.vec[3]/2.0)/sqrt(aa.vec[0]*aa.vec[0]+aa.vec[1]*aa.vec[1]+aa.vec[2]*aa.vec[2]);
	result.x=aa.vec[0]*norm;
	result.y=aa.vec[1]*norm;
	result.z=aa.vec[2]*norm;
	
	return result;
}

template<typename T>
inline Mat33<T> RotFromEigenvectors(Mat33<T> ev)
{
	Mat33<T> rot;
	// check if eigenvectors are orthogonal
#if DEBUG_LEVEL>2
	cout << "Eigenvectors:\n#1 = (" << ev.ColumnVec3(0).V3Str(',') << ")\n#2 = (" << ev.ColumnVec3(1).V3Str(',') << ")\n#3 = (" << ev.ColumnVec3(2).V3Str(',') << ")\n";
#endif
	bool no12=(fabs(ev.ColumnVec3(0)*ev.ColumnVec3(1))>EPS);
	bool no13=(fabs(ev.ColumnVec3(0)*ev.ColumnVec3(2))>EPS);
	bool no23=(fabs(ev.ColumnVec3(1)*ev.ColumnVec3(2))>EPS);
	bool all_equal=false;
	if(no12 || no13 || no23){
#if DEBUG_LEVEL>2
		cout << "first eigenvector times second: " << ev.ColumnVec3(0)*ev.ColumnVec3(1) << " (" << (ev.ColumnVec3(0)-ev.ColumnVec3(1)).V3Norm() << ")\n";
		cout << "first eigenvector times third: " << ev.ColumnVec3(0)*ev.ColumnVec3(2) << " (" << (ev.ColumnVec3(0)-ev.ColumnVec3(2)).V3Norm() << ")\n";
		cout << "second eigenvector times third: " << ev.ColumnVec3(1)*ev.ColumnVec3(2) << " (" << (ev.ColumnVec3(1)-ev.ColumnVec3(2)).V3Norm() << ")\n";
#endif
		bool diff12=((ev.ColumnVec3(0)-ev.ColumnVec3(1)).V3Norm()<EPS);
		bool diff13=((ev.ColumnVec3(0)-ev.ColumnVec3(2)).V3Norm()<EPS);
		bool diff23=((ev.ColumnVec3(1)-ev.ColumnVec3(2)).V3Norm()<EPS);
		if(!(diff12 || diff13 || diff23)){ // numerical issue (calculating the respective crossproduct again should get solved)
			if(no12) diff12=true;
			if(no13) diff13=true;
			if(no23) diff23=true;
		}
		if((unsigned int)(diff12+diff13+diff23)<=1){
			Vec3<T> a,b;
			unsigned int col=0;
#if DEBUG_LEVEL>1
			cout << "WARNING: Eigenvectors ";
#endif
			if(diff12){ // #1=#2
				col=1;
				a=ev.ColumnVec3(0);
				b=ev.ColumnVec3(2);
#if DEBUG_LEVEL>1
				cout << "#1 and #2";
#endif
			} else{
				if(diff13){ // #1=#3
					col=2;
					a=ev.ColumnVec3(0);
					b=ev.ColumnVec3(1);
#if DEBUG_LEVEL>1
					cout << "#1 and #3";
#endif
				} else{ // diff23 (#2=#3)
					col=2;
					a=ev.ColumnVec3(0);
					b=ev.ColumnVec3(1);
#if DEBUG_LEVEL>1
					cout << "#2 and #3";
#endif
				}
			}
			Vec3<T> new_ev=a.V3Cross(b);
			ev.mat[0][col]=new_ev.vec[0];
			ev.mat[1][col]=new_ev.vec[1];
			ev.mat[2][col]=new_ev.vec[2];
#if DEBUG_LEVEL>1
			cout << " were not orthogonal";
#if DEBUG_LEVEL>2
			cout << ", new eigenvector #" << col+1 << " is: " << new_ev.V3Str(',') << "\n";
#else
			cout << ".\n";
#endif
#endif
		} else all_equal=true;
	}
	// now determine rotation matrix
	rot.M3Eye(); // start with identity matrix
	if(!all_equal){
		bool xz=false; bool yz=false; bool zz=false;
		unsigned int zerocount=0;
		Vec3<T> x=ev.ColumnVec3(0);
		T r=x.V3Norm();
		if(r>EPS) x/=r; else{ xz=true; zerocount++; }
		Vec3<T> y=ev.ColumnVec3(1);
		r=y.V3Norm();
		if(r>EPS) y/=r; else{ yz=true; zerocount++; }
		Vec3<T> z=ev.ColumnVec3(2);
		r=z.V3Norm();
		if(r>EPS) z/=r; else{ zz=true; zerocount++; }
		if(zerocount<=1){
			if(xz){
				x=y;
				x.V3Cross(z);
			} else{
				if(yz){
					y=z;
					y.V3Cross(x);
				} else{
					if(zz){
						z=x;
						x.V3Cross(y);
					}
				}
			}
			rot.mat[0][0]=x.vec[0]; rot.mat[0][1]=y.vec[0]; rot.mat[0][2]=z.vec[0];
			rot.mat[1][0]=x.vec[1]; rot.mat[1][1]=y.vec[1]; rot.mat[1][2]=z.vec[1];
			rot.mat[2][0]=x.vec[2]; rot.mat[2][1]=y.vec[2]; rot.mat[2][2]=z.vec[2];
		} else{
			if(zerocount==2){
				if(!xz){
					rot=Vec2Rot(x,0);
				} else{
					if(!yz){
						rot=Vec2Rot(y,1);
					} else{
						if(!zz){
							rot=Vec2Rot(z,2);
						}
					}
				}
			}
		}
	}
	// Need to check determinat being +1, if -1 coordinate system ended up being right-handed ...
	T determinant=rot.M3Det();
	if(fabs(fabs(determinant)-1.0)>EPS){ // looks weird, but is correct b/c det(rot)=-1 is also acceptable here (recoverable)
		cout << "Could not determine proper rotation matrix, det(rot)=" << determinant << ", which is not +1\n";
		exit(2);
	}
	if(determinant<0.0){ // right-handed coordinate system needs to be changed to left-handed
		rot.mat[0][0]*=-1.0; // do this by flipping x around
		rot.mat[1][0]*=-1.0;
		rot.mat[2][0]*=-1.0;
#if DEBUG_LEVEL>2
		cout << "Fixed improper rotation matrix\n";
#endif
	}
#if DEBUG_LEVEL>2
	cout << "rotation matrix:\n" << rot.M3Str() << "\n";
	cout << "det(rot) = " << rot.M3Det() << "\n";
	bool sing=false;
	cout << "inverse rotation matrix:\n" << (rot.M3Inv(sing)).M3Str() << "\n";
	cout << "rotated x: " << (rot*Vec3<T>(1.0,0.0,0.0)).V3Str(',') << "\n";
	cout << "rotated y: " << (rot*Vec3<T>(0.0,1.0,0.0)).V3Str(',') << "\n";
	cout << "rotated z: " << (rot*Vec3<T>(0.0,0.0,1.0)).V3Str(',') << "\n";
#endif
	return rot;
}

template<typename T>
inline void BackSubstitute(Mat33<T> &A, Vec3<T>& solution, unsigned l)
{
	bool diagzero[3], columnzero[3];
	for(unsigned int j=0; j<3; j++){
		diagzero[j]=(fabs(A.mat[j][j])<=EPS);
		columnzero[j]=false;
		for(unsigned int k=0; k<3; k++) columnzero[j]&=(fabs(A.mat[k][j])<=EPS);
	}
	bool goty=false;
	if(diagzero[2]){
		if(!columnzero[2] && (l==2)){
			if(diagzero[1]){
				if(fabs(A.mat[1][2])>EPS){
					solution.vec[2] = solution.vec[2]/A.mat[1][2]; // z-solution is uniquely defined
				} else{ // element 02 is non-zero
					if(fabs(A.mat[0][1])>EPS){ // y-solution is defined here and can compensate any z-solution (may as well be 1 then)
						solution.vec[2]=1.0;
					} else{ // y-solution is arbitrary
						if(diagzero[0]){ // x-solution is arbitary, but z-solution is defined
							solution.vec[2] = solution.vec[2]/A.mat[0][2];
						} else solution.vec[2] = 1.0; // since x-solution exists it can compensate any z-solution (may as well be 1)
					}
				}
			} else{ // unique y-solution exists
				if(fabs(A.mat[1][2])>EPS){ // we care for z-solution existing but not really for y-solution, so let's try to eliminate the y-solution
					// element 01 *must* exist (otherwise Gaussian elimination is not sorted)
					// no matter if x-solution exists or not, y and z-solutions need to be the
					// same for the first and the second row (and we can set y-solution to zero)
					solution.vec[2] = solution.vec[1]/A.mat[1][2];
					solution.vec[1] = 0.0;
					goty=true;
				} else{ // second row has nothing to do with z-solution which defines y-solution (and means z-solution is in first row)
					solution.vec[1] = solution.vec[1]/A.mat[1][1];
					goty=true;
					if(fabs(A.mat[0][1])>EPS){
						if(diagzero[0]){ // x-solution is arbitrary (has nothing to do with first row), hence y-solution solves unique z-solution
							solution.vec[2] = (solution.vec[0] - solution.vec[1]*A.mat[0][1])/A.mat[0][2];
						} else solution.vec[2] = 1.0; // x-solution can compensate any y and z-solution (which may as well be 1 then)
					} else solution.vec[2] = 1.0; // z-solution is independent of y-solution, and no matter what the x-solution is the z-solution is either arbitrary or compensate (aka may as well be 1)
				}
			}
		} else solution.vec[2]=0.0; // it does not matter what number we choose for z-solution, may as well get rid of it ...
	} else solution.vec[2] = solution.vec[2]/A.mat[2][2]; // z-solution is uniquely defined
	if(!goty){
		if(diagzero[1]){
			if(!columnzero[1] && (l==1)){ // y-solution is defined in first row and we care (b/c it's corresponding to the largest eigenvalue)
				if(diagzero[0]){ // x-solution is arbitrary and has nothing to do with anything, hence y-solution is defined
					solution.vec[1] = (solution.vec[0] - solution.vec[2]*A.mat[0][2])/A.mat[0][1];
				} else{ // x-solution exists and can compensate for y-solution (which may as well be 1 then)
					solution.vec[1] = 1.0;
				}
			} else solution.vec[1] = 0.0; // we can't be bothered
		} else solution.vec[1] = (solution.vec[1] - A.mat[1][2]*solution.vec[2])/A.mat[1][1]; // y-solution is uniquely defined
	}
	if(diagzero[0]){ // x-solution can be anything
		if(l==0) solution.vec[0]=1.0; else solution.vec[0]=0.0;
	} else solution.vec[0] = (solution.vec[0] - A.mat[0][1]*solution.vec[1] - A.mat[0][2]*solution.vec[2])/A.mat[0][0];
}

template<typename T>
inline bool M3equal(Mat33<T> A, Mat33<T> B, T error)
{
	for(int i = 0; i < 3; i++){
		if(fabs(A.mat[i][0]-B.mat[i][0])>error) return false;
		if(fabs(A.mat[i][1]-B.mat[i][1])>error) return false;
		if(fabs(A.mat[i][2]-B.mat[i][2])>error) return false;
	}
	return true;
}

template<typename T>
inline T touch_sphere_sigma(Vec3<T> &r, Vec3<T> &saxes, Mat33<T> &rot, T rT)
{
	if(rT<EPS) return EllipsoidRmin(r,saxes,rot);
	// Minimization vectors
	Vec3<T> V, CI, saxes2;
	T t;
	
	// Create lab frame version of A and B matrices (both are symmetric matrices R * L_A/B * R^T)
	// A/B_ij=sum_k L_A/B_k*R_ik*R_jk
	// move into A ellipsoid frame of reference
	// doing so:
	// - replaces 36 multiplication and 12 additions with 36 multiplications and 24 additions once
	// - saves 4 multiplications and 6 additions in loop
	Vec3<T> z; // multply with kk's transposed (inverse) rotation matrix
	z.vec[0] = r.vec[0]*rot.mat[0][0]+r.vec[1]*rot.mat[1][0]+r.vec[2]*rot.mat[2][0];
	z.vec[1] = r.vec[0]*rot.mat[0][1]+r.vec[1]*rot.mat[1][1]+r.vec[2]*rot.mat[2][1];
	z.vec[2] = r.vec[0]*rot.mat[0][2]+r.vec[1]*rot.mat[1][2]+r.vec[2]*rot.mat[2][2];
	T r2=r*r;
	T d=r2/(saxes.vec[0]*saxes.vec[1]*saxes.vec[2]);
	T e=d*d;
	T rT2 = rT*rT*e;
	saxes2.vec[0] = saxes.vec[0]*saxes.vec[0]*e;
	saxes2.vec[1] = saxes.vec[1]*saxes.vec[1]*e;
	saxes2.vec[2] = saxes.vec[2]*saxes.vec[2]*e;
	z*=d;
	
	T lambda=0.5;
	T xlx=1.0; // x=lambda/(1-lambda) -> do manual calculation with value above
	T Var = 1.0; // trying different forms found 6 loops gives 7 figs for Fx
	T VAV, VBV, det, Vz, V22;
	while(Var > 1E-6){ // loop until variance in distance is sufficiently small
		Var = lambda; // keep lambda around but do CI matrix in terms of x (scale CI by 1/(1-lambda))
		//Populate CI matrix
		CI.vec[0] = xlx*rT2+saxes2.vec[0];
		CI.vec[1] = xlx*rT2+saxes2.vec[1];
		CI.vec[2] = xlx*rT2+saxes2.vec[2];
		
		// Solve z = CI*V for V using inverse => V=CI^-1*z
		V.vec[0] = CI.vec[1]*CI.vec[2]*z.vec[0];
		V.vec[1] = CI.vec[0]*CI.vec[2]*z.vec[1];
		V.vec[2] = CI.vec[0]*CI.vec[1]*z.vec[2];
		
		// VAV=V*A_LF*V (uses fact that A_LF is symmetric)
		V22=V.vec[2]*V.vec[2];
		VAV = V.vec[0]*V.vec[0]*saxes2.vec[0]+V.vec[1]*V.vec[1]*saxes2.vec[1]+V22*saxes2.vec[2];
		// denominator=V*B_LF*V (uses fact that B_LF is symmetric)
		VBV = (V.vec[0]*V.vec[0]+V.vec[1]*V.vec[1]+V22)*rT2;
		//Calculate minimization parameter lambda
/*		if(VBV < EPS*EPS){
			cout << "ERROR: Denominator between oids in touch is too close to zero (" << VBV << ").\n";
			exit(3);
		}*/
		xlx = sqrt(VAV/VBV); // independent of z-scaling (and determinant) -> also, interesting note: the sqrt is better than anything else in terms of speed and convergence
		lambda = xlx/(1.0+xlx);
		Var -= lambda;
		Var *= Var;
	}
	
	//Reconstruct CI and run a final iteration once converged
	CI.vec[0] = xlx*rT2+saxes2.vec[0];
	CI.vec[1] = xlx*rT2+saxes2.vec[1];
	CI.vec[2] = xlx*rT2+saxes2.vec[2];
	
	t=CI.vec[1]*CI.vec[2];
	det = CI.vec[0]*t;
	V22=2.0*z.vec[2];
	Vz=z.vec[0]*z.vec[0]*t+CI.vec[0]*CI.vec[2]*z.vec[1]*z.vec[1]+CI.vec[0]*CI.vec[1]*z.vec[2]*z.vec[2];
	// return sigma
	return sqrt(r2*det/(lambda*Vz));
}

#endif

