//============================================================================
// Name        : Libint_Interface.cpp
// Author      : Jon Dullea
// Version     :
// Copyright   : Your copyright notice
// Description : Interface to Libint Integral Engine: not working as of now
//============================================================================

#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <libint2.h>
#include <boys.h>
#include <utility>
#include <algorithm>
#include <numeric>
#include <util.h>
#include <array>
#include <ctime>
#include "libint2_params.h"
#include "libint2_types.h"
#include<Eigen/Dense>
#include <array>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <iterator>
#include<stdlib.h>
#include <fstream>

#define INT_NCART_NN(am) ((((am)+2)*((am)+1))>>1)

using namespace std;

libint2::FmEval_Chebyshev3 fmeval_chebyshev(30);

class shellSharedData{  //class that holds data that is shared by basis functions of the same type
public:
	std::vector <double> exponents_;
	std::vector <double> coeff_;
	short shellAngMo_;
	short contDepth_;

	shellSharedData(vector<double> exp, vector<double> coeff, short angMo)
	:exponents_ (exp),
	coeff_ (coeff),
	shellAngMo_ (angMo)
	{
		contDepth_ = coeff_.size();
	}
};

class atom{
public:
	short Zval;
	std::array<double, 3> center_;
	atom(int z, double center[3])
		:Zval (z)
	{
		for(int i = 0; i< 3; i++) {center_[i] = center[i];}
	}

};

class shell{                              //bad names in this shell class, contraction used instead of primitive a lot.  Will change soon
public:
	atom* atom_;
	shellSharedData* shellData_;
	int shellNum = -1;

	shell(atom* Atom, shellSharedData* shared)
			:atom_ (Atom),
			 shellData_ (shared){}


	short getAngMo(){return shellData_->shellAngMo_;}
	short getContDepth(){return shellData_->contDepth_;}          //num primitives
	//TODO incorporate numbr of basis functions here?
	const std::vector <double>& getExp() const {return shellData_->exponents_;}
	double getExp(short p){return shellData_->exponents_[p];}
	const std::vector <double>& getCoeff() const {return shellData_->coeff_;}
	double getCoeff(short p){return shellData_->coeff_[p];}
	const std::array<double, 3>& getCenter() const {return atom_->center_;}
	double getCenter(short p){return atom_->center_[p];}
	void printShell(){
		cout<<"Shell Number: "<<shellNum<<endl;
		cout<<"\tAm: "<<shellData_->shellAngMo_<<endl;
		cout<<"\tNumPrim: "<<this->getContDepth()<<endl;
		cout<<"\tOrbital Exp: ";
		for(int i = 0; i<this->getContDepth(); i++){
			cout<<shellData_->exponents_[i]<<" ";
		}
		cout<<"\n\tOrbital Coeff: ";
		for(int i = 0; i<this->getContDepth(); i++){
			cout<<shellData_->coeff_[i]<<" ";
		}
		cout<<"\n\tZ coord: "<< atom_->center_[2]<<endl;

	}
};

#define S_Shells 0
#define P_Shells 1
#define D_Shells 2
#define F_Shells 3
#define G_Shells 4
#define H_Shells 5

std::vector <int> numShellType(3);  //next few functions define a basis set for Hydrogen.  ccpvtz data (kind of, added some functions from some deleted some from others)
int numShells = 0;
std::vector< std::vector < shell > > Vec_Shells;
std::vector<shell> shells;

void make_H_P_S(atom* patom, std::vector<shellSharedData>& shared){
	{
		shell s0(patom, &shared[0]);
    	Vec_Shells[S_Shells].push_back(s0);
    	shells.push_back(s0);
	}
    	numShellType[0] += 1;

	{
    	shell s3(patom, &shared[3]);
    	Vec_Shells[P_Shells].push_back(s3);
    	shells.push_back(s3);
	}
	numShellType[1] += 1;

}

void Make_H(atom* patom, std::vector<shellSharedData>& shared){
	{
		shell s0(patom, &shared[0]);
    	Vec_Shells[S_Shells].push_back(s0);
    	shells.push_back(s0);
	}
	{
		shell s1(patom, &shared[1]);
    	Vec_Shells[S_Shells].push_back(s1);
    	shells.push_back(s1);
	}
	{
		shell s2(patom, &shared[2]);
    	Vec_Shells[S_Shells].push_back(s2);
    	shells.push_back(s2);
	}

    	numShellType[0] += 3;

	{
    	shell s3(patom, &shared[3]);
    	Vec_Shells[P_Shells].push_back(s3);
    	shells.push_back(s3);
	}
	{

    	shell s4(patom, &shared[4]);
    	Vec_Shells[P_Shells].push_back(s4);
    	shells.push_back(s4);
	}
	numShellType[1] += 2;
	{
    	shell s5(patom, &shared[5]);
    	Vec_Shells[D_Shells].push_back(s5);
    	shells.push_back(s5);
	}


    	numShellType[2] += 1;
    	numShells+=6;
}

void Make_H_withG(atom* patom, std::vector<shellSharedData>& shared){
	{
		shell s0(patom, &shared[0]);
    	Vec_Shells[S_Shells].push_back(s0);
    	shells.push_back(s0);
	}
	{
		shell s1(patom, &shared[1]);
    	Vec_Shells[S_Shells].push_back(s1);
    	shells.push_back(s1);
	}
	{
		shell s2(patom, &shared[2]);
    	Vec_Shells[S_Shells].push_back(s2);
    	shells.push_back(s2);
	}

    	numShellType[0] += 3;

	{
    	shell s3(patom, &shared[3]);
    	Vec_Shells[P_Shells].push_back(s3);
    	shells.push_back(s3);
	}
	{

    	shell s4(patom, &shared[4]);
    	Vec_Shells[P_Shells].push_back(s4);
    	shells.push_back(s4);
	}
	numShellType[1] += 2;
	{
    	shell s5(patom, &shared[5]);
    	Vec_Shells[D_Shells].push_back(s5);
    	shells.push_back(s5);
	}
	{
		shell s6(patom, &shared[6]);
    	Vec_Shells[F_Shells].push_back(s6);
    	shells.push_back(s6);
	}
	{
		shell s7(patom, &shared[7]);
    	Vec_Shells[G_Shells].push_back(s7);
    	shells.push_back(s7);
	}
    	numShellType[2] += 1;
    	numShellType[3] += 1;
    	numShellType[4] += 1;
    	numShells+=8;
}

void Make_H_withF(atom* patom, std::vector<shellSharedData>& shared){
	{
		shell s0(patom, &shared[0]);
    	Vec_Shells[S_Shells].push_back(s0);
    	shells.push_back(s0);
	}
	{
		shell s1(patom, &shared[1]);
    	Vec_Shells[S_Shells].push_back(s1);
    	shells.push_back(s1);
	}
	{
		shell s2(patom, &shared[2]);
    	Vec_Shells[S_Shells].push_back(s2);
    	shells.push_back(s2);
	}

    	numShellType[0] += 3;

	{
    	shell s3(patom, &shared[3]);
    	Vec_Shells[P_Shells].push_back(s3);
    	shells.push_back(s3);
	}
	{

    	shell s4(patom, &shared[4]);
    	Vec_Shells[P_Shells].push_back(s4);
    	shells.push_back(s4);
	}
	numShellType[1] += 2;
	{
    	shell s5(patom, &shared[5]);
    	Vec_Shells[D_Shells].push_back(s5);
    	shells.push_back(s5);
	}
	{
		shell s6(patom, &shared[6]);
    	Vec_Shells[F_Shells].push_back(s6);
    	shells.push_back(s6);
	}

    	numShellType[2] += 1;
    	numShellType[3] += 1;
}

void createSharedData(std::vector <shellSharedData>& sharedData){
	{
	shellSharedData shared({33.8700000, 5.0950000, 1.1590000}, {0.0060680, 0.0453080,0.2028220}, 0);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared({.3258000}, {1.0}, 0);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared({0.102700}, {1.0}, 0);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared({1.4070000}, {1.0}, 1);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared({0.3880000}, {1.0}, 1);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared( {1.0570000}, {1.0}, 2);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared( {2.35689099}, {1.0}, 3);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared( {1.95092199}, {1.0}, 4);
	sharedData.push_back(shared);
	}
	{
	shellSharedData shared( {1.95092199}, {1.0}, 5);
	sharedData.push_back(shared);
	}

}


int numH;  //creates a chain of hydrogens
void create_H_Chain(std::vector <atom> &atoms ,std::vector <shellSharedData>& sharedData, int numH){
	for(int i = 0; i < numH ; i++){
		double center[3] = {0.0, 0.0, double(i)};
		atom tempAtom(1 , center);
		atoms.push_back(tempAtom);
	}
	for(int i = 0; i< numH; i++){
		Make_H(&atoms[i], sharedData);
	}
}


const short maxAm = 2;
const static int vectorLength = 2;

/*The following functions are used to initialize the SIMD vectors from arrays of data.  they essentially have the functionality of
 * 			return data in and array. SIMD vectors are initialized from the return values */
template <short Vec_Len_,typename T >
class init_center{
public:
	static void result(T* f_vec, std::vector< std::vector < shell* > >& data,
													 short ABCD, short xyz){
		const auto out = f_vec+1;
		*f_vec = data[vectorLength-Vec_Len_][ABCD]->getCenter()[xyz];
		init_center<Vec_Len_-1,T>::result(out, data , ABCD ,xyz);
	}
};

template <typename T>
class init_center<1,T>{
public:
	static void result(T* f_vec, std::vector< std::vector < shell* > >& data,short ABCD, short xyz){
		*f_vec = data[vectorLength-1][ABCD]->getCenter()[xyz];
	}
};

template <short Vec_Len_,typename T >
class init_exp{
public:
	static void result(T* f_vec, std::vector< std::vector < shell* > >& data,
													 short ABCD, short p){
		const auto out = f_vec+1;
		*f_vec = data[vectorLength-Vec_Len_][ABCD]->getExp(p);
		init_exp<Vec_Len_-1,T>::result(out, data , ABCD , p);
	}
};

template <typename T>
class init_exp<1,T>{
public:
	static void result(T* f_vec, std::vector< std::vector < shell* > >& data,short ABCD, short p){
		*f_vec = data[vectorLength-1][ABCD]->getExp(p);
	}
};

template <short Vec_Len_,typename T >
class init_coeff{
public:
	static void result(T* f_vec, std::vector< std::vector < shell* > >& data,
													 short ABCD, short p){
		*f_vec = data[vectorLength-Vec_Len_][ABCD]->getCoeff(p);
		init_coeff<Vec_Len_-1,T>::result(f_vec+1, data , ABCD ,p);
	}
};

template <typename T>
class init_coeff<1,T>{
public:
	static void result(T* f_vec, std::vector< std::vector < shell* > >& data,short ABCD, short p){
		*f_vec = data[vectorLength-1][ABCD]->getCoeff(p);
	}
};

#define x 0
#define y 1
#define z 2

#define A_ 0
#define B_ 1
#define C_ 2
#define D_ 3


template <typename LibintEval>
void inline preplibint2(std::vector<LibintEval>& erievals,
		std::vector < std::vector < shell*> > &IntegralSet, std::array <short, 4> angMo){

		/*permute only needed for verification.  actual loop ensures libint ordering
		 * this was not used in this version of the code.  if conditions never evaluate to true for all shells: am0 >= am1 , am2 >= am3 & am01 <= am23*/

		if(angMo[0] < angMo[1]){
			std::swap(angMo[1], angMo[0]);
			std::swap(IntegralSet[0][0], IntegralSet[0][1]);
			std::swap(IntegralSet[1][0], IntegralSet[1][1]);
		}
		if(angMo[2] < angMo[3]){
			std::swap(angMo[2], angMo[3]);
			std::swap(IntegralSet[0][2], IntegralSet[0][3]);
			std::swap(IntegralSet[1][2], IntegralSet[1][3]);
		}
		if(angMo[0]+angMo[1] > angMo[2]+angMo[3]){
			std::swap(angMo[0], angMo[2]);
			std::swap(angMo[1], angMo[3]);

			std::swap(IntegralSet[0][0], IntegralSet[0][2]);
			std::swap(IntegralSet[1][0], IntegralSet[1][2]);

			std::swap(IntegralSet[0][1], IntegralSet[0][3]);
			std::swap(IntegralSet[1][1], IntegralSet[1][3]);

		}


		int derivOrder = 0;
		double static F[vectorLength][LIBINT_MAX_AM*4 + 6];

		const LIBINT2_REALTYPE two_times_PI_to_2_5 = {34.986836655249725693, 34.986836655249725693};

		int contrdepth_i = IntegralSet[0][0]->getContDepth();
		int contrdepth_j = IntegralSet[0][1]->getContDepth();
		int contrdepth_k = IntegralSet[0][2]->getContDepth();
		int contrdepth_l = IntegralSet[0][3]->getContDepth();


		unsigned short int am0  = angMo[0];
		unsigned short int am1  = angMo[1];
		unsigned short int am2  = angMo[2];
		unsigned short int am3  = angMo[3];
		unsigned short int am01 = am0 + am1;
		unsigned short int am23 = am2 + am3;
		unsigned short int amtot= am01 + am23;

		double temp[vectorLength];

		init_center<vectorLength, double>::result(temp, IntegralSet, A_ , x);       //these are the aforementioned initialization functions. A_ and x are macros
		const LIBINT2_REALTYPE Ax(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, A_ , y);
		const LIBINT2_REALTYPE Ay(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, A_ , z);
		const LIBINT2_REALTYPE Az(temp);

		init_center<vectorLength, double>::result(temp, IntegralSet, B_ , x);
		const LIBINT2_REALTYPE Bx(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, B_ , y);
		const LIBINT2_REALTYPE By(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, B_ , z);
		const LIBINT2_REALTYPE Bz(temp);

		const LIBINT2_REALTYPE AB_x = (Ax-Bx);
		//cout<<AB_x<<endl;
		const LIBINT2_REALTYPE AB_y = (Ay-By);
		const LIBINT2_REALTYPE AB_z = (Az-Bz);

		LIBINT2_REALTYPE AB2 = AB_x*AB_x + AB_y*AB_y + AB_z*AB_z;

		init_center<vectorLength, double>::result(temp, IntegralSet, C_ , x);
		const LIBINT2_REALTYPE Cx(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, C_ , y);
		const LIBINT2_REALTYPE Cy(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, C_ , z);
		const LIBINT2_REALTYPE Cz(temp);

		init_center<vectorLength, double>::result(temp, IntegralSet, D_ , x);
		const LIBINT2_REALTYPE Dx(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, D_ , y);
		const LIBINT2_REALTYPE Dy(temp);
		init_center<vectorLength, double>::result(temp, IntegralSet, D_ , z);
		const LIBINT2_REALTYPE Dz(temp);

		const LIBINT2_REALTYPE CD_x  = (Cx-Dx);
		const LIBINT2_REALTYPE CD_y  = (Cy-Dy);
		const LIBINT2_REALTYPE CD_z  = (Cz-Dz);




		LIBINT2_REALTYPE CD2 = CD_x*CD_x + CD_y*CD_y + CD_z*CD_z;

		LIBINT2_REALTYPE gammaP;
		LIBINT2_REALTYPE one_over_gammaP;
		LIBINT2_REALTYPE rhoP;

		LIBINT2_REALTYPE Px;
		LIBINT2_REALTYPE Py;
		LIBINT2_REALTYPE Pz;

		LIBINT2_REALTYPE PAx;
		LIBINT2_REALTYPE PAy;
		LIBINT2_REALTYPE PAz;

		LIBINT2_REALTYPE PBx;
		LIBINT2_REALTYPE PBy;
		LIBINT2_REALTYPE PBz;

		LIBINT2_REALTYPE gammaQ;
		LIBINT2_REALTYPE one_over_gammaQ;
		LIBINT2_REALTYPE rhoQ;

		LIBINT2_REALTYPE Qx;
		LIBINT2_REALTYPE Qy;
		LIBINT2_REALTYPE Qz;

		LIBINT2_REALTYPE QCx;
		LIBINT2_REALTYPE QCy;
		LIBINT2_REALTYPE QCz;

		LIBINT2_REALTYPE QDx;
		LIBINT2_REALTYPE QDy;
		LIBINT2_REALTYPE QDz;

		LIBINT2_REALTYPE one_over_gamma_plus_gamma;
		LIBINT2_REALTYPE gammaPQ;
		LIBINT2_REALTYPE gammap_over_gammap_gammaq;
		LIBINT2_REALTYPE gammaq_over_gammap_gammaq;

		LIBINT2_REALTYPE PQx;
		LIBINT2_REALTYPE PQy;
		LIBINT2_REALTYPE PQz;
		LIBINT2_REALTYPE PQ2;

		LIBINT2_REALTYPE Wx;
		LIBINT2_REALTYPE Wy;
		LIBINT2_REALTYPE Wz;

		LIBINT2_REALTYPE WPx;
		LIBINT2_REALTYPE WPy;
		LIBINT2_REALTYPE WPz;

		LIBINT2_REALTYPE WQx;
		LIBINT2_REALTYPE WQy;
		LIBINT2_REALTYPE WQz;

		LIBINT2_REALTYPE K1;
		LIBINT2_REALTYPE K2;

		LIBINT2_REALTYPE pfac01;
		LIBINT2_REALTYPE pfac012;
		LIBINT2_REALTYPE pfac;



		int p0123 = 0;
	  for (uint p0 = 0; p0 < contrdepth_i; p0++) {
		  init_exp<vectorLength, double>::result(temp, IntegralSet, A_ , p0);
		  const LIBINT2_REALTYPE alpha0(temp);
		  init_coeff<vectorLength, double>::result(temp, IntegralSet, A_ , p0);
		  LIBINT2_REALTYPE pfac0(temp);
		  pfac0 = pfac0 * two_times_PI_to_2_5;

		  for (uint p1 = 0; p1 < contrdepth_j; p1++) {
			  init_exp<vectorLength, double>::result(temp, IntegralSet, B_ , p1);
			  const LIBINT2_REALTYPE alpha1(temp);
			  init_coeff<vectorLength, double>::result(temp, IntegralSet, B_ , p1);
			  const LIBINT2_REALTYPE c1(temp);

			  gammaP = (alpha0 + alpha1);
			  one_over_gammaP = 1.0 / gammaP;
			  rhoP = alpha0 * alpha1 * one_over_gammaP;
			  Px = (Ax*alpha0 + Bx*alpha1)*one_over_gammaP;
			  Py = (Ay*alpha0 + By*alpha1)*one_over_gammaP;
			  Pz = (Az*alpha0 + Bz*alpha1)*one_over_gammaP;
			  K1 = exp(-rhoP*AB2)*one_over_gammaP;
			  pfac01 = pfac0 * c1 * K1;

			  PAx = Px - Ax;
			  PAy = Py - Ay;
			  PAz = Pz - Az;

			  PBx = Px - Bx;
			  PBy = Py - By;
			  PBz = Pz - Bz;

			  for (uint p2 = 0; p2 < contrdepth_k; p2++) {
				  init_exp<vectorLength, double>::result(temp, IntegralSet, C_, p2);
				  const LIBINT2_REALTYPE alpha2(temp);
				  init_coeff<vectorLength, double>::result(temp, IntegralSet, C_, p2);
				  const LIBINT2_REALTYPE c2(temp);
				  pfac012 = pfac01 * c2;

				  for (uint p3 = 0; p3 < contrdepth_l; p3++, ++p0123) {
						 LibintEval* erieval = &erievals[p0123];

						 erieval->veclen = 1;
						 init_exp<vectorLength, double>::result(temp, IntegralSet, D_, p3);
						 const LIBINT2_REALTYPE alpha3(temp);
						 init_coeff<vectorLength, double>::result(temp, IntegralSet, D_, p3);

						 const LIBINT2_REALTYPE c3(temp);
						 pfac = pfac012 * c3;

						 gammaQ = alpha2 + alpha3;
						 one_over_gammaQ = 1.0 / gammaQ;
						 rhoQ = alpha2 * alpha3 * one_over_gammaQ;

						 one_over_gamma_plus_gamma = 1.0/(gammaP+gammaQ);
						 gammaPQ = (gammaP * gammaQ *one_over_gamma_plus_gamma);
						 gammap_over_gammap_gammaq = (gammaPQ * one_over_gammaQ);
						 gammaq_over_gammap_gammaq = (gammaPQ * one_over_gammaP);

						 Qx = (Cx*alpha2 + Dx*alpha3)*one_over_gammaQ;
						 Qy = (Cy*alpha2 + Dy*alpha3)*one_over_gammaQ;
						 Qz = (Cz*alpha2 + Dz*alpha3)*one_over_gammaQ;


						 QCx = Qx - Cx;
						 QCy = Qy - Cy;
						 QCz = Qz - Cz;

						 QDx = Qx - Dx;
						 QDy = Qy - Dy;
						 QDz = Qz - Dz;

						 PQx = Px - Qx;
						 PQy = Py - Qy;
						 PQz = Pz - Qz;

						 PQ2  = PQx*PQx + PQy*PQy + PQz*PQz;

						 Wx = (Px*gammap_over_gammap_gammaq+
								Qx*gammaq_over_gammap_gammaq);
						 Wy = (Py*gammap_over_gammap_gammaq+
								Qy*gammaq_over_gammap_gammaq);
						 Wz = (Pz*gammap_over_gammap_gammaq+
								Qz*gammaq_over_gammap_gammaq);

						 WPx = Wx - Px;
						 WPy = Wy - Py;
						 WPz = Wz - Pz;

						 WQx = Wx - Qx;
						 WQy = Wy - Qy;
						 WQz = Wz - Qz;

						 K2 = exp(-rhoQ*CD2)*one_over_gammaQ;

						 pfac = K2 * sqrt(one_over_gamma_plus_gamma) * pfac;

						 double PQ2_s[vectorLength];
						 PQ2.convert(PQ2_s);
						 double gammaPQ_s[vectorLength];
						 gammaPQ.convert(gammaPQ_s);
						 fmeval_chebyshev.eval(F[0],PQ2_s[0]*gammaPQ_s[0],amtot);
						 fmeval_chebyshev.eval(F[1],PQ2_s[1]*gammaPQ_s[0],amtot);

						 double pfac_s[vectorLength];
						 pfac.convert(pfac_s);
						 LIBINT2_REALTYPE* ssss_ptr = erieval->LIBINT_T_SS_EREP_SS(0);
						 for(unsigned int l=0; l<=amtot; ++l, ssss_ptr++){
							 *ssss_ptr = {pfac_s[0]*F[0][l], pfac_s[1]*F[1][l]};

								}

				#if LIBINT2_DEFINED(eri,PA_x)
						erieval->PA_x[0] = PAx;
				#endif
				#if LIBINT2_DEFINED(eri,PA_y)
						erieval->PA_y[0] = PAy;
				#endif
				#if LIBINT2_DEFINED(eri,PA_z)
						 erieval->PA_z[0] = PAz;
				#endif
				#if LIBINT2_DEFINED(eri,PB_x)
						 erieval->PB_x[0] = PBx;
				#endif
				#if LIBINT2_DEFINED(eri,PB_y)
						 erieval->PB_y[0] = PBy;
				#endif
				#if LIBINT2_DEFINED(eri,PB_z)
						 erieval->PB_z[0] = PBz;
				#endif
				#if LIBINT2_DEFINED(eri,AB_x)
							erieval->AB_x[0] = AB_x;
				#endif
				#if LIBINT2_DEFINED(eri,AB_y)
							erieval->AB_y[0] = AB_y;
				#endif
				#if LIBINT2_DEFINED(eri,AB_z)
							erieval->AB_z[0] = AB_z;
				#endif
				#if LIBINT2_DEFINED(eri,BA_x)
							erieval->BA_x[0] = -AB_x;
				#endif
				#if LIBINT2_DEFINED(eri,BA_y)
							erieval->BA_y[0] = -AB_y;
				#endif
				#if LIBINT2_DEFINED(eri,BA_z)
							erieval->BA_z[0] = -AB_z;
				#endif
				#if LIBINT2_DEFINED(eri,oo2z)
							erieval->oo2z[0] = 0.5*one_over_gammaP;
				#endif
				#if LIBINT2_DEFINED(eri,QC_x)
							erieval->QC_x[0] = QCx;
				#endif
				#if LIBINT2_DEFINED(eri,QC_y)
							erieval->QC_y[0] =  QCy;
				#endif
				#if LIBINT2_DEFINED(eri,QC_z)
							erieval->QC_z[0] =  QCz;
				#endif
				#if LIBINT2_DEFINED(eri,QD_x)
							erieval->QD_x[0] =  QDx;
				#endif
				#if LIBINT2_DEFINED(eri,QD_y)
							erieval->QD_y[0] =  QDy;
				#endif
				#if LIBINT2_DEFINED(eri,QD_z)
							erieval->QD_z[0] =  QDz;
				#endif

				#if LIBINT2_DEFINED(eri,CD_x)
							erieval->CD_x[0] = CD_x;
				#endif
				#if LIBINT2_DEFINED(eri,CD_y)
							erieval->CD_y[0] = CD_y;
				#endif
				#if LIBINT2_DEFINED(eri,CD_z)
							erieval->CD_z[0] =  CD_z;
				#endif
				#if LIBINT2_DEFINED(eri,DC_x)
							erieval->DC_x[0] = -CD_x;
				#endif
				#if LIBINT2_DEFINED(eri,DC_y)
							erieval->DC_y[0] = -CD_y;
				#endif
				#if LIBINT2_DEFINED(eri,DC_z)
							erieval->DC_z[0] = -CD_z;
				#endif
				#if LIBINT2_DEFINED(eri,oo2e)
							erieval->oo2e[0] = 0.5*one_over_gammaQ;
				#endif
					// Prefactors for interelectron transfer relation
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
							erieval->TwoPRepITR_pfac0_0_0_x[0] = - (alpha1*AB_x + alpha3*CD_x)*oogammaP);
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
							erieval->TwoPRepITR_pfac0_0_0_y[0] = - (alpha1*AB_y + alpha3*CD_y)*oogammaP);
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
							erieval->TwoPRepITR_pfac0_0_0_z[0] =  - (alpha1*AB_z + alpha3*CD_z)*oogammaP);
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
							erieval->TwoPRepITR_pfac0_1_0_x[0] = - (alpha1*AB_x + alpha3*CD_x)*oogammaQ);
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
							erieval->TwoPRepITR_pfac0_1_0_y[0] = - (alpha1*AB_y + alpha3*CD_y)*oogammaQ);
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
							erieval->TwoPRepITR_pfac0_1_0_z[0] = - (alpha1*AB_z + alpha3*CD_z)*oogammaQ);
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
							erieval->TwoPRepITR_pfac0_0_1_x[0] = {(alpha0[0]*AB[0][0] + alpha2[0]*CD[0][0])*oogammaP[0],
																		alpha0[1]*AB[1][0] + alpha2[1]*CD[1][0])*oogammaP[1]};
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
							erieval->TwoPRepITR_pfac0_0_1_y[0] = {(alpha0[0]*AB[0][1] + alpha2[0]*CD[0][1])*oogammaP[0],
																		alpha0[1]*AB[1][1] + alpha2[1]*CD[1][1])*oogammaP[1]};
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
							erieval->TwoPRepITR_pfac0_0_1_z[0] = {(alpha0[0]*AB[0][2] + alpha2[0]*CD[0][2])*oogammaP[0],
																		alpha0[1]*AB[1][2] + alpha2[1]*CD[1][2])*oogammap[1]};
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
							erieval->TwoPRepITR_pfac0_1_1_x[0] = {(alpha0[0]*AB[0][0] + alpha2[0]*CD[0][0])*oogammaQ[0],
																		(alpha0[1]*AB[1][0] + alpha2[1]*CD[1][0])*oogammaQ[1]};
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
							erieval->TwoPRepITR_pfac0_1_1_y[0] = {(alpha0[0]*AB[0][1] + alpha2[0]*CD[0][1])*oogammaQ[0],
																		(alpha0[1]*AB[1][1] + alpha2[1]*CD[1][1])*oogammaQ[1]};
				#endif
				#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
							erieval->TwoPRepITR_pfac0_1_1_z[0] = {(alpha0[0]*AB[0][2] + alpha2[0]*CD[0][2])*oogammaQ[0],
																		(alpha0[1]*AB[1][2] + alpha2[1]*CD[1][2])*oogammaQ[1]};
				#endif
				#if LIBINT2_DEFINED(eri,eoz)
							erieval->eoz[0] = {gammaQ[0]*one_over_gammaP[0], gammaQ[1]*one_over_gammaP[1]};
				#endif
				#if LIBINT2_DEFINED(eri,zoe)
							erieval->zoe[0] = {gammap[0]*one_over_gammaQ[0], gammap[1]*one_over_gammaQ[1]};
				#endif


				#if LIBINT2_DEFINED(eri,WP_x)
							erieval->WP_x[0] = WPx;
				#endif
				#if LIBINT2_DEFINED(eri,WP_y)
							erieval->WP_y[0] = WPy;
				#endif
				#if LIBINT2_DEFINED(eri,WP_z)
							erieval->WP_z[0] = WPz;
				#endif
				#if LIBINT2_DEFINED(eri,WQ_x)
							erieval->WQ_x[0] = WQx;
				#endif
				#if LIBINT2_DEFINED(eri,WQ_y)
							erieval->WQ_y[0] = WQy;
				#endif
				#if LIBINT2_DEFINED(eri,WQ_z)
							erieval->WQ_z[0] = WQz;
				#endif
				#if LIBINT2_DEFINED(eri,oo2ze)
							erieval->oo2ze[0] = 0.5/(gammaP+gammaQ);
				#endif
				#if LIBINT2_DEFINED(eri,roz)
							erieval->roz[0] = gammaPQ*one_over_gammaP;
				#endif
				#if LIBINT2_DEFINED(eri,roe)
							erieval->roe[0] = gammaPQ*one_over_gammaQ;
				#endif
							// prefactors for derivative ERI relations
							if (derivOrder > 0) {
				#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
							erieval->alpha1_rho_over_zeta2[0] = {alpha0[0] * gammaPQ[0] / (gammaP[0] * gammaP[0]),
																	alpha0[1] * gammaPQ[1] / (gammaP[1] * gammaP[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
							erieval->alpha2_rho_over_zeta2[0] = {alpha1[0] * gammaPQ[0] / (gammaP[0] * gammaP[0]),
																	alpha1[1] * gammaPQ[1] / (gammaP[1] * gammaP[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
							erieval->alpha3_rho_over_eta2[0] = {alpha2[0] * gammaPQ[0] / (gammaP[0] * gammaP[0]),
																	alpha2[1] * gammaPQ[1] / (gammaP[1] * gammaP[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
							erieval->alpha4_rho_over_eta2[0] = {alpha3[0] * gammaPQ[0] / (gammaP[0] * gammaP[0]),
																	alpha3[1] * gammaPQ[1] / (gammaP[1] * gammaP[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
							erieval->alpha1_over_zetapluseta[0] = {alpha0[0] / (gammaP[0] + gammaQ[0]),
																	alpha0[1] / (gammaP[1] + gammaQ[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
							erieval->alpha2_over_zetapluseta[0] = {alpha1[0] / (gammaP[0] + gammaQ[0]),
																	alpha1[1] / (gammaP[1] + gammaQ[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
							erieval->alpha3_over_zetapluseta[0] = {alpha2[0] / (gammaP[0] + gammaQ[0]),
																	alpha2[1] / (gammaP[1] + gammaQ[1])};
				#endif
				#if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
							erieval->alpha4_over_zetapluseta[0] = {alpha3[0] / (gammaP[0] + gammaQ[0]),
																	alpha3[1] / (gammaP[1] + gammaQ[1])};
				#endif
				#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
							erieval->rho12_over_alpha1[0] = {rhoP[0] / alpha0[0], rhoP[1] / alpha0[1]};
				#endif
				#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
							erieval->rho12_over_alpha2[0] = {rhoP[0] / alpha1[0], rhoP[1] / alpha1[1]};
				#endif
				#if LIBINT2_DEFINED(eri,rho34_over_alpha3)
							erieval->rho34_over_alpha3[0] = {rhoQ[0] / alpha2[0], rhoQ[1] / alpha2[1]};
				#endif
				#if LIBINT2_DEFINED(eri,rho34_over_alpha4)
							erieval->rho34_over_alpha4[0] = {rhoQ[0] / alpha3[0], rhoQ[1] / alpha3[1]};
				#endif
				#if LIBINT2_DEFINED(eri,two_alpha0_bra)
							erieval->two_alpha0_bra[0] = {2.0 * alpha0[0], 2.0 * alpha0[1]};
				#endif
				#if LIBINT2_DEFINED(eri,two_alpha0_ket)
							erieval->two_alpha0_ket[0] = {2.0 * alpha1[0], 2.0 * alpha1[1]};
				#endif
				#if LIBINT2_DEFINED(eri,two_alpha1_bra)
							erieval->two_alpha1_bra[0] = {2.0 * alpha2[0], 2.0 * alpha2[1]};
				#endif
				#if LIBINT2_DEFINED(eri,two_alpha1_ket)
							erieval->two_alpha1_ket[0] = {2.0 * alpha3[0], 2.0 * alpha3[1]};
				#endif
				}
				} //end of Shell group loop
			}
		  }
		}
	  }

std::vector<int> map_shell_to_basis(const std::vector<shell>& shells){
	std::vector < int > result;
	result.reserve(shells.size());

	int n = 0;
	for(auto shell: shells){
		result.push_back(n);
		n += INT_NCART_NN(shell.getAngMo());
	}
	result.push_back(n);
	return result;
}

/*Normalization of the basis functions  this is not implementd in the current code for the sole reason that I was not getting consistnt answers with the refrence code
 * this may have been the source of the issue.  I will look into this soon.*/

static constexpr std::array<int64_t,31> df = {{1, 1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135,
                       645120, 2027025, 10321920, 34459425, 185794560, 654729075,
                       3715891200, 13749310575, 81749606400, 316234143225, 1961990553600,
                       7905853580625, 51011754393600, 213458046676875, 1428329123020800,
                       6190283353629375}};

void renorm(shell* s){              //was not used, nt sure it is correct
	const auto sqrt_pi_cubed = double{5.56832799683170784528481798212};
	const auto np = s->getContDepth();
	for(int p = 0; p!=np; ++p){
		const auto am = s->shellData_->shellAngMo_;
		const auto two_alpha = 2 * s->getExp(p);
		const auto two_alpha_to_am32 = pow(two_alpha, s->shellData_->shellAngMo_+1)*sqrt(two_alpha);
		const auto norm = sqrt(pow(2, am)*two_alpha_to_am32/(sqrt_pi_cubed * df[2*am] ));

		s->shellData_->coeff_[p] *= norm;
	}

}


int main(int argc, char *argv[]) {
	short numH = 0;
	if(argc == 2)
		numH = atoi(argv[1]);           //pass in number of hydrogens as a command line parameter
	else{
		cout<<"enter number of H: ";	//else ask user
		cin>>numH;

	}

	for(int i = 0; i< 6; i++){
		Vec_Shells.push_back(std::vector < shell>() );
	}
	std::vector <shellSharedData> sharedData;
	std::vector <atom> atoms;
	createSharedData(sharedData);
	create_H_Chain(atoms, sharedData, numH);         //creates a chain of hydrogens to test on

	int maxContDepth = 3;

//	std::array offest[maxAm];
//	for(int i = 0; i< maxAm; i++){
//		offest[i] = Vec_Shells.size()-1;
//	}


	std::vector < std::vector < shell*> > integralSet(vectorLength);        //this is what is passed to prep libint function
	for(int i = 0; i< vectorLength; i++){
		integralSet.push_back(std::vector < shell* >() );
	}

	std::vector < std::vector < std::vector <shell* > > > integralSets;
	const int maxPrimComb = 3*3*3*3;                                         //hard coded to 81 possible combinations of basis functions; way too many. lot of wasted storage
	std::vector < int > numInts(maxPrimComb,0);

	for(int i = 0; i< maxPrimComb; i++){
		integralSets.push_back(std::vector < std::vector < shell* > > () );
		for(int j = 0; j<vectorLength; j++){
			integralSets[i].push_back(std::vector <shell*> () );
			for(int k = 0; k< 4; k++){
				integralSets[i][j].push_back(NULL);
			}
		}
	}

	std::vector <short> amVec(4,0);
    std::vector<Libint_t> Libint_(maxPrimComb);

    LIBINT2_PREFIXED_NAME(libint2_static_init)();
    LIBINT2_PREFIXED_NAME( libint2_init_eri)(&Libint_[0], maxAm, 0);

    std::clock_t start;
    double duration;
    start = std::clock();

	cout<<"Number of Shells: "<<numShells<<endl;

	int numbasisFunc = 0;
	for(int am = 0; am <=  maxAm; am++){
		numbasisFunc += Vec_Shells[am].size() * INT_NCART_NN(am);
	}
	cout<<numbasisFunc<<endl;

	int shell_group_offset[maxAm+1]; //array of numbers of shells with particular am
	shell_group_offset[S_Shells] = 0;
	for(int i = 1; i <= maxAm; ++i){
		shell_group_offset[i] = shell_group_offset[i-1] + Vec_Shells[i-1].size();
	}


	Eigen::MatrixXd Fock = Eigen::MatrixXd::Zero(numbasisFunc, numbasisFunc);
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(numbasisFunc, numbasisFunc);


//	std::vector<shell> shells;
//	for(int am = 0; am <= maxAm; am++){
//		for(int i = 0; i< Vec_Shells[am].size(); i++){
////			Vec_Shells[am][i].shellNum = shellCounter;
////			shellCounter++;
//			shells.push_back(Vec_Shells[am][i]);
//		}
//	}
	for(int i = 0; i<shells.size(); i++){
		shells[i].shellNum = i;
		//renorm(&(shells[i]));
	}

	int num_shells = shells.size();

	// Sort shells

	stable_sort(shells.begin(), shells.end(),[](shell x_, shell y_){return x_.getAngMo() < y_.getAngMo();});    //sorts the shells based on AM

	cout << "ordered am" << endl;

	auto shell_b = shells.begin();
	std::vector<decltype(shell_b)> points;
	points.push_back(shell_b);
	auto it = shells.begin();
	for(auto i = 0; i < maxAm+1; ++i){
		points.push_back(std::partition_point(it, shells.end(), [&](shell x_){return x_.getAngMo() == i;}));   //creates itterators to the points that define the partition between the shells
		it = points.back();
	}
	points.push_back(shells.end());
	cout<<"aftre orderdAm"<<endl;
	std::vector<std::string> am ({"s","p","d","f","g","h","i","j"});
	//for(int i = 0; i< numShells; i++){
	//	cout<<"Shell num: " << i << " with am : "<<shells[i].getAngMo()<<endl;
	//}

	std::vector <std::vector <std::vector <shell*> > > shells_by_cont;                   //used inside AM loops, shells with the same number of primitives are computed at the same time.  this stores them until there are enough to compute
	for(int i = 0; i <maxPrimComb ; i++){
		shells_by_cont.push_back(std::vector< std::vector <shell*> > () );
		for(int j = 0; j< vectorLength+1; j++){
			shells_by_cont[i].push_back(std::vector<shell*> () );
			for(int k = 0; k < 4; k++){
				shells_by_cont[i][j].push_back(NULL);
			}
		}
	}

	std::vector<int> howMany;                         //stores how many integrals of a certain type are ready to be computed.  (i.e. n integrals with m primitive combinations)
	for(int i = 0; i< maxPrimComb; i++){
		howMany.push_back(0);
	}

	std::vector< std::vector <std::vector<int> > > firstBasisFuncts;               //vector of the first basis functions for a given shell
	std::vector < std::vector < int > > s1234_degen_vec;                          
	for(int i = 0; i< maxPrimComb; i++){
		firstBasisFuncts.push_back(std::vector<std::vector<int> > () );
		s1234_degen_vec.push_back( std::vector <int> () );
		for(int j = 0; j< vectorLength; j++){
			s1234_degen_vec[i].push_back(0);
			firstBasisFuncts[i].push_back(std::vector <int> () );
			for(int k = 0; k < 4; k++){
				firstBasisFuncts[i][j].push_back(0);
			}
		}
	}

        /* NOTE:  The termonolgy is wrong from here forth. going to fix this soon, ran out of time.  Contraction depth actually refers to number of primitive gaussians per basis function */

	std::vector <int> contStride = {maxContDepth*maxContDepth*maxContDepth, maxContDepth*maxContDepth, maxContDepth, 1};  //way to index primitive combs.

	double tensorInts[numbasisFunc][numbasisFunc][numbasisFunc][numbasisFunc];
	for(int i = 0; i<numbasisFunc; i++){
		for(int j = 0; j<numbasisFunc; j++){
			for(int k = 0; k<numbasisFunc; k++){
				for(int l = 0; l < numbasisFunc; l++){
					tensorInts[i][j][k][l] = -1.0;
				}
			}
		}
	}

	for(int i = 0; i< numShells; i++){
		shells[i].shellNum = i;
	}

	auto shell_to_bf = map_shell_to_basis(shells);

	std::vector<shell*> tempShell(4);
	unsigned long count = 0;

	short am1 = 0;                                                                                 //loop over AM, ensure Libint ordering.  A little strange with the itterators but it works.
	for(auto it_am1 = points[0]; it_am1 != shells.end() && am1 != maxAm+1; am1++){
		it_am1 = points[am1];
		auto it_am1_end = points[am1+1];
		const int  numPrim_1 = INT_NCART_NN(am1);

		short am2 = 0;
		for(auto it_am2 = points[0]; it_am2 != it_am1_end; am2++){
			it_am2 = points[am2];
			auto it_am2_end = points[am2+1];
			const int numPrim_2 = INT_NCART_NN(am2);
			const short am12 = am1+am2;
			if( am1 >= am2){

				short am3 = 0;
				for( auto it_am3 = points[0]; it_am3 != shells.end() && am3!= maxAm+1; ++am3){
					it_am3 = points[am3];
					auto it_am3_end = points[am3 + 1];
					const int numPrim_3 = INT_NCART_NN(am3);
					const short am123 = am12 +am3;

					short am4 = 0;
					for(auto it_am4 = points[0]; it_am4 != it_am3_end; ++am4){
						it_am4 = points[am4];
						auto it_am4_end = points[am4+1];
						const int  numPrim_4 = INT_NCART_NN(am4);

						if(am3 >= am4 && am4 + am3 >= am12 ){
							const short am_tot = am123+am4;
							cout << am[am1] << " " << am[am2]<< " " << am[am3] << " "<< am[am4] <<endl;

							//cout << numPrim_1 << " " << numPrim_2<< " " << numPrim_3 << " "<< numPrim_4 <<endl;

							int numContrComb = 0;
							int maxContComb = 0;
							for(auto it_shell_1 = it_am1; it_shell_1 != it_am1_end; it_shell_1++){                 //loop over shells ensuring that only unique combinations are hit.  
								tempShell[0] = &(*it_shell_1);
								const int shell_num_1 = it_shell_1->shellNum;
								auto bf1_first = shell_to_bf[shell_num_1];
								const int numCont1 = it_shell_1->getContDepth();                         //primitives
								const int contIndex1 = (numCont1 -1) * contStride[0];                   //primitive stride

								for(auto it_shell_2 = it_am2; it_shell_2 != it_am2_end; ++it_shell_2){
									tempShell[1] = &(*it_shell_2);
									const int shell_num_2 = it_shell_2->shellNum;
									auto bf2_first = shell_to_bf[shell_num_2];
									const int numCont2 = it_shell_2->getContDepth();
									const int contIndex12 = (numCont2 -1) * contStride[1] + contIndex1;
									const int numCont12 = numCont1 * numCont2;
									const int  s12_deg = (shell_num_1 == shell_num_2) ? 1.0 : 2.0;

									for(auto it_shell_3 = it_am3; it_shell_3 != it_am3_end; ++it_shell_3){
										tempShell[2] = &(*it_shell_3);
										const int shell_num_3 = it_shell_3->shellNum;
										auto bf3_first = shell_to_bf[shell_num_3];
										const int numCont3 = it_shell_3->getContDepth();
										const int contIndex123 = contIndex12 + (numCont3-1)*contStride[2];  //make cont stride an int not vector
										const int numCont123 = numCont12 * numCont3;

										for(auto it_shell_4 = points[0]; it_shell_4 != it_am4_end; ++it_shell_4){
											tempShell[3] = &(*it_shell_4);
											const int shell_num_4 = it_shell_4->shellNum;
											auto bf4_first = shell_to_bf[shell_num_4];
											numContrComb = numCont123 * it_shell_4->getContDepth();              //Again this is the wrong terminology
											maxContComb = max(maxContComb, numContrComb);

											int contIndex = contIndex123 + ((it_shell_4->getContDepth())-1);          //compute index based on number of primitives

											const int s34_deg = (shell_num_3 == shell_num_4) ? 1.0 : 2.0;
											const int  s12_34_deg = (shell_num_1 == shell_num_3) ? (shell_num_2 == shell_num_4 ? 1.0 : 2.0) : 2.0;
											int s1234_deg = s12_deg * s34_deg * s12_34_deg;

											const int how_many = howMany[contIndex];
											shells_by_cont[contIndex][how_many] = tempShell;
											firstBasisFuncts[contIndex][how_many] = {bf1_first, bf2_first, bf3_first, bf4_first};  //store shell number based on primitive combs
											s1234_degen_vec[contIndex][how_many] = s1234_deg;
											howMany[contIndex]++;  
											++count;
											LIBINT2_REALTYPE* prim_ints = NULL;
											//cout<<"how many p: "<<howMany[contIndex]<<endl;  //something wrong here

											if(howMany[contIndex] == vectorLength){          //if vecotor length of them are ready to be computed do so, else keep looking
												howMany[contIndex] = 0;
												preplibint2<Libint_t> (Libint_, shells_by_cont[contIndex], {{am1, am2, am3, am4}});   //TODO fix the am entry could declare vector much earlier
												if(am_tot){
													LIBINT2_PREFIXED_NAME(libint2_build_eri)[am1][am2][am3][am4](&Libint_[0]);
													prim_ints = Libint_[0].targets[0];
												}
												else{  // (ss|ss) type integral libint does not compute
													LIBINT2_REALTYPE ssss = 0.0;
													for(int p=0; p<numContrComb; ++p){
														ssss += Libint_[p].LIBINT_T_SS_EREP_SS(0)[0];
													}
													prim_ints = &ssss;
													//TODO there will only be one so add to Fock immediately
												}

												auto ss_1 = shells_by_cont[contIndex][0];
												auto ss_2 = shells_by_cont[contIndex][1];
 
												if(ss_1[0]->getContDepth() != ss_2[0]->getContDepth()){cout<<"broken"<<endl;}
												if(ss_1[1]->getContDepth() != ss_2[1]->getContDepth()){cout<<"broken"<<endl;}
												if(ss_1[2]->getContDepth() != ss_2[2]->getContDepth()){cout<<"broken"<<endl;}
												if(ss_1[3]->getContDepth() != ss_2[3]->getContDepth()){cout<<"broken"<<endl;}

												for(auto f1 = 0, f1234 = 0; f1 != numPrim_1; f1++){        //compute and get data out
													for(auto f2 = 0; f2 != numPrim_2; f2++){
														for(auto f3 = 0; f3 != numPrim_3; f3++){
															for(auto f4 = 0; f4!=numPrim_4; f4++, f1234++){
																auto value = prim_ints[f1234];
																double result[vectorLength];          //TODO may need to align stack buffer for performance
																value.convert(result);

																for(int vec = 0; vec < vectorLength; vec++){
																	const int bf1 = f1 + firstBasisFuncts[contIndex][vec][0];  //TODO fix this
																	const int bf2 = f2 + firstBasisFuncts[contIndex][vec][1];
																	const int bf3 = f3 + firstBasisFuncts[contIndex][vec][2];
																	const int bf4 = f4 + firstBasisFuncts[contIndex][vec][3];
		//															cout<<"first index         "<<firstBasisFuncts[contIndex][vec][0]<<" "
		//																						<<firstBasisFuncts[contIndex][vec][1]<<" "
		//																						<<firstBasisFuncts[contIndex][vec][2]<<" "
		//																						<<firstBasisFuncts[contIndex][vec][3]<<endl;
		//															cout<<"(f1,f2,f3,f4)      (" << f1 <<","<< f2 <<","<< f3 <<","<< f4 <<")"<<endl;
		//															cout<<"(bf1,bf2,bf3,bf4)  (" << bf1 <<","<<bf2 <<","<<bf3 <<","<<bf4 <<")"<<endl;

																	//cout<<"\tvalue: "<<result[vec]<<endl;
																	const auto val_scaled = result[vec] * s1234_degen_vec[contIndex][vec];

																	Fock(bf1 , bf2) += D(bf3 , bf4) * val_scaled;
																	Fock(bf3 , bf4) += D(bf1 , bf2) * val_scaled;
																	Fock(bf1,bf3) -= 0.25 * D(bf2,bf4) * val_scaled;
																	Fock(bf2,bf4) -= 0.25 * D(bf1,bf3) * val_scaled;
																	Fock(bf1,bf4) -= 0.25 * D(bf2,bf3) * val_scaled;
																	Fock(bf2,bf3) -= 0.25 * D(bf1,bf4) * val_scaled;

																}
															}
														}
													}
												}
											}
										}
									}
								}
							}  
							for(int p = 0; p < maxContComb; p++){                 //look at how_many.  if any elements still have integrals to be computed, do so
								//cout<<"howMany"<<howMany[p]<<endl;
								if(howMany[p] != 0){
									cout<<"howMany[p] = "<<howMany[p]<<endl;
									cout << "\t\tcompute with " << vectorLength - howMany[p] << " fake entries." << endl;
									cout << "\t\tnum cont comb " << p <<endl;
									cout<<"\t\tnumber of integrals to compute: "<<howMany[p]<<endl;
									const short numRelevant = howMany[p];
									howMany[p] = 0;
									preplibint2<Libint_t> (Libint_, shells_by_cont[p] , {{am1, am2, am3, am4}});
									LIBINT2_REALTYPE* prim_ints = NULL;
									if(am_tot){
										LIBINT2_PREFIXED_NAME(libint2_build_eri)[am1][am2][am3][am4](&Libint_[0]);
										prim_ints = Libint_[0].targets[0];
										//cout<<prim_ints[0]<<endl;
									}
									else{
										auto s1 = shells_by_cont[p][0][0];
										auto s2 = shells_by_cont[p][0][1];
										auto s3 = shells_by_cont[p][0][2];
										auto s4 = shells_by_cont[p][0][3];
										numContrComb = s1->getContDepth() * s2->getContDepth() * s3->getContDepth() * s4->getContDepth();
										for(int cont = 0; cont < numContrComb; cont++){
											LIBINT2_REALTYPE ssss = 0.0;
											for(int p=0; p<numContrComb; ++p){
												ssss += Libint_[p].LIBINT_T_SS_EREP_SS(0)[0];
										}
											prim_ints = &ssss;
									}

									int bf1_first[vectorLength];
									int bf2_first[vectorLength];
									int bf3_first[vectorLength];
									int bf4_first[vectorLength];

									for(int vec = 0; vec < numRelevant; vec++){
										const auto s1 = shells_by_cont[p][vec][0];
										const auto s2 = shells_by_cont[p][vec][1];
										const auto s3 = shells_by_cont[p][vec][2];
										const auto s4 = shells_by_cont[p][vec][3];
										auto s12_deg = (s1->shellNum == s2->shellNum) ? 1.0 : 2.0;
										auto s34_deg = (s3->shellNum == s4->shellNum) ? 1.0 : 2.0;
										auto s12_34_deg = (s1->shellNum == s3->shellNum) ? (s2->shellNum == s4->shellNum ? 1.0 : 2.0) : 2.0;
										s1234_degen_vec[p][vec] = s12_deg * s34_deg * s12_34_deg;
										bf1_first[vec] = shell_to_bf[s1->shellNum];
										bf2_first[vec] = shell_to_bf[s2->shellNum];
										bf3_first[vec] = shell_to_bf[s3->shellNum];
										bf4_first[vec] = shell_to_bf[s4->shellNum];
									}

									auto ss_1 = shells_by_cont[p][0];
									auto ss_2 = shells_by_cont[p][0];


									for(auto f1 = 0, f1234 = 0; f1 != numPrim_1; f1++){
										for(auto f2 = 0; f2 != numPrim_2; f2++){
											for(auto f3 = 0; f3 != numPrim_3; f3++){
												for(auto f4 = 0; f4!=numPrim_4; f4++, f1234++){
													auto value = prim_ints[f1234];
													double result[vectorLength];  //TODO may need to align stack buffer for performance
													value.convert(result);
													for(int vec = 0; vec< numRelevant; vec++){
														//cout<<"hey "<<s1->shellNum<<" "<<s2->shellNum<<" "<<s3->shellNum<<" "<<s4->shellNum<<endl;
														//cout<<result[vec]<<endl;
														const auto bf1 = f1 + bf1_first[vec];
														const auto bf2 = f2 + bf2_first[vec];
														const auto bf3 = f3 + bf3_first[vec];
														const auto bf4 = f4 + bf4_first[vec];

														const auto val_scaled = result[vectorLength - vec - 1] * s1234_degen_vec[p][vec];

														Fock(bf1 , bf2) += D(bf3 , bf4) * val_scaled;
														Fock(bf3 , bf4) += D(bf1 , bf2) * val_scaled;
														Fock(bf1,bf3) -= 0.25 * D(bf2,bf4) * val_scaled;
														Fock(bf2,bf4) -= 0.25 * D(bf1,bf3) * val_scaled;
														Fock(bf1,bf4) -= 0.25 * D(bf2,bf3) * val_scaled;
														Fock(bf2,bf3) -= 0.25 * D(bf1,bf4) * val_scaled;

													}
												}
											}
										}
									}//end of prim loop
								}
							}
						}
					}// end of if(am1+am2<=am3+am4)
					}
				}
			}


		}
}//end of am loop

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"time for First: " <<duration <<std::endl<<std::endl;

    Eigen::MatrixXd Fock_T = Fock.transpose();
    Fock = 0.5*(Fock+Fock_T);
    cout<<Fock<<endl;
    Eigen::EigenSolver <Eigen::MatrixXd> es(Fock);
    cout<<"eigen Vals"<<es.eigenvalues()<<endl;


    return 0;



}
