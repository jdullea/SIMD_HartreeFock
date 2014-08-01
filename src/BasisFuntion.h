/*
 * BasisFuntion.h
 *
 *  Created on: Jun 23, 2014
 *      Author: jon
 *
 *
 *  */

#ifndef BASISFUNTION_H_
#define BASISFUNTION_H_

#include <vector>
#include <iostream>

/* A class that defines a basis function*/

class BasisFunction
{
public:

	BasisFunction(std::vector <double> &exponents,
			std::vector <double> &contCoefficents,
			std::vector<double> &coordinates, int angularMomentum);

	std::vector <double> exp;          //orbital exponents
	std::vector <double> contCoeff;          //contraction coefficients
	std::vector <double> coord;                   //center on which the function sits
	int angMo;

};



#endif /* BASISFUNTION_H_ */
