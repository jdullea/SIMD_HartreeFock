/*
 * BasisFunction.cpp
 *
 *  Created on: Jun 23, 2014
 *      Author: jon
 */

#include "BasisFuntion.h"
#include <vector>
#include <iostream>

BasisFunction::BasisFunction(std::vector <double> &exponents,
								std::vector <double> &contCoefficents,
								std::vector<double> &coordinates, int angularMomentum){
	exp=exponents;
	contCoeff=contCoefficents;
	coord=coordinates;
	angMo=angularMomentum;

}

