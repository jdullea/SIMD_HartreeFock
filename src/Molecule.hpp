/*
 * Molecule.hpp
 *
 *  Created on: Jun 23, 2014
 *      Author: jon
 */

#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

#include <vector>

class Molecule
{
public:
	Molecule(std::vector < std::vector <double> > coordinates, std::vector atomicNumbers, int numberAtoms);
	Molecule();
	std::vector < std::vector < double> > coord;
	std::vector <int>::atomicNums;
	int numAtoms;

	void getBasisFunctions();
};



#endif /* MOLECULE_HPP_ */
