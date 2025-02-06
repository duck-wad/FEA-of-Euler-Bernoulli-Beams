#include <fstream>

#include "Mesh.h"

Mesh::Mesh(std::string fileName) {
	ReadFile(fileName);


}

void Mesh::ReadFile(std::string fileName) {

	std::string junk;

	std::ifstream infile(fileName);
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
	}

	infile >> junk;
	infile >> junk >> E >> junk >> I;

	//read node coordinates
	infile >> junk;
	infile >> junk >> maxnode;
	globalCoordinates.resize(maxnode, 0.0);

	for (size_t i = 0; i < maxnode; i++) {
		infile >> junk >> junk >> globalCoordinates[i];
	}
	
	std::cout << "test" << std::endl;


}