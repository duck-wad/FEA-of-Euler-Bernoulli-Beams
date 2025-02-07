#include <fstream>

#include "Mesh.h"

Element::Element() {
	std::cout << "FUCK";
}

Mesh::Mesh() {
	ReadFile("INPUT.txt");


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
	
	//read element connectivity
	infile >> junk;
	infile >> junk >> numelem;
	elements.resize(numelem);
	//first entry in each connectivity vector is start, second is end
	connectivity.resize(numelem, std::vector<int>(2));
	
	for (size_t i = 0; i < numelem; i++) {
		infile >> junk >> junk >> connectivity[i][0] >> junk >> connectivity[i][1];
	}

	//read the BCs
	infile >> junk;
	infile >> junk >> numbcs;
	std::string tempstring;
	boundaryConditions.resize(numbcs);
	for (size_t i = 0; i < numbcs; i++) {
		infile >> junk >> junk >> tempstring >> junk >> boundaryConditions[i].node;
		if (tempstring == "pin") {
			boundaryConditions[i].type = BCType::PIN;
		}
		else if (tempstring == "roller") {
			boundaryConditions[i].type = BCType::ROLLER;
		}
		else if (tempstring == "clamp") {
			boundaryConditions[i].type = BCType::CLAMP;
		}
		else throw std::invalid_argument("Not a valid boundary condition");
	}

	//read the loads
	infile >> junk;
	infile >> junk >> numloads;
	loads.resize(numloads);
	for (size_t i = 0; i < numloads; i++) {
		infile >> junk >> junk >> tempstring >> junk >> loads[i].startNode >> junk;
		if (junk == "end:") {
			infile >> loads[i].endNode >> junk >> loads[i].magnitude;
		}
		else if (junk == "magnitude:") {
			infile >> loads[i].magnitude;
		}

		if (tempstring == "force") {
			loads[i].type = LoadType::FORCE;
		}
		else if (tempstring == "moment") {
			loads[i].type = LoadType::MOMENT;
		}
		else if (tempstring == "distributed") {
			loads[i].type = LoadType::DISTRIBUTED;
		}
		else throw std::invalid_argument("Not a valid load condition");
	}
}