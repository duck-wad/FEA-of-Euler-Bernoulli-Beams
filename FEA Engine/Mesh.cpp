#include <fstream>

#include "Mesh.h"
#include "Utils.h"

//constructor makes the element stiffness matrix using 2-point gauss quadrature
Element::Element(double L, double E, double I) {
	
	elementStiffness.resize(4, std::vector<double>(4));
	elementForce.resize(4);

	double Jacobian = L / 2.0;

	//store the B vectors
	std::vector<std::vector<double>> B(gausspoint2.size(), std::vector<double>(4));

	for (size_t i = 0; i < B.size(); i++) {
		B[i][0] = 3.0 * gausspoint2[i] / 2.0;
		B[i][1] = L / 4.0 * (3.0 * gausspoint2[i] - 1.0);
		B[i][2] = -3.0 * gausspoint2[i] / 2.0;
		B[i][3] = L / 4.0 * (3.0 * gausspoint2[i] + 1.0);
		B[i] *= (4.0 / pow(L, 2));
	}

	for (size_t i = 0; i < B.size(); i++) {
		elementStiffness += (outerProduct(B[i], B[i]) * gaussweight2[i]);
	}

	elementStiffness *= (E * I * Jacobian);

	//printMatrix(elementStiffness);
}

void Element::ConstructForce() {
	std::cout << "no";
}

Mesh::Mesh() {
	ReadFile("INPUT.txt");

	Discretize();
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

//builds the elemental stiffness matrices
void Mesh::Discretize() {

	//loop through connectivity list. get the node ids, retrieve the node locations
	//find distance between nodes and set it as a temp length
	//pass the temp length into the Element constructor, along with E and I.

	for (size_t i = 0; i < numelem; i++) {
		//-1 because connectivity indices start at 1
		double startPos = globalCoordinates[connectivity[i][0] - 1];
		double endPos = globalCoordinates[connectivity[i][1] - 1];
		double length = abs(endPos - startPos);
		elements.emplace_back(Element(length, E, I));
	}

	//do something with the force vectors

}

