#include <fstream>

#include "Mesh.h"
#include "Utils.h"

//constructor makes the element stiffness matrix using 2-point gauss quadrature
Element::Element(double L, double E, double I){

	Length = L;
	
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
}

void Element::ConstructForce(double w1, double w2) {
	double temp1 = w1 * Length / 60.0;
	double temp2 = w2 * Length / 60.0;

	std::vector<double> vec1 = { 21.0, 3.0 * Length, 9.0, -2.0 * Length };
	std::vector<double> vec2 = { 9.0, 2.0 * Length, 21.0, -3.0 * Length };
	elementForce = vec1 * temp1 + vec2 * temp2;

	printVector(elementForce);
}

Mesh::Mesh() {
	ReadFile("INPUT.txt");

	Discretize();
	Assemble();
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
			infile >> loads[i].endNode >> junk >> loads[i].startMagnitude >> junk >> loads[i].endMagnitude;
		}
		else if (junk == "magnitude:") {
			infile >> loads[i].startMagnitude;
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

	//loop through force list. only apply the distributed loads to elemental force. apply point loads to global force
	//check if the startNode and endNode match with any element in connectivity list
	//if yes pass in start and end magnitude
	for (const auto& load : loads) {
		if (load.type == LoadType::DISTRIBUTED) {
			for (size_t i = 0; i < connectivity.size(); i++) {
				int eStart = connectivity[i][0];
				int eEnd = connectivity[i][1];

				if (eStart == load.startNode && eEnd == load.endNode) {
					elements[i].ConstructForce(load.startMagnitude, load.endMagnitude);
					break;
				}
				else if (eStart == load.endNode && eEnd == load.startNode) {
					elements[i].ConstructForce(load.endMagnitude, load.startMagnitude);
					break;
				}
			}
		}
	}

}

void Mesh::Assemble() {
	globalStiffness.resize(2 * maxnode, std::vector<double>(2 * maxnode));
	for (size_t i = 0; i < numelem; i++) {
		int node1 = connectivity[i][0];
		int node2 = connectivity[i][1];

		std::vector<std::vector<double>> ke = elements[i].GetStiffness();

		std::vector<int> globalDOFs = {
			2 * (node1-1), 2 * (node1-1) + 1,
			2 * (node2-1), 2 * (node2-1) + 1
		};

		for (size_t j = 0; j < 4; j++) {
			for (size_t k = 0; k < 4; k++) {
				globalStiffness[globalDOFs[j]][globalDOFs[k]] += ke[j][k];
			}
		}
	}
	
	//assemble global force vector from just distributed loads

}