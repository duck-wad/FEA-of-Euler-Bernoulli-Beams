#include <fstream>

#include "Mesh.h"
#include "Utils.h"
#include "MatrixSolver.h"

//constructor makes the element stiffness matrix using 2-point gauss quadrature
Element::Element(double L, double E, double I, int start, int end){

	Length = L;
	startNode = start;
	endNode = end;
	elemE = E;
	elemI = I;
	
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

	elementStiffness *= (elemE * elemI * Jacobian);
}

void Element::ConstructForce(double w1, double w2) {
	double temp1 = w1 * Length / 60.0;
	double temp2 = w2 * Length / 60.0;

	std::vector<double> vec1 = { 21.0, 3.0 * Length, 9.0, -2.0 * Length };
	std::vector<double> vec2 = { 9.0, 2.0 * Length, 21.0, -3.0 * Length };
	elementForce = vec1 * temp1 + vec2 * temp2;

	hasLoad = true;
}

Mesh::Mesh(std::string infile) {
	ReadFile(infile);

	Discretize();
	Assemble();
	ApplyBCs();

	displacements = GaussSeidel(globalStiffnessBC, globalForceBC);
	printVector(displacements);

	SolveReactions();

	CalculateMoment();
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

	//temporary variables
	int startNode;
	int endNode;
	double startMag;
	double endMag;

	for (size_t i = 0; i < numloads; i++) {

		infile >> junk >> junk >> tempstring;

		if (tempstring == "distributed") {
			infile >> junk >> startNode >> junk >> endNode >> junk >> startMag >> junk >> endMag;

			distributedLoads.emplace_back(Load{ LoadType::DISTRIBUTED, startNode, endNode, startMag, endMag });
		}
		else if (tempstring == "force" || tempstring == "moment") {
			infile >> junk >> startNode >> junk >> startMag;
			if (tempstring == "force") {
				pointLoads.emplace_back(Load{ LoadType::FORCE, startNode, -1, startMag, 0 });
			}
			else {
				pointLoads.emplace_back(Load{ LoadType::MOMENT, startNode, -1, startMag, 0 });
			}
		}
		else {
			throw std::invalid_argument("Not a valid load type");
		}
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
		elements.emplace_back(Element(length, E, I, connectivity[i][0], connectivity[i][1]));
	}

	//loop through distributed force list. only apply the distributed loads to elemental force. apply point loads to global force
	//check if the startNode and endNode match with any element in connectivity list
	//if yes pass in start and end magnitude
	for (const auto& load : distributedLoads) {

		for (size_t i = 0; i < connectivity.size(); i++) {
			int eStart = connectivity[i][0];
			int eEnd = connectivity[i][1];

			if (eStart == load.startNode && eEnd == load.endNode) {
				elements[i].ConstructForce(load.startMagnitude, load.endMagnitude);
				break;
			}
			else if (eStart == load.endNode && eEnd == load.startNode) {
				elements[i].ConstructForce(load.endMagnitude,load.startMagnitude);
				break;
			}
		}
		
	}

}

//assembles the global stiffness matrix and force vector
//point loads are applied directly to global force vector
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
	globalForce.resize(2 * maxnode);
	for (size_t i = 0; i < numelem; i++) {
		if (elements[i].hasLoad == true) {
			int globalDOF1 = 2*(connectivity[i][0] - 1);
			int globalDOF2 = 2*(connectivity[i][1] - 1);

			globalForce[globalDOF1] += (elements[i].GetForce())[0];
			globalForce[globalDOF1 + 1] += (elements[i].GetForce())[1];
			globalForce[globalDOF2] += (elements[i].GetForce())[2];
			globalForce[globalDOF2 + 1] += (elements[i].GetForce())[3];
		}
	}

	//now put point loads into the force vector
	for (size_t i = 0; i < pointLoads.size(); i++) {
		int globalDOF = 2*(pointLoads[i].startNode - 1);
		if (pointLoads[i].type == LoadType::MOMENT) {
			globalDOF += 1;
		}
		globalForce[globalDOF] += pointLoads[i].startMagnitude;
	}

	//printVector(globalForce);
	//writeMatrixToCSV(globalStiffness, "GLOBAL_STIFFNESS.csv");
	//writeVectorToCSV(globalForce, "GLOBAL_FORCE.csv");
}

//for each node that has a support, set the row in stiffness matrix to zero
//except the diagonal which is set to 1
//set the entry in force vector equal to zero
//pin and roller have the vertical displacement set to zero
//clamp has displacement zero and slope zero
void Mesh::ApplyBCs() {

	globalStiffnessBC = globalStiffness;
	globalForceBC = globalForce;

	for (size_t i = 0; i < boundaryConditions.size(); i++) {
		int globalDOF = 2 * (boundaryConditions[i].node - 1);
		for (size_t i = 0; i < globalStiffnessBC.size(); i++) {
			if (i == globalDOF) {
				globalStiffnessBC[globalDOF][i] = 1.0;
			}
			else {
				globalStiffnessBC[globalDOF][i] = 0.0;
			}
		}
		globalForceBC[globalDOF] = 0.0;

		if (boundaryConditions[i].type == BCType::CLAMP) {
			for (size_t i = 0; i < globalStiffnessBC.size(); i++) {
				if (i == (globalDOF + 1)) {
					globalStiffnessBC[globalDOF + 1][i] = 1.0;
				}
				else {
					globalStiffnessBC[globalDOF + 1][i] = 0.0;
				}
			}
			globalForceBC[globalDOF + 1] = 0.0;
		}
	}
}

void Mesh::SolveReactions() {
	std::vector<double> tempReactions = globalStiffness * displacements - globalForce;

	//depending on fineness of mesh, some non-support nodes may have non-zero reactions
	//loop through the boundary conditions, set all non-support node reactions to zero
	reactions.resize(tempReactions.size(), 0.0);

	for (size_t i = 0; i < boundaryConditions.size(); i++) {
		int globalDOF = 2 * (boundaryConditions[i].node - 1);
		
		reactions[globalDOF] = tempReactions[globalDOF];
		
		if (boundaryConditions[i].type == BCType::CLAMP) {
			reactions[globalDOF+1] = tempReactions[globalDOF+1];
		}
	}
	printVector(reactions);
}

//moments are related to displacement by M=EI(d^2w/dx2)
//displacement within element approx by w(x)=N1(x)w1+N2(x)theta1... (cubic polynomial)
//take derivative of the shape functions wrt x
void Mesh::CalculateMoment() {
	
	//each element contributes some moment to the global node moment
	//keep track of how many elements contribute in momentCount then take average
	//since continuity is not enforced for moment across elements this is required to have a smooth BMD
	moments.resize(maxnode, 0.0);
	std::vector<int> momentCount(maxnode, 0);

	for (size_t i = 0; i < elements.size(); i++) {
		int startNode = elements[i].GetStart();
		int endNode = elements[i].GetEnd();
		double L = elements[i].GetLength();
		double EI = elements[i].GetE() * elements[i].GetI();

		int globalDOF1 = 2 * (startNode - 1);
		int globalDOF2 = 2 * (endNode - 1);

		double w1 = displacements[globalDOF1];
		double theta1 = displacements[globalDOF1 + 1];
		double w2 = displacements[globalDOF2];
		double theta2 = displacements[globalDOF2 + 1];

		//moment at left node x=0, moment at right node x=L
		moments[startNode-1] += EI * (-6.0 / (L * L) * w1 - 4.0 / L * theta1 + 6.0 / (L * L) * w2 - 2.0 / L * theta2);
		moments[endNode-1] += EI * (6.0 / (L * L) * w1 + 2.0 / L * theta1 - 6.0 / (L * L) * w2 + 4.0 / L * theta2);

		momentCount[startNode - 1]++;
		momentCount[endNode - 1]++;
	}

	//now average each moment by the corresponding momentCount
	for (size_t i = 0; i < maxnode; i++) {
		if (momentCount[i] > 0) {
			moments[i] /= momentCount[i];
		}
	}
	printVector(moments);
	printVector(momentCount);
}