#include <fstream>

#include "Mesh.h"
#include "Utils.h"
#include "MatrixSolver.h"

//constructor makes the element stiffness matrix using 2-point gauss quadrature
Element::Element(double L, double E, double I, double A, int start, int end){

	Length = L;
	startNode = start;
	endNode = end;
	elemE = E;
	elemI = I;
	elemA = A;
	
	elementStiffness.resize(6, std::vector<double>(6));
	elementForce.resize(6);

	//have temp matrices for the components of element stiffness for transverse and axial
	std::vector<std::vector<double>> transverseStiffness(4, std::vector<double>(4, 0.0));
	std::vector<std::vector<double>> axialStiffness(2, std::vector<double>(2, 0.0));

	//construct the transverse part of the element stiffness
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
		transverseStiffness += (outerProduct(B[i], B[i]) * gaussweight2[i]);
	}

	transverseStiffness *= (elemE * elemI * Jacobian);

	//construct the axial part of element stiffness. Just define it easily as AE/L
	axialStiffness[0][0] = elemE * elemA / L;
	axialStiffness[0][1] = -elemE * elemA / L;
	axialStiffness[1][0] = -elemE * elemA / L;
	axialStiffness[1][1] = elemE * elemA / L;

	//combine into elementStiffness
	//assemble the transverse part first
	int skiprow = 1;
	for (size_t i = 0; i < 4; i++) {
		//skip the first and third column/row which should only contain axial terms
		if (i == 2) {
			skiprow = 2;
		}
		int skipcol = 1;
		for (size_t j = 0; j < 4; j++) {

			if (j == 2) {
				skipcol = 2;
			}
			elementStiffness[i+skiprow][j+skipcol] = transverseStiffness[i][j];
		}
	}
	//now assemble the axial part, filling only first/third column/row
	elementStiffness[0][0] = axialStiffness[0][0];
	elementStiffness[0][3] = axialStiffness[0][1];
	elementStiffness[3][0] = axialStiffness[1][0];
	elementStiffness[3][3] = axialStiffness[1][1];
}

void Element::ConstructForce(double w1, double w2) {
	double temp1 = w1 * Length / 60.0;
	double temp2 = w2 * Length / 60.0;

	//first and fourth DOF are for axial and therefore are zero for a vertical dist. load
	std::vector<double> vec1 = { 0.0, 21.0, 3.0 * Length, 0.0, 9.0, -2.0 * Length };
	std::vector<double> vec2 = { 0.0, 9.0, 2.0 * Length, 0.0, 21.0, -3.0 * Length };
	elementForce = vec1 * temp1 + vec2 * temp2;

	hasLoad = true;
}

Mesh::Mesh(std::string infile) {
	ReadFile(infile);

	Discretize();
	Assemble();
	ApplyBCs();

	//displacements = GaussSeidel(globalStiffnessBC, globalForceBC, 0.001, 5000);
	displacements = GaussianElimination(globalStiffnessBC, globalForceBC);

	SolveReactions();
	//CalculateMoment();

	WriteResults();
}

void Mesh::ReadFile(std::string fileName) {

	std::string junk;

	std::ifstream infile(fileName);
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
	}

	infile >> junk;
	infile >> junk >> E >> junk >> I >> junk >> A;

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
		else if (tempstring == "vforce" || tempstring == "hforce" || tempstring == "moment") {
			infile >> junk >> startNode >> junk >> startMag;
			if (tempstring == "vforce") {
				pointLoads.emplace_back(Load{ LoadType::VFORCE, startNode, -1, startMag, 0 });
			}
			else if (tempstring == "hforce") {
				pointLoads.emplace_back(Load{ LoadType::HFORCE, startNode, -1, startMag, 0 });
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
		elements.emplace_back(Element(length, E, I, A, connectivity[i][0], connectivity[i][1]));
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
	//3 DOF per node for axial, transverse, rotation
	globalStiffness.resize(3 * maxnode, std::vector<double>(3 * maxnode));
	for (size_t i = 0; i < numelem; i++) {
		int node1 = connectivity[i][0];
		int node2 = connectivity[i][1];

		std::vector<std::vector<double>> ke = elements[i].GetStiffness();

		std::vector<int> globalDOFs = {
			3 * (node1-1), 3 * (node1-1) + 1, 3 * (node1-1)+2,
			3 * (node2-1), 3 * (node2-1) + 1, 3 * (node2-1)+2
		};

		for (size_t j = 0; j < 6; j++) {
			for (size_t k = 0; k < 6; k++) {
				globalStiffness[globalDOFs[j]][globalDOFs[k]] += ke[j][k];
			}
		}
	}
	
	//assemble global force vector from just distributed loads
	globalForce.resize(3 * maxnode);
	for (size_t i = 0; i < numelem; i++) {
		if (elements[i].hasLoad == true) {
			int globalDOF1 = 3*(connectivity[i][0] - 1);
			int globalDOF2 = 3*(connectivity[i][1] - 1);

			globalForce[globalDOF1] += (elements[i].GetForce())[0];
			globalForce[globalDOF1 + 1] += (elements[i].GetForce())[1];
			globalForce[globalDOF1 + 2] += (elements[i].GetForce())[2];
			globalForce[globalDOF2] += (elements[i].GetForce())[3];
			globalForce[globalDOF2 + 1] += (elements[i].GetForce())[4];
			globalForce[globalDOF2 + 2] += (elements[i].GetForce())[5];
		}
	}

	//now put point loads into the force vector
	for (size_t i = 0; i < pointLoads.size(); i++) {
		int globalDOF = 3*(pointLoads[i].startNode - 1);
		//first DOF for axial, second for transverse, third rotational
		if (pointLoads[i].type == LoadType::VFORCE) {
			globalDOF += 1;
		}
		else if (pointLoads[i].type == LoadType::MOMENT) {
			globalDOF += 2;
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
		int globalDOF = 3 * (boundaryConditions[i].node - 1);

		//for all BC types (pin, roller, clamp) the vertical is fixed
		for (size_t i = 0; i < globalStiffnessBC.size(); i++) {
			if (i == (globalDOF + 1)) {
				globalStiffnessBC[globalDOF+1][i] = 1.0;
			}
			else {
				globalStiffnessBC[globalDOF+1][i] = 0.0;
			}
		}
		globalForceBC[globalDOF+1] = 0.0;

		//pin BC fixes the horizontal movement
		if (boundaryConditions[i].type == BCType::PIN) {
			for (size_t i = 0; i < globalStiffnessBC.size(); i++) {
				if (i == globalDOF) {
					globalStiffnessBC[globalDOF][i] = 1.0;
				}
				else {
					globalStiffnessBC[globalDOF][i] = 0.0;
				}
			}
			globalForceBC[globalDOF] = 0.0;
		}

		//clamp fixes vertical, horizontal and rotation
		else if (boundaryConditions[i].type == BCType::CLAMP) {
			for (size_t i = 0; i < globalStiffnessBC.size(); i++) {
				if (i == globalDOF) {
					globalStiffnessBC[globalDOF][i] = 1.0;
				}
				else {
					globalStiffnessBC[globalDOF][i] = 0.0;
				}
			}
			globalForceBC[globalDOF] = 0.0;
			for (size_t i = 0; i < globalStiffnessBC.size(); i++) {
				if (i == (globalDOF + 2)) {
					globalStiffnessBC[globalDOF + 2][i] = 1.0;
				}
				else {
					globalStiffnessBC[globalDOF + 2][i] = 0.0;
				}
			}
			globalForceBC[globalDOF + 2] = 0.0;
		}
	}
}

void Mesh::SolveReactions() {
	std::vector<double> tempReactions = globalStiffness * displacements - globalForce;

	//depending on fineness of mesh, some non-support nodes may have non-zero reactions
	//loop through the boundary conditions, set all non-support node reactions to zero
	reactions.resize(tempReactions.size(), 0.0);

	for (size_t i = 0; i < boundaryConditions.size(); i++) {
		int globalDOF = 3 * (boundaryConditions[i].node - 1);
		//there will always be vertical reaction regardless of BC type
		reactions[globalDOF+1] = tempReactions[globalDOF+1];
		if (boundaryConditions[i].type == BCType::PIN) {
			reactions[globalDOF] = tempReactions[globalDOF];
		}
		
		else if (boundaryConditions[i].type == BCType::CLAMP) {
			reactions[globalDOF] = tempReactions[globalDOF];
			reactions[globalDOF+2] = tempReactions[globalDOF+2];
		}
	}
}

void Mesh::WriteResults() {
	std::vector<std::vector<double>> results;
	results.emplace_back(globalCoordinates);

	std::vector<double> axialTranslations;
	std::vector<double> transverseTranslations;
	std::vector<double> rotations;
	std::vector<double> horizontalReactions;
	std::vector<double> verticalReactions;
	std::vector<double> momentReactions;

	//split up dof related to horizontal translation, vertical translation, and rotation
	for (size_t i = 0; i < displacements.size(); i++) {
		if (i % 3 == 0) {
			axialTranslations.emplace_back(displacements[i]);
			horizontalReactions.emplace_back(reactions[i]);
		}
		else if ((i + 2) % 3 == 0) {
			transverseTranslations.emplace_back(displacements[i]);
			verticalReactions.emplace_back(reactions[i]);
		}
		else {
			rotations.emplace_back(displacements[i]);
			momentReactions.emplace_back(reactions[i]);
		}
	}

	results.emplace_back(axialTranslations);
	results.emplace_back(transverseTranslations);
	results.emplace_back(rotations);
	results.emplace_back(horizontalReactions);
	results.emplace_back(verticalReactions);
	results.emplace_back(momentReactions);

	results = transpose(results);

	std::vector<std::string> titles;
	titles.emplace_back("Node location (m)");
	titles.emplace_back("Nodal axial displacement (m)");
	titles.emplace_back("Nodal transverse displacement (m)");
	titles.emplace_back("Nodal rotation (rad)");
	titles.emplace_back("Horizontal force reactions (N)");
	titles.emplace_back("Vertical force reactions (N)");
	titles.emplace_back("Moment reactions (Nm)");

	writeMatrixToCSV(results, "./Output/RESULTS.csv", titles);
}