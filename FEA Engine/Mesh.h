#pragma once

#include <iostream>
#include <vector>

const std::vector<double> gausspoint2 = { -0.5773502692, 0.5773502692 };
const std::vector<double> gaussweight2 = { 1, 1 };

//store the BCs in a list of structs
enum class BCType {PIN, ROLLER, CLAMP};
struct BC {
	BCType type;
	int node;
};

//store the loads in a list of structs
enum class LoadType {FORCE, MOMENT, DISTRIBUTED};
struct Load {
	LoadType type;
	//for point loads, only startNode is used
	int startNode; 
	int endNode;
	double startMagnitude;
	double endMagnitude;
};

class Element {
public:
	Element(double L, double E, double I, int start, int end);
	void ConstructForce(double w1, double w2);
	std::vector<std::vector<double>> GetStiffness() { return elementStiffness; }
	std::vector<double> GetForce() { return elementForce; }

	bool hasLoad = false;

	double GetLength() { return Length; }
	int GetStart() { return startNode; }
	int GetEnd() { return endNode; }
	double GetE() { 
		return elemE; 
	}
	double GetI() { return elemI; }

protected:
	std::vector<std::vector<double>> elementStiffness;
	//if distributed load is applied to an element it will go into the force vector
	std::vector<double> elementForce;
	double Length;
	int startNode;
	int endNode;
	double elemE;
	double elemI;
};

class Mesh
{
public: 
	Mesh(std::string infile);

	void ReadFile(std::string fileName);
	void Discretize();
	void Assemble();
	void ApplyBCs();
	void SolveReactions();

	//post-processing functions
	void CalculateMoment();
	void CalculateShear();

protected:
	double E;
	double I;

	//coordinates is 1D vector since only x values
	std::vector<double> globalCoordinates;
	std::vector<std::vector<int>> connectivity;

	int numelem;
	int maxnode;
	int numbcs;
	int numloads;

	//vector of Elements. each Element stores its own stiffness and force
	std::vector<Element> elements;

	//vector of BCs and loads
	std::vector<BC> boundaryConditions;
	std::vector<Load> loads;

	std::vector<Load> pointLoads;
	std::vector<Load> distributedLoads;

	//global stiffness and force
	std::vector<std::vector<double>> globalStiffness;
	std::vector<double> globalForce;
	//since reactions are computed in post-processing, only alter the below matrix and vector for boundary conditions
	std::vector<std::vector<double>> globalStiffnessBC;
	std::vector<double> globalForceBC;

	//vector containing the resulting displacements/rotations
	std::vector<double> displacements;
	//vector with support reactions calculated in post-processing Kd - f
	std::vector<double> reactions;

	std::vector<double> moments;
	std::vector<double> shears;
};

