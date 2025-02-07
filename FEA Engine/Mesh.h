#pragma once

#include <iostream>
#include <vector>

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
	double magnitude;
};

class Element {
public:
	Element();

protected:
	std::vector<std::vector<double>> elementStiffness;
	//if distributed load is applied to an element it will go into the force vector
	std::vector<double> elementForce;
};

class Mesh
{
public: 
	Mesh();

	void ReadFile(std::string fileName);
	void Discretize();

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

	//global stiffness and force
	std::vector<std::vector<double>> globalStiffness;
	std::vector<double> globalForce;
};

