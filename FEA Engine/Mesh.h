#pragma once

#include <iostream>
#include <vector>

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
	Mesh(std::string fileName);

	void ReadFile(std::string fileName);

protected:
	double E;
	double I;

	//coordinates is 1D vector since only x values
	std::vector<double> globalCoordinates;
	std::vector<std::vector<int>> connectivity;

	int numelem;
	int maxnode;

	//vector of Elements. each Element stores its own stiffness and force
	std::vector<Element> elements;

	//global stiffness and force
	std::vector<std::vector<double>> globalStiffness;
	std::vector<double> globalForce;
};

