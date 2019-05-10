#pragma once
#include "Node.h"

class Element
{
public:
	int id;
	long double matrixH[4][4];
	int node[4];
	bool boundary[4];
	Element();
	Element(int, int ,int ,int , int);
	~Element();
};

