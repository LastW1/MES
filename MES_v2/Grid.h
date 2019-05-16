#pragma once
#include "Grid.h"
#include <vector>
#include "Element.h"
#include "Node.h"
#include "Data.h"

using namespace std;

class Grid
{
vector<Node> gridNodes;
vector<Element> gridElements;


public:

	void show_coordinate(long int);
	void show(Data*);
	Grid();
	~Grid();
	void universal(Data*);
	void matrixH(Data*);
	void zeros(Data*);
	void calculate_HB(long double [][4], long double, long double, long double, long double,long double, long double,long double[],long double);
	long double*gauss(int, long double **, long double *);
	const double epsilon = 1e-11;
};

