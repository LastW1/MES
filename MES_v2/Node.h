#pragma once
class Node
{
public:
	long double x;  //wsp x
	long double y;   //wsp y
	long int id; 
	bool boundary;
	long double LocalT;
	Node();
	Node(long double, long double,int,bool);
	~Node();
};

