#include "Element.h"



Element::Element()
{
}

Element::Element(int no1, int no2, int no3, int no4, int i) {
	node[0] = no1;
	node[1] = no2;
	node[2] = no3;
	node[3] = no4;
	id = i;


	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrixH[i][j] = 0;
	
}


Element::~Element()
{
}
