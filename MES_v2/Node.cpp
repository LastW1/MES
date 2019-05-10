#include "Node.h"



Node::Node()
{
}

Node::Node(long double x1, long double y1, int i, bool b) {

	x = x1;
	y = y1;
	id = i;
	boundary = b;

}

Node::~Node()
{
}
