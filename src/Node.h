#pragma once
#include "Point.h"

class Node
{
public:
	int globalID;
	Point p;

	Node(int globalID, Point p);
	Node(Point p);
	Node();

	bool Node::operator<(const Node& val) const;

};

