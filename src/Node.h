#pragma once
#include "Point.h"
#include <functional>

class Node
{
public:
	int globalID;
	Point p;
	bool hasBcond1=false;
	std::function<double(Point)> bcond1;

	Node(int globalID, Point p);
	Node(Point p);
	Node();

	bool Node::operator<(const Node& val) const;

	void setBcond1(std::function<double(Point)>);
};

