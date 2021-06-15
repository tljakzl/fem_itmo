#include "Node.h"



Node::Node(int globalID, Point p) {
	this->globalID = globalID;
	this->p = p;
}

Node::Node(Point p) {
	this->p = p;
}

Node::Node() {

}

void Node::setBcond1(std::function<double(Point)> f) {
	this->hasBcond1 = true;
	this->bcond1 = f;
}

bool Node::operator<(const Node& val) const {
	return p < val.p && priority < val.priority;
}