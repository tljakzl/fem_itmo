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


bool Node::operator<(const Node& val) const {
	return p < val.p;
}