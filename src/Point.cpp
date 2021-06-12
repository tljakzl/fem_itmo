#include "Point.h"

Point::Point(double x, double y) {
	this->x = x;
	this->y = y;
}

Point::Point() {
	this->x = 0;
	this->y = 0;
}

bool Point::operator<(const Point& val) const {
	if (y != val.y) {
		return y < val.y;
	}
	return x < val.x;
}