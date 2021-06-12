#pragma once
class Point
{

public:
	double x;
	double y;
	Point();
	Point(double x, double y);
	bool Point::operator<(const Point& val) const;

};

