#pragma once
#include "Node.h"
#include <vector>
#include <functional>
#include <map>

class FE
{
public:
	int id;
	int nodeCount;
	std::vector<Node*> nodes;

	double approximateFunc(std::function<double(Point)>);
	std::vector<std::vector<double>> _A, _RP, _G, _M, _b, _As3x, _As3y, _bs3x, _bs3y, _bs2x, _bs2y;
	std::vector<std::vector<double>> calcA(); // локальная матрица А (A = G + M + As3x + As3y)
	std::vector<std::vector<double>> calcRP();  // локальный вектор правой части (RP = b + bs3x + bs3y + bs2x + bs2y)

	std::vector<std::vector<double>> calcG();
	std::vector<std::vector<double>> calcM();
	std::vector<std::vector<double>> calcb();
	std::vector<std::vector<double>> calcAs3x();
	std::vector<std::vector<double>> calcAs3y();
	std::vector<std::vector<double>> calcbs3x();
	std::vector<std::vector<double>> calcbs3y();
	std::vector<std::vector<double>> calcbs2x();
	std::vector<std::vector<double>> calcbs2y();

	void makeImpactToGlobalMatrix(std::map<int, std::vector<std::pair<int, double>>>&);
	void makeImpactToGlobalRP(std::vector<double>&);

	
	double hx, hy;
	std::function<double(Point)> _lambda, _gamma, _rp, _bcond2x_rp, _bcond2y_rp, _bcond3x, _bcond3y, _bcond3x_rp, _bcond3y_rp;

	FE(std::vector<Node*> nodes,
		int id,
		std::function<double(Point)> lambda,
		std::function<double(Point)> gamma,
		std::function<double(Point)> rp,
		std::function<double(Point)> bcond2x_rp,
		std::function<double(Point)> bcond2y_rp,
		std::function<double(Point)> bcond3x,
		std::function<double(Point)> bcond3y,
		std::function<double(Point)> bcond3x_rp,
		std::function<double(Point)> bcond3y_rp
	);
};

