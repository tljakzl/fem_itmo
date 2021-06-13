#pragma once
#include <vector>
#include <map>
#include <algorithm>
#include "Node.h"
#include "FE.h"
#include "DMat.h"

class Grid
{
public:
	std::vector<Node*> allNodes;
	std::vector<FE*> finiteElements;


	/// <summary>
	/// Передаем координаты сетки, привязку и гу создается сетка, создаются КЭ
	/// </summary>
	/// <param name="points">Массив координат узлов (геометрически предсортированный)</param>
	/// <param name="rectangles">Привязка КЭ к номеру узла</param>
	Grid(std::vector<Point> points, 
		std::map<int, std::vector<int>> rectangles, 
		std::map<int, int> bcond1, 
		std::map<int, int> rp, 
		std::map<int, int> lambda, 
		std::map<int, int> gamma, 
		std::map<std::pair<int, int>, int> bcond2, 
		std::map<std::pair<int, int>, int> bcond3);
	~Grid();
	std::vector<Point> scaffoldFEPoints(std::vector<Point> points);
	

	void calculateLocalMatrices();
	DMat<double> buildGlobalMatrix();
	std::vector<double> buildGlobalRP();
};

