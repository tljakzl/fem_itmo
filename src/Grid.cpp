#include "Grid.h"
#include <set>
#include <Config.h>


Grid::Grid(std::vector<Point> points,
	std::map<int, std::vector<int>> rectangles,
	std::map<int, int> bcond1,
	std::map<int, int> rp,
	std::map<int, int> lambda,
	std::map<int, int> gamma,
	std::map<std::pair<int, int>, int> bcond2,
	std::map<std::pair<int, int>, int> bcond3) {


	// нужно для удаления повторяющихся нод в глобальном смысле
	std::map<std::pair<double, double>, Node*> uniqueNodes;

	for (auto rectangle_kvp : rectangles) {
		std::vector<Node*> FENodes;
		std::vector<Point> FEPoints;
		std::vector<Point> initialPoints;
		for (int p_i : rectangle_kvp.second) {
			initialPoints.push_back(points[p_i - 1]);
		}
		FEPoints = scaffoldFEPoints(initialPoints);
		
		for (auto p : FEPoints) {

			if (uniqueNodes.count({ p.x, p.y }) == 0) {
				Node* node = new Node(p);
				uniqueNodes[{p.x, p.y}] = node;
				
			}
			FENodes.push_back(uniqueNodes[{ p.x, p.y }]);

		}

		std::function<double(Point)> bcond3x = Config::GET_ZERO_FUNC(), 
									 bcond3y = Config::GET_ZERO_FUNC(), 
									 bcond3x_rp = Config::GET_ZERO_FUNC(),
									 bcond3y_rp = Config::GET_ZERO_FUNC(),
									 bcond2x_rp = Config::GET_ZERO_FUNC(),
									 bcond2y_rp = Config::GET_ZERO_FUNC();

		for (int p_i : rectangle_kvp.second) {
			for (int p_j : rectangle_kvp.second) {
				if (p_i != p_j) {
					if (bcond2.count({ p_i, p_j })) {
						if (points[p_i - 1].y - points[p_j - 1].y == 0) {  // горизонтальное ГУ
							bcond2x_rp = Config::GET_BCOND2_FUNC(bcond2[{p_i, p_j}]);
						}
						else {
							bcond2y_rp = Config::GET_BCOND2_FUNC(bcond2[{p_i, p_j}]);
						}
					}
					if (bcond3.count({ p_i, p_j })) {
						if (points[p_i - 1].y - points[p_j - 1].y == 0) {  // горизонтальное ГУ
							bcond3x_rp = Config::GET_BCOND3_RP_FUNC(bcond3[{p_i, p_j}]);
							bcond3x = Config::GET_BCOND3_FUNC(bcond3[{p_i, p_j}]);
						}
						else {
							bcond3y_rp = Config::GET_BCOND3_RP_FUNC(bcond3[{p_i, p_j}]);
							bcond3y = Config::GET_BCOND3_FUNC(bcond3[{p_i, p_j}]);
						}
					}
				}
			}
		}


		int _n_idx = rectangle_kvp.first;
		FE* fe = new FE(FENodes, 
						_n_idx,
						Config::GET_LAMBDA_FUNC(lambda[_n_idx]),
						Config::GET_GAMMA_FUNC(gamma[_n_idx]),
						Config::GET_RP_FUNC(rp[_n_idx]),
						bcond2x_rp,
						bcond2y_rp,
						bcond3x,
						bcond3y,
						bcond3x_rp,
						bcond3y_rp
						);
		finiteElements.push_back(fe);
	}

	for (auto unode_kvp : uniqueNodes) {
		allNodes.push_back(unode_kvp.second);
	}

	sort(allNodes.begin(), allNodes.end(), [](const Node* n1, const Node* n2) {
		return *n1 < *n2;
	});
	
	for (int g_i = 1; g_i <= allNodes.size(); g_i++) {
		allNodes[g_i - 1]->globalID = g_i;
	}

}

std::vector<Point> Grid::scaffoldFEPoints(std::vector<Point> points) {
	std::sort(points.begin(), points.end());

	double  max_x = DBL_MIN, min_x = DBL_MAX,
			max_y = DBL_MIN, min_y = DBL_MAX;

	for (auto& point : points) {
		max_x = std::max(max_x, point.x);
		min_x = std::min(min_x, point.x);
		max_y = std::max(max_y, point.y);
		min_y = std::min(min_y, point.y);
	}

	std::vector<Point> auxilaryPoints;

	for (auto p : points) auxilaryPoints.push_back(p);


	double  mid_x = min_x + (max_x - min_x) / 2,
			mid_y = min_y + (max_y - min_y) / 2;


	auxilaryPoints.push_back(
		*new Point(mid_x, min_y)
	);

	auxilaryPoints.push_back(
		*new Point(min_x, mid_y)
	);
	auxilaryPoints.push_back(
		*new Point(mid_x, mid_y)
	);
	auxilaryPoints.push_back(
		*new Point(max_x, mid_y)
	);

	auxilaryPoints.push_back(
		*new Point(mid_x, max_y)
	);

	return auxilaryPoints;

}

void Grid::calculateLocalMatrices() {
	for (auto& fe : finiteElements) {
		fe->calcA(); // матрица А
		fe->calcRP(); // вектор правой части
	}
}