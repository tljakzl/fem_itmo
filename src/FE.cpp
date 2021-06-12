#include "FE.h"
#include <algorithm>


FE::FE(std::vector<Node*> nodes,
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
) : _lambda(lambda), _gamma(gamma), _rp(rp), _bcond2x_rp(bcond2x_rp), _bcond2y_rp(bcond2y_rp),
	_bcond3x(bcond3x), _bcond3y(bcond3y), _bcond3x_rp(bcond3x_rp), _bcond3y_rp(bcond3y_rp), id(id) {

	sort(nodes.begin(), nodes.end(), [](const Node* n1, const Node* n2) {
		return *n1 < *n2;
	});
	this->nodes = nodes;
	this->nodeCount = nodes.size();

	double  max_x = DBL_MIN, min_x = DBL_MAX,
		max_y = DBL_MIN, min_y = DBL_MAX;

	for (auto& node : nodes) {
		max_x = std::max(max_x, node->p.x);
		min_x = std::min(min_x, node->p.x);
		max_y = std::max(max_y, node->p.y);
		min_y = std::min(min_y, node->p.y);
	}

	hx = max_x - min_x;
	hy = max_y - min_y;
}

double FE::approximateFunc(std::function<double(Point)> func) {
	double res = 0;
	for (int i = 0; i < nodeCount; i++) {
		res += func(nodes[i]->p);
	}
	return res / nodeCount;
}


std::vector<std::vector<double>> FE::calcG() {
	std::vector<std::vector<double>> __g_1 = {
		{77 * _lambda(nodes[0]->p) + 21 * _lambda(nodes[2]->p) + 11 * _lambda(nodes[6]->p) + 3 * _lambda(nodes[8]->p), -84 * _lambda(nodes[0]->p) - 28 * _lambda(nodes[2]->p) - 12 * _lambda(nodes[6]->p) - 4 * _lambda(nodes[8]->p), 7 * _lambda(nodes[0]->p) + 7 * _lambda(nodes[2]->p) + _lambda(nodes[6]->p) + _lambda(nodes[8]->p), 44 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p), -48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p), -11 * _lambda(nodes[0]->p) - 3 * _lambda(nodes[2]->p) - 11 * _lambda(nodes[6]->p) - 3 * _lambda(nodes[8]->p), 12 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p)},
		{-84 * _lambda(nodes[0]->p) - 28 * _lambda(nodes[2]->p) - 12 * _lambda(nodes[6]->p) - 4 * _lambda(nodes[8]->p), 112 * _lambda(nodes[0]->p) + 112 * _lambda(nodes[2]->p) + 16 * _lambda(nodes[6]->p) + 16 * _lambda(nodes[8]->p), -28 * _lambda(nodes[0]->p) - 84 * _lambda(nodes[2]->p) - 4 * _lambda(nodes[6]->p) - 12 * _lambda(nodes[8]->p), -48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p), 64 * _lambda(nodes[0]->p) + 64 * _lambda(nodes[2]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[2]->p), 12 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p)},
		{7 * _lambda(nodes[0]->p) + 7 * _lambda(nodes[2]->p) + _lambda(nodes[6]->p) + _lambda(nodes[8]->p), -28 * _lambda(nodes[0]->p) - 84 * _lambda(nodes[2]->p) - 4 * _lambda(nodes[6]->p) - 12 * _lambda(nodes[8]->p), 21 * _lambda(nodes[0]->p) + 77 * _lambda(nodes[2]->p) + 3 * _lambda(nodes[6]->p) + 11 * _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[2]->p), 12 * _lambda(nodes[0]->p) + 44 * _lambda(nodes[2]->p), -_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -3 * _lambda(nodes[0]->p) - 11 * _lambda(nodes[2]->p) - 3 * _lambda(nodes[6]->p) - 11 * _lambda(nodes[8]->p)},
		{44 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p), -48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p), 176 * _lambda(nodes[0]->p) + 48 * _lambda(nodes[2]->p) + 176 * _lambda(nodes[6]->p) + 48 * _lambda(nodes[8]->p), -192 * _lambda(nodes[0]->p) - 64 * _lambda(nodes[2]->p) - 192 * _lambda(nodes[6]->p) - 64 * _lambda(nodes[8]->p), 16 * _lambda(nodes[0]->p) + 16 * _lambda(nodes[2]->p) + 16 * _lambda(nodes[6]->p) + 16 * _lambda(nodes[8]->p), 44 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -48 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p)},
		{-48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p), 64 * _lambda(nodes[0]->p) + 64 * _lambda(nodes[2]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[2]->p), -192 * _lambda(nodes[0]->p) - 64 * _lambda(nodes[2]->p) - 192 * _lambda(nodes[6]->p) - 64 * _lambda(nodes[8]->p), 256 * _lambda(nodes[0]->p) + 256 * _lambda(nodes[2]->p) + 256 * _lambda(nodes[6]->p) + 256 * _lambda(nodes[8]->p), -64 * _lambda(nodes[0]->p) - 192 * _lambda(nodes[2]->p) - 64 * _lambda(nodes[6]->p) - 192 * _lambda(nodes[8]->p), -48 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 64 * _lambda(nodes[6]->p) + 64 * _lambda(nodes[8]->p), -16 * _lambda(nodes[6]->p) - 48 * _lambda(nodes[8]->p)},
		{4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[2]->p), 12 * _lambda(nodes[0]->p) + 44 * _lambda(nodes[2]->p), 16 * _lambda(nodes[0]->p) + 16 * _lambda(nodes[2]->p) + 16 * _lambda(nodes[6]->p) + 16 * _lambda(nodes[8]->p), -64 * _lambda(nodes[0]->p) - 192 * _lambda(nodes[2]->p) - 64 * _lambda(nodes[6]->p) - 192 * _lambda(nodes[8]->p), 48 * _lambda(nodes[0]->p) + 176 * _lambda(nodes[2]->p) + 48 * _lambda(nodes[6]->p) + 176 * _lambda(nodes[8]->p), 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -16 * _lambda(nodes[6]->p) - 48 * _lambda(nodes[8]->p), 12 * _lambda(nodes[6]->p) + 44 * _lambda(nodes[8]->p)},
		{-11 * _lambda(nodes[0]->p) - 3 * _lambda(nodes[2]->p) - 11 * _lambda(nodes[6]->p) - 3 * _lambda(nodes[8]->p), 12 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p), 44 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -48 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), 11 * _lambda(nodes[0]->p) + 3 * _lambda(nodes[2]->p) + 77 * _lambda(nodes[6]->p) + 21 * _lambda(nodes[8]->p), -12 * _lambda(nodes[0]->p) - 4 * _lambda(nodes[2]->p) - 84 * _lambda(nodes[6]->p) - 28 * _lambda(nodes[8]->p), _lambda(nodes[0]->p) + _lambda(nodes[2]->p) + 7 * _lambda(nodes[6]->p) + 7 * _lambda(nodes[8]->p)},
		{12 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -48 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 64 * _lambda(nodes[6]->p) + 64 * _lambda(nodes[8]->p), -16 * _lambda(nodes[6]->p) - 48 * _lambda(nodes[8]->p), -12 * _lambda(nodes[0]->p) - 4 * _lambda(nodes[2]->p) - 84 * _lambda(nodes[6]->p) - 28 * _lambda(nodes[8]->p), 16 * _lambda(nodes[0]->p) + 16 * _lambda(nodes[2]->p) + 112 * _lambda(nodes[6]->p) + 112 * _lambda(nodes[8]->p), -4 * _lambda(nodes[0]->p) - 12 * _lambda(nodes[2]->p) - 28 * _lambda(nodes[6]->p) - 84 * _lambda(nodes[8]->p)},
		{-_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -3 * _lambda(nodes[0]->p) - 11 * _lambda(nodes[2]->p) - 3 * _lambda(nodes[6]->p) - 11 * _lambda(nodes[8]->p), 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -16 * _lambda(nodes[6]->p) - 48 * _lambda(nodes[8]->p), 12 * _lambda(nodes[6]->p) + 44 * _lambda(nodes[8]->p), _lambda(nodes[0]->p) + _lambda(nodes[2]->p) + 7 * _lambda(nodes[6]->p) + 7 * _lambda(nodes[8]->p), -4 * _lambda(nodes[0]->p) - 12 * _lambda(nodes[2]->p) - 28 * _lambda(nodes[6]->p) - 84 * _lambda(nodes[8]->p), 3 * _lambda(nodes[0]->p) + 11 * _lambda(nodes[2]->p) + 21 * _lambda(nodes[6]->p) + 77 * _lambda(nodes[8]->p)},
	};

	std::vector<std::vector<double>> __g_2 = {
		{77 * _lambda(nodes[0]->p) + 11 * _lambda(nodes[2]->p) + 21 * _lambda(nodes[6]->p) + 3 * _lambda(nodes[8]->p), 44 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[6]->p), -11 * _lambda(nodes[0]->p) - 11 * _lambda(nodes[2]->p) - 3 * _lambda(nodes[6]->p) - 3 * _lambda(nodes[8]->p), -84 * _lambda(nodes[0]->p) - 12 * _lambda(nodes[2]->p) - 28 * _lambda(nodes[6]->p) - 4 * _lambda(nodes[8]->p), -48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[6]->p), 12 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), 7 * _lambda(nodes[0]->p) + _lambda(nodes[2]->p) + 7 * _lambda(nodes[6]->p) + _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[6]->p), -_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p)},
		{44 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[6]->p), 176 * _lambda(nodes[0]->p) + 176 * _lambda(nodes[2]->p) + 48 * _lambda(nodes[6]->p) + 48 * _lambda(nodes[8]->p), 44 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[8]->p), -48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[6]->p), -192 * _lambda(nodes[0]->p) - 192 * _lambda(nodes[2]->p) - 64 * _lambda(nodes[6]->p) - 64 * _lambda(nodes[8]->p), -48 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[6]->p), 16 * _lambda(nodes[0]->p) + 16 * _lambda(nodes[2]->p) + 16 * _lambda(nodes[6]->p) + 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[8]->p)},
		{-11 * _lambda(nodes[0]->p) - 11 * _lambda(nodes[2]->p) - 3 * _lambda(nodes[6]->p) - 3 * _lambda(nodes[8]->p), 44 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[8]->p), 11 * _lambda(nodes[0]->p) + 77 * _lambda(nodes[2]->p) + 3 * _lambda(nodes[6]->p) + 21 * _lambda(nodes[8]->p), 12 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -48 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[8]->p), -12 * _lambda(nodes[0]->p) - 84 * _lambda(nodes[2]->p) - 4 * _lambda(nodes[6]->p) - 28 * _lambda(nodes[8]->p), -_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p), 4 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[8]->p), _lambda(nodes[0]->p) + 7 * _lambda(nodes[2]->p) + _lambda(nodes[6]->p) + 7 * _lambda(nodes[8]->p)},
		{-84 * _lambda(nodes[0]->p) - 12 * _lambda(nodes[2]->p) - 28 * _lambda(nodes[6]->p) - 4 * _lambda(nodes[8]->p), -48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[6]->p), 12 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), 112 * _lambda(nodes[0]->p) + 16 * _lambda(nodes[2]->p) + 112 * _lambda(nodes[6]->p) + 16 * _lambda(nodes[8]->p), 64 * _lambda(nodes[0]->p) + 64 * _lambda(nodes[6]->p), -16 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), -28 * _lambda(nodes[0]->p) - 4 * _lambda(nodes[2]->p) - 84 * _lambda(nodes[6]->p) - 12 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[6]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p)},
		{-48 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[6]->p), -192 * _lambda(nodes[0]->p) - 192 * _lambda(nodes[2]->p) - 64 * _lambda(nodes[6]->p) - 64 * _lambda(nodes[8]->p), -48 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[8]->p), 64 * _lambda(nodes[0]->p) + 64 * _lambda(nodes[6]->p), 256 * _lambda(nodes[0]->p) + 256 * _lambda(nodes[2]->p) + 256 * _lambda(nodes[6]->p) + 256 * _lambda(nodes[8]->p), 64 * _lambda(nodes[2]->p) + 64 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[6]->p), -64 * _lambda(nodes[0]->p) - 64 * _lambda(nodes[2]->p) - 192 * _lambda(nodes[6]->p) - 192 * _lambda(nodes[8]->p), -16 * _lambda(nodes[2]->p) - 48 * _lambda(nodes[8]->p)},
		{12 * _lambda(nodes[0]->p) + 12 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[6]->p) + 4 * _lambda(nodes[8]->p), -48 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[8]->p), -12 * _lambda(nodes[0]->p) - 84 * _lambda(nodes[2]->p) - 4 * _lambda(nodes[6]->p) - 28 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 16 * _lambda(nodes[2]->p) - 16 * _lambda(nodes[6]->p) - 16 * _lambda(nodes[8]->p), 64 * _lambda(nodes[2]->p) + 64 * _lambda(nodes[8]->p), 16 * _lambda(nodes[0]->p) + 112 * _lambda(nodes[2]->p) + 16 * _lambda(nodes[6]->p) + 112 * _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -16 * _lambda(nodes[2]->p) - 48 * _lambda(nodes[8]->p), -4 * _lambda(nodes[0]->p) - 28 * _lambda(nodes[2]->p) - 12 * _lambda(nodes[6]->p) - 84 * _lambda(nodes[8]->p)},
		{7 * _lambda(nodes[0]->p) + _lambda(nodes[2]->p) + 7 * _lambda(nodes[6]->p) + _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[6]->p), -_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p), -28 * _lambda(nodes[0]->p) - 4 * _lambda(nodes[2]->p) - 84 * _lambda(nodes[6]->p) - 12 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[6]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), 21 * _lambda(nodes[0]->p) + 3 * _lambda(nodes[2]->p) + 77 * _lambda(nodes[6]->p) + 11 * _lambda(nodes[8]->p), 12 * _lambda(nodes[0]->p) + 44 * _lambda(nodes[6]->p), -3 * _lambda(nodes[0]->p) - 3 * _lambda(nodes[2]->p) - 11 * _lambda(nodes[6]->p) - 11 * _lambda(nodes[8]->p)},
		{4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[6]->p), 16 * _lambda(nodes[0]->p) + 16 * _lambda(nodes[2]->p) + 16 * _lambda(nodes[6]->p) + 16 * _lambda(nodes[8]->p), 4 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[8]->p), -16 * _lambda(nodes[0]->p) - 48 * _lambda(nodes[6]->p), -64 * _lambda(nodes[0]->p) - 64 * _lambda(nodes[2]->p) - 192 * _lambda(nodes[6]->p) - 192 * _lambda(nodes[8]->p), -16 * _lambda(nodes[2]->p) - 48 * _lambda(nodes[8]->p), 12 * _lambda(nodes[0]->p) + 44 * _lambda(nodes[6]->p), 48 * _lambda(nodes[0]->p) + 48 * _lambda(nodes[2]->p) + 176 * _lambda(nodes[6]->p) + 176 * _lambda(nodes[8]->p), 12 * _lambda(nodes[2]->p) + 44 * _lambda(nodes[8]->p)},
		{-_lambda(nodes[0]->p) - _lambda(nodes[2]->p) - _lambda(nodes[6]->p) - _lambda(nodes[8]->p), 4 * _lambda(nodes[2]->p) + 4 * _lambda(nodes[8]->p), _lambda(nodes[0]->p) + 7 * _lambda(nodes[2]->p) + _lambda(nodes[6]->p) + 7 * _lambda(nodes[8]->p), 4 * _lambda(nodes[0]->p) + 4 * _lambda(nodes[2]->p) + 12 * _lambda(nodes[6]->p) + 12 * _lambda(nodes[8]->p), -16 * _lambda(nodes[2]->p) - 48 * _lambda(nodes[8]->p), -4 * _lambda(nodes[0]->p) - 28 * _lambda(nodes[2]->p) - 12 * _lambda(nodes[6]->p) - 84 * _lambda(nodes[8]->p), -3 * _lambda(nodes[0]->p) - 3 * _lambda(nodes[2]->p) - 11 * _lambda(nodes[6]->p) - 11 * _lambda(nodes[8]->p), 12 * _lambda(nodes[2]->p) + 44 * _lambda(nodes[8]->p), 3 * _lambda(nodes[0]->p) + 21 * _lambda(nodes[2]->p) + 11 * _lambda(nodes[6]->p) + 77 * _lambda(nodes[8]->p)},
	};

	_G.assign(nodeCount, std::vector<double>(nodeCount, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_G[i][j] = (hy / (hx * 360)) * __g_1[i][j] + (hx / (hy * 360)) * __g_2[i][j];
		}
	}
	
	return _G;
}


std::vector<std::vector<double>> FE::calcM() {
	std::vector<std::vector<double>> __m = {
		{16.0, 8.0, -4.0, 8.0, 4.0, -2.0, -4.0, -2.0, 1.0},
		{8.0, 64.0, 8.0, 4.0, 32.0, 4.0, -2.0, -16.0, -2.0},
		{-4.0, 8.0, 16.0, -2.0, 4.0, 8.0, 1.0, -2.0, -4.0},
		{8.0, 4.0, -2.0, 64.0, 32.0, -16.0, 8.0, 4.0, -2.0},
		{4.0, 32.0, 4.0, 32.0, 256.0, 32.0, 4.0, 32.0, 4.0},
		{-2.0, 4.0, 8.0, -16.0, 32.0, 64.0, -2.0, 4.0, 8.0},
		{-4.0, -2.0, 1.0, 8.0, 4.0, -2.0, 16.0, 8.0, -4.0},
		{-2.0, -16.0, -2.0, 4.0, 32.0, 4.0, 8.0, 64.0, 8.0},
		{1.0, -2.0, -4.0, -2.0, 4.0, 8.0, -4.0, 8.0, 16.0},
	};

	_M.assign(nodeCount, std::vector<double>(nodeCount, 0));
	
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_M[i][j] = approximateFunc(_gamma) * hx * hy / 900 * __m[i][j];
		}
	}

	return _M;
}


std::vector<std::vector<double>> FE::calcb() {
	std::vector<std::vector<double>> __b = {
		{16 * _rp(nodes[0]->p) + 8 * _rp(nodes[1]->p) - 4 * _rp(nodes[2]->p) + 8 * _rp(nodes[3]->p) + 4 * _rp(nodes[4]->p) - 2 * _rp(nodes[5]->p) - 4 * _rp(nodes[6]->p) - 2 * _rp(nodes[7]->p) + _rp(nodes[8]->p)},
		{8 * _rp(nodes[0]->p) + 64 * _rp(nodes[1]->p) + 8 * _rp(nodes[2]->p) + 4 * _rp(nodes[3]->p) + 32 * _rp(nodes[4]->p) + 4 * _rp(nodes[5]->p) - 2 * _rp(nodes[6]->p) - 16 * _rp(nodes[7]->p) - 2 * _rp(nodes[8]->p)},
		{-4 * _rp(nodes[0]->p) + 8 * _rp(nodes[1]->p) + 16 * _rp(nodes[2]->p) - 2 * _rp(nodes[3]->p) + 4 * _rp(nodes[4]->p) + 8 * _rp(nodes[5]->p) + _rp(nodes[6]->p) - 2 * _rp(nodes[7]->p) - 4 * _rp(nodes[8]->p)},
		{8 * _rp(nodes[0]->p) + 4 * _rp(nodes[1]->p) - 2 * _rp(nodes[2]->p) + 64 * _rp(nodes[3]->p) + 32 * _rp(nodes[4]->p) - 16 * _rp(nodes[5]->p) + 8 * _rp(nodes[6]->p) + 4 * _rp(nodes[7]->p) - 2 * _rp(nodes[8]->p)},
		{4 * _rp(nodes[0]->p) + 32 * _rp(nodes[1]->p) + 4 * _rp(nodes[2]->p) + 32 * _rp(nodes[3]->p) + 256 * _rp(nodes[4]->p) + 32 * _rp(nodes[5]->p) + 4 * _rp(nodes[6]->p) + 32 * _rp(nodes[7]->p) + 4 * _rp(nodes[8]->p)},
		{-2 * _rp(nodes[0]->p) + 4 * _rp(nodes[1]->p) + 8 * _rp(nodes[2]->p) - 16 * _rp(nodes[3]->p) + 32 * _rp(nodes[4]->p) + 64 * _rp(nodes[5]->p) - 2 * _rp(nodes[6]->p) + 4 * _rp(nodes[7]->p) + 8 * _rp(nodes[8]->p)},
		{-4 * _rp(nodes[0]->p) - 2 * _rp(nodes[1]->p) + _rp(nodes[2]->p) + 8 * _rp(nodes[3]->p) + 4 * _rp(nodes[4]->p) - 2 * _rp(nodes[5]->p) + 16 * _rp(nodes[6]->p) + 8 * _rp(nodes[7]->p) - 4 * _rp(nodes[8]->p)},
		{-2 * _rp(nodes[0]->p) - 16 * _rp(nodes[1]->p) - 2 * _rp(nodes[2]->p) + 4 * _rp(nodes[3]->p) + 32 * _rp(nodes[4]->p) + 4 * _rp(nodes[5]->p) + 8 * _rp(nodes[6]->p) + 64 * _rp(nodes[7]->p) + 8 * _rp(nodes[8]->p)},
		{_rp(nodes[0]->p) - 2 * _rp(nodes[1]->p) - 4 * _rp(nodes[2]->p) - 2 * _rp(nodes[3]->p) + 4 * _rp(nodes[4]->p) + 8 * _rp(nodes[5]->p) - 4 * _rp(nodes[6]->p) + 8 * _rp(nodes[7]->p) + 16 * _rp(nodes[8]->p)}
	};

	_b.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_b[i][j] = hx * hy * 900 * __b[i][j];
		}
	}

	return _b;
}

std::vector<std::vector<double>> FE::calcAs3x() {
	std::vector<std::vector<double>> __as3x = {
		{4.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2.0, 16.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	};

	_As3x.assign(nodeCount, std::vector<double>(nodeCount, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_As3x[i][j] = approximateFunc(_bcond3x) * hx / 30 * __as3x[i][j];
		}
	}

	return _As3x;

}


std::vector<std::vector<double>> FE::calcAs3y() {
	std::vector<std::vector<double>> __as3y = {
		{4.0, 0.0, 0.0, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2.0, 0.0, 0.0, 16.0, 0.0, 0.0, 2.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{-1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	};

	_As3y.assign(nodeCount, std::vector<double>(nodeCount, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_As3y[i][j] = approximateFunc(_bcond3y) * hy / 30 * __as3y[i][j];
		}
	}

	return _As3y;
}


std::vector<std::vector<double>> FE::calcbs3x() {

	std::vector<std::vector<double>> __bs3x = {
		{4 * _bcond3x_rp(nodes[0]->p) + 2 * _bcond3x_rp(nodes[1]->p) - _bcond3x_rp(nodes[2]->p)},
		{2 * _bcond3x_rp(nodes[0]->p) + 16 * _bcond3x_rp(nodes[1]->p) + 2 * _bcond3x_rp(nodes[2]->p)},
		{-_bcond3x_rp(nodes[0]->p) + 2 * _bcond3x_rp(nodes[1]->p) + 4 * _bcond3x_rp(nodes[2]->p)},
		{0},
		{0},
		{0},
		{0},
		{0},
		{0}
	};

	_bs3x.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_bs3x[i][j] = approximateFunc(_bcond3x) * hx / 30 * __bs3x[i][j];
		}
	}

	return _bs3x;
}

std::vector<std::vector<double>> FE::calcbs3y() {
	std::vector<std::vector<double>> __bs3y = {
		{4 * _bcond3y_rp(nodes[0]->p) + 2 * _bcond3y_rp(nodes[3]->p) - _bcond3y_rp(nodes[6]->p)},
		{0},
		{0},
		{2 * _bcond3y_rp(nodes[0]->p) + 16 * _bcond3y_rp(nodes[3]->p) + 2 * _bcond3y_rp(nodes[6]->p)},
		{0},
		{0},
		{-_bcond3y_rp(nodes[0]->p) + 2 * _bcond3y_rp(nodes[3]->p) + 4 * _bcond3y_rp(nodes[6]->p)},
		{0},
		{0},
	};

	_bs3y.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_bs3y[i][j] = approximateFunc(_bcond3y) * hy / 30 * __bs3y[i][j];
		}
	}

	return _bs3y;
}


std::vector<std::vector<double>> FE::calcbs2x() {

	std::vector<std::vector<double>> __bs2x = {
		{4 * _bcond2x_rp(nodes[0]->p) + 2 * _bcond2x_rp(nodes[1]->p) - _bcond2x_rp(nodes[2]->p)},
		{2 * _bcond2x_rp(nodes[0]->p) + 16 * _bcond2x_rp(nodes[1]->p) + 2 * _bcond2x_rp(nodes[2]->p)},
		{-_bcond2x_rp(nodes[0]->p) + 2 * _bcond2x_rp(nodes[1]->p) + 4 * _bcond2x_rp(nodes[2]->p)},
		{0},
		{0},
		{0},
		{0},
		{0},
		{0}
	};

	_bs2x.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_bs2x[i][j] = hx / 30 * __bs2x[i][j];
		}
	}

	return _bs2x;
}

std::vector<std::vector<double>> FE::calcbs2y() {
	std::vector<std::vector<double>> __bs2y = {
		{4 * _bcond2y_rp(nodes[0]->p) + 2 * _bcond2y_rp(nodes[3]->p) - _bcond2y_rp(nodes[6]->p)},
		{0},
		{0},
		{2 * _bcond2y_rp(nodes[0]->p) + 16 * _bcond2y_rp(nodes[3]->p) + 2 * _bcond2y_rp(nodes[6]->p)},
		{0},
		{0},
		{-_bcond2y_rp(nodes[0]->p) + 2 * _bcond2y_rp(nodes[3]->p) + 4 * _bcond2y_rp(nodes[6]->p)},
		{0},
		{0},
	};

	_bs2y.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_bs2y[i][j] = hy / 30 * __bs2y[i][j];
		}
	}

	return _bs2y;
}

std::vector<std::vector<double>> FE::calcA() {
	calcG();
	calcM();
	calcAs3x();
	calcAs3y();

	_A.assign(nodeCount, std::vector<double>(nodeCount, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_A[i][j] = _G[i][j] + _M[i][j] + _As3x[i][j] + _As3y[i][j];
		}
	}

	return _A;
}

std::vector<std::vector<double>> FE::calcRP() {
	calcb();
	calcbs2x();
	calcbs2y();
	calcbs3x();
	calcbs3y();

	_RP.assign(nodeCount, std::vector<double>(1, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_RP[i][j] = _b[i][j] + _bs2x[i][j] + _bs2y[i][j] + _bs3x[i][j] + _bs3y[i][j];
		}
	}

	return _RP;
}