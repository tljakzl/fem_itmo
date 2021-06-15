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
	std::vector<std::vector<double>> __g_1 =     {
            {432*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 432*_lambda(nodes[4]->p)/7 + 648*_lambda(nodes[8]->p)/35, -216*_lambda(nodes[0]->p)/7 - 108*_lambda(nodes[12]->p)/7 - 360*_lambda(nodes[4]->p)/7 - 324*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 288*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 144*_lambda(nodes[8]->p)/35, -432*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 432*_lambda(nodes[4]->p)/7 - 648*_lambda(nodes[8]->p)/35, -144*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 432*_lambda(nodes[4]->p)/7 - 216*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 288*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 96*_lambda(nodes[8]->p)/35, 486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, -243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 81*_lambda(nodes[4]->p)/7 - 243*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 + 3*_lambda(nodes[4]->p)/7 - 27*_lambda(nodes[8]->p)/35, -486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, -162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 162*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35},
            {-216*_lambda(nodes[0]->p)/7 - 108*_lambda(nodes[12]->p)/7 - 360*_lambda(nodes[4]->p)/7 - 324*_lambda(nodes[8]->p)/35, 216*_lambda(nodes[0]->p)/7 + 1332*_lambda(nodes[12]->p)/35 + 888*_lambda(nodes[4]->p)/7 + 324*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 144*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 592*_lambda(nodes[12]->p)/35 + 111*_lambda(nodes[4]->p)/5 + 144*_lambda(nodes[8]->p)/35, 216*_lambda(nodes[0]->p)/7 + 108*_lambda(nodes[12]->p)/7 + 360*_lambda(nodes[4]->p)/7 + 324*_lambda(nodes[8]->p)/35, 24*_lambda(nodes[0]->p) + 252*_lambda(nodes[12]->p)/5 + 168*_lambda(nodes[4]->p) + 36*_lambda(nodes[8]->p)/5, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 144*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 16*_lambda(nodes[8]->p)/5, -243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 81*_lambda(nodes[4]->p)/7 - 243*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 999*_lambda(nodes[12]->p)/35 + 999*_lambda(nodes[4]->p)/35 + 243*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 + 3*_lambda(nodes[4]->p)/7 - 27*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 111*_lambda(nodes[12]->p)/35 - 37*_lambda(nodes[4]->p)/35 + 27*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 81*_lambda(nodes[4]->p)/7 + 243*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 189*_lambda(nodes[4]->p)/5 + 27*_lambda(nodes[8]->p)/5, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 - 3*_lambda(nodes[4]->p)/7 + 27*_lambda(nodes[8]->p)/35, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 - 7*_lambda(nodes[4]->p)/5 + 3*_lambda(nodes[8]->p)/5},
            {54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 288*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 144*_lambda(nodes[8]->p)/35, 234*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 234*_lambda(nodes[4]->p)/35 + 1062*_lambda(nodes[8]->p)/35, -117*_lambda(nodes[0]->p)/35 - 177*_lambda(nodes[12]->p)/7 - 39*_lambda(nodes[4]->p)/7 - 531*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 288*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 96*_lambda(nodes[8]->p)/35, -234*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 234*_lambda(nodes[4]->p)/35 - 1062*_lambda(nodes[8]->p)/35, -78*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 234*_lambda(nodes[4]->p)/35 - 354*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 1476*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 738*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 828*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 15*_lambda(nodes[4]->p)/7 - 414*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 1476*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 492*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 828*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 276*_lambda(nodes[8]->p)/35},
            {-27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 144*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 592*_lambda(nodes[12]->p)/35 + 111*_lambda(nodes[4]->p)/5 + 144*_lambda(nodes[8]->p)/35, -117*_lambda(nodes[0]->p)/35 - 177*_lambda(nodes[12]->p)/7 - 39*_lambda(nodes[4]->p)/7 - 531*_lambda(nodes[8]->p)/35, 117*_lambda(nodes[0]->p)/35 + 2183*_lambda(nodes[12]->p)/35 + 481*_lambda(nodes[4]->p)/35 + 531*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 144*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 16*_lambda(nodes[8]->p)/5, 117*_lambda(nodes[0]->p)/35 + 177*_lambda(nodes[12]->p)/7 + 39*_lambda(nodes[4]->p)/7 + 531*_lambda(nodes[8]->p)/35, 13*_lambda(nodes[0]->p)/5 + 413*_lambda(nodes[12]->p)/5 + 91*_lambda(nodes[4]->p)/5 + 59*_lambda(nodes[8]->p)/5, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 738*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 3034*_lambda(nodes[12]->p)/35 + 111*_lambda(nodes[4]->p)/5 + 738*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 15*_lambda(nodes[4]->p)/7 - 414*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/7 + 1702*_lambda(nodes[12]->p)/35 + 37*_lambda(nodes[4]->p)/7 + 414*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 738*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 82*_lambda(nodes[8]->p)/5, 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 15*_lambda(nodes[4]->p)/7 + 414*_lambda(nodes[8]->p)/35, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 7*_lambda(nodes[4]->p) + 46*_lambda(nodes[8]->p)/5},
            {-432*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 432*_lambda(nodes[4]->p)/7 - 648*_lambda(nodes[8]->p)/35, 216*_lambda(nodes[0]->p)/7 + 108*_lambda(nodes[12]->p)/7 + 360*_lambda(nodes[4]->p)/7 + 324*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 288*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 144*_lambda(nodes[8]->p)/35, 432*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 432*_lambda(nodes[4]->p)/7 + 648*_lambda(nodes[8]->p)/35, 144*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 432*_lambda(nodes[4]->p)/7 + 216*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 288*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 96*_lambda(nodes[8]->p)/35, -486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 81*_lambda(nodes[4]->p)/7 + 243*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 - 3*_lambda(nodes[4]->p)/7 + 27*_lambda(nodes[8]->p)/35, 486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 162*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35},
            {-144*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 432*_lambda(nodes[4]->p)/7 - 216*_lambda(nodes[8]->p)/35, 24*_lambda(nodes[0]->p) + 252*_lambda(nodes[12]->p)/5 + 168*_lambda(nodes[4]->p) + 36*_lambda(nodes[8]->p)/5, -18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 96*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 16*_lambda(nodes[8]->p)/5, 144*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 432*_lambda(nodes[4]->p)/7 + 216*_lambda(nodes[8]->p)/35, 240*_lambda(nodes[0]->p)/7 + 2376*_lambda(nodes[12]->p)/35 + 1584*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/7, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 96*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p) + 1056*_lambda(nodes[12]->p)/35 + 198*_lambda(nodes[4]->p)/5 + 32*_lambda(nodes[8]->p)/7, -162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 162*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 189*_lambda(nodes[4]->p)/5 + 27*_lambda(nodes[8]->p)/5, 6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 - 7*_lambda(nodes[4]->p)/5 + 3*_lambda(nodes[8]->p)/5, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 162*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/7 + 1782*_lambda(nodes[12]->p)/35 + 1782*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/7, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -2*_lambda(nodes[0]->p)/7 + 198*_lambda(nodes[12]->p)/35 - 66*_lambda(nodes[4]->p)/35 + 6*_lambda(nodes[8]->p)/7},
            {-54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 288*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 144*_lambda(nodes[8]->p)/35, -234*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 234*_lambda(nodes[4]->p)/35 - 1062*_lambda(nodes[8]->p)/35, 117*_lambda(nodes[0]->p)/35 + 177*_lambda(nodes[12]->p)/7 + 39*_lambda(nodes[4]->p)/7 + 531*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 288*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 96*_lambda(nodes[8]->p)/35, 234*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 234*_lambda(nodes[4]->p)/35 + 1062*_lambda(nodes[8]->p)/35, 78*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 234*_lambda(nodes[4]->p)/35 + 354*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 1476*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 738*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 828*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 15*_lambda(nodes[4]->p)/7 + 414*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 1476*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 492*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 828*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 276*_lambda(nodes[8]->p)/35},
            {-18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 96*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 16*_lambda(nodes[8]->p)/5, -78*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 234*_lambda(nodes[4]->p)/35 - 354*_lambda(nodes[8]->p)/35, 13*_lambda(nodes[0]->p)/5 + 413*_lambda(nodes[12]->p)/5 + 91*_lambda(nodes[4]->p)/5 + 59*_lambda(nodes[8]->p)/5, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 96*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p) + 1056*_lambda(nodes[12]->p)/35 + 198*_lambda(nodes[4]->p)/5 + 32*_lambda(nodes[8]->p)/7, 78*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 234*_lambda(nodes[4]->p)/35 + 354*_lambda(nodes[8]->p)/35, 26*_lambda(nodes[0]->p)/7 + 3894*_lambda(nodes[12]->p)/35 + 858*_lambda(nodes[4]->p)/35 + 118*_lambda(nodes[8]->p)/7, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 492*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 82*_lambda(nodes[8]->p)/5, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 276*_lambda(nodes[8]->p)/35, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 7*_lambda(nodes[4]->p) + 46*_lambda(nodes[8]->p)/5, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 492*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p) + 5412*_lambda(nodes[12]->p)/35 + 198*_lambda(nodes[4]->p)/5 + 164*_lambda(nodes[8]->p)/7, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 276*_lambda(nodes[8]->p)/35, 10*_lambda(nodes[0]->p)/7 + 3036*_lambda(nodes[12]->p)/35 + 66*_lambda(nodes[4]->p)/7 + 92*_lambda(nodes[8]->p)/7},
            {486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, -243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 81*_lambda(nodes[4]->p)/7 - 243*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 1476*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 738*_lambda(nodes[8]->p)/35, -486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, -162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 162*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 1476*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 492*_lambda(nodes[8]->p)/35, 648*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 648*_lambda(nodes[4]->p)/35 + 432*_lambda(nodes[8]->p)/7, -324*_lambda(nodes[0]->p)/35 - 360*_lambda(nodes[12]->p)/7 - 108*_lambda(nodes[4]->p)/7 - 216*_lambda(nodes[8]->p)/7, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 216*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 24*_lambda(nodes[4]->p)/7 - 108*_lambda(nodes[8]->p)/7, -648*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 648*_lambda(nodes[4]->p)/35 - 432*_lambda(nodes[8]->p)/7, -216*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 648*_lambda(nodes[4]->p)/35 - 144*_lambda(nodes[8]->p)/7, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 216*_lambda(nodes[8]->p)/7, -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 72*_lambda(nodes[8]->p)/7},
            {-243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 81*_lambda(nodes[4]->p)/7 - 243*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 999*_lambda(nodes[12]->p)/35 + 999*_lambda(nodes[4]->p)/35 + 243*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 9*_lambda(nodes[4]->p) - 738*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 3034*_lambda(nodes[12]->p)/35 + 111*_lambda(nodes[4]->p)/5 + 738*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 81*_lambda(nodes[4]->p)/7 + 243*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 189*_lambda(nodes[4]->p)/5 + 27*_lambda(nodes[8]->p)/5, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 738*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 82*_lambda(nodes[8]->p)/5, -324*_lambda(nodes[0]->p)/35 - 360*_lambda(nodes[12]->p)/7 - 108*_lambda(nodes[4]->p)/7 - 216*_lambda(nodes[8]->p)/7, 324*_lambda(nodes[0]->p)/35 + 888*_lambda(nodes[12]->p)/7 + 1332*_lambda(nodes[4]->p)/35 + 216*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 24*_lambda(nodes[4]->p)/7 - 108*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 444*_lambda(nodes[12]->p)/7 + 296*_lambda(nodes[4]->p)/35 + 108*_lambda(nodes[8]->p)/7, 324*_lambda(nodes[0]->p)/35 + 360*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 216*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/5 + 168*_lambda(nodes[12]->p) + 252*_lambda(nodes[4]->p)/5 + 24*_lambda(nodes[8]->p), 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 24*_lambda(nodes[4]->p)/7 + 108*_lambda(nodes[8]->p)/7, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 56*_lambda(nodes[4]->p)/5 + 12*_lambda(nodes[8]->p)},
            {-18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 + 3*_lambda(nodes[4]->p)/7 - 27*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 828*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 15*_lambda(nodes[4]->p)/7 - 414*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 828*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 276*_lambda(nodes[8]->p)/35, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 216*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 24*_lambda(nodes[4]->p)/7 - 108*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/35 + 144*_lambda(nodes[8]->p)/7, -36*_lambda(nodes[0]->p)/35 - 120*_lambda(nodes[12]->p)/7 - 12*_lambda(nodes[4]->p)/7 - 72*_lambda(nodes[8]->p)/7, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 216*_lambda(nodes[8]->p)/7, -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 72*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/35 - 144*_lambda(nodes[8]->p)/7, -24*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/35 - 48*_lambda(nodes[8]->p)/7},
            {9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 + 3*_lambda(nodes[4]->p)/7 - 27*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 111*_lambda(nodes[12]->p)/35 - 37*_lambda(nodes[4]->p)/35 + 27*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 15*_lambda(nodes[4]->p)/7 - 414*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/7 + 1702*_lambda(nodes[12]->p)/35 + 37*_lambda(nodes[4]->p)/7 + 414*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 - 3*_lambda(nodes[4]->p)/7 + 27*_lambda(nodes[8]->p)/35, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 - 7*_lambda(nodes[4]->p)/5 + 3*_lambda(nodes[8]->p)/5, 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 15*_lambda(nodes[4]->p)/7 + 414*_lambda(nodes[8]->p)/35, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 7*_lambda(nodes[4]->p) + 46*_lambda(nodes[8]->p)/5, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 24*_lambda(nodes[4]->p)/7 - 108*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 444*_lambda(nodes[12]->p)/7 + 296*_lambda(nodes[4]->p)/35 + 108*_lambda(nodes[8]->p)/7, -36*_lambda(nodes[0]->p)/35 - 120*_lambda(nodes[12]->p)/7 - 12*_lambda(nodes[4]->p)/7 - 72*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/35 + 296*_lambda(nodes[12]->p)/7 + 148*_lambda(nodes[4]->p)/35 + 72*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 24*_lambda(nodes[4]->p)/7 + 108*_lambda(nodes[8]->p)/7, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 56*_lambda(nodes[4]->p)/5 + 12*_lambda(nodes[8]->p), 36*_lambda(nodes[0]->p)/35 + 120*_lambda(nodes[12]->p)/7 + 12*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/7, 4*_lambda(nodes[0]->p)/5 + 56*_lambda(nodes[12]->p) + 28*_lambda(nodes[4]->p)/5 + 8*_lambda(nodes[8]->p)},
            {-486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 81*_lambda(nodes[4]->p)/7 + 243*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 1476*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 9*_lambda(nodes[4]->p) + 738*_lambda(nodes[8]->p)/35, 486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 162*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 1476*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 492*_lambda(nodes[8]->p)/35, -648*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 648*_lambda(nodes[4]->p)/35 - 432*_lambda(nodes[8]->p)/7, 324*_lambda(nodes[0]->p)/35 + 360*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 216*_lambda(nodes[8]->p)/7, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 216*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 24*_lambda(nodes[4]->p)/7 + 108*_lambda(nodes[8]->p)/7, 648*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 648*_lambda(nodes[4]->p)/35 + 432*_lambda(nodes[8]->p)/7, 216*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 648*_lambda(nodes[4]->p)/35 + 144*_lambda(nodes[8]->p)/7, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 216*_lambda(nodes[8]->p)/7, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 72*_lambda(nodes[8]->p)/7},
            {-162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 162*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 189*_lambda(nodes[4]->p)/5 + 27*_lambda(nodes[8]->p)/5, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/5 - 492*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 147*_lambda(nodes[4]->p)/5 + 82*_lambda(nodes[8]->p)/5, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 162*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/7 + 1782*_lambda(nodes[12]->p)/35 + 1782*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/7, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/5 + 492*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p) + 5412*_lambda(nodes[12]->p)/35 + 198*_lambda(nodes[4]->p)/5 + 164*_lambda(nodes[8]->p)/7, -216*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 648*_lambda(nodes[4]->p)/35 - 144*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/5 + 168*_lambda(nodes[12]->p) + 252*_lambda(nodes[4]->p)/5 + 24*_lambda(nodes[8]->p), -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 72*_lambda(nodes[8]->p)/7, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 56*_lambda(nodes[4]->p)/5 + 12*_lambda(nodes[8]->p), 216*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 648*_lambda(nodes[4]->p)/35 + 144*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/7 + 1584*_lambda(nodes[12]->p)/7 + 2376*_lambda(nodes[4]->p)/35 + 240*_lambda(nodes[8]->p)/7, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 72*_lambda(nodes[8]->p)/7, 16*_lambda(nodes[0]->p)/7 + 792*_lambda(nodes[12]->p)/7 + 528*_lambda(nodes[4]->p)/35 + 120*_lambda(nodes[8]->p)/7},
            {18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 - 3*_lambda(nodes[4]->p)/7 + 27*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 828*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 15*_lambda(nodes[4]->p)/7 + 414*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 828*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 276*_lambda(nodes[8]->p)/35, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 216*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 24*_lambda(nodes[4]->p)/7 + 108*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/35 - 144*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/35 + 120*_lambda(nodes[12]->p)/7 + 12*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/7, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 216*_lambda(nodes[8]->p)/7, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 72*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/35 + 144*_lambda(nodes[8]->p)/7, 24*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/35 + 48*_lambda(nodes[8]->p)/7},
            {6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 - 7*_lambda(nodes[4]->p)/5 + 3*_lambda(nodes[8]->p)/5, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/7 - 276*_lambda(nodes[8]->p)/35, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 7*_lambda(nodes[4]->p) + 46*_lambda(nodes[8]->p)/5, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -2*_lambda(nodes[0]->p)/7 + 198*_lambda(nodes[12]->p)/35 - 66*_lambda(nodes[4]->p)/35 + 6*_lambda(nodes[8]->p)/7, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/7 + 276*_lambda(nodes[8]->p)/35, 10*_lambda(nodes[0]->p)/7 + 3036*_lambda(nodes[12]->p)/35 + 66*_lambda(nodes[4]->p)/7 + 92*_lambda(nodes[8]->p)/7, -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 72*_lambda(nodes[8]->p)/7, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 56*_lambda(nodes[4]->p)/5 + 12*_lambda(nodes[8]->p), -24*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/35 - 48*_lambda(nodes[8]->p)/7, 4*_lambda(nodes[0]->p)/5 + 56*_lambda(nodes[12]->p) + 28*_lambda(nodes[4]->p)/5 + 8*_lambda(nodes[8]->p), 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 72*_lambda(nodes[8]->p)/7, 16*_lambda(nodes[0]->p)/7 + 792*_lambda(nodes[12]->p)/7 + 528*_lambda(nodes[4]->p)/35 + 120*_lambda(nodes[8]->p)/7, 24*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/35 + 48*_lambda(nodes[8]->p)/7, 8*_lambda(nodes[0]->p)/7 + 528*_lambda(nodes[12]->p)/7 + 264*_lambda(nodes[4]->p)/35 + 80*_lambda(nodes[8]->p)/7},
    };

	std::vector<std::vector<double>> __g_2 = {
            {432*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 648*_lambda(nodes[4]->p)/35 + 432*_lambda(nodes[8]->p)/7, 54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 288*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, -216*_lambda(nodes[0]->p)/7 - 108*_lambda(nodes[12]->p)/7 - 324*_lambda(nodes[4]->p)/35 - 360*_lambda(nodes[8]->p)/7, -27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), 486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, -243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 243*_lambda(nodes[4]->p)/35 - 81*_lambda(nodes[8]->p)/7, 9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 - 27*_lambda(nodes[4]->p)/35 + 3*_lambda(nodes[8]->p)/7, -432*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 648*_lambda(nodes[4]->p)/35 - 432*_lambda(nodes[8]->p)/7, -54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 288*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -144*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 216*_lambda(nodes[4]->p)/35 - 432*_lambda(nodes[8]->p)/7, -18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 96*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 162*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35},
            {54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 288*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 234*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 1062*_lambda(nodes[4]->p)/35 + 234*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), -117*_lambda(nodes[0]->p)/35 - 177*_lambda(nodes[12]->p)/7 - 531*_lambda(nodes[4]->p)/35 - 39*_lambda(nodes[8]->p)/7, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 1476*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 828*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 738*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 414*_lambda(nodes[4]->p)/35 - 15*_lambda(nodes[8]->p)/7, -54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 288*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -234*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 1062*_lambda(nodes[4]->p)/35 - 234*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 96*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -78*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 354*_lambda(nodes[4]->p)/35 - 234*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 1476*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 828*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 492*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 276*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7},
            {-216*_lambda(nodes[0]->p)/7 - 108*_lambda(nodes[12]->p)/7 - 324*_lambda(nodes[4]->p)/35 - 360*_lambda(nodes[8]->p)/7, -27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), 216*_lambda(nodes[0]->p)/7 + 1332*_lambda(nodes[12]->p)/35 + 324*_lambda(nodes[4]->p)/35 + 888*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 592*_lambda(nodes[12]->p)/35 + 144*_lambda(nodes[4]->p)/35 + 111*_lambda(nodes[8]->p)/5, -243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 243*_lambda(nodes[4]->p)/35 - 81*_lambda(nodes[8]->p)/7, 9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 - 27*_lambda(nodes[4]->p)/35 + 3*_lambda(nodes[8]->p)/7, 243*_lambda(nodes[0]->p)/35 + 999*_lambda(nodes[12]->p)/35 + 243*_lambda(nodes[4]->p)/35 + 999*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 111*_lambda(nodes[12]->p)/35 + 27*_lambda(nodes[4]->p)/35 - 37*_lambda(nodes[8]->p)/35, 216*_lambda(nodes[0]->p)/7 + 108*_lambda(nodes[12]->p)/7 + 324*_lambda(nodes[4]->p)/35 + 360*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), 24*_lambda(nodes[0]->p) + 252*_lambda(nodes[12]->p)/5 + 36*_lambda(nodes[4]->p)/5 + 168*_lambda(nodes[8]->p), 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 16*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 243*_lambda(nodes[4]->p)/35 + 81*_lambda(nodes[8]->p)/7, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 + 27*_lambda(nodes[4]->p)/35 - 3*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 27*_lambda(nodes[4]->p)/5 + 189*_lambda(nodes[8]->p)/5, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 + 3*_lambda(nodes[4]->p)/5 - 7*_lambda(nodes[8]->p)/5},
            {-27*_lambda(nodes[0]->p)/5 - 48*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), -117*_lambda(nodes[0]->p)/35 - 177*_lambda(nodes[12]->p)/7 - 531*_lambda(nodes[4]->p)/35 - 39*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 592*_lambda(nodes[12]->p)/35 + 144*_lambda(nodes[4]->p)/35 + 111*_lambda(nodes[8]->p)/5, 117*_lambda(nodes[0]->p)/35 + 2183*_lambda(nodes[12]->p)/35 + 531*_lambda(nodes[4]->p)/35 + 481*_lambda(nodes[8]->p)/35, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 738*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 414*_lambda(nodes[4]->p)/35 - 15*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 3034*_lambda(nodes[12]->p)/35 + 738*_lambda(nodes[4]->p)/35 + 111*_lambda(nodes[8]->p)/5, 9*_lambda(nodes[0]->p)/7 + 1702*_lambda(nodes[12]->p)/35 + 414*_lambda(nodes[4]->p)/35 + 37*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), 117*_lambda(nodes[0]->p)/35 + 177*_lambda(nodes[12]->p)/7 + 531*_lambda(nodes[4]->p)/35 + 39*_lambda(nodes[8]->p)/7, 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 16*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, 13*_lambda(nodes[0]->p)/5 + 413*_lambda(nodes[12]->p)/5 + 59*_lambda(nodes[4]->p)/5 + 91*_lambda(nodes[8]->p)/5, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 738*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 414*_lambda(nodes[4]->p)/35 + 15*_lambda(nodes[8]->p)/7, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 82*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 46*_lambda(nodes[4]->p)/5 + 7*_lambda(nodes[8]->p)},
            {486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 1476*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, -243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 243*_lambda(nodes[4]->p)/35 - 81*_lambda(nodes[8]->p)/7, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 738*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), 648*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 432*_lambda(nodes[4]->p)/7 + 648*_lambda(nodes[8]->p)/35, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, -324*_lambda(nodes[0]->p)/35 - 360*_lambda(nodes[12]->p)/7 - 216*_lambda(nodes[4]->p)/7 - 108*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 108*_lambda(nodes[4]->p)/7 - 24*_lambda(nodes[8]->p)/7, -486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 1476*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 162*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 492*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -648*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 432*_lambda(nodes[4]->p)/7 - 648*_lambda(nodes[8]->p)/35, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 216*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, -216*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/7 - 648*_lambda(nodes[8]->p)/35, -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35},
            {-18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 828*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, 9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 - 27*_lambda(nodes[4]->p)/35 + 3*_lambda(nodes[8]->p)/7, -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 414*_lambda(nodes[4]->p)/35 - 15*_lambda(nodes[8]->p)/7, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, 72*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/35, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 108*_lambda(nodes[4]->p)/7 - 24*_lambda(nodes[8]->p)/7, -36*_lambda(nodes[0]->p)/35 - 120*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/7 - 12*_lambda(nodes[8]->p)/7, 18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 828*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, 6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 276*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 216*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, -72*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/7 - 72*_lambda(nodes[8]->p)/35, -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, -24*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 48*_lambda(nodes[4]->p)/7 - 72*_lambda(nodes[8]->p)/35},
            {-243*_lambda(nodes[0]->p)/35 - 81*_lambda(nodes[12]->p)/7 - 243*_lambda(nodes[4]->p)/35 - 81*_lambda(nodes[8]->p)/7, -27*_lambda(nodes[0]->p)/5 - 246*_lambda(nodes[12]->p)/7 - 738*_lambda(nodes[4]->p)/35 - 9*_lambda(nodes[8]->p), 243*_lambda(nodes[0]->p)/35 + 999*_lambda(nodes[12]->p)/35 + 243*_lambda(nodes[4]->p)/35 + 999*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 3034*_lambda(nodes[12]->p)/35 + 738*_lambda(nodes[4]->p)/35 + 111*_lambda(nodes[8]->p)/5, -324*_lambda(nodes[0]->p)/35 - 360*_lambda(nodes[12]->p)/7 - 216*_lambda(nodes[4]->p)/7 - 108*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 108*_lambda(nodes[4]->p)/7 - 24*_lambda(nodes[8]->p)/7, 324*_lambda(nodes[0]->p)/35 + 888*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 1332*_lambda(nodes[8]->p)/35, 72*_lambda(nodes[0]->p)/35 + 444*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 296*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 243*_lambda(nodes[4]->p)/35 + 81*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 738*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 27*_lambda(nodes[4]->p)/5 + 189*_lambda(nodes[8]->p)/5, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 82*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, 324*_lambda(nodes[0]->p)/35 + 360*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 108*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 24*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/5 + 168*_lambda(nodes[12]->p) + 24*_lambda(nodes[4]->p) + 252*_lambda(nodes[8]->p)/5, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 12*_lambda(nodes[4]->p) + 56*_lambda(nodes[8]->p)/5},
            {9*_lambda(nodes[0]->p)/35 - 9*_lambda(nodes[12]->p)/7 - 27*_lambda(nodes[4]->p)/35 + 3*_lambda(nodes[8]->p)/7, -9*_lambda(nodes[0]->p)/7 - 138*_lambda(nodes[12]->p)/7 - 414*_lambda(nodes[4]->p)/35 - 15*_lambda(nodes[8]->p)/7, -9*_lambda(nodes[0]->p)/35 + 111*_lambda(nodes[12]->p)/35 + 27*_lambda(nodes[4]->p)/35 - 37*_lambda(nodes[8]->p)/35, 9*_lambda(nodes[0]->p)/7 + 1702*_lambda(nodes[12]->p)/35 + 414*_lambda(nodes[4]->p)/35 + 37*_lambda(nodes[8]->p)/7, -72*_lambda(nodes[0]->p)/35 - 180*_lambda(nodes[12]->p)/7 - 108*_lambda(nodes[4]->p)/7 - 24*_lambda(nodes[8]->p)/7, -36*_lambda(nodes[0]->p)/35 - 120*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/7 - 12*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 444*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 296*_lambda(nodes[8]->p)/35, 36*_lambda(nodes[0]->p)/35 + 296*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 148*_lambda(nodes[8]->p)/35, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 + 27*_lambda(nodes[4]->p)/35 - 3*_lambda(nodes[8]->p)/7, 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 414*_lambda(nodes[4]->p)/35 + 15*_lambda(nodes[8]->p)/7, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 + 3*_lambda(nodes[4]->p)/5 - 7*_lambda(nodes[8]->p)/5, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 46*_lambda(nodes[4]->p)/5 + 7*_lambda(nodes[8]->p), 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 24*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/35 + 120*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 12*_lambda(nodes[8]->p)/7, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 12*_lambda(nodes[4]->p) + 56*_lambda(nodes[8]->p)/5, 4*_lambda(nodes[0]->p)/5 + 56*_lambda(nodes[12]->p) + 8*_lambda(nodes[4]->p) + 28*_lambda(nodes[8]->p)/5},
            {-432*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 648*_lambda(nodes[4]->p)/35 - 432*_lambda(nodes[8]->p)/7, -54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 288*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, 216*_lambda(nodes[0]->p)/7 + 108*_lambda(nodes[12]->p)/7 + 324*_lambda(nodes[4]->p)/35 + 360*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), -486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 243*_lambda(nodes[4]->p)/35 + 81*_lambda(nodes[8]->p)/7, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 + 27*_lambda(nodes[4]->p)/35 - 3*_lambda(nodes[8]->p)/7, 432*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 648*_lambda(nodes[4]->p)/35 + 432*_lambda(nodes[8]->p)/7, 54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 288*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 144*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 216*_lambda(nodes[4]->p)/35 + 432*_lambda(nodes[8]->p)/7, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 96*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 162*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35},
            {-54*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 288*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -234*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 1062*_lambda(nodes[4]->p)/35 - 234*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 48*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), 117*_lambda(nodes[0]->p)/35 + 177*_lambda(nodes[12]->p)/7 + 531*_lambda(nodes[4]->p)/35 + 39*_lambda(nodes[8]->p)/7, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 1476*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 828*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 738*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 414*_lambda(nodes[4]->p)/35 + 15*_lambda(nodes[8]->p)/7, 54*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 288*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 234*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 1062*_lambda(nodes[4]->p)/35 + 234*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 96*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 78*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 354*_lambda(nodes[4]->p)/35 + 234*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 1476*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 828*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 492*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 276*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7},
            {-144*_lambda(nodes[0]->p)/7 - 648*_lambda(nodes[12]->p)/35 - 216*_lambda(nodes[4]->p)/35 - 432*_lambda(nodes[8]->p)/7, -18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 96*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, 24*_lambda(nodes[0]->p) + 252*_lambda(nodes[12]->p)/5 + 36*_lambda(nodes[4]->p)/5 + 168*_lambda(nodes[8]->p), 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 16*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, -162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 162*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 27*_lambda(nodes[4]->p)/5 + 189*_lambda(nodes[8]->p)/5, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 + 3*_lambda(nodes[4]->p)/5 - 7*_lambda(nodes[8]->p)/5, 144*_lambda(nodes[0]->p)/7 + 648*_lambda(nodes[12]->p)/35 + 216*_lambda(nodes[4]->p)/35 + 432*_lambda(nodes[8]->p)/7, 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 96*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 240*_lambda(nodes[0]->p)/7 + 2376*_lambda(nodes[12]->p)/35 + 72*_lambda(nodes[4]->p)/7 + 1584*_lambda(nodes[8]->p)/7, 6*_lambda(nodes[0]->p) + 1056*_lambda(nodes[12]->p)/35 + 32*_lambda(nodes[4]->p)/7 + 198*_lambda(nodes[8]->p)/5, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 162*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/7 + 1782*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/7 + 1782*_lambda(nodes[8]->p)/35, -2*_lambda(nodes[0]->p)/7 + 198*_lambda(nodes[12]->p)/35 + 6*_lambda(nodes[4]->p)/7 - 66*_lambda(nodes[8]->p)/35},
            {-18*_lambda(nodes[0]->p)/5 - 288*_lambda(nodes[12]->p)/35 - 96*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -78*_lambda(nodes[0]->p)/35 - 1062*_lambda(nodes[12]->p)/35 - 354*_lambda(nodes[4]->p)/35 - 234*_lambda(nodes[8]->p)/35, 21*_lambda(nodes[0]->p)/5 + 112*_lambda(nodes[12]->p)/5 + 16*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, 13*_lambda(nodes[0]->p)/5 + 413*_lambda(nodes[12]->p)/5 + 59*_lambda(nodes[4]->p)/5 + 91*_lambda(nodes[8]->p)/5, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 492*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 276*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 82*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 46*_lambda(nodes[4]->p)/5 + 7*_lambda(nodes[8]->p), 18*_lambda(nodes[0]->p)/5 + 288*_lambda(nodes[12]->p)/35 + 96*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 78*_lambda(nodes[0]->p)/35 + 1062*_lambda(nodes[12]->p)/35 + 354*_lambda(nodes[4]->p)/35 + 234*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p) + 1056*_lambda(nodes[12]->p)/35 + 32*_lambda(nodes[4]->p)/7 + 198*_lambda(nodes[8]->p)/5, 26*_lambda(nodes[0]->p)/7 + 3894*_lambda(nodes[12]->p)/35 + 118*_lambda(nodes[4]->p)/7 + 858*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 492*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 276*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, 6*_lambda(nodes[0]->p) + 5412*_lambda(nodes[12]->p)/35 + 164*_lambda(nodes[4]->p)/7 + 198*_lambda(nodes[8]->p)/5, 10*_lambda(nodes[0]->p)/7 + 3036*_lambda(nodes[12]->p)/35 + 92*_lambda(nodes[4]->p)/7 + 66*_lambda(nodes[8]->p)/7},
            {-486*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 486*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, -54*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 1476*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, 243*_lambda(nodes[0]->p)/35 + 81*_lambda(nodes[12]->p)/7 + 243*_lambda(nodes[4]->p)/35 + 81*_lambda(nodes[8]->p)/7, 27*_lambda(nodes[0]->p)/5 + 246*_lambda(nodes[12]->p)/7 + 738*_lambda(nodes[4]->p)/35 + 9*_lambda(nodes[8]->p), -648*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 432*_lambda(nodes[4]->p)/7 - 648*_lambda(nodes[8]->p)/35, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 216*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, 324*_lambda(nodes[0]->p)/35 + 360*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 108*_lambda(nodes[8]->p)/7, 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 24*_lambda(nodes[8]->p)/7, 486*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 486*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, 54*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 1476*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 162*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 492*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 648*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 432*_lambda(nodes[4]->p)/7 + 648*_lambda(nodes[8]->p)/35, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, 216*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/7 + 648*_lambda(nodes[8]->p)/35, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35},
            {18*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 54*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 828*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, -9*_lambda(nodes[0]->p)/35 + 9*_lambda(nodes[12]->p)/7 + 27*_lambda(nodes[4]->p)/35 - 3*_lambda(nodes[8]->p)/7, 9*_lambda(nodes[0]->p)/7 + 138*_lambda(nodes[12]->p)/7 + 414*_lambda(nodes[4]->p)/35 + 15*_lambda(nodes[8]->p)/7, -144*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 216*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, -72*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/7 - 72*_lambda(nodes[8]->p)/35, 72*_lambda(nodes[0]->p)/35 + 180*_lambda(nodes[12]->p)/7 + 108*_lambda(nodes[4]->p)/7 + 24*_lambda(nodes[8]->p)/7, 36*_lambda(nodes[0]->p)/35 + 120*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 12*_lambda(nodes[8]->p)/7, -18*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 828*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 276*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, 144*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 216*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, 72*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/35, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, 24*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 48*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/35},
            {-162*_lambda(nodes[0]->p)/35 - 486*_lambda(nodes[12]->p)/35 - 162*_lambda(nodes[4]->p)/35 - 486*_lambda(nodes[8]->p)/35, -18*_lambda(nodes[0]->p)/5 - 1476*_lambda(nodes[12]->p)/35 - 492*_lambda(nodes[4]->p)/35 - 54*_lambda(nodes[8]->p)/5, 27*_lambda(nodes[0]->p)/5 + 189*_lambda(nodes[12]->p)/5 + 27*_lambda(nodes[4]->p)/5 + 189*_lambda(nodes[8]->p)/5, 21*_lambda(nodes[0]->p)/5 + 574*_lambda(nodes[12]->p)/5 + 82*_lambda(nodes[4]->p)/5 + 147*_lambda(nodes[8]->p)/5, -216*_lambda(nodes[0]->p)/35 - 432*_lambda(nodes[12]->p)/7 - 144*_lambda(nodes[4]->p)/7 - 648*_lambda(nodes[8]->p)/35, -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, 36*_lambda(nodes[0]->p)/5 + 168*_lambda(nodes[12]->p) + 24*_lambda(nodes[4]->p) + 252*_lambda(nodes[8]->p)/5, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 12*_lambda(nodes[4]->p) + 56*_lambda(nodes[8]->p)/5, 162*_lambda(nodes[0]->p)/35 + 486*_lambda(nodes[12]->p)/35 + 162*_lambda(nodes[4]->p)/35 + 486*_lambda(nodes[8]->p)/35, 18*_lambda(nodes[0]->p)/5 + 1476*_lambda(nodes[12]->p)/35 + 492*_lambda(nodes[4]->p)/35 + 54*_lambda(nodes[8]->p)/5, 54*_lambda(nodes[0]->p)/7 + 1782*_lambda(nodes[12]->p)/35 + 54*_lambda(nodes[4]->p)/7 + 1782*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p) + 5412*_lambda(nodes[12]->p)/35 + 164*_lambda(nodes[4]->p)/7 + 198*_lambda(nodes[8]->p)/5, 216*_lambda(nodes[0]->p)/35 + 432*_lambda(nodes[12]->p)/7 + 144*_lambda(nodes[4]->p)/7 + 648*_lambda(nodes[8]->p)/35, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, 72*_lambda(nodes[0]->p)/7 + 1584*_lambda(nodes[12]->p)/7 + 240*_lambda(nodes[4]->p)/7 + 2376*_lambda(nodes[8]->p)/35, 16*_lambda(nodes[0]->p)/7 + 792*_lambda(nodes[12]->p)/7 + 120*_lambda(nodes[4]->p)/7 + 528*_lambda(nodes[8]->p)/35},
            {6*_lambda(nodes[0]->p)/35 - 54*_lambda(nodes[12]->p)/35 - 18*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/35, -6*_lambda(nodes[0]->p)/7 - 828*_lambda(nodes[12]->p)/35 - 276*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/7, -_lambda(nodes[0]->p)/5 + 21*_lambda(nodes[12]->p)/5 + 3*_lambda(nodes[4]->p)/5 - 7*_lambda(nodes[8]->p)/5, _lambda(nodes[0]->p) + 322*_lambda(nodes[12]->p)/5 + 46*_lambda(nodes[4]->p)/5 + 7*_lambda(nodes[8]->p), -48*_lambda(nodes[0]->p)/35 - 216*_lambda(nodes[12]->p)/7 - 72*_lambda(nodes[4]->p)/7 - 144*_lambda(nodes[8]->p)/35, -24*_lambda(nodes[0]->p)/35 - 144*_lambda(nodes[12]->p)/7 - 48*_lambda(nodes[4]->p)/7 - 72*_lambda(nodes[8]->p)/35, 8*_lambda(nodes[0]->p)/5 + 84*_lambda(nodes[12]->p) + 12*_lambda(nodes[4]->p) + 56*_lambda(nodes[8]->p)/5, 4*_lambda(nodes[0]->p)/5 + 56*_lambda(nodes[12]->p) + 8*_lambda(nodes[4]->p) + 28*_lambda(nodes[8]->p)/5, -6*_lambda(nodes[0]->p)/35 + 54*_lambda(nodes[12]->p)/35 + 18*_lambda(nodes[4]->p)/35 - 18*_lambda(nodes[8]->p)/35, 6*_lambda(nodes[0]->p)/7 + 828*_lambda(nodes[12]->p)/35 + 276*_lambda(nodes[4]->p)/35 + 18*_lambda(nodes[8]->p)/7, -2*_lambda(nodes[0]->p)/7 + 198*_lambda(nodes[12]->p)/35 + 6*_lambda(nodes[4]->p)/7 - 66*_lambda(nodes[8]->p)/35, 10*_lambda(nodes[0]->p)/7 + 3036*_lambda(nodes[12]->p)/35 + 92*_lambda(nodes[4]->p)/7 + 66*_lambda(nodes[8]->p)/7, 48*_lambda(nodes[0]->p)/35 + 216*_lambda(nodes[12]->p)/7 + 72*_lambda(nodes[4]->p)/7 + 144*_lambda(nodes[8]->p)/35, 24*_lambda(nodes[0]->p)/35 + 144*_lambda(nodes[12]->p)/7 + 48*_lambda(nodes[4]->p)/7 + 72*_lambda(nodes[8]->p)/35, 16*_lambda(nodes[0]->p)/7 + 792*_lambda(nodes[12]->p)/7 + 120*_lambda(nodes[4]->p)/7 + 528*_lambda(nodes[8]->p)/35, 8*_lambda(nodes[0]->p)/7 + 528*_lambda(nodes[12]->p)/7 + 80*_lambda(nodes[4]->p)/7 + 264*_lambda(nodes[8]->p)/35},
    };

	_G.assign(nodeCount, std::vector<double>(nodeCount, 0));

    auto && dop = [this](int t){
        if(t % 4 == 1)
            return hx;
        if(t % 4 == 2)
            return hy;
        if(t % 4 == 3)
            return hx*hy;
        return 1.0;
    };

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_G[i][j] = dop(i)*dop(j)*((hy/hx) * __g_1[i][j] + (hx/hy) * __g_2[i][j]); // домножить на hx hy как в теории
		}
	}
	
	return _G;
}


std::vector<std::vector<double>> FE::calcM() {
	std::vector<std::vector<double>> __m = {
            {124.16326530612245, 29.448979591836736, 29.448979591836736, 6.98469387755102, 42.97959183673469, 1.5918367346938775, 10.193877551020408, 0.37755102040816324, 42.97959183673469, 10.193877551020408, 1.5918367346938775, 0.37755102040816324, 14.877551020408163, 0.5510204081632653, 0.5510204081632653, 0.02040816326530612},
            {29.448979591836736, 57.30612244897959, 6.98469387755102, 13.591836734693878, 81.9795918367347, 40.59183673469388, 19.443877551020407, 9.627551020408163, 10.193877551020408, 19.836734693877553, 0.37755102040816324, 0.7346938775510204, 28.377551020408163, 14.051020408163266, 1.0510204081632653, 0.5204081632653061},
            {29.448979591836736, 6.98469387755102, 57.30612244897959, 13.591836734693878, 10.193877551020408, 0.37755102040816324, 19.836734693877553, 0.7346938775510204, 81.9795918367347, 19.443877551020407, 40.59183673469388, 9.627551020408163, 28.377551020408163, 1.0510204081632653, 14.051020408163266, 0.5204081632653061},
            {6.98469387755102, 13.591836734693878, 13.591836734693878, 26.448979591836736, 19.443877551020407, 9.627551020408163, 37.83673469387755, 18.73469387755102, 19.443877551020407, 37.83673469387755, 9.627551020408163, 18.73469387755102, 54.12755102040816, 26.801020408163264, 26.801020408163264, 13.270408163265307},
            {42.97959183673469, 81.9795918367347, 10.193877551020408, 19.443877551020407, 124.16326530612245, 54.12244897959184, 29.448979591836736, 12.83673469387755, 14.877551020408163, 28.377551020408163, 0.5510204081632653, 1.0510204081632653, 42.97959183673469, 18.73469387755102, 1.5918367346938775, 0.6938775510204082},
            {1.5918367346938775, 40.59183673469388, 0.37755102040816324, 9.627551020408163, 54.12244897959184, 35.02040816326531, 12.83673469387755, 8.306122448979592, 0.5510204081632653, 14.051020408163266, 0.02040816326530612, 0.5204081632653061, 18.73469387755102, 12.122448979591837, 0.6938775510204082, 0.4489795918367347},
            {10.193877551020408, 19.443877551020407, 19.836734693877553, 37.83673469387755, 29.448979591836736, 12.83673469387755, 57.30612244897959, 24.979591836734695, 28.377551020408163, 54.12755102040816, 14.051020408163266, 26.801020408163264, 81.9795918367347, 35.734693877551024, 40.59183673469388, 17.693877551020407},
            {0.37755102040816324, 9.627551020408163, 0.7346938775510204, 18.73469387755102, 12.83673469387755, 8.306122448979592, 24.979591836734695, 16.163265306122447, 1.0510204081632653, 26.801020408163264, 0.5204081632653061, 13.270408163265307, 35.734693877551024, 23.122448979591837, 17.693877551020407, 11.448979591836734},
            {42.97959183673469, 10.193877551020408, 81.9795918367347, 19.443877551020407, 14.877551020408163, 0.5510204081632653, 28.377551020408163, 1.0510204081632653, 124.16326530612245, 29.448979591836736, 54.12244897959184, 12.83673469387755, 42.97959183673469, 1.5918367346938775, 18.73469387755102, 0.6938775510204082},
            {10.193877551020408, 19.836734693877553, 19.443877551020407, 37.83673469387755, 28.377551020408163, 14.051020408163266, 54.12755102040816, 26.801020408163264, 29.448979591836736, 57.30612244897959, 12.83673469387755, 24.979591836734695, 81.9795918367347, 40.59183673469388, 35.734693877551024, 17.693877551020407},
            {1.5918367346938775, 0.37755102040816324, 40.59183673469388, 9.627551020408163, 0.5510204081632653, 0.02040816326530612, 14.051020408163266, 0.5204081632653061, 54.12244897959184, 12.83673469387755, 35.02040816326531, 8.306122448979592, 18.73469387755102, 0.6938775510204082, 12.122448979591837, 0.4489795918367347},
            {0.37755102040816324, 0.7346938775510204, 9.627551020408163, 18.73469387755102, 1.0510204081632653, 0.5204081632653061, 26.801020408163264, 13.270408163265307, 12.83673469387755, 24.979591836734695, 8.306122448979592, 16.163265306122447, 35.734693877551024, 17.693877551020407, 23.122448979591837, 11.448979591836734},
            {14.877551020408163, 28.377551020408163, 28.377551020408163, 54.12755102040816, 42.97959183673469, 18.73469387755102, 81.9795918367347, 35.734693877551024, 42.97959183673469, 81.9795918367347, 18.73469387755102, 35.734693877551024, 124.16326530612245, 54.12244897959184, 54.12244897959184, 23.591836734693878},
            {0.5510204081632653, 14.051020408163266, 1.0510204081632653, 26.801020408163264, 18.73469387755102, 12.122448979591837, 35.734693877551024, 23.122448979591837, 1.5918367346938775, 40.59183673469388, 0.6938775510204082, 17.693877551020407, 54.12244897959184, 35.02040816326531, 23.591836734693878, 15.26530612244898},
            {0.5510204081632653, 1.0510204081632653, 14.051020408163266, 26.801020408163264, 1.5918367346938775, 0.6938775510204082, 40.59183673469388, 17.693877551020407, 18.73469387755102, 35.734693877551024, 12.122448979591837, 23.122448979591837, 54.12244897959184, 23.591836734693878, 35.02040816326531, 15.26530612244898},
            {0.02040816326530612, 0.5204081632653061, 0.5204081632653061, 13.270408163265307, 0.6938775510204082, 0.4489795918367347, 17.693877551020407, 11.448979591836734, 0.6938775510204082, 17.693877551020407, 0.4489795918367347, 11.448979591836734, 23.591836734693878, 15.26530612244898, 15.26530612244898, 9.877551020408163},
    };

	_M.assign(nodeCount, std::vector<double>(nodeCount, 0));
	
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_M[i][j] = approximateFunc(_gamma) * hx * hy * __m[i][j];
		}
	}

	return _M;
}

std::vector<std::vector<double>> FE::calcb() {
	std::vector<std::vector<double>> __b = {
            {6084 * _rp(nodes[0]->p) / 49 + 78 * _rp(nodes[10]->p) / 49 + 37 * _rp(nodes[11]->p) / 98 +
             729 * _rp(nodes[12]->p) / 49 + 27 * _rp(nodes[13]->p) / 49 + 27 * _rp(nodes[14]->p) / 49 +
             _rp(nodes[15]->p) / 49 + 1443 * _rp(nodes[1]->p) / 49 + 1443 * _rp(nodes[2]->p) / 49 +
             1369 * _rp(nodes[3]->p) / 196 + 2106 * _rp(nodes[4]->p) / 49 + 78 * _rp(nodes[5]->p) / 49 +
             999 * _rp(nodes[6]->p) / 98 + 37 * _rp(nodes[7]->p) / 98 + 2106 * _rp(nodes[8]->p) / 49 +
             999 * _rp(nodes[9]->p) / 98},
            {1443 * _rp(nodes[0]->p) / 49 + 37 * _rp(nodes[10]->p) / 98 + 36 * _rp(nodes[11]->p) / 49 +
             2781 * _rp(nodes[12]->p) / 98 + 1377 * _rp(nodes[13]->p) / 98 + 103 * _rp(nodes[14]->p) / 98 +
             51 * _rp(nodes[15]->p) / 98 + 2808 * _rp(nodes[1]->p) / 49 + 1369 * _rp(nodes[2]->p) / 196 +
             666 * _rp(nodes[3]->p) / 49 + 4017 * _rp(nodes[4]->p) / 49 + 1989 * _rp(nodes[5]->p) / 49 +
             3811 * _rp(nodes[6]->p) / 196 + 1887 * _rp(nodes[7]->p) / 196 + 999 * _rp(nodes[8]->p) / 98 +
             972 * _rp(nodes[9]->p) / 49},
            {1443 * _rp(nodes[0]->p) / 49 + 1989 * _rp(nodes[10]->p) / 49 + 1887 * _rp(nodes[11]->p) / 196 +
             2781 * _rp(nodes[12]->p) / 98 + 103 * _rp(nodes[13]->p) / 98 + 1377 * _rp(nodes[14]->p) / 98 +
             51 * _rp(nodes[15]->p) / 98 + 1369 * _rp(nodes[1]->p) / 196 + 2808 * _rp(nodes[2]->p) / 49 +
             666 * _rp(nodes[3]->p) / 49 + 999 * _rp(nodes[4]->p) / 98 + 37 * _rp(nodes[5]->p) / 98 +
             972 * _rp(nodes[6]->p) / 49 + 36 * _rp(nodes[7]->p) / 49 + 4017 * _rp(nodes[8]->p) / 49 +
             3811 * _rp(nodes[9]->p) / 196},
            {1369 * _rp(nodes[0]->p) / 196 + 1887 * _rp(nodes[10]->p) / 196 + 918 * _rp(nodes[11]->p) / 49 +
             10609 * _rp(nodes[12]->p) / 196 + 5253 * _rp(nodes[13]->p) / 196 + 5253 * _rp(nodes[14]->p) / 196 +
             2601 * _rp(nodes[15]->p) / 196 + 666 * _rp(nodes[1]->p) / 49 + 666 * _rp(nodes[2]->p) / 49 +
             1296 * _rp(nodes[3]->p) / 49 + 3811 * _rp(nodes[4]->p) / 196 + 1887 * _rp(nodes[5]->p) / 196 +
             1854 * _rp(nodes[6]->p) / 49 + 918 * _rp(nodes[7]->p) / 49 + 3811 * _rp(nodes[8]->p) / 196 +
             1854 * _rp(nodes[9]->p) / 49},
            {2106 * _rp(nodes[0]->p) / 49 + 27 * _rp(nodes[10]->p) / 49 + 103 * _rp(nodes[11]->p) / 98 +
             2106 * _rp(nodes[12]->p) / 49 + 918 * _rp(nodes[13]->p) / 49 + 78 * _rp(nodes[14]->p) / 49 +
             34 * _rp(nodes[15]->p) / 49 + 4017 * _rp(nodes[1]->p) / 49 + 999 * _rp(nodes[2]->p) / 98 +
             3811 * _rp(nodes[3]->p) / 196 + 6084 * _rp(nodes[4]->p) / 49 + 2652 * _rp(nodes[5]->p) / 49 +
             1443 * _rp(nodes[6]->p) / 49 + 629 * _rp(nodes[7]->p) / 49 + 729 * _rp(nodes[8]->p) / 49 +
             2781 * _rp(nodes[9]->p) / 98},
            {78 * _rp(nodes[0]->p) / 49 + _rp(nodes[10]->p) / 49 + 51 * _rp(nodes[11]->p) / 98 +
             918 * _rp(nodes[12]->p) / 49 + 594 * _rp(nodes[13]->p) / 49 + 34 * _rp(nodes[14]->p) / 49 +
             22 * _rp(nodes[15]->p) / 49 + 1989 * _rp(nodes[1]->p) / 49 + 37 * _rp(nodes[2]->p) / 98 +
             1887 * _rp(nodes[3]->p) / 196 + 2652 * _rp(nodes[4]->p) / 49 + 1716 * _rp(nodes[5]->p) / 49 +
             629 * _rp(nodes[6]->p) / 49 + 407 * _rp(nodes[7]->p) / 49 + 27 * _rp(nodes[8]->p) / 49 +
             1377 * _rp(nodes[9]->p) / 98},
            {999 * _rp(nodes[0]->p) / 98 + 1377 * _rp(nodes[10]->p) / 98 + 5253 * _rp(nodes[11]->p) / 196 +
             4017 * _rp(nodes[12]->p) / 49 + 1751 * _rp(nodes[13]->p) / 49 + 1989 * _rp(nodes[14]->p) / 49 +
             867 * _rp(nodes[15]->p) / 49 + 3811 * _rp(nodes[1]->p) / 196 + 972 * _rp(nodes[2]->p) / 49 +
             1854 * _rp(nodes[3]->p) / 49 + 1443 * _rp(nodes[4]->p) / 49 + 629 * _rp(nodes[5]->p) / 49 +
             2808 * _rp(nodes[6]->p) / 49 + 1224 * _rp(nodes[7]->p) / 49 + 2781 * _rp(nodes[8]->p) / 98 +
             10609 * _rp(nodes[9]->p) / 196},
            {37 * _rp(nodes[0]->p) / 98 + 51 * _rp(nodes[10]->p) / 98 + 2601 * _rp(nodes[11]->p) / 196 +
             1751 * _rp(nodes[12]->p) / 49 + 1133 * _rp(nodes[13]->p) / 49 + 867 * _rp(nodes[14]->p) / 49 +
             561 * _rp(nodes[15]->p) / 49 + 1887 * _rp(nodes[1]->p) / 196 + 36 * _rp(nodes[2]->p) / 49 +
             918 * _rp(nodes[3]->p) / 49 + 629 * _rp(nodes[4]->p) / 49 + 407 * _rp(nodes[5]->p) / 49 +
             1224 * _rp(nodes[6]->p) / 49 + 792 * _rp(nodes[7]->p) / 49 + 103 * _rp(nodes[8]->p) / 98 +
             5253 * _rp(nodes[9]->p) / 196},
            {2106 * _rp(nodes[0]->p) / 49 + 2652 * _rp(nodes[10]->p) / 49 + 629 * _rp(nodes[11]->p) / 49 +
             2106 * _rp(nodes[12]->p) / 49 + 78 * _rp(nodes[13]->p) / 49 + 918 * _rp(nodes[14]->p) / 49 +
             34 * _rp(nodes[15]->p) / 49 + 999 * _rp(nodes[1]->p) / 98 + 4017 * _rp(nodes[2]->p) / 49 +
             3811 * _rp(nodes[3]->p) / 196 + 729 * _rp(nodes[4]->p) / 49 + 27 * _rp(nodes[5]->p) / 49 +
             2781 * _rp(nodes[6]->p) / 98 + 103 * _rp(nodes[7]->p) / 98 + 6084 * _rp(nodes[8]->p) / 49 +
             1443 * _rp(nodes[9]->p) / 49},
            {999 * _rp(nodes[0]->p) / 98 + 629 * _rp(nodes[10]->p) / 49 + 1224 * _rp(nodes[11]->p) / 49 +
             4017 * _rp(nodes[12]->p) / 49 + 1989 * _rp(nodes[13]->p) / 49 + 1751 * _rp(nodes[14]->p) / 49 +
             867 * _rp(nodes[15]->p) / 49 + 972 * _rp(nodes[1]->p) / 49 + 3811 * _rp(nodes[2]->p) / 196 +
             1854 * _rp(nodes[3]->p) / 49 + 2781 * _rp(nodes[4]->p) / 98 + 1377 * _rp(nodes[5]->p) / 98 +
             10609 * _rp(nodes[6]->p) / 196 + 5253 * _rp(nodes[7]->p) / 196 + 1443 * _rp(nodes[8]->p) / 49 +
             2808 * _rp(nodes[9]->p) / 49},
            {78 * _rp(nodes[0]->p) / 49 + 1716 * _rp(nodes[10]->p) / 49 + 407 * _rp(nodes[11]->p) / 49 +
             918 * _rp(nodes[12]->p) / 49 + 34 * _rp(nodes[13]->p) / 49 + 594 * _rp(nodes[14]->p) / 49 +
             22 * _rp(nodes[15]->p) / 49 + 37 * _rp(nodes[1]->p) / 98 + 1989 * _rp(nodes[2]->p) / 49 +
             1887 * _rp(nodes[3]->p) / 196 + 27 * _rp(nodes[4]->p) / 49 + _rp(nodes[5]->p) / 49 +
             1377 * _rp(nodes[6]->p) / 98 + 51 * _rp(nodes[7]->p) / 98 + 2652 * _rp(nodes[8]->p) / 49 +
             629 * _rp(nodes[9]->p) / 49},
            {37 * _rp(nodes[0]->p) / 98 + 407 * _rp(nodes[10]->p) / 49 + 792 * _rp(nodes[11]->p) / 49 +
             1751 * _rp(nodes[12]->p) / 49 + 867 * _rp(nodes[13]->p) / 49 + 1133 * _rp(nodes[14]->p) / 49 +
             561 * _rp(nodes[15]->p) / 49 + 36 * _rp(nodes[1]->p) / 49 + 1887 * _rp(nodes[2]->p) / 196 +
             918 * _rp(nodes[3]->p) / 49 + 103 * _rp(nodes[4]->p) / 98 + 51 * _rp(nodes[5]->p) / 98 +
             5253 * _rp(nodes[6]->p) / 196 + 2601 * _rp(nodes[7]->p) / 196 + 629 * _rp(nodes[8]->p) / 49 +
             1224 * _rp(nodes[9]->p) / 49},
            {729 * _rp(nodes[0]->p) / 49 + 918 * _rp(nodes[10]->p) / 49 + 1751 * _rp(nodes[11]->p) / 49 +
             6084 * _rp(nodes[12]->p) / 49 + 2652 * _rp(nodes[13]->p) / 49 + 2652 * _rp(nodes[14]->p) / 49 +
             1156 * _rp(nodes[15]->p) / 49 + 2781 * _rp(nodes[1]->p) / 98 + 2781 * _rp(nodes[2]->p) / 98 +
             10609 * _rp(nodes[3]->p) / 196 + 2106 * _rp(nodes[4]->p) / 49 + 918 * _rp(nodes[5]->p) / 49 +
             4017 * _rp(nodes[6]->p) / 49 + 1751 * _rp(nodes[7]->p) / 49 + 2106 * _rp(nodes[8]->p) / 49 +
             4017 * _rp(nodes[9]->p) / 49},
            {27 * _rp(nodes[0]->p) / 49 + 34 * _rp(nodes[10]->p) / 49 + 867 * _rp(nodes[11]->p) / 49 +
             2652 * _rp(nodes[12]->p) / 49 + 1716 * _rp(nodes[13]->p) / 49 + 1156 * _rp(nodes[14]->p) / 49 +
             748 * _rp(nodes[15]->p) / 49 + 1377 * _rp(nodes[1]->p) / 98 + 103 * _rp(nodes[2]->p) / 98 +
             5253 * _rp(nodes[3]->p) / 196 + 918 * _rp(nodes[4]->p) / 49 + 594 * _rp(nodes[5]->p) / 49 +
             1751 * _rp(nodes[6]->p) / 49 + 1133 * _rp(nodes[7]->p) / 49 + 78 * _rp(nodes[8]->p) / 49 +
             1989 * _rp(nodes[9]->p) / 49},
            {27 * _rp(nodes[0]->p) / 49 + 594 * _rp(nodes[10]->p) / 49 + 1133 * _rp(nodes[11]->p) / 49 +
             2652 * _rp(nodes[12]->p) / 49 + 1156 * _rp(nodes[13]->p) / 49 + 1716 * _rp(nodes[14]->p) / 49 +
             748 * _rp(nodes[15]->p) / 49 + 103 * _rp(nodes[1]->p) / 98 + 1377 * _rp(nodes[2]->p) / 98 +
             5253 * _rp(nodes[3]->p) / 196 + 78 * _rp(nodes[4]->p) / 49 + 34 * _rp(nodes[5]->p) / 49 +
             1989 * _rp(nodes[6]->p) / 49 + 867 * _rp(nodes[7]->p) / 49 + 918 * _rp(nodes[8]->p) / 49 +
             1751 * _rp(nodes[9]->p) / 49},
            {_rp(nodes[0]->p) / 49 + 22 * _rp(nodes[10]->p) / 49 + 561 * _rp(nodes[11]->p) / 49 +
             1156 * _rp(nodes[12]->p) / 49 + 748 * _rp(nodes[13]->p) / 49 + 748 * _rp(nodes[14]->p) / 49 +
             484 * _rp(nodes[15]->p) / 49 + 51 * _rp(nodes[1]->p) / 98 + 51 * _rp(nodes[2]->p) / 98 +
             2601 * _rp(nodes[3]->p) / 196 + 34 * _rp(nodes[4]->p) / 49 + 22 * _rp(nodes[5]->p) / 49 +
             867 * _rp(nodes[6]->p) / 49 + 561 * _rp(nodes[7]->p) / 49 + 34 * _rp(nodes[8]->p) / 49 +
             867 * _rp(nodes[9]->p) / 49}
    };

	_b.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_b[i][j] = hx * hy * __b[i][j];
		}
	}

	return _b;
}

std::vector<std::vector<double>> FE::calcAs3x() {
	std::vector<std::vector<double>> __as3x = {
            {11.142857142857142, 2.642857142857143, 11.142857142857142, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {2.642857142857143, 5.142857142857143, 2.642857142857143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {11.142857142857142, 2.642857142857143, 11.142857142857142, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    };

	_As3x.assign(nodeCount, std::vector<double>(nodeCount, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_As3x[i][j] = approximateFunc(_bcond3x) * hx * __as3x[i][j];
		}
	}

	return _As3x;

}


std::vector<std::vector<double>> FE::calcAs3y() {
	std::vector<std::vector<double>> __as3y = {
            {11.142857142857142, 0.0, 0.0, 2.642857142857143, 0.0, 0.0, 2.642857142857143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {2.642857142857143, 0.0, 0.0, 5.142857142857143, 0.0, 0.0, 5.142857142857143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {2.642857142857143, 0.0, 0.0, 5.142857142857143, 0.0, 0.0, 5.142857142857143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    };

	_As3y.assign(nodeCount, std::vector<double>(nodeCount, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_As3y[i][j] = approximateFunc(_bcond3y) * hy * __as3y[i][j];
		}
	}

	return _As3y;
}


std::vector<std::vector<double>> FE::calcbs3x() {

	std::vector<std::vector<double>> __bs3x = {
            {78 * _bcond3x_rp(nodes[0]->p) / 7 + 37 * _bcond3x_rp(nodes[1]->p) / 14 +
             78 * _bcond3x_rp(nodes[2]->p) / 7},
            {37 * _bcond3x_rp(nodes[0]->p) / 14 + 36 * _bcond3x_rp(nodes[1]->p) / 7 +
             37 * _bcond3x_rp(nodes[2]->p) / 14},
            {78 * _bcond3x_rp(nodes[0]->p) / 7 + 37 * _bcond3x_rp(nodes[1]->p) / 14 +
             78 * _bcond3x_rp(nodes[2]->p) / 7},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
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
			_bs3x[i][j] = approximateFunc(_bcond3x) * hx * __bs3x[i][j];
		}
	}

	return _bs3x;
}

std::vector<std::vector<double>> FE::calcbs3y() {
	std::vector<std::vector<double>> __bs3y = {
            {78 * _bcond3x_rp(nodes[0]->p) / 7 + 37 * _bcond3x_rp(nodes[3]->p) / 14 +
             37 * _bcond3x_rp(nodes[6]->p) / 14},
            {0},
            {0},
            {37 * _bcond3x_rp(nodes[0]->p) / 14 + 36 * _bcond3x_rp(nodes[3]->p) / 7 +
             36 * _bcond3x_rp(nodes[6]->p) / 7},
            {0},
            {0},
            {37 * _bcond3x_rp(nodes[0]->p) / 14 + 36 * _bcond3x_rp(nodes[3]->p) / 7 +
             36 * _bcond3x_rp(nodes[6]->p) / 7},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0}
    };

	_bs3y.assign(nodeCount, std::vector<double>(1, 0));

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_bs3y[i][j] = approximateFunc(_bcond3y) * hy * __bs3y[i][j];
		}
	}

	return _bs3y;
}


std::vector<std::vector<double>> FE::calcbs2x() {

	std::vector<std::vector<double>> __bs2x = {
            {78 * _bcond2y_rp(nodes[0]->p) / 7 + 37 * _bcond2y_rp(nodes[1]->p) / 14 +
             78 * _bcond2y_rp(nodes[2]->p) / 7},
            {37 * _bcond2y_rp(nodes[0]->p) / 14 + 36 * _bcond2y_rp(nodes[1]->p) / 7 +
             37 * _bcond2y_rp(nodes[2]->p) / 14},
            {78 * _bcond2y_rp(nodes[0]->p) / 7 + 37 * _bcond2y_rp(nodes[1]->p) / 14 +
             78 * _bcond2y_rp(nodes[2]->p) / 7},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
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
			_bs2x[i][j] = hx * __bs2x[i][j];
		}
	}

	return _bs2x;
}

std::vector<std::vector<double>> FE::calcbs2y() {
	std::vector<std::vector<double>> __bs2y = {
            {78 * _bcond2y_rp(nodes[0]->p) / 7 + 37 * _bcond2y_rp(nodes[3]->p) / 14 +
             37 * _bcond3x_rp(nodes[6]->p) / 14},
            {0},
            {0},
            {37 * _bcond2y_rp(nodes[0]->p) / 14 + 36 * _bcond2y_rp(nodes[3]->p) / 7 +
             36 * _bcond2y_rp(nodes[6]->p) / 7},
            {0},
            {0},
            {37 * _bcond2y_rp(nodes[0]->p) / 14 + 36 * _bcond2y_rp(nodes[3]->p) / 7 +
             36 * _bcond2y_rp(nodes[6]->p) / 7},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0}
	};

	_bs2y.assign(nodeCount, std::vector<double>(1, 0));


	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < 1; j++) {
			_bs2y[i][j] = hy * __bs2y[i][j];
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


void FE::makeImpactToGlobalMatrix(std::map<int, std::vector<std::pair<int, double>>>& _m) {
	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < nodeCount; j++) {
			_m[nodes[i]->globalID].push_back({ nodes[j]->globalID, _A[i][j] });
		}
	}
}

void FE::makeImpactToGlobalRP(std::vector<double>& _rp) {
	for (int i = 0; i < nodeCount; i++) {
		_rp[nodes[i]->globalID] += _RP[i][0];
	}
}