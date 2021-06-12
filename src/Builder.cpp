#include "Builder.h"
#include <fstream>



Grid* Builder::createGrid() {
	std::vector<Point> points;
	std::map<int, std::vector<int>> rectangles;
	std::map<int, int> bcond1, rp, gamma, lambda;
	std::map<std::pair<int, int>, int> bcond2, bcond3;
	double x, y;
	int k, k1, n1, n2, n3, n4, n;
	std::ifstream ifs;

	ifs.open("xy.txt");
	while (ifs >> x >> y) {
		Point* p = new Point(x, y);
		points.push_back(*p);
	}
	ifs.close();


	ifs.open("rectangles.txt");
	while (ifs >> k >> n1 >> n2 >> n3 >> n4) {
		rectangles[k] = { n1, n2, n3, n4 };
	}
	ifs.close();

	ifs.open("bcond1.txt");
	while (ifs >> k >> n) {
		bcond1[k] = n;
	}
	ifs.close();

	ifs.open("rp.txt");
	while (ifs >> k >> n) {
		rp[k] = n;
	}
	ifs.close();

	ifs.open("lambda.txt");
	while (ifs >> k >> n) {
		lambda[k] = n;
	}
	ifs.close();

	ifs.open("gamma.txt");
	while (ifs >> k >> n) {
		gamma[k] = n;
	}
	ifs.close();

	ifs.open("bcond2.txt");
	while (ifs >> k >> k1 >> n) {
		bcond2[{k, k1}] = n;
	}
	ifs.close();

	ifs.open("bcond3.txt");
	while (ifs >> k >> k1 >> n) {
		bcond3[{k, k1}] = n;
	}
	ifs.close();

	Grid* grid = new Grid(points, rectangles, bcond1, rp, lambda, gamma, bcond2, bcond3);
	return grid;
}