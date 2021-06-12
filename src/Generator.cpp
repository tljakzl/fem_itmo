#include "Generator.h"
#include <fstream>


/// <summary>
/// Сгенерируем прямоугольную сетку
/// </summary>
void Generator::generateGrid(int n, int m, double dx, double dy) {
	std::ofstream xy;
	xy.open("xy.txt");
	double start_x = 0, start_y = 0;
	double x = start_x, y = start_y;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			xy << x << " " << y << '\n';
			x += dx;
		}
		y += dy;
		x = start_x;
	}
	xy.close();

	xy.open("rectangles.txt");
	int el = 1;
	for (int i = 1; i <= n - 1; i++) {
		for (int j = 1; j <= m - 1; j++) {
			xy << el++ << " " << ((i - 1) * m) + j << " " << ((i - 1) * m) + j + 1 << " " << (i * m) + j << " " << (i * m) + j + 1 << '\n';
		}
	}
	xy.close();

	xy.open("bcond1.txt");
	for (int i = 1; i <= m; i++) {
		xy << i << " " << 0 << '\n';
	}
	xy.close();

	xy.open("bcond2.txt");
	for (int i = 1; i <= n - 1; i++) {
		xy << (i - 1) * m + 1 << " " << i * m + 1 << " " << 0 << '\n';
	}
	xy.close();

	xy.open("bcond3.txt");
	for (int i = 1; i <= n - 1; i++) {
		xy << (i - 1) * m + m << " " << i * m + m << " " << 0 << '\n';
	}
	xy.close();

	xy.open("rp.txt");
	for (int i = 1; i < el; i++) {
		xy << i << " " << 0 << '\n';
	}
	xy.close();

	xy.open("labmda.txt");
	for (int i = 1; i < el; i++) {
		xy << i << " " << 0 << '\n';
	}
	xy.close();

	xy.open("gamma.txt");
	for (int i = 1; i < el; i++) {
		xy << i << " " << 0 << '\n';
	}
	xy.close();

}