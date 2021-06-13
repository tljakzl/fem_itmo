#pragma once
#include <vector>
#include "DMat.h"

class MSG {
public:
	static std::vector<double> solve(DMat<double>& A, std::vector<double> f, double eps, int maxiter);
};