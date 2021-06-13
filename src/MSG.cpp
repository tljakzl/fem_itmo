#include "MSG.h"
#include <iostream>


std::vector<double> plus(const std::vector<double>& v1, const std::vector<double>& v2) {
	std::vector<double> ans(v1.size());
	for (int i = 1; i < v1.size(); i++) {
		ans[i] = v1[i] + v2[i];
	}
	return ans;
}

std::vector<double> minus(const std::vector<double>& v1, const std::vector<double>& v2) {
	std::vector<double> ans(v1.size());
	for (int i = 1; i < v1.size(); i++) {
		ans[i] = v1[i] - v2[i];
	}
	return ans;
}

std::vector<double> times(double a, const std::vector<double>& v) {
	std::vector<double> ans(v.size());
	for (int i = 1; i < v.size(); i++) {
		ans[i] = a * v[i];
	}
	return ans;
}

double dot(const std::vector<double>& v1, const std::vector<double>& v2) {
	double ans = 0;
	for (int i = 1; i < v1.size(); i++) {
		ans += v1[i] * v2[i];
	}
	return ans;
}


double norm(const std::vector<double>& v1) {
	double ans = 0;
	for (int i = 1; i < v1.size(); i++) {
		ans += v1[i] * v1[i];
	}
	return sqrt(ans);
}


std::vector<double> MSG::solve(DMat<double>& A, std::vector<double> f, double eps, int maxiter) {
	int n = f.size();
	
	DMat<double> M = DMat<double>(n - 1);
	for (int i = 1; i < n; i++) {
		M.setDiag(i, 1 / A.get(i, i));
	}

	std::vector<double> px(n, 1);
	std::vector<double> pr = minus(f, A.product(px));
	std::vector<double> pz = M.product(pr);

	for (int k = 1; k <= maxiter; k++) {
		std::cout << "Iter " << k << std::endl;
		double ak = dot(M.product(pr), pr) / dot(A.product(pz), pz);
		std::vector<double> x = plus(px, times(ak, pz));
		std::vector<double> r = minus(pr, times(ak, A.product(pz)));
		double bk = dot(M.product(r), r) / dot(M.product(pr), pr);
		std::vector<double> z = plus(M.product(r), times(bk, pz));

		px = x;
		pr = r;
		pz = z;

		double npr = norm(pr);
		//double nf = norm(f);
		// будем считать абсолютную погрешность невязки вместо относительной, 
		// потому что норма вектора f очень велика из-за учета граничных условий первого рода
		// еще одно возможное решение - считать норму Чебышева
		// еще одно возможное решение - на больших матрицах таки использовать евклидову, но с очень низким эпсилон
		if (npr < eps) return px;
	}

	return px;
}


