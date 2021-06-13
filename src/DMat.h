#pragma once
#include <vector>
#include <map>

template <typename T>
class DMat {
public:
	int dim;
	std::vector<T> di;

	std::vector<int> l_ia;
	std::vector<int> l_ja;
	std::vector<T> al;

	std::vector<int> u_ia;
	std::vector<int> u_ja;
	std::vector<T> au;
	

	DMat(int _dim) {
		dim = _dim;
		di.resize(dim + 1);

		for (int i = 0; i <= _dim; i++) {
			l_ia.push_back(1);
			u_ia.push_back(1);
		}
	}

	DMat(std::map<int, std::vector<std::pair<int, T>>> nze, int _dim) {
		dim = _dim;
		di.resize(dim + 1);

		std::map<int, std::vector<std::pair<int, T>>> c_wise;

		l_ia.push_back(1);
		for (int i = 1; i <= _dim; i++) {
			std::vector<std::pair<int, T>> _tmp = nze[i];
			sort(_tmp.begin(), _tmp.end());
			int m = _tmp.size();
			int cnt = 0;
			int prev_c = -1;
			for (int j = 0; j < m; j++) {
				if (_tmp[j].first > i) {
					c_wise[_tmp[j].first].push_back({ i, _tmp[j].second });
					continue;
				}
				if (_tmp[j].first == i) {
					di[i] = _tmp[j].second;
					continue;
				}

				if (_tmp[j].first != prev_c) {
					l_ja.push_back(_tmp[j].first);
					al.push_back(_tmp[j].second);
					prev_c = _tmp[j].first;
					cnt++;
				}
				else {
					al.back() += _tmp[j].second;
				}
			}
			l_ia.push_back(l_ia.back() + cnt);
		}



		u_ia.push_back(1);
		u_ia.push_back(1);
		for (int i = 2; i <= _dim; i++) {
			std::vector<std::pair<int, T>> _tmp = c_wise[i];
			sort(_tmp.begin(), _tmp.end());
			int m = _tmp.size();
			int cnt = 0;
			int prev_r = -1;
			for (int j = 0; j < m; j++) {
				if (_tmp[j].first != prev_r) {
					u_ja.push_back(_tmp[j].first);
					au.push_back(_tmp[j].second);
					prev_r = _tmp[j].first;
					cnt++;
				}
				else {
					au.back() += _tmp[j].second;
				}
			}
			u_ia.push_back(u_ia.back() + cnt);
		}

	}

	T get(int i, int j) {
		if (i == j) return di[i];

		auto _get = [](std::vector<int>& ia, std::vector<int>& ja, std::vector<T>& a, int i, int j) -> T {
			int first = ia[i - 1] - 1;
			int cnt = ia[i] - ia[i - 1];
			if (!cnt) return (T)0;
			std::vector<int>::iterator beg = ja.begin() + first;
			std::vector<int>::iterator en = ja.begin() + first + cnt - 1;
			auto lb = lower_bound(beg, en, j);
			if (*lb == j) {
				int _idx = lb - ja.begin();
				return a[_idx];
			}
			else {
				return (T)0;
			}
		};

		if (i > j) {
			return _get(l_ia, l_ja, al, i, j);
		}
		else {
			return _get(u_ia, u_ja, au, j, i);
		}
	}

	void setDiag(int d, T val) {
		di[d] = val;
	}

	void print() {
		for (int i = 1; i <= dim; i++) {
			for (int j = 1; j <= dim; j++) {
				std::cout << get(i, j) << " ";
			}
			std::cout << std::endl;
		}
	}

	// 1-based
	std::vector<T> product(std::vector<T> &v) {
		std::vector<T> ans(v.size());
		
		for (int row = 1; row < v.size(); row++) {
			// диагональ
			ans[row] += di[row] * v[row];

			// кусок строки в НТМ
			int first = l_ia[row - 1] - 1;
			int cnt = l_ia[row] - l_ia[row - 1];
			for (int j = 0; j < cnt; j++) {
				int pos = first + j;
				ans[row] += al[pos] * v[l_ja[pos]];
			}
		}

		// по столбцам ВТМ
		for (int col = 1; col < v.size(); col++) {
			int first = u_ia[col - 1] - 1;
			int cnt = u_ia[col] - u_ia[col - 1];
			for (int j = 0; j < cnt; j++) {
				int pos = first + j;
				ans[u_ja[pos]] += au[pos] * v[col];
			}
		}
		return ans;
	}
};