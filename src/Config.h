#pragma once
#include <vector>
#include <functional>
#include <Point.h>



class Config
{

private:
	static std::vector<std::function<double(Point)>> BCOND1_FUNCS;
	static std::vector<std::function<double(Point)>> BCOND2_FUNCS;
	static std::vector<std::function<double(Point)>> BCOND3_FUNCS;
	static std::vector<std::function<double(Point)>> BCOND3_RP_FUNCS;
	static std::vector<std::function<double(Point)>> RP_FUNCS;
	static std::vector<std::function<double(Point)>> LAMBDA_FUNCS;
	static std::vector<std::function<double(Point)>> GAMMA_FUNCS;
	static std::function<double(Point)> ZERO_FUNC;

public: 

	static std::function<double(Point)> GET_BCOND1_FUNC(size_t _B1_IDX) {
		return BCOND1_FUNCS[_B1_IDX];
	}



	static std::function<double(Point)> GET_BCOND2_FUNC(size_t _B2_IDX) {
		return BCOND2_FUNCS[_B2_IDX];
	}


	static std::function<double(Point)> GET_BCOND3_FUNC(size_t _B3_IDX) {
		return BCOND3_FUNCS[_B3_IDX];
	}

	static std::function<double(Point)> GET_BCOND3_RP_FUNC(size_t _B3_IDX) {
		return BCOND3_RP_FUNCS[_B3_IDX];
	}


	static std::function<double(Point)> GET_RP_FUNC(size_t _RP_IDX) {
		return RP_FUNCS[_RP_IDX];
	}


	static std::function<double(Point)> GET_LAMBDA_FUNC(size_t _LAMBDA_IDX) {
		return LAMBDA_FUNCS[_LAMBDA_IDX];
	}


	static std::function<double(Point)> GET_GAMMA_FUNC(size_t _GAMMA_IDX) {
		return GAMMA_FUNCS[_GAMMA_IDX];
	}

	static std::function<double(Point)> GET_ZERO_FUNC() {
		return ZERO_FUNC;
	}
};

