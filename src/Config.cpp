#include "Config.h"

std::vector<std::function<double(Point)>> Config::BCOND1_FUNCS = {
	[](Point p) {
		return p.x + p.y + 1;
	}
};

std::vector<std::function<double(Point)>> Config::BCOND2_FUNCS = {
	[](Point p) {
		return p.x + p.y;
	}
};

std::vector<std::function<double(Point)>> Config::BCOND3_FUNCS = {
	[](Point p) {
		return p.x + p.y;
	}
};

std::vector<std::function<double(Point)>> Config::BCOND3_RP_FUNCS = {
	[](Point p) {
		return p.x + p.y;
	}
};

std::vector<std::function<double(Point)>> Config::RP_FUNCS = {
	[](Point p) {
		return p.x + p.y;
	}
};

std::vector<std::function<double(Point)>> Config::LAMBDA_FUNCS = {
	[](Point p) {
		return p.x + p.y;
	}
};

std::vector<std::function<double(Point)>> Config::GAMMA_FUNCS = {
	[](Point p) {
		return p.x + p.y;
	}
};

std::function<double(Point)> Config::ZERO_FUNC = {
	[](Point p) {
		return 0;
	}
};