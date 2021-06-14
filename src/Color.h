#pragma once
#include <algorithm>
#include <cmath>

struct Color
{
    constexpr Color(double r, double g, double b) noexcept : data{ r, g, b } {}

    double data[3];

    constexpr double r() const noexcept { return data[0]; }
    constexpr double g() const noexcept { return data[1]; }
    constexpr double b() const noexcept { return data[2]; }

    double& operator[](std::size_t n) noexcept { return data[n]; }
    constexpr double operator[](std::size_t n) const noexcept { return data[n]; }
    double& operator()(std::size_t n) noexcept { return data[n]; }
    constexpr double operator()(std::size_t n) const noexcept { return data[n]; }

    friend constexpr Color operator+(const Color& c0, const Color& c1) noexcept
    {
        return { c0.r() + c1.r(), c0.g() + c1.g(), c0.b() + c1.b() };
    }

    friend constexpr Color operator*(double s, const Color& c) noexcept
    {
        return { s * c.r(), s * c.g(), s * c.b() };
    }
};

namespace internal {
    inline constexpr double Clamp01(double x) noexcept {
        return (x < 0.0) ? 0.0 : (x > 1.0) ? 1.0 : x;
    }

    template<std::size_t N>
    Color CalcLerp(double x, const Color (&data)[N]) {
        const double a = Clamp01(x) * (N - 1);
        const double i = std::floor(a);
        const double t = a - i;
        const Color &c0 = data[static_cast<std::size_t>(i)];
        const Color &c1 = data[static_cast<std::size_t>(std::ceil(a))];

        return (1.0 - t) * c0 + t * c1;
    }
}

Color GetParulaColor(double x);
Color GetColor(double x);
