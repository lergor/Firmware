#pragma once

#include <iostream>
#include <cmath>
#include <tuple>
#include <cassert>
#include <functional>

#include "polynomial.h"

#include <matrix/matrix/math.hpp>

namespace ratio_solver {

    static const double G = 9.80665; // m / s 2

    struct drone_params {
        double Kf = 6.41E-6; // Ns 2 /rad 2
        double MASS = 0.5; // kg
        double Kt = 1.69E-2;// Nm/ N
        double GAMMA = 2.75E-3; // N m s / rad
        double L = 0.17; // m
        double Ixx_T = 3.2E-3; // kg m 2
        double Izz_T = 5.5E-3; // kg m 2
        double Ixx_P = 0; // kg m 2
        double Izz_P = 1.5E-9; // kg mm 2
        double Fmin = 0.2; // N
        double Fmax = 3.8; // N

        matrix::Vector3f acc{0, 0, 0};
        matrix::Vector3f g{0, 0, -9.80665};

        inline double accel() const {
            return static_cast<double>((acc - g).norm());
        }
    };

    struct params {

        double p = 0, q = 0, r = 0;
        double nx = 0, ny = 0, nz = 0;
        double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

        inline double f_sum() const {
            return f1 + f2 + f3 + f4;
        }

        inline double w_sum(double Kf) const {
            return (std::sqrt(f1) + std::sqrt(f2) + std::sqrt(f3) + std::sqrt(f4)) / std::sqrt(Kf);
        }

        friend std::ostream&operator<<(std::ostream& os, const params& par);
    };

    params calculate_params_for_ratio(double rho, const drone_params &drone = drone_params(), double EPS = 1E-06);

    double calculate_f1_for_ratio(double rho, const drone_params &drone, double EPS = 1E-06);

    double calculate_optimal_ratio(const drone_params &drone = drone_params(), double accuracy = 0.001, double EPS = 1E-06);

    void plot_dependencies(const drone_params &drone = drone_params(), double rho_max = 10., double EPS = 1E-06);

} // namespace ratio_solver
