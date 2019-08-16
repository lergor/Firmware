#include "solver.h"

namespace ratio_solver {

    std::ostream &operator<<(std::ostream &os, const params &par) {
        os << "w: {" << par.p << ", " << par.q << ", " << par.r << "}" << std::endl;
        os << "n: {" << par.nx << ", " << par.ny << ", " << par.nz << "}" << std::endl;
        os << "f: {" << par.f1 << ", " << par.f2 << ", " << par.f3 << ", " << par.f4 << "}" << std::endl;
        return os;
    }

    double calculate_f1_for_ratio(double rho, const drone_params &drone, double EPS) {
        if (math::is_equal(rho, 2.))
            return drone.Fmax * 10000;

        if (math::is_equal(rho, 0.))
            return drone.MASS * drone.accel() / 2;

        const double A = rho * drone.L;
        const double B = drone.Izz_T - drone.Ixx_T;
        const double C = drone.Izz_P * (2 + std::sqrt(rho)) / std::sqrt(drone.Kf);
        const double M = drone.MASS * drone.accel() / (2 + rho);
        const double R = drone.Kt * (2 - rho) / drone.GAMMA;

        const double Z = B * R * R / M;
        const double X = C * R / M;
        const double M2 = M * M;

        const double Z2 = Z * Z;
        const double XZ = C / (B * R);
        const double XZ2 = XZ * XZ;

        std::vector<double> coeffs(9, 0);
        coeffs[0] = -A * A / Z2;
        coeffs[2] = -XZ2 * M2;
        coeffs[3] = -2 * XZ * M2;
        coeffs[4] = -M2;
        coeffs[6] = XZ2;
        coeffs[7] = 2 * XZ;
        coeffs[8] = 1;

        math::polynomial<double> poly(coeffs);
        std::vector<double> roots = poly.solve(EPS);

        std::function<double(double)> equation = [=](double sqrt_f) {
            const double tmp = std::sqrt(std::pow(sqrt_f, 4) - M * M);
            return A - Z * tmp * sqrt_f * sqrt_f - X * tmp * sqrt_f;
        };

        if (roots.empty())
            return drone.Fmax * 10000;

        for (auto &root : roots) {
            if (math::is_less(0.0, root)) {
                double zero = poly.evaluate(root);
                if (math::is_equal(zero, 0., EPS) && math::is_equal(equation(root), 0., EPS)) {
                    return root * root;
                }
            }
        }
        return 0.;
    }

    params calculate_params_for_ratio(double rho, const drone_params &drone, double EPS) {
        params result;
        const double M = drone.MASS * drone.accel() / (2 + rho);
        const double R = drone.Kt * (2 - rho) / drone.GAMMA;

        result.f1 = result.f3 = calculate_f1_for_ratio(rho, drone, EPS);

        if (math::is_equal(result.f1, 0., EPS))
            return result;

        result.f2 = rho * result.f1;
        result.nz = M / result.f1;
        result.ny = std::sqrt(1 - result.nz * result.nz);
        result.r = R * result.f1;
        result.q = result.ny * result.r / result.nz;

        return result;
    }

    double calculate_optimal_ratio(const drone_params &drone, double accuracy, double EPS) {
        const double rho_min = 0.4;
        const double rho_max = 2 - 2 * accuracy; // not the case of two broken propellers
        const double step = accuracy;

        double w_min = 1E09;
        double optimal_rho_w = 0;
        double f_min = 1E09;
        double optimal_rho_f = 0;

        double ratio = rho_min;
        while (math::is_less(rho_max, ratio)) {
            auto params = calculate_params_for_ratio(ratio, drone, EPS);
            if (math::is_less(params.f1, f_min)) {
                if (math::is_less(drone.Fmin, params.f1) && math::is_less(params.f1, drone.Fmax)) {
                    f_min = params.f1;
                    optimal_rho_f = ratio;
                } else {
                    continue;
                }
            }

            double w_b = std::sqrt(params.q * params.q + params.r * params.r);
            if (math::is_less(w_b, w_min)) {
                w_min = w_b;
                optimal_rho_w = ratio;
            }
            ratio += step;
        }

        if (std::abs(optimal_rho_f - optimal_rho_w) <= 0.2 * rho_max) {
            // little trick : take the value away from  the saturation
            return 0.7 * (optimal_rho_f + optimal_rho_w) / 2;
        }
        // TODO
        return optimal_rho_f;
    }

//    void plot_dependencies(const drone_params &drone, double rho_max, double EPS) {
//        const double step = 0.1;
//        const int iters = static_cast<int>(rho_max / step);
//        std::vector<double> rho(iters, 0);
//        std::vector<double> zeros(iters, 0);
//        std::vector<double> f1(iters, 0);
//        std::vector<double> f2(iters, 0);
//        std::vector<double> ny(iters, 0);
//        std::vector<double> nz(iters, 0);
//        std::vector<double> q(iters, 0);
//        std::vector<double> r(iters, 0);
//        std::vector<double> w_len(iters, 0);
//
//        double ratio = 0;
//        for (int i = 0; i <= iters; ++i) {
//            auto params = ratio_solver::calculate_params_for_ratio(ratio, drone, EPS);
//            rho[i] = ratio;
//            f1[i] = params.f1;
//            f2[i] = params.f2;
//            ny[i] = params.ny;
//            nz[i] = params.nz;
//            q[i] = params.q;
//            r[i] = params.r;
//            w_len[i] = std::sqrt(q[i] * q[i] + r[i] * r[i]);
//            ratio += step;
//        }
//
//        int x_max = static_cast<int>(rho_max + 0.5);
//        {
//            plt::figure();
//            plt::xlim(0, x_max);
//            plt::ylim(0, 6);
//            plt::named_plot("f1, f3", rho, f1);
//            plt::named_plot("f2", rho, f2);
//            plt::grid(true);
//            plt::xlabel("rho");
//            plt::show();
//        }
//
//        {
//            plt::figure();
//            plt::xlim(0, x_max);
//            plt::ylim(0, 1);
//            plt::named_plot("ny", rho, ny);
//            plt::named_plot("nz", rho, nz);
//            plt::grid(true);
//            plt::xlabel("rho");
//            plt::show();
//        }
//
//        {
//            plt::figure();
//            plt::xlim(0, x_max);
//            plt::ylim(-30, 30);
//            plt::named_plot("q", rho, q);
//            plt::named_plot("r", rho, r);
//            plt::named_plot("||w_b||", rho, w_len);
//            plt::grid(true);
//            plt::xlabel("rho");
//            plt::show();
//        }
//    }

} // namespace ratio_solver
